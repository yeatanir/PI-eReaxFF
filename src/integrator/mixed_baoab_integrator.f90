module ereaxff_mixed_baoab_integrator
  use ereaxff_kinds, only: dp
  use ereaxff_units, only: amu_angstrom2_per_fs2_to_eV
  use ereaxff_ring_polymer, only: ring_polymer_state
  use ereaxff_normal_modes, only: bead_to_nm, nm_to_bead
  use ereaxff_free_ring_propagator, only: exact_free_propagation
  use ereaxff_charge_coupling_policy, only: charge_memory_state
  use ereaxff_bead_force_manager, only: oracle_bead_summary, evaluate_oracle_beads
  use ereaxff_rng, only: rng_state
  use ereaxff_thermostat_pile_l, only: pile_l_thermostat_state, pile_l_step
  use ereaxff_oracle_adapter, only: oracle_adapter_state
  use ereaxff_system_state, only: system_state, wrap_position_into_cell
  implicit none
  private

  public :: mixed_baoab_quantum_step

contains

  subroutine mixed_baoab_quantum_step(positions_A, q_nm, v_nm, quantum_indices, quantum_masses_amu, &
      rp, quantum_bead_counts, dt_fs, adapter, charge_memory, solver_guess_by_bead, charge_policy, summary, &
      quantum_thermostats, rng, state_in, step_index)
    real(dp), intent(inout) :: positions_A(:, :, :)
    real(dp), intent(inout) :: q_nm(:, :, :)
    real(dp), intent(inout) :: v_nm(:, :, :)
    integer, intent(in) :: quantum_indices(:)
    real(dp), intent(in) :: quantum_masses_amu(:), dt_fs
    type(ring_polymer_state), intent(in) :: rp(:)
    integer, intent(in) :: quantum_bead_counts(:)
    type(oracle_adapter_state), intent(inout) :: adapter
    type(charge_memory_state), intent(inout) :: charge_memory
    real(dp), intent(inout) :: solver_guess_by_bead(:, :)
    character(len=*), intent(in) :: charge_policy
    type(oracle_bead_summary), intent(inout) :: summary
    type(pile_l_thermostat_state), intent(inout) :: quantum_thermostats(:)
    type(rng_state), intent(inout) :: rng
    type(system_state), intent(in) :: state_in
    integer, intent(in), optional :: step_index
    integer :: q
    integer :: step_value
    if (size(quantum_indices) < 1) error stop "mixed_baoab_quantum_step: quantum index list is empty"
    if (size(quantum_indices) /= size(quantum_masses_amu)) error stop "mixed_baoab_quantum_step: mass count mismatch"
    if (size(q_nm, 1) /= size(positions_A, 1)) error stop "mixed_baoab_quantum_step: dimension mismatch"
    if (size(q_nm, 2) /= size(quantum_indices)) error stop "mixed_baoab_quantum_step: quantum count mismatch"
    if (size(v_nm) /= size(q_nm)) error stop "mixed_baoab_quantum_step: velocity mismatch"
    if (size(q_nm, 3) /= size(positions_A, 3)) error stop "mixed_baoab_quantum_step: bead mismatch"
    if (size(rp) /= size(quantum_indices)) error stop "mixed_baoab_quantum_step: ring-polymer count mismatch"
    if (size(quantum_bead_counts) /= size(quantum_indices)) error stop "mixed_baoab_quantum_step: bead-count mismatch"
    if (size(quantum_thermostats) /= size(quantum_indices)) error stop "mixed_baoab_quantum_step: thermostat count mismatch"
    if (dt_fs <= 0.0_dp) error stop "mixed_baoab_quantum_step: timestep must be positive"
    if (any(quantum_masses_amu <= 0.0_dp)) error stop "mixed_baoab_quantum_step: masses must be positive"
    if (any(quantum_bead_counts < 1)) error stop "mixed_baoab_quantum_step: bead counts must be positive"
    if (any(quantum_bead_counts > size(q_nm, 3))) error stop "mixed_baoab_quantum_step: bead count exceeds storage"

    call update_quantum_positions(positions_A, q_nm, rp, quantum_indices, quantum_bead_counts, state_in)

    if (present(step_index)) then
      step_value = step_index
    else
      step_value = -1
    end if

    call evaluate_oracle_beads(adapter, positions_A, charge_memory, solver_guess_by_bead, &
        charge_policy, summary, step_value, "quantum_pre_force")

    do q = 1, size(quantum_indices)
      call half_kick_from_summary(summary, quantum_indices(q), quantum_masses_amu(q), dt_fs, rp(q), &
          quantum_bead_counts(q), v_nm(:, q, :))
    end do

    do q = 1, size(quantum_indices)
      call exact_free_propagation(q_nm(:, q, 1:quantum_bead_counts(q)), &
          v_nm(:, q, 1:quantum_bead_counts(q)), 0.5_dp * dt_fs, rp(q)%omega_k)
      if (allocated(quantum_thermostats(q)%c1)) then
        call pile_l_step(quantum_thermostats(q), v_nm(:, q, 1:quantum_bead_counts(q)), &
            [quantum_masses_amu(q), quantum_masses_amu(q), quantum_masses_amu(q)], rng)
      end if
      call exact_free_propagation(q_nm(:, q, 1:quantum_bead_counts(q)), &
          v_nm(:, q, 1:quantum_bead_counts(q)), 0.5_dp * dt_fs, rp(q)%omega_k)
    end do

    call update_quantum_positions(positions_A, q_nm, rp, quantum_indices, quantum_bead_counts, state_in)

    call evaluate_oracle_beads(adapter, positions_A, charge_memory, solver_guess_by_bead, &
        charge_policy, summary, step_value, "quantum_post_force")

    do q = 1, size(quantum_indices)
      call half_kick_from_summary(summary, quantum_indices(q), quantum_masses_amu(q), dt_fs, rp(q), &
          quantum_bead_counts(q), v_nm(:, q, :))
    end do
  end subroutine mixed_baoab_quantum_step

  subroutine update_quantum_positions(positions_A, q_nm, rp, quantum_indices, quantum_bead_counts, state_in)
    real(dp), intent(inout) :: positions_A(:, :, :)
    real(dp), intent(in) :: q_nm(:, :, :)
    type(ring_polymer_state), intent(in) :: rp(:)
    integer, intent(in) :: quantum_indices(:)
    integer, intent(in) :: quantum_bead_counts(:)
    type(system_state), intent(in) :: state_in
    real(dp) :: q_bead(size(q_nm, 3))
    real(dp) :: centroid(size(q_nm, 1))
    integer :: coord
    integer :: q
    integer :: bead

    do q = 1, size(quantum_indices)
      if (rp(q)%beads /= quantum_bead_counts(q)) error stop "update_quantum_positions: ring-polymer bead mismatch"
      centroid = 0.0_dp
      do coord = 1, size(q_nm, 1)
        call nm_to_bead(rp(q)%cmat, q_nm(coord, q, 1:quantum_bead_counts(q)), q_bead(1:quantum_bead_counts(q)))
        positions_A(coord, quantum_indices(q), 1:quantum_bead_counts(q)) = q_bead(1:quantum_bead_counts(q))
        centroid(coord) = sum(q_bead(1:quantum_bead_counts(q))) / real(quantum_bead_counts(q), dp)
      end do
      do bead = 1, quantum_bead_counts(q)
        call wrap_position_into_cell(state_in, positions_A(:, quantum_indices(q), bead))
      end do
      call wrap_position_into_cell(state_in, centroid)
      if (quantum_bead_counts(q) < size(positions_A, 3)) then
        positions_A(:, quantum_indices(q), quantum_bead_counts(q) + 1:size(positions_A, 3)) = &
            spread(centroid, dim=2, ncopies=size(positions_A, 3) - quantum_bead_counts(q))
      end if
    end do
  end subroutine update_quantum_positions

  subroutine half_kick_from_summary(summary, quantum_index, quantum_mass_amu, dt_fs, rp, active_beads, v_nm)
    type(oracle_bead_summary), intent(in) :: summary
    integer, intent(in) :: quantum_index
    real(dp), intent(in) :: quantum_mass_amu, dt_fs
    type(ring_polymer_state), intent(in) :: rp
    integer, intent(in) :: active_beads
    real(dp), intent(inout) :: v_nm(:, :)
    real(dp) :: force_bead(size(v_nm, 1), active_beads)
    real(dp) :: force_nm(size(v_nm, 1), active_beads)
    real(dp) :: force_rescale
    integer :: coord

    if (active_beads /= rp%beads) error stop "half_kick_from_summary: ring-polymer bead mismatch"
    if (active_beads > size(v_nm, 2)) error stop "half_kick_from_summary: active beads exceed storage"
    if (active_beads > size(summary%scaled_forces_eV_per_A, 3)) error stop "half_kick_from_summary: force bead mismatch"

    force_rescale = real(size(summary%scaled_forces_eV_per_A, 3), dp) / real(active_beads, dp)
    force_bead = summary%scaled_forces_eV_per_A(:, quantum_index, 1:active_beads) * force_rescale
    do coord = 1, size(force_bead, 1)
      call bead_to_nm(rp%cmat, force_bead(coord, :), force_nm(coord, :))
      v_nm(coord, 1:active_beads) = v_nm(coord, 1:active_beads) + 0.5_dp * dt_fs * force_nm(coord, :) / &
          (quantum_mass_amu * amu_angstrom2_per_fs2_to_eV)
    end do
  end subroutine half_kick_from_summary

end module ereaxff_mixed_baoab_integrator
