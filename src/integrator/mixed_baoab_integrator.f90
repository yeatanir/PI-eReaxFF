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
  implicit none
  private

  public :: mixed_baoab_electron_step

contains

  subroutine mixed_baoab_electron_step(positions_A, q_nm, v_nm, electron_index, electron_mass_amu, &
      rp, dt_fs, adapter, charge_memory, solver_guess_by_bead, charge_policy, summary, &
      electron_thermostat, rng)
    real(dp), intent(inout) :: positions_A(:, :, :)
    real(dp), intent(inout) :: q_nm(:, :)
    real(dp), intent(inout) :: v_nm(:, :)
    integer, intent(in) :: electron_index
    real(dp), intent(in) :: electron_mass_amu, dt_fs
    type(ring_polymer_state), intent(in) :: rp
    type(oracle_adapter_state), intent(inout) :: adapter
    type(charge_memory_state), intent(inout) :: charge_memory
    real(dp), intent(inout) :: solver_guess_by_bead(:, :)
    character(len=*), intent(in) :: charge_policy
    type(oracle_bead_summary), intent(inout) :: summary
    type(pile_l_thermostat_state), intent(inout) :: electron_thermostat
    type(rng_state), intent(inout) :: rng
    real(dp) :: force_bead(size(q_nm, 1), size(q_nm, 2))
    real(dp) :: force_nm(size(q_nm, 1), size(q_nm, 2))
    if (electron_index < 1 .or. electron_index > size(positions_A, 2)) then
      error stop "mixed_baoab_electron_step: electron index out of range"
    end if
    if (size(q_nm, 1) /= size(positions_A, 1)) error stop "mixed_baoab_electron_step: dimension mismatch"
    if (size(q_nm, 2) /= size(positions_A, 3)) error stop "mixed_baoab_electron_step: bead mismatch"
    if (size(v_nm) /= size(q_nm)) error stop "mixed_baoab_electron_step: velocity mismatch"
    if (electron_mass_amu <= 0.0_dp) error stop "mixed_baoab_electron_step: mass must be positive"
    if (dt_fs <= 0.0_dp) error stop "mixed_baoab_electron_step: timestep must be positive"

    call update_electron_positions(positions_A, q_nm, rp)

    call evaluate_oracle_beads(adapter, positions_A, charge_memory, solver_guess_by_bead, &
        charge_policy, summary)

    call half_kick_from_summary(summary, electron_index, electron_mass_amu, dt_fs, rp, force_bead, force_nm, v_nm)

    call exact_free_propagation(q_nm, v_nm, 0.5_dp * dt_fs, rp%omega_k)
    if (allocated(electron_thermostat%c1)) then
      call pile_l_step(electron_thermostat, v_nm, [electron_mass_amu, electron_mass_amu, electron_mass_amu], rng)
    end if
    call exact_free_propagation(q_nm, v_nm, 0.5_dp * dt_fs, rp%omega_k)

    call update_electron_positions(positions_A, q_nm, rp)

    call evaluate_oracle_beads(adapter, positions_A, charge_memory, solver_guess_by_bead, &
        charge_policy, summary)

    call half_kick_from_summary(summary, electron_index, electron_mass_amu, dt_fs, rp, force_bead, force_nm, v_nm)
  end subroutine mixed_baoab_electron_step

  subroutine update_electron_positions(positions_A, q_nm, rp)
    real(dp), intent(inout) :: positions_A(:, :, :)
    real(dp), intent(in) :: q_nm(:, :)
    type(ring_polymer_state), intent(in) :: rp
    real(dp) :: q_bead(size(q_nm, 2))
    integer :: coord

    do coord = 1, size(q_nm, 1)
      call nm_to_bead(rp%cmat, q_nm(coord, :), q_bead)
      positions_A(coord, 1, :) = q_bead
    end do
  end subroutine update_electron_positions

  subroutine half_kick_from_summary(summary, electron_index, electron_mass_amu, dt_fs, rp, force_bead, force_nm, v_nm)
    type(oracle_bead_summary), intent(in) :: summary
    integer, intent(in) :: electron_index
    real(dp), intent(in) :: electron_mass_amu, dt_fs
    type(ring_polymer_state), intent(in) :: rp
    real(dp), intent(out) :: force_bead(:, :)
    real(dp), intent(out) :: force_nm(:, :)
    real(dp), intent(inout) :: v_nm(:, :)
    integer :: coord

    force_bead = summary%scaled_forces_eV_per_A(:, electron_index, :)
    do coord = 1, size(force_bead, 1)
      call bead_to_nm(rp%cmat, force_bead(coord, :), force_nm(coord, :))
      v_nm(coord, :) = v_nm(coord, :) + 0.5_dp * dt_fs * force_nm(coord, :) / &
          (electron_mass_amu * amu_angstrom2_per_fs2_to_eV)
    end do
  end subroutine half_kick_from_summary

end module ereaxff_mixed_baoab_integrator
