module ereaxff_baoab_integrator
  use ereaxff_kinds, only: dp
  use ereaxff_units, only: amu_angstrom2_per_fs2_to_eV
  use ereaxff_rng, only: rng_state
  use ereaxff_normal_modes, only: bead_to_nm, nm_to_bead
  use ereaxff_free_ring_propagator, only: exact_free_propagation
  use ereaxff_thermostat_pile_l, only: pile_l_thermostat_state, pile_l_step
  implicit none
  private

  public :: baoab_harmonic_pimd_step_1d

contains

  subroutine baoab_harmonic_pimd_step_1d(q_nm, v_nm, cmat, omega_k, mass_amu, k_eV_per_A2, &
      dt_fs, thermostat, rng, physical_force_scale)
    real(dp), intent(inout) :: q_nm(:)
    real(dp), intent(inout) :: v_nm(:)
    real(dp), intent(in) :: cmat(:, :)
    real(dp), intent(in) :: omega_k(:)
    real(dp), intent(in) :: mass_amu, k_eV_per_A2, dt_fs
    type(pile_l_thermostat_state), intent(inout) :: thermostat
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in), optional :: physical_force_scale
    real(dp) :: q2d(1, size(q_nm)), v2d(1, size(v_nm)), masses(1)
    real(dp) :: force_nm(size(q_nm))
    real(dp) :: scale

    if (size(q_nm) /= size(v_nm)) error stop "baoab_harmonic_pimd_step_1d: q/v size mismatch"
    if (size(omega_k) /= size(q_nm)) error stop "baoab_harmonic_pimd_step_1d: omega size mismatch"
    if (mass_amu <= 0.0_dp) error stop "baoab_harmonic_pimd_step_1d: mass must be positive"
    if (k_eV_per_A2 <= 0.0_dp) error stop "baoab_harmonic_pimd_step_1d: k must be positive"

    scale = 1.0_dp
    if (present(physical_force_scale)) scale = physical_force_scale

    call harmonic_force_nm(q_nm, cmat, k_eV_per_A2, scale, force_nm)
    v_nm = v_nm + 0.5_dp * dt_fs * force_nm / (mass_amu * amu_angstrom2_per_fs2_to_eV)

    q2d(1, :) = q_nm
    v2d(1, :) = v_nm
    call exact_free_propagation(q2d, v2d, 0.5_dp * dt_fs, omega_k)

    masses = [mass_amu]
    call pile_l_step(thermostat, v2d, masses, rng)

    call exact_free_propagation(q2d, v2d, 0.5_dp * dt_fs, omega_k)
    q_nm = q2d(1, :)
    v_nm = v2d(1, :)

    call harmonic_force_nm(q_nm, cmat, k_eV_per_A2, scale, force_nm)
    v_nm = v_nm + 0.5_dp * dt_fs * force_nm / (mass_amu * amu_angstrom2_per_fs2_to_eV)
  end subroutine baoab_harmonic_pimd_step_1d

  subroutine harmonic_force_nm(q_nm, cmat, k_eV_per_A2, physical_force_scale, force_nm)
    real(dp), intent(in) :: q_nm(:)
    real(dp), intent(in) :: cmat(:, :)
    real(dp), intent(in) :: k_eV_per_A2
    real(dp), intent(in) :: physical_force_scale
    real(dp), intent(out) :: force_nm(:)
    real(dp) :: q_bead(size(q_nm)), force_bead(size(q_nm))

    call nm_to_bead(cmat, q_nm, q_bead)
    force_bead = -k_eV_per_A2 * q_bead * physical_force_scale
    call bead_to_nm(cmat, force_bead, force_nm)
  end subroutine harmonic_force_nm

end module ereaxff_baoab_integrator
