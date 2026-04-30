module ereaxff_thermostat_langevin
  use ereaxff_kinds, only: dp
  use ereaxff_constants, only: boltzmann_eV_per_K
  use ereaxff_rng, only: rng_state, rng_normal
  use ereaxff_thermostat_base, only: kinetic_energy_1d_eV, maxwell_velocity_sigma
  implicit none
  private

  type, public :: langevin_thermostat_state
    real(dp) :: temperature_K = 0.0_dp
    real(dp) :: gamma_fs_inv = 0.0_dp
    real(dp) :: dt_fs = 0.0_dp
    real(dp) :: c1 = 1.0_dp
    real(dp) :: heat_eV = 0.0_dp
  end type langevin_thermostat_state

  public :: langevin_init
  public :: langevin_step

contains

  subroutine langevin_init(state, temperature_K, gamma_fs_inv, dt_fs)
    type(langevin_thermostat_state), intent(inout) :: state
    real(dp), intent(in) :: temperature_K, gamma_fs_inv, dt_fs

    if (temperature_K < 0.0_dp) error stop "langevin_init: temperature must be non-negative"
    if (gamma_fs_inv < 0.0_dp) error stop "langevin_init: gamma must be non-negative"
    if (dt_fs < 0.0_dp) error stop "langevin_init: timestep must be non-negative"

    state%temperature_K = temperature_K
    state%gamma_fs_inv = gamma_fs_inv
    state%dt_fs = dt_fs
    state%c1 = exp(-gamma_fs_inv * dt_fs)
    state%heat_eV = 0.0_dp
  end subroutine langevin_init

  subroutine langevin_step(state, v, masses_amu, rng)
    type(langevin_thermostat_state), intent(inout) :: state
    real(dp), intent(inout) :: v(:)
    real(dp), intent(in) :: masses_amu(:)
    type(rng_state), intent(inout) :: rng
    integer :: i
    real(dp) :: ke_before, ke_after, kBT_eV, c2_unit

    if (size(v) /= size(masses_amu)) error stop "langevin_step: mass vector size mismatch"

    ke_before = kinetic_energy_1d_eV(v, masses_amu)
    kBT_eV = boltzmann_eV_per_K * state%temperature_K

    do i = 1, size(v)
      c2_unit = maxwell_velocity_sigma(kBT_eV, masses_amu(i)) * &
          sqrt(max(0.0_dp, 1.0_dp - state%c1**2))
      v(i) = state%c1 * v(i) + c2_unit * rng_normal(rng)
    end do

    ke_after = kinetic_energy_1d_eV(v, masses_amu)
    state%heat_eV = state%heat_eV + (ke_before - ke_after)
  end subroutine langevin_step

end module ereaxff_thermostat_langevin

