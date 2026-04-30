module ereaxff_thermostat_pile_l
  use ereaxff_kinds, only: dp
  use ereaxff_constants, only: boltzmann_eV_per_K
  use ereaxff_rng, only: rng_state, rng_normal
  use ereaxff_thermostat_base, only: kinetic_energy_2d_eV, maxwell_velocity_sigma
  implicit none
  private

  type, public :: pile_l_thermostat_state
    real(dp) :: temperature_K = 0.0_dp
    real(dp) :: centroid_gamma_fs_inv = 0.0_dp
    real(dp) :: lambda = 1.0_dp
    real(dp) :: dt_fs = 0.0_dp
    real(dp) :: heat_eV = 0.0_dp
    real(dp), allocatable :: gamma_mode_fs_inv(:)
    real(dp), allocatable :: c1(:)
  end type pile_l_thermostat_state

  public :: pile_l_init
  public :: pile_l_step

contains

  subroutine pile_l_init(state, temperature_K, centroid_gamma_fs_inv, lambda, dt_fs, omega_k)
    type(pile_l_thermostat_state), intent(inout) :: state
    real(dp), intent(in) :: temperature_K, centroid_gamma_fs_inv, lambda, dt_fs
    real(dp), intent(in) :: omega_k(:)
    integer :: mode

    if (temperature_K < 0.0_dp) error stop "pile_l_init: temperature must be non-negative"
    if (centroid_gamma_fs_inv < 0.0_dp) error stop "pile_l_init: centroid gamma must be non-negative"
    if (lambda < 0.0_dp) error stop "pile_l_init: lambda must be non-negative"
    if (dt_fs < 0.0_dp) error stop "pile_l_init: timestep must be non-negative"
    if (size(omega_k) < 1) error stop "pile_l_init: omega_k must not be empty"

    state%temperature_K = temperature_K
    state%centroid_gamma_fs_inv = centroid_gamma_fs_inv
    state%lambda = lambda
    state%dt_fs = dt_fs
    state%heat_eV = 0.0_dp

    allocate(state%gamma_mode_fs_inv(size(omega_k)))
    allocate(state%c1(size(omega_k)))
    state%gamma_mode_fs_inv(1) = centroid_gamma_fs_inv
    do mode = 2, size(omega_k)
      state%gamma_mode_fs_inv(mode) = 2.0_dp * lambda * omega_k(mode)
    end do
    state%c1 = exp(-state%gamma_mode_fs_inv * dt_fs)
  end subroutine pile_l_init

  subroutine pile_l_step(state, v_nm, masses_amu, rng)
    type(pile_l_thermostat_state), intent(inout) :: state
    real(dp), intent(inout) :: v_nm(:, :)
    real(dp), intent(in) :: masses_amu(:)
    type(rng_state), intent(inout) :: rng
    integer :: coord, mode
    real(dp) :: ke_before, ke_after, kBT_eV, sigma

    if (size(v_nm, 1) /= size(masses_amu)) error stop "pile_l_step: mass vector size mismatch"
    if (size(v_nm, 2) /= size(state%c1)) error stop "pile_l_step: mode count mismatch"

    ke_before = kinetic_energy_2d_eV(v_nm, masses_amu)
    kBT_eV = boltzmann_eV_per_K * state%temperature_K

    do mode = 1, size(v_nm, 2)
      do coord = 1, size(v_nm, 1)
        sigma = maxwell_velocity_sigma(kBT_eV, masses_amu(coord)) * &
            sqrt(max(0.0_dp, 1.0_dp - state%c1(mode)**2))
        v_nm(coord, mode) = state%c1(mode) * v_nm(coord, mode) + sigma * rng_normal(rng)
      end do
    end do

    ke_after = kinetic_energy_2d_eV(v_nm, masses_amu)
    state%heat_eV = state%heat_eV + (ke_before - ke_after)
  end subroutine pile_l_step

end module ereaxff_thermostat_pile_l

