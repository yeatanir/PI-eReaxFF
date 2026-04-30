program test_thermostat_ou_statistics
  use ereaxff_kinds, only: dp, i8
  use ereaxff_constants, only: boltzmann_eV_per_K
  use ereaxff_units, only: amu_angstrom2_per_fs2_to_eV
  use ereaxff_rng, only: rng_state, rng_seed
  use ereaxff_thermostat_langevin, only: langevin_thermostat_state, langevin_init, langevin_step
  use ereaxff_thermostat_none, only: thermostat_none_state, thermostat_none_init, thermostat_none_step_1d
  use ereaxff_thermostat_registry, only: thermostat_is_supported, thermostat_requires_stub_failure
  implicit none

  integer, parameter :: nsteps = 250000
  real(dp), parameter :: temperature_K = 300.0_dp
  real(dp), parameter :: mass_amu(1) = [12.0_dp]
  real(dp), parameter :: gamma_fs_inv = 0.2_dp
  real(dp), parameter :: dt_fs = 0.5_dp
  type(rng_state) :: rng
  type(langevin_thermostat_state) :: langevin
  type(thermostat_none_state) :: none_state
  real(dp) :: v(1), v_none(1), sum_v2, variance, expected_variance, rel_error
  integer :: step

  call rng_seed(rng, 424242_i8)
  call langevin_init(langevin, temperature_K, gamma_fs_inv, dt_fs)

  v = 0.0_dp
  sum_v2 = 0.0_dp
  do step = 1, nsteps
    call langevin_step(langevin, v, mass_amu, rng)
    if (step > nsteps / 5) sum_v2 = sum_v2 + v(1)**2
  end do

  variance = sum_v2 / real(nsteps - nsteps / 5, dp)
  expected_variance = boltzmann_eV_per_K * temperature_K / &
      (mass_amu(1) * amu_angstrom2_per_fs2_to_eV)
  rel_error = abs(variance - expected_variance) / expected_variance
  call assert_less("Langevin OU velocity variance relative error", rel_error, 0.025_dp)
  call assert_true("Langevin heat accounting changed", abs(langevin%heat_eV) > 1.0e-8_dp)

  call thermostat_none_init(none_state)
  v_none = [0.123_dp]
  call thermostat_none_step_1d(none_state, v_none)
  call assert_less("none thermostat does not modify velocity", abs(v_none(1) - 0.123_dp), 0.0_dp)
  call assert_less("none thermostat heat is zero", abs(none_state%heat_eV), 0.0_dp)

  call assert_true("registry supports langevin", thermostat_is_supported("langevin"))
  call assert_true("registry rejects invalid thermostat", .not. thermostat_is_supported("invalid"))
  call assert_true("registry marks pile_g as stub", thermostat_requires_stub_failure("pile_g"))

  write(*,'(a)') "test_04_thermostat_ou_statistics: PASS"
  write(*,'(a,es12.4)') "  sampled_variance_A2_fs2 = ", variance
  write(*,'(a,es12.4)') "  expected_variance       = ", expected_variance
  write(*,'(a,es12.4)') "  relative_error          = ", rel_error
  write(*,'(a,es12.4)') "  accumulated_heat_eV     = ", langevin%heat_eV

contains

  subroutine assert_less(label, value, threshold)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: value, threshold

    if (value > threshold) then
      write(*,'(a)') "FAIL: " // trim(label)
      write(*,'(a,es24.16)') "  value = ", value
      write(*,'(a,es24.16)') "  max   = ", threshold
      error stop "thermostat OU regression failed"
    end if
  end subroutine assert_less

  subroutine assert_true(label, condition)
    character(len=*), intent(in) :: label
    logical, intent(in) :: condition

    if (.not. condition) then
      write(*,'(a)') "FAIL: " // trim(label)
      error stop "thermostat registry/heat regression failed"
    end if
  end subroutine assert_true

end program test_thermostat_ou_statistics

