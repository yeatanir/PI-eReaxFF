program test_pile_mode_temperatures
  use ereaxff_kinds, only: dp, i8
  use ereaxff_constants, only: boltzmann_eV_per_K
  use ereaxff_units, only: amu_angstrom2_per_fs2_to_eV
  use ereaxff_rng, only: rng_state, rng_seed
  use ereaxff_ring_polymer, only: ring_polymer_state, init_ring_polymer
  use ereaxff_thermostat_pile_l, only: pile_l_thermostat_state, pile_l_init, pile_l_step
  implicit none

  integer, parameter :: ndim = 3
  integer, parameter :: p = 16
  integer, parameter :: nsteps = 180000
  real(dp), parameter :: temperature_K = 300.0_dp
  real(dp), parameter :: dt_fs = 0.05_dp
  real(dp), parameter :: centroid_gamma_fs_inv = 0.02_dp
  real(dp), parameter :: lambda = 1.0_dp
  real(dp), parameter :: masses_amu(ndim) = [1.0_dp, 4.0_dp, 12.0_dp]
  type(ring_polymer_state) :: rp
  type(pile_l_thermostat_state) :: pile_l
  type(rng_state) :: rng
  real(dp) :: v_nm(ndim, p), mode_ke_sum(p), mode_temp(p), expected_temp
  real(dp) :: max_internal_rel_error
  integer :: step, mode, coord, nsample

  call init_ring_polymer(rp, p, temperature_K, "standard_ipi")
  call pile_l_init(pile_l, temperature_K, centroid_gamma_fs_inv, lambda, dt_fs, rp%omega_k)
  call rng_seed(rng, 777_i8)

  v_nm = 0.0_dp
  mode_ke_sum = 0.0_dp
  nsample = 0

  do step = 1, nsteps
    call pile_l_step(pile_l, v_nm, masses_amu, rng)
    if (step > nsteps / 5) then
      nsample = nsample + 1
      do mode = 1, p
        do coord = 1, ndim
          mode_ke_sum(mode) = mode_ke_sum(mode) + 0.5_dp * masses_amu(coord) * &
              v_nm(coord, mode)**2 * amu_angstrom2_per_fs2_to_eV
        end do
      end do
    end if
  end do

  expected_temp = temperature_K
  do mode = 1, p
    mode_temp(mode) = 2.0_dp * (mode_ke_sum(mode) / real(nsample, dp)) / &
        (real(ndim, dp) * boltzmann_eV_per_K)
  end do

  max_internal_rel_error = maxval(abs(mode_temp(2:p) - expected_temp) / expected_temp)
  call assert_less("PILE-L internal mode temperature relative error", max_internal_rel_error, 0.08_dp)
  call assert_true("PILE-L internal damping uses omega_k", all(pile_l%gamma_mode_fs_inv(2:p) > 0.0_dp))
  call assert_true("PILE-L heat accounting changed", abs(pile_l%heat_eV) > 1.0e-8_dp)

  write(*,'(a)') "test_05_pile_mode_temperatures: PASS"
  write(*,'(a,es12.4)') "  centroid_temp_K              = ", mode_temp(1)
  write(*,'(a,es12.4)') "  max_internal_rel_error       = ", max_internal_rel_error
  write(*,'(a,es12.4)') "  accumulated_heat_eV          = ", pile_l%heat_eV

contains

  subroutine assert_less(label, value, threshold)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: value, threshold

    if (value > threshold) then
      write(*,'(a)') "FAIL: " // trim(label)
      write(*,'(a,es24.16)') "  value = ", value
      write(*,'(a,es24.16)') "  max   = ", threshold
      error stop "PILE-L temperature regression failed"
    end if
  end subroutine assert_less

  subroutine assert_true(label, condition)
    character(len=*), intent(in) :: label
    logical, intent(in) :: condition

    if (.not. condition) then
      write(*,'(a)') "FAIL: " // trim(label)
      error stop "PILE-L thermostat regression failed"
    end if
  end subroutine assert_true

end program test_pile_mode_temperatures

