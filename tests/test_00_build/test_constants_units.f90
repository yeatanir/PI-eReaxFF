program test_constants_units
  use ereaxff_kinds, only: dp, i8
  use ereaxff_constants, only: boltzmann_eV_per_K, hbar_eV_fs
  use ereaxff_units, only: amu_angstrom2_per_fs2_to_eV, kcalmol_to_ev_value, ev_to_kcalmol_value
  use ereaxff_rng, only: rng_state, rng_seed, rng_uniform
  implicit none

  type(rng_state) :: rng1, rng2
  real(dp) :: u1, u2

  call assert_close("kB eV/K", boltzmann_eV_per_K, 8.617333262145179e-5_dp, 1.0e-18_dp)
  call assert_close("hbar eV fs", hbar_eV_fs, 0.6582119569509067_dp, 1.0e-16_dp)
  call assert_close("amu A2/fs2 to eV", amu_angstrom2_per_fs2_to_eV, 103.642696671608_dp, 1.0e-12_dp)
  call assert_close("kcal/mol roundtrip", ev_to_kcalmol_value(kcalmol_to_ev_value(10.0_dp)), 10.0_dp, 1.0e-12_dp)

  call rng_seed(rng1, 12345_i8)
  call rng_seed(rng2, 12345_i8)
  u1 = rng_uniform(rng1)
  u2 = rng_uniform(rng2)
  call assert_close("rng reproducibility", u1, u2, 0.0_dp)

  write(*,'(a)') "test_00_build: PASS"
  write(*,'(a,es12.4)') "  amu_A2_fs2_to_eV = ", amu_angstrom2_per_fs2_to_eV
  write(*,'(a,es12.4)') "  first_rng_uniform = ", u1

contains

  subroutine assert_close(label, value, expected, threshold)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: value, expected, threshold

    if (abs(value - expected) > threshold) then
      write(*,'(a)') "FAIL: " // trim(label)
      write(*,'(a,es24.16)') "  value    = ", value
      write(*,'(a,es24.16)') "  expected = ", expected
      write(*,'(a,es24.16)') "  tol      = ", threshold
      error stop "constants/units regression failed"
    end if
  end subroutine assert_close

end program test_constants_units
