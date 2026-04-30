program test_normal_modes
  use ereaxff_kinds, only: dp
  use ereaxff_constants, only: pi, boltzmann_eV_per_K, hbar_eV_fs
  use ereaxff_normal_modes, only: build_cmat, bead_to_nm, nm_to_bead, &
      compute_omega_p, compute_omega_k, centroid_from_nm
  implicit none

  integer, parameter :: p = 16
  real(dp), parameter :: temperature_K = 300.0_dp
  real(dp), parameter :: tol = 1.0e-12_dp
  real(dp) :: cmat(p, p), ident(p, p), x_bead(p), x_nm(p), x_roundtrip(p)
  real(dp) :: omega_k(p)
  real(dp) :: max_orth, max_roundtrip, centroid, bead_average
  real(dp) :: omega_p, omega_p_expected, omega_p_sqrt_legacy
  real(dp) :: omega_expected, omega_sqrt_legacy
  integer :: i

  call build_cmat(p, cmat)
  ident = matmul(cmat, transpose(cmat))
  do i = 1, p
    ident(i, i) = ident(i, i) - 1.0_dp
  end do
  max_orth = maxval(abs(ident))
  call assert_close("C*C^T orthonormality", max_orth, 0.0_dp, tol)

  do i = 1, p
    x_bead(i) = sin(0.37_dp * real(i, dp)) + 0.1_dp * real(i * i, dp)
  end do
  call bead_to_nm(cmat, x_bead, x_nm)
  call nm_to_bead(cmat, x_nm, x_roundtrip)
  max_roundtrip = maxval(abs(x_roundtrip - x_bead))
  call assert_close("bead->nm->bead roundtrip", max_roundtrip, 0.0_dp, tol)

  centroid = centroid_from_nm(x_nm)
  bead_average = sum(x_bead) / real(p, dp)
  call assert_close("centroid equals bead average", centroid, bead_average, tol)

  omega_p = compute_omega_p(p, temperature_K)
  omega_p_expected = real(p, dp) * boltzmann_eV_per_K * temperature_K / hbar_eV_fs
  omega_p_sqrt_legacy = sqrt(real(p, dp)) * boltzmann_eV_per_K * temperature_K / hbar_eV_fs
  call assert_close("omega_P standard prefactor", omega_p, omega_p_expected, 1.0e-14_dp)
  if (abs(omega_p - omega_p_sqrt_legacy) < 1.0e-14_dp) then
    error stop "omega_P regression: legacy sqrt(P) prefactor detected"
  end if

  call compute_omega_k(p, temperature_K, omega_k)
  call assert_close("centroid omega_k", omega_k(1), 0.0_dp, 0.0_dp)
  do i = 2, p
    omega_expected = 2.0_dp * omega_p_expected * sin(pi * real(i - 1, dp) / real(p, dp))
    omega_sqrt_legacy = 2.0_dp * omega_p_sqrt_legacy * sin(pi * real(i - 1, dp) / real(p, dp))
    call assert_close("omega_k standard formula", omega_k(i), omega_expected, 1.0e-14_dp)
    if (abs(omega_k(i) - omega_sqrt_legacy) < 1.0e-14_dp) then
      error stop "omega_k regression: legacy sqrt(P) prefactor detected"
    end if
  end do

  write(*,'(a)') "test_01_normal_modes: PASS"
  write(*,'(a,es12.4)') "  max_orthogonality_error = ", max_orth
  write(*,'(a,es12.4)') "  max_roundtrip_error     = ", max_roundtrip
  write(*,'(a,es12.4)') "  omega_P_fs_inv          = ", omega_p

contains

  subroutine assert_close(label, value, expected, threshold)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: value, expected, threshold

    if (abs(value - expected) > threshold) then
      write(*,'(a)') "FAIL: " // trim(label)
      write(*,'(a,es24.16)') "  value    = ", value
      write(*,'(a,es24.16)') "  expected = ", expected
      write(*,'(a,es24.16)') "  tol      = ", threshold
      error stop "normal-mode regression failed"
    end if
  end subroutine assert_close

end program test_normal_modes

