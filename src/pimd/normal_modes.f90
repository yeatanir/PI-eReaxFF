module ereaxff_normal_modes
  use ereaxff_kinds, only: dp
  use ereaxff_constants, only: pi, boltzmann_eV_per_K, hbar_eV_fs
  implicit none
  private

  character(len=*), parameter, public :: ring_polymer_convention_standard_ipi = "standard_ipi"

  public :: build_cmat
  public :: bead_to_nm
  public :: nm_to_bead
  public :: compute_omega_p
  public :: compute_omega_k
  public :: centroid_from_nm

contains

  subroutine build_cmat(p, cmat)
    integer, intent(in) :: p
    real(dp), intent(out) :: cmat(p, p)
    integer :: bead, mode, n, max_pair

    if (p < 1) error stop "build_cmat: bead count must be positive"

    cmat = 0.0_dp
    cmat(1, :) = 1.0_dp / sqrt(real(p, dp))

    max_pair = (p - 1) / 2
    do mode = 1, max_pair
      do bead = 1, p
        n = bead - 1
        cmat(2 * mode, bead) = sqrt(2.0_dp / real(p, dp)) * &
            cos(2.0_dp * pi * real(mode * n, dp) / real(p, dp))
        cmat(2 * mode + 1, bead) = sqrt(2.0_dp / real(p, dp)) * &
            sin(2.0_dp * pi * real(mode * n, dp) / real(p, dp))
      end do
    end do

    if (mod(p, 2) == 0) then
      do bead = 1, p
        if (mod(bead - 1, 2) == 0) then
          cmat(p, bead) = 1.0_dp / sqrt(real(p, dp))
        else
          cmat(p, bead) = -1.0_dp / sqrt(real(p, dp))
        end if
      end do
    end if
  end subroutine build_cmat

  subroutine bead_to_nm(cmat, x_bead, x_nm)
    real(dp), intent(in) :: cmat(:, :)
    real(dp), intent(in) :: x_bead(:)
    real(dp), intent(out) :: x_nm(:)

    call assert_transform_shape(cmat, size(x_bead), size(x_nm), "bead_to_nm")
    x_nm = matmul(cmat, x_bead)
  end subroutine bead_to_nm

  subroutine nm_to_bead(cmat, x_nm, x_bead)
    real(dp), intent(in) :: cmat(:, :)
    real(dp), intent(in) :: x_nm(:)
    real(dp), intent(out) :: x_bead(:)

    call assert_transform_shape(cmat, size(x_bead), size(x_nm), "nm_to_bead")
    x_bead = matmul(transpose(cmat), x_nm)
  end subroutine nm_to_bead

  pure real(dp) function compute_omega_p(p, temperature_K) result(omega_p)
    integer, intent(in) :: p
    real(dp), intent(in) :: temperature_K

    if (p < 1) error stop "compute_omega_p: bead count must be positive"
    if (temperature_K < 0.0_dp) error stop "compute_omega_p: temperature must be non-negative"

    omega_p = real(p, dp) * boltzmann_eV_per_K * temperature_K / hbar_eV_fs
  end function compute_omega_p

  subroutine compute_omega_k(p, temperature_K, omega_k)
    integer, intent(in) :: p
    real(dp), intent(in) :: temperature_K
    real(dp), intent(out) :: omega_k(p)
    integer :: k
    real(dp) :: omega_p

    omega_p = compute_omega_p(p, temperature_K)
    omega_k(1) = 0.0_dp
    do k = 2, p
      omega_k(k) = 2.0_dp * omega_p * sin(pi * real(k - 1, dp) / real(p, dp))
    end do
  end subroutine compute_omega_k

  pure real(dp) function centroid_from_nm(x_nm) result(centroid)
    real(dp), intent(in) :: x_nm(:)

    if (size(x_nm) < 1) error stop "centroid_from_nm: empty normal-mode vector"
    centroid = x_nm(1) / sqrt(real(size(x_nm), dp))
  end function centroid_from_nm

  subroutine assert_transform_shape(cmat, nbead, nnm, caller)
    real(dp), intent(in) :: cmat(:, :)
    integer, intent(in) :: nbead, nnm
    character(len=*), intent(in) :: caller

    if (size(cmat, 1) /= size(cmat, 2)) error stop caller // ": transform matrix must be square"
    if (size(cmat, 1) /= nbead) error stop caller // ": bead vector size mismatch"
    if (size(cmat, 2) /= nnm) error stop caller // ": normal-mode vector size mismatch"
  end subroutine assert_transform_shape

end module ereaxff_normal_modes

