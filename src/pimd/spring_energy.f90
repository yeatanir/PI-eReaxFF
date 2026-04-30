module ereaxff_spring_energy
  use ereaxff_kinds, only: dp
  use ereaxff_units, only: amu_angstrom2_per_fs2_to_eV
  implicit none
  private

  public :: compute_spring_energy_nm
  public :: compute_free_ring_energy_nm

contains

  real(dp) function compute_spring_energy_nm(q_nm, masses_amu, omega_k) result(energy_eV)
    real(dp), intent(in) :: q_nm(:, :)
    real(dp), intent(in) :: masses_amu(:)
    real(dp), intent(in) :: omega_k(:)
    integer :: coord, mode

    call assert_energy_shapes(q_nm, masses_amu, omega_k, "compute_spring_energy_nm")

    energy_eV = 0.0_dp
    do mode = 2, size(q_nm, 2)
      do coord = 1, size(q_nm, 1)
        energy_eV = energy_eV + 0.5_dp * masses_amu(coord) * omega_k(mode)**2 * &
            q_nm(coord, mode)**2 * amu_angstrom2_per_fs2_to_eV
      end do
    end do
  end function compute_spring_energy_nm

  real(dp) function compute_free_ring_energy_nm(q_nm, v_nm, masses_amu, omega_k) result(energy_eV)
    real(dp), intent(in) :: q_nm(:, :)
    real(dp), intent(in) :: v_nm(:, :)
    real(dp), intent(in) :: masses_amu(:)
    real(dp), intent(in) :: omega_k(:)
    integer :: coord, mode

    call assert_energy_shapes(q_nm, masses_amu, omega_k, "compute_free_ring_energy_nm")
    if (any(shape(v_nm) /= shape(q_nm))) then
      error stop "compute_free_ring_energy_nm: velocity shape mismatch"
    end if

    energy_eV = 0.0_dp
    do mode = 1, size(q_nm, 2)
      do coord = 1, size(q_nm, 1)
        energy_eV = energy_eV + 0.5_dp * masses_amu(coord) * v_nm(coord, mode)**2 * &
            amu_angstrom2_per_fs2_to_eV
      end do
    end do
    energy_eV = energy_eV + compute_spring_energy_nm(q_nm, masses_amu, omega_k)
  end function compute_free_ring_energy_nm

  subroutine assert_energy_shapes(q_nm, masses_amu, omega_k, caller)
    real(dp), intent(in) :: q_nm(:, :)
    real(dp), intent(in) :: masses_amu(:)
    real(dp), intent(in) :: omega_k(:)
    character(len=*), intent(in) :: caller

    if (size(q_nm, 1) /= size(masses_amu)) error stop caller // ": mass vector size mismatch"
    if (size(q_nm, 2) /= size(omega_k)) error stop caller // ": omega_k size mismatch"
  end subroutine assert_energy_shapes

end module ereaxff_spring_energy

