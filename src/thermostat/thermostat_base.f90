module ereaxff_thermostat_base
  use ereaxff_kinds, only: dp
  use ereaxff_units, only: amu_angstrom2_per_fs2_to_eV
  implicit none
  private

  public :: kinetic_energy_1d_eV
  public :: kinetic_energy_2d_eV
  public :: maxwell_velocity_sigma

contains

  pure real(dp) function kinetic_energy_1d_eV(v, masses_amu) result(ke_eV)
    real(dp), intent(in) :: v(:)
    real(dp), intent(in) :: masses_amu(:)
    integer :: i

    if (size(v) /= size(masses_amu)) error stop "kinetic_energy_1d_eV: mass vector size mismatch"

    ke_eV = 0.0_dp
    do i = 1, size(v)
      ke_eV = ke_eV + 0.5_dp * masses_amu(i) * v(i)**2 * amu_angstrom2_per_fs2_to_eV
    end do
  end function kinetic_energy_1d_eV

  pure real(dp) function kinetic_energy_2d_eV(v, masses_amu) result(ke_eV)
    real(dp), intent(in) :: v(:, :)
    real(dp), intent(in) :: masses_amu(:)
    integer :: coord, mode

    if (size(v, 1) /= size(masses_amu)) error stop "kinetic_energy_2d_eV: mass vector size mismatch"

    ke_eV = 0.0_dp
    do mode = 1, size(v, 2)
      do coord = 1, size(v, 1)
        ke_eV = ke_eV + 0.5_dp * masses_amu(coord) * v(coord, mode)**2 * &
            amu_angstrom2_per_fs2_to_eV
      end do
    end do
  end function kinetic_energy_2d_eV

  pure real(dp) function maxwell_velocity_sigma(kBT_eV, mass_amu) result(sigma_A_per_fs)
    real(dp), intent(in) :: kBT_eV
    real(dp), intent(in) :: mass_amu

    if (mass_amu <= 0.0_dp) error stop "maxwell_velocity_sigma: mass must be positive"
    if (kBT_eV < 0.0_dp) error stop "maxwell_velocity_sigma: kBT must be non-negative"

    sigma_A_per_fs = sqrt(kBT_eV / (mass_amu * amu_angstrom2_per_fs2_to_eV))
  end function maxwell_velocity_sigma

end module ereaxff_thermostat_base

