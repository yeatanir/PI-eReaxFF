module ereaxff_units
  use ereaxff_kinds, only: dp
  use ereaxff_constants, only: kcalmol_to_eV
  implicit none
  private

  real(dp), parameter, public :: joule_per_eV = 1.602176634e-19_dp
  real(dp), parameter, public :: kg_per_amu = 1.66053906892e-27_dp
  real(dp), parameter, public :: meter_per_angstrom = 1.0e-10_dp
  real(dp), parameter, public :: second_per_fs = 1.0e-15_dp

  real(dp), parameter, public :: amu_angstrom2_per_fs2_to_eV = &
      kg_per_amu * (meter_per_angstrom / second_per_fs)**2 / joule_per_eV

  public :: kcalmol_to_ev_value
  public :: ev_to_kcalmol_value
  public :: kinetic_energy_eV

contains

  pure real(dp) function kcalmol_to_ev_value(value_kcalmol) result(value_eV)
    real(dp), intent(in) :: value_kcalmol

    value_eV = value_kcalmol * kcalmol_to_eV
  end function kcalmol_to_ev_value

  pure real(dp) function ev_to_kcalmol_value(value_eV) result(value_kcalmol)
    real(dp), intent(in) :: value_eV

    value_kcalmol = value_eV / kcalmol_to_eV
  end function ev_to_kcalmol_value

  pure real(dp) function kinetic_energy_eV(mass_amu, velocity_A_per_fs) result(ke_eV)
    real(dp), intent(in) :: mass_amu
    real(dp), intent(in) :: velocity_A_per_fs

    ke_eV = 0.5_dp * mass_amu * velocity_A_per_fs**2 * amu_angstrom2_per_fs2_to_eV
  end function kinetic_energy_eV

end module ereaxff_units

