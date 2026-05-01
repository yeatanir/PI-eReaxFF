module ereaxff_diagnostics
  use ereaxff_kinds, only: dp
  implicit none
  private

  public :: nvt_conserved_like

contains

  pure real(dp) function nvt_conserved_like(kinetic_eV, potential_eV, spring_eV, thermostat_heat_eV) result(value)
    real(dp), intent(in) :: kinetic_eV, potential_eV, spring_eV, thermostat_heat_eV

    value = kinetic_eV + potential_eV + spring_eV + thermostat_heat_eV
  end function nvt_conserved_like

end module ereaxff_diagnostics

