module ereaxff_thermostat_pile_g
  implicit none
  private

  public :: pile_g_not_implemented

contains

  subroutine pile_g_not_implemented()
    error stop "PILE-G requested but not implemented; use PILE-L or implement global centroid SVR."
  end subroutine pile_g_not_implemented

end module ereaxff_thermostat_pile_g

