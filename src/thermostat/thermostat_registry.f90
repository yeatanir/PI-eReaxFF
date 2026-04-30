module ereaxff_thermostat_registry
  implicit none
  private

  public :: thermostat_is_supported
  public :: thermostat_requires_stub_failure

contains

  logical function thermostat_is_supported(name) result(is_supported)
    character(len=*), intent(in) :: name

    select case (trim(name))
    case ("none", "langevin", "pile_l", "pile_g")
      is_supported = .true.
    case default
      is_supported = .false.
    end select
  end function thermostat_is_supported

  logical function thermostat_requires_stub_failure(name) result(requires_failure)
    character(len=*), intent(in) :: name

    requires_failure = trim(name) == "pile_g"
  end function thermostat_requires_stub_failure

end module ereaxff_thermostat_registry

