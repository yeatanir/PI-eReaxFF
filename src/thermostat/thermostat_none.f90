module ereaxff_thermostat_none
  use ereaxff_kinds, only: dp
  implicit none
  private

  type, public :: thermostat_none_state
    real(dp) :: heat_eV = 0.0_dp
  end type thermostat_none_state

  public :: thermostat_none_init
  public :: thermostat_none_step_1d
  public :: thermostat_none_step_2d

contains

  subroutine thermostat_none_init(state)
    type(thermostat_none_state), intent(inout) :: state

    state%heat_eV = 0.0_dp
  end subroutine thermostat_none_init

  subroutine thermostat_none_step_1d(state, v)
    type(thermostat_none_state), intent(inout) :: state
    real(dp), intent(inout) :: v(:)

    state%heat_eV = state%heat_eV + 0.0_dp * sum(v)
  end subroutine thermostat_none_step_1d

  subroutine thermostat_none_step_2d(state, v)
    type(thermostat_none_state), intent(inout) :: state
    real(dp), intent(inout) :: v(:, :)

    state%heat_eV = state%heat_eV + 0.0_dp * sum(v)
  end subroutine thermostat_none_step_2d

end module ereaxff_thermostat_none
