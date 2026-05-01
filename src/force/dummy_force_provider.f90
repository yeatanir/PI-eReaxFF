module ereaxff_dummy_force_provider
  use ereaxff_kinds, only: dp
  use ereaxff_force_provider_base, only: force_result, allocate_force_result
  implicit none
  private

  type, public :: dummy_force_provider_state
    logical :: initialized = .false.
  end type dummy_force_provider_state

  public :: dummy_force_provider_init
  public :: dummy_force_provider_evaluate

contains

  subroutine dummy_force_provider_init(provider)
    type(dummy_force_provider_state), intent(inout) :: provider

    provider%initialized = .true.
  end subroutine dummy_force_provider_init

  subroutine dummy_force_provider_evaluate(provider, positions_A, result)
    type(dummy_force_provider_state), intent(in) :: provider
    real(dp), intent(in) :: positions_A(:, :)
    type(force_result), intent(inout) :: result

    if (.not. provider%initialized) error stop "dummy_force_provider_evaluate: provider not initialized"

    call allocate_force_result(result, size(positions_A, 1), size(positions_A, 2), .false.)
    result%energy_eV = 0.0_dp
    result%forces_eV_per_A = 0.0_dp
    result%status = 0
    result%diagnostic = "dummy force provider"
  end subroutine dummy_force_provider_evaluate

end module ereaxff_dummy_force_provider

