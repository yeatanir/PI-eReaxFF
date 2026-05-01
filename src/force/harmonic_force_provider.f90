module ereaxff_harmonic_force_provider
  use ereaxff_kinds, only: dp
  use ereaxff_force_provider_base, only: force_result, allocate_force_result
  implicit none
  private

  type, public :: harmonic_force_provider_state
    real(dp) :: k_eV_per_A2 = 0.0_dp
    real(dp), allocatable :: center_A(:, :)
    logical :: initialized = .false.
  end type harmonic_force_provider_state

  public :: harmonic_force_provider_init
  public :: harmonic_force_provider_evaluate

contains

  subroutine harmonic_force_provider_init(provider, k_eV_per_A2, center_A)
    type(harmonic_force_provider_state), intent(inout) :: provider
    real(dp), intent(in) :: k_eV_per_A2
    real(dp), intent(in), optional :: center_A(:, :)

    if (k_eV_per_A2 < 0.0_dp) error stop "harmonic_force_provider_init: k must be non-negative"

    provider%k_eV_per_A2 = k_eV_per_A2
    if (allocated(provider%center_A)) deallocate(provider%center_A)
    if (present(center_A)) then
      allocate(provider%center_A(size(center_A, 1), size(center_A, 2)))
      provider%center_A = center_A
    end if
    provider%initialized = .true.
  end subroutine harmonic_force_provider_init

  subroutine harmonic_force_provider_evaluate(provider, positions_A, result)
    type(harmonic_force_provider_state), intent(in) :: provider
    real(dp), intent(in) :: positions_A(:, :)
    type(force_result), intent(inout) :: result
    real(dp), allocatable :: displacement(:, :)

    if (.not. provider%initialized) error stop "harmonic_force_provider_evaluate: provider not initialized"

    if (allocated(provider%center_A)) then
      if (any(shape(provider%center_A) /= shape(positions_A))) then
        error stop "harmonic_force_provider_evaluate: center shape mismatch"
      end if
      displacement = positions_A - provider%center_A
    else
      displacement = positions_A
    end if

    call allocate_force_result(result, size(positions_A, 1), size(positions_A, 2), .false.)
    result%energy_eV = 0.5_dp * provider%k_eV_per_A2 * sum(displacement**2)
    result%forces_eV_per_A = -provider%k_eV_per_A2 * displacement
    result%status = 0
    result%diagnostic = "harmonic force provider"
  end subroutine harmonic_force_provider_evaluate

end module ereaxff_harmonic_force_provider

