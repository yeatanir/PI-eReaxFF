module ereaxff_force_provider_base
  use ereaxff_kinds, only: dp
  implicit none
  private

  type, public :: force_result
    real(dp) :: energy_eV = 0.0_dp
    real(dp), allocatable :: forces_eV_per_A(:, :)
    real(dp), allocatable :: charges(:)
    integer :: status = 0
    character(len=256) :: diagnostic = ""
  end type force_result

  type, public :: batched_force_result
    integer :: nbeads = 0
    real(dp), allocatable :: energies_eV(:)
    real(dp), allocatable :: forces_eV_per_A(:, :, :)
    real(dp), allocatable :: charges(:, :)
    integer, allocatable :: status(:)
    character(len=256), allocatable :: diagnostics(:)
  end type batched_force_result

  public :: allocate_force_result
  public :: allocate_batched_force_result
  public :: assert_position_force_shape

contains

  subroutine allocate_force_result(result, ndim, nparticles, with_charges)
    type(force_result), intent(inout) :: result
    integer, intent(in) :: ndim, nparticles
    logical, intent(in), optional :: with_charges
    logical :: allocate_charges

    if (ndim < 1) error stop "allocate_force_result: ndim must be positive"
    if (nparticles < 1) error stop "allocate_force_result: nparticles must be positive"

    allocate_charges = .false.
    if (present(with_charges)) allocate_charges = with_charges

    if (allocated(result%forces_eV_per_A)) deallocate(result%forces_eV_per_A)
    if (allocated(result%charges)) deallocate(result%charges)

    allocate(result%forces_eV_per_A(ndim, nparticles))
    result%forces_eV_per_A = 0.0_dp
    if (allocate_charges) then
      allocate(result%charges(nparticles))
      result%charges = 0.0_dp
    end if

    result%energy_eV = 0.0_dp
    result%status = 0
    result%diagnostic = ""
  end subroutine allocate_force_result

  subroutine allocate_batched_force_result(result, ndim, nparticles, nbeads, with_charges)
    type(batched_force_result), intent(inout) :: result
    integer, intent(in) :: ndim, nparticles, nbeads
    logical, intent(in), optional :: with_charges
    logical :: allocate_charges

    if (ndim < 1) error stop "allocate_batched_force_result: ndim must be positive"
    if (nparticles < 1) error stop "allocate_batched_force_result: nparticles must be positive"
    if (nbeads < 1) error stop "allocate_batched_force_result: nbeads must be positive"

    allocate_charges = .false.
    if (present(with_charges)) allocate_charges = with_charges

    if (allocated(result%energies_eV)) deallocate(result%energies_eV)
    if (allocated(result%forces_eV_per_A)) deallocate(result%forces_eV_per_A)
    if (allocated(result%charges)) deallocate(result%charges)
    if (allocated(result%status)) deallocate(result%status)
    if (allocated(result%diagnostics)) deallocate(result%diagnostics)

    result%nbeads = nbeads
    allocate(result%energies_eV(nbeads))
    allocate(result%forces_eV_per_A(ndim, nparticles, nbeads))
    allocate(result%status(nbeads))
    allocate(result%diagnostics(nbeads))
    result%energies_eV = 0.0_dp
    result%forces_eV_per_A = 0.0_dp
    result%status = 0
    result%diagnostics = ""
    if (allocate_charges) then
      allocate(result%charges(nparticles, nbeads))
      result%charges = 0.0_dp
    end if
  end subroutine allocate_batched_force_result

  subroutine assert_position_force_shape(positions_A, result, caller)
    real(dp), intent(in) :: positions_A(:, :)
    type(force_result), intent(in) :: result
    character(len=*), intent(in) :: caller

    if (.not. allocated(result%forces_eV_per_A)) error stop caller // ": result forces not allocated"
    if (any(shape(result%forces_eV_per_A) /= shape(positions_A))) then
      error stop caller // ": force result shape mismatch"
    end if
  end subroutine assert_position_force_shape

end module ereaxff_force_provider_base
