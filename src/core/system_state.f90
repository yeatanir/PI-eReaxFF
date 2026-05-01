module ereaxff_system_state
  use ereaxff_kinds, only: dp
  implicit none
  private

  integer, parameter, public :: element_len = 8

  type, public :: system_state
    integer :: ndim = 0
    integer :: nparticles = 0
    integer, allocatable :: id(:)
    character(len=element_len), allocatable :: element(:)
    real(dp), allocatable :: mass_amu(:)
    real(dp), allocatable :: position_A(:, :)
    real(dp), allocatable :: velocity_A_per_fs(:, :)
  end type system_state

  public :: init_system_state
  public :: set_particle

contains

  subroutine init_system_state(state, ndim, nparticles)
    type(system_state), intent(inout) :: state
    integer, intent(in) :: ndim, nparticles
    integer :: i

    if (ndim < 1) error stop "init_system_state: ndim must be positive"
    if (nparticles < 1) error stop "init_system_state: nparticles must be positive"

    state%ndim = ndim
    state%nparticles = nparticles

    if (allocated(state%id)) deallocate(state%id)
    if (allocated(state%element)) deallocate(state%element)
    if (allocated(state%mass_amu)) deallocate(state%mass_amu)
    if (allocated(state%position_A)) deallocate(state%position_A)
    if (allocated(state%velocity_A_per_fs)) deallocate(state%velocity_A_per_fs)

    allocate(state%id(nparticles))
    allocate(state%element(nparticles))
    allocate(state%mass_amu(nparticles))
    allocate(state%position_A(ndim, nparticles))
    allocate(state%velocity_A_per_fs(ndim, nparticles))

    state%id = [(i, i = 1, nparticles)]
    state%element = ""
    state%mass_amu = 0.0_dp
    state%position_A = 0.0_dp
    state%velocity_A_per_fs = 0.0_dp
  end subroutine init_system_state

  subroutine set_particle(state, index, id, element, mass_amu, position_A, velocity_A_per_fs)
    type(system_state), intent(inout) :: state
    integer, intent(in) :: index, id
    character(len=*), intent(in) :: element
    real(dp), intent(in) :: mass_amu
    real(dp), intent(in) :: position_A(:)
    real(dp), intent(in), optional :: velocity_A_per_fs(:)

    if (index < 1 .or. index > state%nparticles) error stop "set_particle: index out of range"
    if (size(position_A) /= state%ndim) error stop "set_particle: position dimension mismatch"
    if (present(velocity_A_per_fs)) then
      if (size(velocity_A_per_fs) /= state%ndim) error stop "set_particle: velocity dimension mismatch"
    end if
    if (mass_amu <= 0.0_dp) error stop "set_particle: mass must be positive"

    state%id(index) = id
    state%element(index) = adjustl(element)
    state%mass_amu(index) = mass_amu
    state%position_A(:, index) = position_A
    if (present(velocity_A_per_fs)) state%velocity_A_per_fs(:, index) = velocity_A_per_fs
  end subroutine set_particle

end module ereaxff_system_state
