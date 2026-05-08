module ereaxff_system_state
  use ereaxff_kinds, only: dp
  implicit none
  private

  integer, parameter, public :: element_len = 8

  type, public :: cell_state
    logical :: defined = .false.
    logical :: periodic = .false.
    real(dp) :: crystx(6) = 0.0_dp
    real(dp) :: matrix_A(3, 3) = 0.0_dp
    real(dp) :: inverse_A(3, 3) = 0.0_dp
  end type cell_state

  type, public :: system_state
    integer :: ndim = 0
    integer :: nparticles = 0
    integer, allocatable :: id(:)
    character(len=element_len), allocatable :: element(:)
    real(dp), allocatable :: mass_amu(:)
    real(dp), allocatable :: position_A(:, :)
    real(dp), allocatable :: velocity_A_per_fs(:, :)
    type(cell_state) :: cell
  end type system_state

  public :: init_system_state
  public :: set_particle
  public :: set_cell_from_crystx
  public :: set_orthorhombic_cell
  public :: cell_is_defined
  public :: cell_lengths_A
  public :: minimum_image_displacement
  public :: wrap_position_into_cell

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
    state%cell%defined = .false.
    state%cell%periodic = .false.
    state%cell%crystx = 0.0_dp
    state%cell%matrix_A = 0.0_dp
    state%cell%inverse_A = 0.0_dp
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

  subroutine set_cell_from_crystx(state, crystx, periodic)
    type(system_state), intent(inout) :: state
    real(dp), intent(in) :: crystx(6)
    logical, intent(in), optional :: periodic
    real(dp) :: a, b, c, alpha, beta, gamma
    real(dp) :: cos_alpha, cos_beta, cos_gamma, sin_gamma, volume_term

    a = crystx(1)
    b = crystx(2)
    c = crystx(3)
    alpha = crystx(4) * acos(-1.0_dp) / 180.0_dp
    beta = crystx(5) * acos(-1.0_dp) / 180.0_dp
    gamma = crystx(6) * acos(-1.0_dp) / 180.0_dp

    if (a <= 0.0_dp .or. b <= 0.0_dp .or. c <= 0.0_dp) error stop "set_cell_from_crystx: lengths must be positive"
    if (abs(sin(gamma)) <= 1.0e-12_dp) error stop "set_cell_from_crystx: gamma angle is singular"

    cos_alpha = cos(alpha)
    cos_beta = cos(beta)
    cos_gamma = cos(gamma)
    sin_gamma = sin(gamma)
    volume_term = 1.0_dp - cos_alpha**2 - cos_beta**2 - cos_gamma**2 + 2.0_dp*cos_alpha*cos_beta*cos_gamma
    if (volume_term < 0.0_dp) then
      if (volume_term > -1.0e-12_dp) then
        volume_term = 0.0_dp
      else
        error stop "set_cell_from_crystx: invalid crystal geometry"
      end if
    end if

    state%cell%defined = .true.
    if (present(periodic)) then
      state%cell%periodic = periodic
    else
      state%cell%periodic = .true.
    end if
    state%cell%crystx = crystx
    state%cell%matrix_A = 0.0_dp
    state%cell%matrix_A(:, 1) = [a, 0.0_dp, 0.0_dp]
    state%cell%matrix_A(1, 2) = b * cos_gamma
    state%cell%matrix_A(2, 2) = b * sin_gamma
    state%cell%matrix_A(1, 3) = c * cos_beta
    state%cell%matrix_A(2, 3) = c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma
    state%cell%matrix_A(3, 3) = c * sqrt(volume_term) / sin_gamma
    call invert_3x3(state%cell%matrix_A, state%cell%inverse_A)
  end subroutine set_cell_from_crystx

  subroutine set_orthorhombic_cell(state, lengths_A, periodic)
    type(system_state), intent(inout) :: state
    real(dp), intent(in) :: lengths_A(3)
    logical, intent(in), optional :: periodic

    call set_cell_from_crystx(state, [lengths_A(1), lengths_A(2), lengths_A(3), 90.0_dp, 90.0_dp, 90.0_dp], periodic)
  end subroutine set_orthorhombic_cell

  logical function cell_is_defined(state) result(is_defined)
    type(system_state), intent(in) :: state

    is_defined = state%cell%defined
  end function cell_is_defined

  pure function cell_lengths_A(state) result(lengths_A)
    type(system_state), intent(in) :: state
    real(dp) :: lengths_A(3)

    if (.not. state%cell%defined) then
      lengths_A = 0.0_dp
      return
    end if
    lengths_A(1) = sqrt(sum(state%cell%matrix_A(:, 1)**2))
    lengths_A(2) = sqrt(sum(state%cell%matrix_A(:, 2)**2))
    lengths_A(3) = sqrt(sum(state%cell%matrix_A(:, 3)**2))
  end function cell_lengths_A

  pure subroutine minimum_image_displacement(state, a, b, displacement_A)
    type(system_state), intent(in) :: state
    real(dp), intent(in) :: a(:), b(:)
    real(dp), intent(out) :: displacement_A(:)
    real(dp) :: fractional(3)

    displacement_A = a - b
    if (.not. state%cell%defined .or. .not. state%cell%periodic) return
    fractional = matmul(state%cell%inverse_A, displacement_A)
    fractional = fractional - nint(fractional)
    displacement_A = matmul(state%cell%matrix_A, fractional)
  end subroutine minimum_image_displacement

  pure subroutine wrap_position_into_cell(state, position_A)
    type(system_state), intent(in) :: state
    real(dp), intent(inout) :: position_A(:)
    real(dp) :: fractional(3)

    if (.not. state%cell%defined .or. .not. state%cell%periodic) return
    if (size(position_A) /= 3) error stop "wrap_position_into_cell: position must be 3D"
    fractional = matmul(state%cell%inverse_A, position_A)
    fractional = fractional - floor(fractional)
    position_A = matmul(state%cell%matrix_A, fractional)
  end subroutine wrap_position_into_cell

  subroutine invert_3x3(matrix_A, inverse_A)
    real(dp), intent(in) :: matrix_A(3, 3)
    real(dp), intent(out) :: inverse_A(3, 3)
    real(dp) :: det

    det = matrix_A(1,1) * (matrix_A(2,2) * matrix_A(3,3) - matrix_A(2,3) * matrix_A(3,2)) - &
          matrix_A(1,2) * (matrix_A(2,1) * matrix_A(3,3) - matrix_A(2,3) * matrix_A(3,1)) + &
          matrix_A(1,3) * (matrix_A(2,1) * matrix_A(3,2) - matrix_A(2,2) * matrix_A(3,1))
    if (abs(det) <= 1.0e-18_dp) error stop "invert_3x3: singular cell matrix"

    inverse_A(1,1) =  (matrix_A(2,2) * matrix_A(3,3) - matrix_A(2,3) * matrix_A(3,2)) / det
    inverse_A(1,2) = -(matrix_A(1,2) * matrix_A(3,3) - matrix_A(1,3) * matrix_A(3,2)) / det
    inverse_A(1,3) =  (matrix_A(1,2) * matrix_A(2,3) - matrix_A(1,3) * matrix_A(2,2)) / det
    inverse_A(2,1) = -(matrix_A(2,1) * matrix_A(3,3) - matrix_A(2,3) * matrix_A(3,1)) / det
    inverse_A(2,2) =  (matrix_A(1,1) * matrix_A(3,3) - matrix_A(1,3) * matrix_A(3,1)) / det
    inverse_A(2,3) = -(matrix_A(1,1) * matrix_A(2,3) - matrix_A(1,3) * matrix_A(2,1)) / det
    inverse_A(3,1) =  (matrix_A(2,1) * matrix_A(3,2) - matrix_A(2,2) * matrix_A(3,1)) / det
    inverse_A(3,2) = -(matrix_A(1,1) * matrix_A(3,2) - matrix_A(1,2) * matrix_A(3,1)) / det
    inverse_A(3,3) =  (matrix_A(1,1) * matrix_A(2,2) - matrix_A(1,2) * matrix_A(2,1)) / det
  end subroutine invert_3x3

end module ereaxff_system_state
