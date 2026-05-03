module ereaxff_particle_groups
  use ereaxff_kinds, only: dp
  use ereaxff_system_state, only: system_state
  use ereaxff_selectors, only: selector_matches, selector_is_supported
  use ereaxff_thermostat_registry, only: thermostat_is_supported
  implicit none
  private

  character(len=*), parameter, public :: dynamics_classical = "classical"
  character(len=*), parameter, public :: dynamics_pimd = "pimd"
  character(len=*), parameter, public :: dynamics_frozen = "frozen"

  type, public :: particle_group_spec
    character(len=64) :: name = ""
    character(len=160) :: selector = ""
    character(len=16) :: dynamics = dynamics_classical
    logical :: mobile = .true.
    integer :: beads = 1
    character(len=32) :: thermostat = "none"
    real(dp) :: temperature_K = 0.0_dp
    real(dp) :: gamma_fs_inv = 0.0_dp
    real(dp) :: centroid_gamma_fs_inv = 0.0_dp
  end type particle_group_spec

  type, public :: particle_group_assignment
    integer, allocatable :: group_index(:)
    character(len=16), allocatable :: dynamics(:)
    logical, allocatable :: mobile(:)
    integer, allocatable :: beads(:)
    character(len=32), allocatable :: thermostat(:)
    real(dp), allocatable :: temperature_K(:)
  end type particle_group_assignment

  public :: validate_group_spec
  public :: assign_particle_groups
  public :: count_particles_with_dynamics

contains

  subroutine validate_group_spec(group)
    type(particle_group_spec), intent(in) :: group

    if (len_trim(group%name) == 0) error stop "validate_group_spec: group name is empty"
    if (.not. selector_is_supported(group%selector)) error stop "validate_group_spec: unsupported selector"
    if (.not. dynamics_is_valid(group%dynamics)) error stop "validate_group_spec: invalid dynamics"
    if (group%beads < 1) error stop "validate_group_spec: beads must be positive"
    if (trim(group%dynamics) /= dynamics_pimd .and. group%beads /= 1) then
      error stop "validate_group_spec: non-PIMD groups must use one bead"
    end if
    if (group%temperature_K < 0.0_dp) error stop "validate_group_spec: temperature must be non-negative"
    if (group%gamma_fs_inv < 0.0_dp) error stop "validate_group_spec: gamma must be non-negative"
    if (group%centroid_gamma_fs_inv < 0.0_dp) error stop "validate_group_spec: centroid gamma must be non-negative"
    if (.not. thermostat_is_supported(group%thermostat)) error stop "validate_group_spec: unsupported thermostat"
  end subroutine validate_group_spec

  subroutine assign_particle_groups(state, groups, assignment)
    type(system_state), intent(in) :: state
    type(particle_group_spec), intent(in) :: groups(:)
    type(particle_group_assignment), intent(inout) :: assignment
    integer :: i, g, matched

    if (size(groups) < 1) error stop "assign_particle_groups: no groups supplied"
    do g = 1, size(groups)
      call validate_group_spec(groups(g))
    end do

    call allocate_assignment(assignment, state%nparticles)

    do i = 1, state%nparticles
      matched = 0
      do g = 1, size(groups)
        if (selector_matches(state, i, groups(g)%selector)) then
          if (matched /= 0) error stop "assign_particle_groups: particle matches multiple groups"
          matched = g
        end if
      end do

      if (matched == 0) error stop "assign_particle_groups: particle did not match any group"

      assignment%group_index(i) = matched
      assignment%dynamics(i) = groups(matched)%dynamics
      assignment%mobile(i) = groups(matched)%mobile .and. trim(groups(matched)%dynamics) /= dynamics_frozen
      assignment%beads(i) = groups(matched)%beads
      assignment%thermostat(i) = groups(matched)%thermostat
      assignment%temperature_K(i) = groups(matched)%temperature_K
    end do
  end subroutine assign_particle_groups

  integer function count_particles_with_dynamics(assignment, dynamics) result(count)
    type(particle_group_assignment), intent(in) :: assignment
    character(len=*), intent(in) :: dynamics
    integer :: i

    count = 0
    do i = 1, size(assignment%dynamics)
      if (trim(assignment%dynamics(i)) == trim(dynamics)) count = count + 1
    end do
  end function count_particles_with_dynamics

  subroutine allocate_assignment(assignment, nparticles)
    type(particle_group_assignment), intent(inout) :: assignment
    integer, intent(in) :: nparticles

    if (allocated(assignment%group_index)) deallocate(assignment%group_index)
    if (allocated(assignment%dynamics)) deallocate(assignment%dynamics)
    if (allocated(assignment%mobile)) deallocate(assignment%mobile)
    if (allocated(assignment%beads)) deallocate(assignment%beads)
    if (allocated(assignment%thermostat)) deallocate(assignment%thermostat)
    if (allocated(assignment%temperature_K)) deallocate(assignment%temperature_K)

    allocate(assignment%group_index(nparticles))
    allocate(assignment%dynamics(nparticles))
    allocate(assignment%mobile(nparticles))
    allocate(assignment%beads(nparticles))
    allocate(assignment%thermostat(nparticles))
    allocate(assignment%temperature_K(nparticles))

    assignment%group_index = 0
    assignment%dynamics = ""
    assignment%mobile = .false.
    assignment%beads = 1
    assignment%thermostat = "none"
    assignment%temperature_K = 0.0_dp
  end subroutine allocate_assignment

  logical function dynamics_is_valid(dynamics) result(is_valid)
    character(len=*), intent(in) :: dynamics

    select case (trim(dynamics))
    case (dynamics_classical, dynamics_pimd, dynamics_frozen)
      is_valid = .true.
    case default
      is_valid = .false.
    end select
  end function dynamics_is_valid

end module ereaxff_particle_groups
