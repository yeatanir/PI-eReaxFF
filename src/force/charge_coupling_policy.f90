module ereaxff_charge_coupling_policy
  use ereaxff_kinds, only: dp
  implicit none
  private

  character(len=*), parameter, public :: charge_policy_bead_local_acks2 = "bead_local_acks2"
  character(len=*), parameter, public :: charge_policy_centroid_acks2 = "centroid_acks2"
  character(len=*), parameter, public :: charge_policy_frozen_charges = "frozen_charges"
  character(len=*), parameter, public :: charge_policy_no_charge_coupling = "no_charge_coupling"

  type, public :: charge_memory_state
    integer :: nparticles = 0
    integer :: beads = 0
    character(len=32) :: policy = charge_policy_no_charge_coupling
    real(dp), allocatable :: charges_by_bead(:, :)
  end type charge_memory_state

  public :: charge_policy_is_valid
  public :: charge_policy_uses_bead_local_memory
  public :: init_charge_memory
  public :: set_bead_charges
  public :: charge_sum_for_bead
  public :: max_charge_difference_across_beads

contains

  logical function charge_policy_is_valid(policy) result(is_valid)
    character(len=*), intent(in) :: policy

    select case (trim(policy))
    case (charge_policy_bead_local_acks2, charge_policy_centroid_acks2, &
        charge_policy_frozen_charges, charge_policy_no_charge_coupling)
      is_valid = .true.
    case default
      is_valid = .false.
    end select
  end function charge_policy_is_valid

  logical function charge_policy_uses_bead_local_memory(policy) result(uses_memory)
    character(len=*), intent(in) :: policy

    uses_memory = trim(policy) == charge_policy_bead_local_acks2
  end function charge_policy_uses_bead_local_memory

  subroutine init_charge_memory(memory, nparticles, beads, policy)
    type(charge_memory_state), intent(inout) :: memory
    integer, intent(in) :: nparticles, beads
    character(len=*), intent(in) :: policy

    if (nparticles < 1) error stop "init_charge_memory: nparticles must be positive"
    if (beads < 1) error stop "init_charge_memory: beads must be positive"
    if (.not. charge_policy_is_valid(policy)) error stop "init_charge_memory: invalid charge policy"

    memory%nparticles = nparticles
    memory%beads = beads
    memory%policy = trim(policy)
    if (allocated(memory%charges_by_bead)) deallocate(memory%charges_by_bead)
    allocate(memory%charges_by_bead(nparticles, beads))
    memory%charges_by_bead = 0.0_dp
  end subroutine init_charge_memory

  subroutine set_bead_charges(memory, bead, charges)
    type(charge_memory_state), intent(inout) :: memory
    integer, intent(in) :: bead
    real(dp), intent(in) :: charges(:)

    if (.not. allocated(memory%charges_by_bead)) error stop "set_bead_charges: memory not initialized"
    if (bead < 1 .or. bead > memory%beads) error stop "set_bead_charges: bead out of range"
    if (size(charges) /= memory%nparticles) error stop "set_bead_charges: charge size mismatch"

    memory%charges_by_bead(:, bead) = charges
  end subroutine set_bead_charges

  real(dp) function charge_sum_for_bead(memory, bead) result(charge_sum)
    type(charge_memory_state), intent(in) :: memory
    integer, intent(in) :: bead

    if (.not. allocated(memory%charges_by_bead)) error stop "charge_sum_for_bead: memory not initialized"
    if (bead < 1 .or. bead > memory%beads) error stop "charge_sum_for_bead: bead out of range"

    charge_sum = sum(memory%charges_by_bead(:, bead))
  end function charge_sum_for_bead

  real(dp) function max_charge_difference_across_beads(memory) result(max_diff)
    type(charge_memory_state), intent(in) :: memory
    integer :: bead_i, bead_j

    if (.not. allocated(memory%charges_by_bead)) error stop "max_charge_difference_across_beads: memory not initialized"

    max_diff = 0.0_dp
    do bead_i = 1, memory%beads
      do bead_j = bead_i + 1, memory%beads
        max_diff = max(max_diff, maxval(abs(memory%charges_by_bead(:, bead_i) - &
            memory%charges_by_bead(:, bead_j))))
      end do
    end do
  end function max_charge_difference_across_beads

end module ereaxff_charge_coupling_policy

