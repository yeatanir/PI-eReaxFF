module ereaxff_ring_polymer
  use ereaxff_kinds, only: dp
  use ereaxff_normal_modes, only: ring_polymer_convention_standard_ipi, build_cmat, compute_omega_k
  implicit none
  private

  type, public :: ring_polymer_state
    integer :: beads = 0
    real(dp) :: temperature_K = 0.0_dp
    character(len=32) :: convention = ring_polymer_convention_standard_ipi
    real(dp), allocatable :: cmat(:, :)
    real(dp), allocatable :: omega_k(:)
  end type ring_polymer_state

  public :: init_ring_polymer

contains

  subroutine init_ring_polymer(rp, beads, temperature_K, convention)
    type(ring_polymer_state), intent(inout) :: rp
    integer, intent(in) :: beads
    real(dp), intent(in) :: temperature_K
    character(len=*), intent(in), optional :: convention
    character(len=32) :: selected

    selected = ring_polymer_convention_standard_ipi
    if (present(convention)) selected = trim(convention)

    if (beads < 1) error stop "init_ring_polymer: bead count must be positive"
    if (temperature_K < 0.0_dp) error stop "init_ring_polymer: temperature must be non-negative"
    if (trim(selected) /= ring_polymer_convention_standard_ipi) then
      error stop "init_ring_polymer: only standard_ipi convention is implemented"
    end if

    rp%beads = beads
    rp%temperature_K = temperature_K
    rp%convention = selected

    allocate(rp%cmat(beads, beads))
    allocate(rp%omega_k(beads))
    call build_cmat(beads, rp%cmat)
    call compute_omega_k(beads, temperature_K, rp%omega_k)
  end subroutine init_ring_polymer

end module ereaxff_ring_polymer

