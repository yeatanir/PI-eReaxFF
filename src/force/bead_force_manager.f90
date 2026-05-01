module ereaxff_bead_force_manager
  use ereaxff_kinds, only: dp
  use ereaxff_force_provider_base, only: force_result
  use ereaxff_harmonic_force_provider, only: harmonic_force_provider_state, harmonic_force_provider_evaluate
  use ereaxff_dummy_force_provider, only: dummy_force_provider_state, dummy_force_provider_evaluate
  implicit none
  private

  type, public :: bead_force_summary
    real(dp) :: pimd_potential_eV = 0.0_dp
    real(dp), allocatable :: bead_energies_eV(:)
    real(dp), allocatable :: scaled_forces_eV_per_A(:, :, :)
  end type bead_force_summary

  public :: evaluate_harmonic_beads
  public :: evaluate_dummy_beads
  public :: average_forces_over_beads

contains

  subroutine evaluate_harmonic_beads(provider, positions_A, summary)
    type(harmonic_force_provider_state), intent(in) :: provider
    real(dp), intent(in) :: positions_A(:, :, :)
    type(bead_force_summary), intent(inout) :: summary

    call evaluate_bead_shapes(positions_A, summary)
    call evaluate_harmonic_loop(provider, positions_A, summary)
  end subroutine evaluate_harmonic_beads

  subroutine evaluate_dummy_beads(provider, positions_A, summary)
    type(dummy_force_provider_state), intent(in) :: provider
    real(dp), intent(in) :: positions_A(:, :, :)
    type(bead_force_summary), intent(inout) :: summary
    type(force_result) :: result
    integer :: bead, p

    result%energy_eV = 0.0_dp
    result%status = 0
    result%diagnostic = ""

    call evaluate_bead_shapes(positions_A, summary)
    p = size(positions_A, 3)

    do bead = 1, p
      call dummy_force_provider_evaluate(provider, positions_A(:, :, bead), result)
      summary%bead_energies_eV(bead) = result%energy_eV
      summary%scaled_forces_eV_per_A(:, :, bead) = result%forces_eV_per_A / real(p, dp)
    end do
    summary%pimd_potential_eV = sum(summary%bead_energies_eV) / real(p, dp)
  end subroutine evaluate_dummy_beads

  subroutine average_forces_over_beads(bead_forces, mean_forces)
    real(dp), intent(in) :: bead_forces(:, :, :)
    real(dp), intent(out) :: mean_forces(:, :)

    if (size(mean_forces, 1) /= size(bead_forces, 1) .or. &
        size(mean_forces, 2) /= size(bead_forces, 2)) then
      error stop "average_forces_over_beads: output shape mismatch"
    end if

    mean_forces = sum(bead_forces, dim=3)
  end subroutine average_forces_over_beads

  subroutine evaluate_bead_shapes(positions_A, summary)
    real(dp), intent(in) :: positions_A(:, :, :)
    type(bead_force_summary), intent(inout) :: summary
    integer :: ndim, nparticles, p

    ndim = size(positions_A, 1)
    nparticles = size(positions_A, 2)
    p = size(positions_A, 3)
    if (ndim < 1) error stop "evaluate_bead_shapes: ndim must be positive"
    if (nparticles < 1) error stop "evaluate_bead_shapes: nparticles must be positive"
    if (p < 1) error stop "evaluate_bead_shapes: bead count must be positive"

    if (allocated(summary%bead_energies_eV)) deallocate(summary%bead_energies_eV)
    if (allocated(summary%scaled_forces_eV_per_A)) deallocate(summary%scaled_forces_eV_per_A)
    allocate(summary%bead_energies_eV(p))
    allocate(summary%scaled_forces_eV_per_A(ndim, nparticles, p))
    summary%bead_energies_eV = 0.0_dp
    summary%scaled_forces_eV_per_A = 0.0_dp
    summary%pimd_potential_eV = 0.0_dp
  end subroutine evaluate_bead_shapes

  subroutine evaluate_harmonic_loop(provider, positions_A, summary)
    type(harmonic_force_provider_state), intent(in) :: provider
    real(dp), intent(in) :: positions_A(:, :, :)
    type(bead_force_summary), intent(inout) :: summary
    type(force_result) :: result
    integer :: bead, p

    result%energy_eV = 0.0_dp
    result%status = 0
    result%diagnostic = ""

    p = size(positions_A, 3)
    do bead = 1, p
      call harmonic_force_provider_evaluate(provider, positions_A(:, :, bead), result)
      summary%bead_energies_eV(bead) = result%energy_eV
      summary%scaled_forces_eV_per_A(:, :, bead) = result%forces_eV_per_A / real(p, dp)
    end do
    summary%pimd_potential_eV = sum(summary%bead_energies_eV) / real(p, dp)
  end subroutine evaluate_harmonic_loop

end module ereaxff_bead_force_manager
