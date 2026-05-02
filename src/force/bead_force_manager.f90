module ereaxff_bead_force_manager
  use ereaxff_kinds, only: dp
  use ereaxff_force_provider_base, only: force_result
  use ereaxff_charge_coupling_policy, only: charge_memory_state, &
      charge_policy_is_valid, charge_policy_uses_bead_local_memory, set_bead_charges
  use ereaxff_harmonic_force_provider, only: harmonic_force_provider_state, harmonic_force_provider_evaluate
  use ereaxff_dummy_force_provider, only: dummy_force_provider_state, dummy_force_provider_evaluate
  use ereaxff_oracle_adapter, only: oracle_adapter_state, oracle_adapter_result, oracle_adapter_eval
  implicit none
  private

  type, public :: bead_force_summary
    real(dp) :: pimd_potential_eV = 0.0_dp
    real(dp), allocatable :: bead_energies_eV(:)
    real(dp), allocatable :: scaled_forces_eV_per_A(:, :, :)
  end type bead_force_summary

  type, public :: oracle_bead_summary
    real(dp) :: pimd_potential_eV = 0.0_dp
    real(dp), allocatable :: bead_energies_eV(:)
    real(dp), allocatable :: scaled_forces_eV_per_A(:, :, :)
    real(dp), allocatable :: bead_charges(:, :)
    real(dp), allocatable :: bead_solver_guess(:, :)
    integer, allocatable :: bead_status(:)
  end type oracle_bead_summary

  public :: evaluate_harmonic_beads
  public :: evaluate_dummy_beads
  public :: evaluate_oracle_beads
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

  subroutine evaluate_oracle_beads(provider, positions_A, charge_memory, solver_guess_by_bead, &
      charge_policy, summary)
    type(oracle_adapter_state), intent(inout) :: provider
    real(dp), intent(in) :: positions_A(:, :, :)
    type(charge_memory_state), intent(inout) :: charge_memory
    real(dp), intent(inout) :: solver_guess_by_bead(:, :)
    character(len=*), intent(in) :: charge_policy
    type(oracle_bead_summary), intent(inout) :: summary
    type(oracle_adapter_result) :: result
    integer :: bead, p

    result%force%energy_eV = 0.0_dp
    result%force%status = 0
    result%force%diagnostic = ""
    if (allocated(result%force%forces_eV_per_A)) deallocate(result%force%forces_eV_per_A)
    if (allocated(result%force%charges)) deallocate(result%force%charges)
    if (allocated(result%solver_guess)) deallocate(result%solver_guess)
    if (.not. charge_policy_is_valid(charge_policy)) then
      error stop "evaluate_oracle_beads: invalid charge policy"
    end if
    if (.not. charge_policy_uses_bead_local_memory(charge_policy)) then
      error stop "evaluate_oracle_beads: only bead_local_acks2 is implemented in the manager path"
    end if
    call evaluate_oracle_shapes(positions_A, solver_guess_by_bead, charge_memory, summary)
    p = size(positions_A, 3)

    do bead = 1, p
      call oracle_adapter_eval(provider, positions_A(:, :, bead), &
          charge_memory%charges_by_bead(:, bead), solver_guess_by_bead(:, bead), bead_id=bead, &
          charge_policy=charge_policy, result=result)
      summary%bead_energies_eV(bead) = result%force%energy_eV
      summary%scaled_forces_eV_per_A(:, :, bead) = result%force%forces_eV_per_A / real(p, dp)
      summary%bead_charges(:, bead) = result%force%charges
      summary%bead_solver_guess(:, bead) = result%solver_guess
      summary%bead_status(bead) = result%force%status
      call set_bead_charges(charge_memory, bead, result%force%charges)
    end do
    summary%pimd_potential_eV = sum(summary%bead_energies_eV) / real(p, dp)
  end subroutine evaluate_oracle_beads

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

  subroutine evaluate_oracle_shapes(positions_A, solver_guess_by_bead, charge_memory, summary)
    real(dp), intent(in) :: positions_A(:, :, :)
    real(dp), intent(inout) :: solver_guess_by_bead(:, :)
    type(charge_memory_state), intent(inout) :: charge_memory
    type(oracle_bead_summary), intent(inout) :: summary
    integer :: ndim, nparticles, p

    ndim = size(positions_A, 1)
    nparticles = size(positions_A, 2)
    p = size(positions_A, 3)
    if (ndim < 1) error stop "evaluate_oracle_shapes: ndim must be positive"
    if (nparticles < 1) error stop "evaluate_oracle_shapes: nparticles must be positive"
    if (p < 1) error stop "evaluate_oracle_shapes: bead count must be positive"
    if (size(solver_guess_by_bead, 1) /= nparticles .or. size(solver_guess_by_bead, 2) /= p) then
      error stop "evaluate_oracle_shapes: solver memory shape mismatch"
    end if
    if (.not. allocated(charge_memory%charges_by_bead)) then
      error stop "evaluate_oracle_shapes: charge memory not initialized"
    end if
    if (size(charge_memory%charges_by_bead, 1) /= nparticles .or. &
        size(charge_memory%charges_by_bead, 2) /= p) then
      error stop "evaluate_oracle_shapes: charge memory shape mismatch"
    end if

    if (allocated(summary%bead_energies_eV)) deallocate(summary%bead_energies_eV)
    if (allocated(summary%scaled_forces_eV_per_A)) deallocate(summary%scaled_forces_eV_per_A)
    if (allocated(summary%bead_charges)) deallocate(summary%bead_charges)
    if (allocated(summary%bead_solver_guess)) deallocate(summary%bead_solver_guess)
    if (allocated(summary%bead_status)) deallocate(summary%bead_status)
    allocate(summary%bead_energies_eV(p))
    allocate(summary%scaled_forces_eV_per_A(ndim, nparticles, p))
    allocate(summary%bead_charges(nparticles, p))
    allocate(summary%bead_solver_guess(nparticles, p))
    allocate(summary%bead_status(p))
    summary%bead_energies_eV = 0.0_dp
    summary%scaled_forces_eV_per_A = 0.0_dp
    summary%bead_charges = 0.0_dp
    summary%bead_solver_guess = 0.0_dp
    summary%bead_status = 0
    summary%pimd_potential_eV = 0.0_dp
  end subroutine evaluate_oracle_shapes

end module ereaxff_bead_force_manager
