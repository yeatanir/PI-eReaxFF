module ereaxff_bead_force_manager
  use ereaxff_kinds, only: dp
  use ereaxff_force_provider_base, only: force_result
  use ereaxff_charge_coupling_policy, only: charge_memory_state, &
      charge_policy_is_valid, charge_policy_supports_oracle_manager, &
      charge_policy_uses_centroid_contraction, charge_policy_uses_bead_local_memory, set_bead_charges
  use ereaxff_harmonic_force_provider, only: harmonic_force_provider_state, harmonic_force_provider_evaluate
  use ereaxff_dummy_force_provider, only: dummy_force_provider_state, dummy_force_provider_evaluate
  use ereaxff_oracle_adapter, only: oracle_adapter_state, oracle_adapter_result, oracle_adapter_eval, &
      oracle_adapter_ensure_workers, oracle_adapter_worker_state
  use ereaxff_fast_jax_adapter, only: fast_jax_adapter_state, fast_jax_batch_result, fast_jax_adapter_eval_batch
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

  type :: oracle_call_timing_record
    integer :: step_index = -1
    integer :: batch_index = 0
    integer :: job_index = 0
    integer :: bead_index = 0
    integer :: worker_id = 0
    real(dp) :: start_s = 0.0_dp
    real(dp) :: end_s = 0.0_dp
    real(dp) :: elapsed_s = 0.0_dp
    integer :: status = 0
    character(len=32) :: dispatch_mode = ""
    character(len=32) :: stage_name = ""
    character(len=64) :: stage_label = ""
    character(len=16) :: worker_state = ""
    character(len=512) :: work_dir = ""
  end type oracle_call_timing_record

  public :: evaluate_harmonic_beads
  public :: evaluate_dummy_beads
  public :: evaluate_zero_oracle_beads
  public :: evaluate_oracle_beads
  public :: evaluate_fast_jax_beads
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

  subroutine evaluate_zero_oracle_beads(positions_A, charge_memory, solver_guess_by_bead, summary)
    real(dp), intent(in) :: positions_A(:, :, :)
    type(charge_memory_state), intent(inout) :: charge_memory
    real(dp), intent(inout) :: solver_guess_by_bead(:, :)
    type(oracle_bead_summary), intent(inout) :: summary
    integer :: bead

    call evaluate_oracle_shapes(positions_A, solver_guess_by_bead, charge_memory, summary)
    do bead = 1, size(positions_A, 3)
      call set_bead_charges(charge_memory, bead, summary%bead_charges(:, bead))
    end do
  end subroutine evaluate_zero_oracle_beads

  subroutine evaluate_oracle_beads(provider, positions_A, charge_memory, solver_guess_by_bead, &
      charge_policy, summary, step_index, stage_label)
    type(oracle_adapter_state), intent(inout) :: provider
    real(dp), intent(in) :: positions_A(:, :, :)
    type(charge_memory_state), intent(inout) :: charge_memory
    real(dp), intent(inout) :: solver_guess_by_bead(:, :)
    character(len=*), intent(in) :: charge_policy
    type(oracle_bead_summary), intent(inout) :: summary
    integer, intent(in), optional :: step_index
    character(len=*), intent(in), optional :: stage_label
    type(oracle_adapter_result) :: result
    type(oracle_adapter_result), allocatable :: bead_results(:)
    type(oracle_adapter_state), allocatable :: bead_providers(:)
    type(oracle_call_timing_record), allocatable :: call_records(:)
    real(dp), allocatable :: batch_elapsed(:)
    integer :: call_base, clock_base, clock_rate
    integer :: bead, p, batch_index, batch_start, batch_end, batch_size
    integer :: step_value, failed_jobs
    character(len=32) :: dispatch_mode
    character(len=64) :: stage_value

    result%force%energy_eV = 0.0_dp
    result%force%status = 0
    result%force%diagnostic = ""
    if (allocated(result%force%forces_eV_per_A)) deallocate(result%force%forces_eV_per_A)
    if (allocated(result%force%charges)) deallocate(result%force%charges)
    if (allocated(result%solver_guess)) deallocate(result%solver_guess)
    if (.not. charge_policy_is_valid(charge_policy)) then
      error stop "evaluate_oracle_beads: invalid charge policy"
    end if
    call evaluate_oracle_shapes(positions_A, solver_guess_by_bead, charge_memory, summary)
    p = size(positions_A, 3)
    if (.not. charge_policy_supports_oracle_manager(charge_policy)) then
      error stop "evaluate_oracle_beads: oracle manager supports only bead_local_acks2 and centroid_acks2"
    end if

    if (present(step_index)) then
      step_value = step_index
    else
      step_value = -1
    end if
    if (present(stage_label)) then
      stage_value = trim(stage_label)
    else
      stage_value = "unspecified"
    end if
    dispatch_mode = trim(adapter_dispatch_mode(provider))
    call system_clock(clock_base, clock_rate)
    if (clock_rate <= 0) error stop "evaluate_oracle_beads: system clock rate unavailable"

    if (charge_policy_uses_bead_local_memory(charge_policy)) then
      allocate(bead_results(p))
      allocate(bead_providers(p))
      allocate(call_records(p))
      allocate(batch_elapsed(p))
      call_base = provider%call_count
      provider%call_count = provider%call_count + p
      if (trim(provider%backend_mode) == "persistent_worker_pool") call oracle_adapter_ensure_workers(provider)
      do bead = 1, p
        bead_providers(bead) = provider
      end do
      select case (dispatch_mode)
      case ("serial")
        batch_index = 0
        do bead = 1, p
          batch_index = batch_index + 1
          if (provider%stagger_launch_ms > 0 .and. bead > 1) call wait_milliseconds(provider%stagger_launch_ms)
          call run_oracle_call(provider, positions_A(:, :, bead), &
              charge_memory%charges_by_bead(:, bead), solver_guess_by_bead(:, bead), bead, call_base + bead, &
              charge_policy, "bead_local", dispatch_mode, batch_index, step_value, stage_value, call_records(bead), &
              bead_results(bead), 1, clock_base, clock_rate)
          batch_elapsed(1) = call_records(bead)%elapsed_s
          call finalize_oracle_batch(provider, step_value, batch_index, "bead_local", dispatch_mode, stage_value, 1, &
              bead_results(bead:bead), call_records(bead:bead), batch_elapsed(1:1), failed_jobs)
          if (provider%fail_fast .and. failed_jobs > 0) error stop "evaluate_oracle_beads: oracle call failed"
        end do
      case ("openmp_all")
        batch_index = 1
        batch_start = 1
        batch_end = p
        batch_size = p
!$omp parallel do default(none) private(bead) shared(batch_start, batch_end, bead_providers, positions_A, &
!$omp& charge_memory, solver_guess_by_bead, charge_policy, bead_results, call_records, call_base, step_value, clock_base, &
!$omp& clock_rate, dispatch_mode, batch_index, stage_value)
        do bead = batch_start, batch_end
          call run_oracle_call(bead_providers(bead), positions_A(:, :, bead), &
              charge_memory%charges_by_bead(:, bead), solver_guess_by_bead(:, bead), bead, call_base + bead, &
              charge_policy, "bead_local", dispatch_mode, batch_index, step_value, stage_value, call_records(bead), &
              bead_results(bead), bead - batch_start + 1, &
              clock_base, clock_rate)
        end do
!$omp end parallel do
        do bead = batch_start, batch_end
          batch_elapsed(bead) = call_records(bead)%elapsed_s
        end do
        call finalize_oracle_batch(provider, step_value, batch_index, "bead_local", dispatch_mode, stage_value, batch_size, &
            bead_results(batch_start:batch_end), call_records(batch_start:batch_end), batch_elapsed(batch_start:batch_end), &
            failed_jobs)
      case ("bounded_pool")
        batch_index = 0
        batch_start = 1
        do while (batch_start <= p)
          batch_index = batch_index + 1
          batch_size = min(provider%max_workers, p - batch_start + 1)
          batch_end = batch_start + batch_size - 1
          if (provider%stagger_launch_ms > 0 .and. batch_index > 1) call wait_milliseconds(provider%stagger_launch_ms)
          if (batch_size > 1) then
!$omp parallel do default(none) private(bead) shared(batch_start, batch_end, bead_providers, positions_A, &
!$omp& charge_memory, solver_guess_by_bead, charge_policy, bead_results, call_records, call_base, step_value, clock_base, &
!$omp& clock_rate, dispatch_mode, batch_index, stage_value)
            do bead = batch_start, batch_end
              call run_oracle_call(bead_providers(bead), positions_A(:, :, bead), &
                  charge_memory%charges_by_bead(:, bead), solver_guess_by_bead(:, bead), bead, call_base + bead, &
              charge_policy, "bead_local", dispatch_mode, batch_index, step_value, stage_value, call_records(bead), &
              bead_results(bead), &
              bead - batch_start + 1, &
              clock_base, clock_rate)
            end do
!$omp end parallel do
          else
            bead = batch_start
              call run_oracle_call(bead_providers(bead), positions_A(:, :, bead), charge_memory%charges_by_bead(:, bead), &
                  solver_guess_by_bead(:, bead), bead, call_base + bead, charge_policy, "bead_local", dispatch_mode, batch_index, &
                step_value, stage_value, call_records(bead), bead_results(bead), 1, clock_base, clock_rate)
          end if
          do bead = batch_start, batch_end
            batch_elapsed(bead) = call_records(bead)%elapsed_s
          end do
          call finalize_oracle_batch(provider, step_value, batch_index, "bead_local", dispatch_mode, stage_value, batch_size, &
              bead_results(batch_start:batch_end), call_records(batch_start:batch_end), &
              batch_elapsed(batch_start:batch_end), failed_jobs)
          if (provider%fail_fast .and. failed_jobs > 0) then
            error stop "evaluate_oracle_beads: oracle call failed"
          end if
          batch_start = batch_end + 1
        end do
      case default
        error stop "evaluate_oracle_beads: unsupported dispatch mode"
      end select
      do bead = 1, p
        summary%bead_energies_eV(bead) = bead_results(bead)%force%energy_eV
        summary%scaled_forces_eV_per_A(:, :, bead) = bead_results(bead)%force%forces_eV_per_A / real(p, dp)
        summary%bead_charges(:, bead) = bead_results(bead)%force%charges
        summary%bead_solver_guess(:, bead) = bead_results(bead)%solver_guess
        summary%bead_status(bead) = bead_results(bead)%force%status
        call set_bead_charges(charge_memory, bead, bead_results(bead)%force%charges)
      end do
      summary%pimd_potential_eV = sum(summary%bead_energies_eV) / real(p, dp)
    else if (charge_policy_uses_centroid_contraction(charge_policy)) then
      call_base = provider%call_count
      provider%call_count = provider%call_count + 1
      call evaluate_oracle_centroid(provider, positions_A, charge_memory, solver_guess_by_bead, charge_policy, summary, &
          result, step_value, call_base, dispatch_mode, stage_value, clock_base, clock_rate)
    end if
  end subroutine evaluate_oracle_beads

  subroutine evaluate_fast_jax_beads(provider, positions_A, charge_memory, solver_guess_by_bead, &
      charge_policy, summary)
    type(fast_jax_adapter_state), intent(inout) :: provider
    real(dp), intent(in) :: positions_A(:, :, :)
    type(charge_memory_state), intent(inout) :: charge_memory
    real(dp), intent(inout) :: solver_guess_by_bead(:, :)
    character(len=*), intent(in) :: charge_policy
    type(oracle_bead_summary), intent(inout) :: summary
    type(fast_jax_batch_result) :: result
    integer :: bead, p

    if (.not. charge_policy_is_valid(charge_policy)) then
      error stop "evaluate_fast_jax_beads: invalid charge policy"
    end if
    if (.not. charge_policy_uses_bead_local_memory(charge_policy)) then
      error stop "evaluate_fast_jax_beads: fast backend currently supports only bead-local charge memory"
    end if
    call evaluate_oracle_shapes(positions_A, solver_guess_by_bead, charge_memory, summary)
    p = size(positions_A, 3)
    call fast_jax_adapter_eval_batch(provider, positions_A, charge_memory%charges_by_bead, solver_guess_by_bead, result)
    do bead = 1, p
      summary%bead_energies_eV(bead) = result%bead_energies_eV(bead)
      summary%scaled_forces_eV_per_A(:, :, bead) = result%forces_eV_per_A(:, :, bead) / real(p, dp)
      summary%bead_charges(:, bead) = result%charges(:, bead)
      summary%bead_solver_guess(:, bead) = result%solver_guess(:, bead)
      summary%bead_status(bead) = result%bead_status(bead)
      call set_bead_charges(charge_memory, bead, result%charges(:, bead))
      solver_guess_by_bead(:, bead) = result%solver_guess(:, bead)
    end do
    summary%pimd_potential_eV = sum(summary%bead_energies_eV) / real(p, dp)
  end subroutine evaluate_fast_jax_beads

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

  pure function adapter_dispatch_mode(provider) result(mode)
    type(oracle_adapter_state), intent(in) :: provider
    character(len=32) :: mode

    if (allocated(provider%backend_mode) .and. trim(provider%backend_mode) == "persistent_worker") then
      mode = "serial"
      return
    else if (allocated(provider%backend_mode) .and. trim(provider%backend_mode) == "persistent_worker_pool") then
      mode = "bounded_pool"
      return
    end if
    if (allocated(provider%dispatch_mode) .and. len_trim(provider%dispatch_mode) > 0) then
      mode = trim(provider%dispatch_mode)
    else
      mode = "openmp_all"
    end if
  end function adapter_dispatch_mode

  subroutine run_oracle_call(provider, positions_A, charge_in, solver_guess_in, bead_index, call_slot, charge_policy, &
      stage_name, dispatch_mode, batch_index, step_index, stage_label, record, result, worker_id, clock_base, clock_rate)
    type(oracle_adapter_state), intent(inout) :: provider
    real(dp), intent(in) :: positions_A(:, :)
    real(dp), intent(in) :: charge_in(:)
    real(dp), intent(in) :: solver_guess_in(:)
    integer, intent(in) :: bead_index, call_slot, batch_index, step_index, worker_id, clock_base, clock_rate
    character(len=*), intent(in) :: charge_policy, stage_name, dispatch_mode, stage_label
    type(oracle_call_timing_record), intent(out) :: record
    type(oracle_adapter_result), intent(inout) :: result
    integer :: start_count, end_count

    call system_clock(start_count)
    call oracle_adapter_eval(provider, positions_A, charge_in, solver_guess_in, bead_id=bead_index, &
        charge_policy=charge_policy, result=result, call_slot=call_slot)
    call system_clock(end_count)

    record%step_index = step_index
    record%batch_index = batch_index
    record%job_index = bead_index
    record%bead_index = bead_index
    record%worker_id = worker_id
    record%start_s = real(start_count - clock_base, dp) / real(clock_rate, dp)
    record%end_s = real(end_count - clock_base, dp) / real(clock_rate, dp)
    record%elapsed_s = max(0.0_dp, record%end_s - record%start_s)
    record%status = result%force%status
    record%dispatch_mode = dispatch_mode
    record%stage_name = stage_name
    record%stage_label = stage_label
    record%worker_state = oracle_adapter_worker_state(provider)
    record%work_dir = trim(provider%work_root) // "/adapter_call_" // pad_int(call_slot, 4)
  end subroutine run_oracle_call

  subroutine finalize_oracle_batch(provider, step_index, batch_index, stage_name, dispatch_mode, stage_label, n_jobs, &
      results, records, elapsed, failed_jobs)
    type(oracle_adapter_state), intent(inout) :: provider
    integer, intent(in) :: step_index, batch_index, n_jobs
    character(len=*), intent(in) :: stage_name, dispatch_mode, stage_label
    type(oracle_adapter_result), intent(in) :: results(:)
    type(oracle_call_timing_record), intent(in) :: records(:)
    real(dp), intent(in) :: elapsed(:)
    integer, intent(out) :: failed_jobs
    integer :: i
    real(dp) :: min_elapsed, max_elapsed, mean_elapsed, batch_start_s, batch_end_s, total_elapsed

    failed_jobs = 0
    if (size(records) /= n_jobs .or. size(results) /= n_jobs .or. size(elapsed) /= n_jobs) then
      error stop "finalize_oracle_batch: size mismatch"
    end if
    if (n_jobs < 1) error stop "finalize_oracle_batch: batch must be non-empty"
    min_elapsed = elapsed(1)
    max_elapsed = elapsed(1)
    mean_elapsed = 0.0_dp
    batch_start_s = records(1)%start_s
    batch_end_s = records(1)%end_s
    do i = 1, n_jobs
      min_elapsed = min(min_elapsed, elapsed(i))
      max_elapsed = max(max_elapsed, elapsed(i))
      mean_elapsed = mean_elapsed + elapsed(i)
      batch_start_s = min(batch_start_s, records(i)%start_s)
      batch_end_s = max(batch_end_s, records(i)%end_s)
      if (results(i)%force%status /= 0) failed_jobs = failed_jobs + 1
    end do
    mean_elapsed = mean_elapsed / real(n_jobs, dp)
    total_elapsed = batch_end_s - batch_start_s
    call append_oracle_timing_log(provider, step_index, batch_index, stage_name, dispatch_mode, stage_label, n_jobs, &
        records, min_elapsed, mean_elapsed, max_elapsed, failed_jobs)
    if (provider%timing_log) then
      write(*,*) "oracle batch", trim(dispatch_mode), "batch", batch_index, "jobs", n_jobs, "step", step_index, &
          "stage", trim(stage_label), "elapsed_s", total_elapsed, "min", min_elapsed, "mean", mean_elapsed, "max", max_elapsed, &
          "failed", failed_jobs
    end if
  end subroutine finalize_oracle_batch

  subroutine append_oracle_timing_log(provider, step_index, batch_index, stage_name, dispatch_mode, stage_label, n_jobs, &
      records, min_elapsed, mean_elapsed, max_elapsed, failed_jobs)
    type(oracle_adapter_state), intent(in) :: provider
    integer, intent(in) :: step_index, batch_index, n_jobs, failed_jobs
    character(len=*), intent(in) :: stage_name, dispatch_mode, stage_label
    type(oracle_call_timing_record), intent(in) :: records(:)
    real(dp), intent(in) :: min_elapsed, mean_elapsed, max_elapsed
    integer :: unit, ios, i
    logical :: exists
    character(len=512) :: path
    real(dp) :: batch_start_s, batch_end_s

    if (.not. provider%timing_log) return
    path = trim(provider%timing_log_path)
    if (len_trim(path) == 0) return
    inquire(file=trim(path), exist=exists)
    if (.not. exists) then
      open(newunit=unit, file=trim(path), status="replace", action="write", iostat=ios)
      if (ios /= 0) error stop "append_oracle_timing_log: unable to create timing log"
      write(unit,'(a)') "# kind step batch job bead worker dispatch_mode stage_kind stage_label worker_state"
      write(unit,'(a)') "# n_jobs start_s end_s elapsed_s status workdir"
      write(unit,'(a)') "# BATCH records summarize n_jobs and batch timing; CALL records capture each oracle call."
      close(unit)
    end if
    open(newunit=unit, file=trim(path), status="unknown", action="write", position="append", iostat=ios)
    if (ios /= 0) error stop "append_oracle_timing_log: unable to append timing log"
    batch_start_s = records(1)%start_s
    batch_end_s = records(1)%end_s
    do i = 1, size(records)
      batch_start_s = min(batch_start_s, records(i)%start_s)
      batch_end_s = max(batch_end_s, records(i)%end_s)
    end do
    write(unit,'(a,1x,i0,1x,i0,1x,a,1x,a,1x,i0,1x,es24.15,1x,es24.15,1x,es24.15,1x,es24.15,1x,es24.15,1x,es24.15,1x,i0)') &
        "BATCH", step_index, batch_index, trim(dispatch_mode), trim(stage_name) // ":" // trim(stage_label), n_jobs, &
        batch_start_s, batch_end_s, batch_end_s - batch_start_s, min_elapsed, mean_elapsed, max_elapsed, failed_jobs
    do i = 1, size(records)
      write(unit,'(a,1x,i0,1x,i0,1x,i0,1x,i0,1x,i0,1x,a,1x,a,1x,a,1x,a,1x,es24.15,1x,es24.15,1x,es24.15,1x,i0,1x,a)') &
          "CALL", records(i)%step_index, records(i)%batch_index, records(i)%job_index, records(i)%bead_index, &
          records(i)%worker_id, trim(records(i)%dispatch_mode), trim(records(i)%stage_name), trim(records(i)%stage_label), &
          trim(records(i)%worker_state), &
          records(i)%start_s, records(i)%end_s, records(i)%elapsed_s, records(i)%status, trim(log_token(records(i)%work_dir))
    end do
    close(unit)
  end subroutine append_oracle_timing_log

  subroutine wait_milliseconds(ms)
    integer, intent(in) :: ms
    integer :: start_count, now_count, rate, target_count

    if (ms <= 0) return
    call system_clock(start_count, rate)
    if (rate <= 0) return
    target_count = start_count + ceiling(real(ms, dp) * real(rate, dp) / 1000.0_dp)
    do
      call system_clock(now_count)
      if (now_count >= target_count) exit
    end do
  end subroutine wait_milliseconds

  subroutine evaluate_oracle_centroid(provider, positions_A, charge_memory, solver_guess_by_bead, &
      charge_policy, summary, result, step_index, call_base, dispatch_mode, stage_label, clock_base, clock_rate)
    type(oracle_adapter_state), intent(inout) :: provider
    real(dp), intent(in) :: positions_A(:, :, :)
    type(charge_memory_state), intent(inout) :: charge_memory
    real(dp), intent(inout) :: solver_guess_by_bead(:, :)
    character(len=*), intent(in) :: charge_policy
    type(oracle_bead_summary), intent(inout) :: summary
    type(oracle_adapter_result), intent(inout) :: result
    integer, intent(in) :: step_index, call_base, clock_base, clock_rate
    character(len=*), intent(in) :: dispatch_mode, stage_label
    real(dp), allocatable :: centroid_positions(:, :)
    real(dp), allocatable :: centroid_charges(:)
    real(dp), allocatable :: centroid_solver_guess(:)
    type(oracle_adapter_result), allocatable :: centroid_results(:)
    type(oracle_call_timing_record), allocatable :: centroid_records(:)
    real(dp), allocatable :: centroid_elapsed(:)
    integer :: bead, p
    type(oracle_call_timing_record) :: call_record
    integer :: failed_jobs

    p = size(positions_A, 3)
    allocate(centroid_positions(size(positions_A, 1), size(positions_A, 2)))
    allocate(centroid_charges(size(charge_memory%charges_by_bead, 1)))
    allocate(centroid_solver_guess(size(solver_guess_by_bead, 1)))
    centroid_positions = sum(positions_A, dim=3) / real(p, dp)
    centroid_charges = sum(charge_memory%charges_by_bead, dim=2) / real(p, dp)
    centroid_solver_guess = sum(solver_guess_by_bead, dim=2) / real(p, dp)
    allocate(centroid_results(1))
    allocate(centroid_records(1))
    allocate(centroid_elapsed(1))

    call run_oracle_call(provider, centroid_positions, centroid_charges, centroid_solver_guess, 1, call_base + 1, &
        charge_policy, "centroid", dispatch_mode, 1, step_index, stage_label, call_record, result, 1, clock_base, clock_rate)
    centroid_results(1) = result
    centroid_records(1) = call_record
    centroid_elapsed(1) = call_record%elapsed_s
    call finalize_oracle_batch(provider, step_index, 1, "centroid", dispatch_mode, stage_label, 1, centroid_results, &
        centroid_records, centroid_elapsed, failed_jobs)

    do bead = 1, p
      summary%bead_energies_eV(bead) = result%force%energy_eV
      summary%scaled_forces_eV_per_A(:, :, bead) = result%force%forces_eV_per_A / real(p, dp)
      summary%bead_charges(:, bead) = result%force%charges
      summary%bead_solver_guess(:, bead) = result%solver_guess
      summary%bead_status(bead) = result%force%status
      call set_bead_charges(charge_memory, bead, result%force%charges)
      solver_guess_by_bead(:, bead) = result%solver_guess
    end do
    summary%pimd_potential_eV = result%force%energy_eV
  end subroutine evaluate_oracle_centroid

  pure function pad_int(value, width) result(text)
    integer, intent(in) :: value, width
    character(len=32) :: text

    write(text,'(i0)') value
    text = adjustl(text)
    if (len_trim(text) < width) then
      text = repeat("0", width - len_trim(text)) // trim(text)
    end if
  end function pad_int

  pure function log_token(text_in) result(text_out)
    character(len=*), intent(in) :: text_in
    character(len=len_trim(text_in)) :: text_out
    integer :: i

    text_out = trim(text_in)
    do i = 1, len(text_out)
      if (text_out(i:i) == " ") text_out(i:i) = "_"
    end do
  end function log_token

end module ereaxff_bead_force_manager
