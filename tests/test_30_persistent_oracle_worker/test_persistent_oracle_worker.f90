program test_persistent_oracle_worker
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use ereaxff_kinds, only: dp
  use ereaxff_system_state, only: system_state, init_system_state, set_particle
  use ereaxff_charge_coupling_policy, only: charge_memory_state, init_charge_memory, &
      charge_policy_bead_local_acks2
  use ereaxff_oracle_adapter, only: oracle_adapter_state, oracle_adapter_init, oracle_adapter_shutdown
  use ereaxff_bead_force_manager, only: oracle_bead_summary, evaluate_oracle_beads
  implicit none

  integer, parameter :: ndim = 3
  integer, parameter :: nparticles = 61
  integer, parameter :: p = 2
  type(system_state) :: state
  real(dp) :: positions(ndim, nparticles, p)
  type(oracle_bead_summary) :: external_summary, persistent_summary, pool_summary
  type(charge_memory_state) :: external_memory, persistent_memory, pool_memory
  character(len=512) :: v1_root, v2_root, work_root, smoke_binary, worker_binary, geo_path
  real(dp) :: external_solver(nparticles, p)
  real(dp) :: persistent_solver(nparticles, p)
  real(dp) :: pool_solver(nparticles, p)

  call get_environment_variable("EREAXFF_V1_ROOT", v1_root)
  call get_environment_variable("EREAXFF_V2_ROOT", v2_root)
  call get_environment_variable("EREAXFF_ADAPTER_WORK_ROOT", work_root)
  call get_environment_variable("EREAXFF_ORACLE_SMOKE_BINARY", smoke_binary)
  call get_environment_variable("EREAXFF_PERSISTENT_WORKER_BINARY", worker_binary)
  if (len_trim(v1_root) == 0) v1_root = "/Users/yeatanir/PhD/PIMD/ereaxff_pimd/V1"
  if (len_trim(v2_root) == 0) v2_root = "/Users/yeatanir/PhD/PIMD/ereaxff_pimd/V2"
  if (len_trim(work_root) == 0) work_root = "build/test_30_persistent_worker"
  if (len_trim(worker_binary) == 0) worker_binary = "build/persistent_oracle_worker"

  call build_phase4_template(state)
  call build_positions(state, positions)

  call run_mode("external_process", "serial", 1, trim(work_root) // "/external", state, positions, v1_root, v2_root, &
      smoke_binary, worker_binary, external_memory, external_solver, external_summary)
  geo_path = trim(work_root) // "/external/adapter_call_0001/geo"
  call run_mode("persistent_worker", "serial", 1, trim(work_root) // "/persistent", state, positions, v1_root, v2_root, &
      smoke_binary, worker_binary, persistent_memory, persistent_solver, persistent_summary, trim(geo_path))
  call run_mode("persistent_worker_pool", "bounded_pool", 2, trim(work_root) // "/pool", state, positions, v1_root, v2_root, &
      smoke_binary, worker_binary, pool_memory, pool_solver, pool_summary, trim(geo_path))

  call assert_close("energy", external_summary%pimd_potential_eV, persistent_summary%pimd_potential_eV, 1.0e-4_dp)
  call assert_close("forces", maxval(abs(external_summary%scaled_forces_eV_per_A - persistent_summary%scaled_forces_eV_per_A)), &
      0.0_dp, 1.0e-5_dp)
  call assert_close("charges", maxval(abs(external_summary%bead_charges - persistent_summary%bead_charges)), 0.0_dp, 1.0e-5_dp)
  call assert_close("solver", maxval(abs(external_summary%bead_solver_guess - persistent_summary%bead_solver_guess)), &
      0.0_dp, 1.0e-5_dp)
  call assert_close("memory", maxval(abs(external_memory%charges_by_bead - persistent_memory%charges_by_bead)), &
      0.0_dp, 1.0e-5_dp)
  call assert_close("pool energy", persistent_summary%pimd_potential_eV, pool_summary%pimd_potential_eV, 1.0e-4_dp)
  call assert_close("pool forces", maxval(abs(persistent_summary%scaled_forces_eV_per_A - pool_summary%scaled_forces_eV_per_A)), &
      0.0_dp, 1.0e-5_dp)
  call assert_close("pool charges", maxval(abs(persistent_summary%bead_charges - pool_summary%bead_charges)), 0.0_dp, 1.0e-5_dp)
  call assert_close("pool solver", maxval(abs(persistent_summary%bead_solver_guess - pool_summary%bead_solver_guess)), &
      0.0_dp, 1.0e-5_dp)
  call assert_close("pool memory", maxval(abs(persistent_memory%charges_by_bead - pool_memory%charges_by_bead)), &
      0.0_dp, 1.0e-5_dp)
  call assert_true("persistent energy finite", ieee_is_finite(persistent_summary%pimd_potential_eV))
  call assert_true("persistent forces finite", all(ieee_is_finite(persistent_summary%scaled_forces_eV_per_A)))
  call assert_true("persistent timing log", file_exists(trim(work_root) // "/persistent/persistent_oracle_timing.log"))
  call assert_true("pool timing log", file_exists(trim(work_root) // "/pool/persistent_oracle_timing.log"))
  call assert_true("external energy terms", &
      file_exists(trim(work_root) // "/external/adapter_call_0001/single_point_energy_terms.out"))
  call assert_true("external bond orders", &
      file_exists(trim(work_root) // "/external/adapter_call_0001/single_point_bond_orders.out"))
  call assert_true("persistent energy terms", &
      file_exists(trim(work_root) // "/persistent/persistent_worker/adapter_call_0001/single_point_energy_terms.out"))
  call assert_true("persistent bond orders", &
      file_exists(trim(work_root) // "/persistent/persistent_worker/adapter_call_0001/single_point_bond_orders.out"))
  call assert_true("pool energy terms", &
      file_exists(trim(work_root) // "/pool/persistent_worker_001/adapter_call_0001/single_point_energy_terms.out"))
  call assert_true("pool bond orders", &
      file_exists(trim(work_root) // "/pool/persistent_worker_001/adapter_call_0001/single_point_bond_orders.out"))

  write(*,'(a)') "test_30_persistent_oracle_worker: PASS"
  write(*,'(a,es24.15)') "  energy_eV = ", persistent_summary%pimd_potential_eV

contains

  subroutine run_mode(backend_mode, dispatch_mode, max_workers, mode_root, state_in, positions_in, v1_root, v2_root, &
      smoke_binary, worker_binary, memory_out, solver_guess_out, summary_out, geo_path)
    character(len=*), intent(in) :: backend_mode, dispatch_mode, mode_root, v1_root, v2_root, smoke_binary, worker_binary
    integer, intent(in) :: max_workers
    type(system_state), intent(in) :: state_in
    real(dp), intent(in) :: positions_in(:, :, :)
    type(charge_memory_state), intent(out) :: memory_out
    real(dp), intent(out) :: solver_guess_out(:, :)
    type(oracle_bead_summary), intent(out) :: summary_out
    character(len=*), intent(in), optional :: geo_path
    type(oracle_adapter_state) :: adapter

    call execute_command_line("mkdir -p '" // trim(mode_root) // "'")
    if (present(geo_path)) then
      call oracle_adapter_init(adapter, state_in, trim(v1_root), trim(v2_root), trim(mode_root), &
          geo_path=trim(geo_path), backend_mode=trim(backend_mode), dispatch_mode=trim(dispatch_mode), max_workers=max_workers, &
          timing_log=.true., &
          stagger_launch_ms=0, fail_fast=.true., timeout_s=90, repeatability_check=.false., smoke_binary=trim(smoke_binary), &
          persistent_worker_start_timeout_s=300, &
          persistent_worker_binary=trim(worker_binary))
    else
      call oracle_adapter_init(adapter, state_in, trim(v1_root), trim(v2_root), trim(mode_root), &
          backend_mode=trim(backend_mode), dispatch_mode=trim(dispatch_mode), max_workers=max_workers, &
          timing_log=.true., &
          stagger_launch_ms=0, fail_fast=.true., timeout_s=90, repeatability_check=.false., smoke_binary=trim(smoke_binary), &
          persistent_worker_start_timeout_s=300, &
          persistent_worker_binary=trim(worker_binary))
    end if
    call init_charge_memory(memory_out, nparticles, p, charge_policy_bead_local_acks2)
    solver_guess_out = 0.0_dp
    call evaluate_oracle_beads(adapter, positions_in, memory_out, solver_guess_out, &
        charge_policy_bead_local_acks2, summary_out, 7, "initial_force")
    call oracle_adapter_shutdown(adapter)
    call assert_true(trim(backend_mode) // " status zero", all(summary_out%bead_status == 0))
    call assert_true(trim(backend_mode) // " energy finite", ieee_is_finite(summary_out%pimd_potential_eV))
  end subroutine run_mode

  subroutine build_phase4_template(state_out)
    type(system_state), intent(inout) :: state_out
    integer :: ix, iy, unit_index, idx
    real(dp), parameter :: a_cc = 1.4210_dp
    real(dp), parameter :: z_c = 10.0_dp
    real(dp), parameter :: el_x = 7.38300_dp
    real(dp), parameter :: el_y = 5.68433_dp
    real(dp), parameter :: el_z = 10.001_dp
    real(dp) :: ax, ay, x, y
    real(dp), dimension(4, 2) :: unit

    unit = 0.0_dp
    unit(1, :) = [0.0_dp, 0.0_dp]
    unit(2, :) = [0.5_dp, 1.0_dp / 6.0_dp]
    unit(3, :) = [0.5_dp, 3.0_dp / 6.0_dp]
    unit(4, :) = [0.0_dp, 4.0_dp / 6.0_dp]
    ax = sqrt(3.0_dp) * a_cc
    ay = 3.0_dp * a_cc

    call init_system_state(state_out, ndim, nparticles)
    call set_particle(state_out, 1, 1, "El", 0.00054858_dp, [el_x, el_y, el_z])
    idx = 2
    do ix = 0, 4
      do iy = 0, 2
        do unit_index = 1, 4
          x = (real(ix, dp) + unit(unit_index, 1)) * ax
          y = (real(iy, dp) + unit(unit_index, 2)) * ay
          call set_particle(state_out, idx, idx, "C", 12.011_dp, [x, y, z_c])
          idx = idx + 1
        end do
      end do
    end do
  end subroutine build_phase4_template

  subroutine build_positions(state_in, positions_out)
    type(system_state), intent(in) :: state_in
    real(dp), intent(out) :: positions_out(:, :, :)
    integer :: bead

    do bead = 1, size(positions_out, 3)
      positions_out(:, :, bead) = state_in%position_A(:, :)
      positions_out(1, 1, bead) = positions_out(1, 1, bead) + 1.0_dp
      positions_out(2, 1, bead) = positions_out(2, 1, bead) + 1.0_dp
    end do
  end subroutine build_positions

  subroutine assert_true(label, condition)
    character(len=*), intent(in) :: label
    logical, intent(in) :: condition

    if (.not. condition) then
      write(*,'(a)') "FAIL: " // trim(label)
      error stop "persistent worker regression failed"
    end if
  end subroutine assert_true

  subroutine assert_close(label, value, expected, tol)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: value, expected, tol

    if (abs(value - expected) > tol) then
      write(*,'(a)') "FAIL: " // trim(label)
      write(*,'(a,es24.15)') "  value    = ", value
      write(*,'(a,es24.15)') "  expected = ", expected
      write(*,'(a,es24.15)') "  tol      = ", tol
      error stop "persistent worker regression failed"
    end if
  end subroutine assert_close

  logical function file_exists(path) result(exists)
    character(len=*), intent(in) :: path
    inquire(file=trim(path), exist=exists)
  end function file_exists

end program test_persistent_oracle_worker
