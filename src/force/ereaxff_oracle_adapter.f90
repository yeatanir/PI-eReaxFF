module ereaxff_oracle_adapter
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use ereaxff_kinds, only: dp
  use ereaxff_force_provider_base, only: force_result, allocate_force_result
  use ereaxff_system_state, only: system_state
  implicit none
  private

  integer, parameter :: shell_path_len = 16384
  real(dp), parameter :: kcal2eV = 0.043364_dp

  type, public :: oracle_adapter_state
    logical :: initialized = .false.
    integer :: call_count = 0
    type(system_state) :: template
    character(len=:), allocatable :: v1_root
    character(len=:), allocatable :: v2_root
    character(len=:), allocatable :: work_root
    character(len=:), allocatable :: smoke_binary
    character(len=:), allocatable :: repeat_binary
    character(len=:), allocatable :: make_geo_script
    character(len=:), allocatable :: timeout_runner
    character(len=:), allocatable :: geo_path
    character(len=:), allocatable :: control_path
    character(len=:), allocatable :: ffield_path
    character(len=:), allocatable :: iopt_path
    character(len=:), allocatable :: backend_mode
    character(len=:), allocatable :: dispatch_mode
    character(len=:), allocatable :: timing_log_path
    character(len=:), allocatable :: persistent_worker_binary
    character(len=:), allocatable :: persistent_worker_root
    character(len=:), allocatable :: persistent_request_path
    character(len=:), allocatable :: persistent_response_path
    character(len=:), allocatable :: persistent_ready_path
    character(len=:), allocatable :: persistent_shutdown_path
    character(len=:), allocatable :: persistent_done_path
    character(len=:), allocatable :: persistent_pid_path
    character(len=:), allocatable :: persistent_log_path
    integer :: max_workers = 1
    integer :: stagger_launch_ms = 0
    integer :: timeout_s = 90
    integer :: persistent_worker_poll_ms = 25
    integer :: persistent_worker_start_timeout_s = 60
    integer :: persistent_worker_request_count = 0
    logical :: timing_log = .false.
    logical :: fail_fast = .true.
    logical :: repeatability_check = .true.
    logical :: persistent_worker_started = .false.
  end type oracle_adapter_state

  type, public :: oracle_adapter_result
    type(force_result) :: force
    real(dp), allocatable :: solver_guess(:)
    character(len=512) :: legacy_summary = ""
    character(len=512) :: diagnostic = ""
  end type oracle_adapter_result

  public :: oracle_adapter_init
  public :: oracle_adapter_eval
  public :: oracle_adapter_shutdown
  public :: oracle_adapter_worker_state
  public :: oracle_adapter_smoke

contains

  subroutine oracle_adapter_init(adapter, template, v1_root, v2_root, work_root, control_path, ffield_path, iopt_path, &
      smoke_binary, repeat_binary, make_geo_script, timeout_runner, geo_path, backend_mode, dispatch_mode, max_workers, &
      timing_log, stagger_launch_ms, fail_fast, timeout_s, repeatability_check, persistent_worker_binary, &
      persistent_worker_poll_ms, persistent_worker_start_timeout_s)
    type(oracle_adapter_state), intent(inout) :: adapter
    type(system_state), intent(in) :: template
    character(len=*), intent(in) :: v1_root, v2_root, work_root
    character(len=*), intent(in), optional :: control_path, ffield_path, iopt_path
    character(len=*), intent(in), optional :: smoke_binary, repeat_binary, make_geo_script, timeout_runner, geo_path
    character(len=*), intent(in), optional :: backend_mode, persistent_worker_binary
    character(len=*), intent(in), optional :: dispatch_mode
    integer, intent(in), optional :: max_workers, stagger_launch_ms, timeout_s
    integer, intent(in), optional :: persistent_worker_poll_ms, persistent_worker_start_timeout_s
    logical, intent(in), optional :: timing_log, fail_fast, repeatability_check

    call assert_nonempty(v1_root, "oracle_adapter_init: v1_root")
    call assert_nonempty(v2_root, "oracle_adapter_init: v2_root")
    call assert_nonempty(work_root, "oracle_adapter_init: work_root")
    if (template%nparticles < 1) error stop "oracle_adapter_init: template state is empty"

    adapter%template = template
    adapter%v1_root = trim(v1_root)
    adapter%v2_root = trim(v2_root)
    adapter%work_root = trim(work_root)
    call ensure_workdir(trim(work_root))
    if (present(smoke_binary) .and. len_trim(smoke_binary) > 0) then
      adapter%smoke_binary = trim(smoke_binary)
    else
      adapter%smoke_binary = trim(v2_root) // "/phase4_single_point_oracle"
    end if
    if (present(repeat_binary) .and. len_trim(repeat_binary) > 0) then
      adapter%repeat_binary = trim(repeat_binary)
    else
      adapter%repeat_binary = trim(v2_root) // "/single_point_repeatability"
    end if
    if (present(make_geo_script) .and. len_trim(make_geo_script) > 0) then
      adapter%make_geo_script = trim(make_geo_script)
    else
      adapter%make_geo_script = trim(v2_root) // "/tests/regression/test_phase4_fixed_graphene/make_geo.py"
    end if
    if (present(timeout_runner) .and. len_trim(timeout_runner) > 0) then
      adapter%timeout_runner = trim(timeout_runner)
    else
      adapter%timeout_runner = trim(v2_root) // "/tests/regression/test_force_repeatability/run_with_timeout.py"
    end if
    if (present(geo_path) .and. len_trim(geo_path) > 0) then
      adapter%geo_path = trim(geo_path)
    else
      adapter%geo_path = ""
    end if
    if (present(backend_mode) .and. len_trim(backend_mode) > 0) then
      adapter%backend_mode = trim(backend_mode)
    else
      adapter%backend_mode = "external_process"
    end if
    if (present(repeatability_check)) then
      adapter%repeatability_check = repeatability_check
    else
      adapter%repeatability_check = .true.
    end if
    if (present(control_path) .and. len_trim(control_path) > 0) then
      adapter%control_path = trim(control_path)
    else
      adapter%control_path = trim(v1_root) // "/tests/T4_60el/control"
    end if
    if (present(ffield_path) .and. len_trim(ffield_path) > 0) then
      adapter%ffield_path = trim(ffield_path)
    else
      adapter%ffield_path = trim(v1_root) // "/tests/T4_60el/ffield"
    end if
    if (present(iopt_path) .and. len_trim(iopt_path) > 0) then
      adapter%iopt_path = trim(iopt_path)
    else
      adapter%iopt_path = trim(v1_root) // "/tests/reference/iopt"
    end if
    if (present(dispatch_mode) .and. len_trim(dispatch_mode) > 0) then
      adapter%dispatch_mode = trim(dispatch_mode)
    else
      adapter%dispatch_mode = "openmp_all"
    end if
    if (present(persistent_worker_binary) .and. len_trim(persistent_worker_binary) > 0) then
      adapter%persistent_worker_binary = trim(persistent_worker_binary)
    else
      adapter%persistent_worker_binary = "build/persistent_oracle_worker"
    end if
    if (present(max_workers)) then
      adapter%max_workers = max_workers
    else
      adapter%max_workers = 1
    end if
    if (present(timing_log)) then
      adapter%timing_log = timing_log
    else
      adapter%timing_log = .false.
    end if
    if (present(stagger_launch_ms)) then
      adapter%stagger_launch_ms = stagger_launch_ms
    else
      adapter%stagger_launch_ms = 0
    end if
    if (present(fail_fast)) then
      adapter%fail_fast = fail_fast
    else
      adapter%fail_fast = .true.
    end if
    if (present(timeout_s)) then
      adapter%timeout_s = timeout_s
    else
      adapter%timeout_s = 90
    end if
    if (present(persistent_worker_poll_ms)) then
      adapter%persistent_worker_poll_ms = persistent_worker_poll_ms
    else
      adapter%persistent_worker_poll_ms = 25
    end if
    if (present(persistent_worker_start_timeout_s)) then
      adapter%persistent_worker_start_timeout_s = persistent_worker_start_timeout_s
    else
      adapter%persistent_worker_start_timeout_s = 60
    end if
    if (trim(adapter%backend_mode) == "persistent_worker") then
      adapter%timing_log_path = trim(work_root) // "/persistent_oracle_timing.log"
      adapter%persistent_worker_root = trim(work_root) // "/persistent_worker"
      adapter%persistent_request_path = trim(adapter%persistent_worker_root) // "/persistent_request.in"
      adapter%persistent_response_path = trim(adapter%persistent_worker_root) // "/persistent_response.ready"
      adapter%persistent_ready_path = trim(adapter%persistent_worker_root) // "/persistent_worker.ready"
      adapter%persistent_shutdown_path = trim(adapter%persistent_worker_root) // "/persistent_shutdown.flag"
      adapter%persistent_done_path = trim(adapter%persistent_worker_root) // "/persistent_worker.done"
      adapter%persistent_pid_path = trim(adapter%persistent_worker_root) // "/persistent_worker.pid"
      adapter%persistent_log_path = trim(adapter%persistent_worker_root) // "/persistent_worker.log"
    else
      adapter%timing_log_path = trim(work_root) // "/oracle_timing.log"
      adapter%persistent_worker_root = ""
      adapter%persistent_request_path = ""
      adapter%persistent_response_path = ""
      adapter%persistent_ready_path = ""
      adapter%persistent_shutdown_path = ""
      adapter%persistent_done_path = ""
      adapter%persistent_pid_path = ""
      adapter%persistent_log_path = ""
    end if
    if (trim(adapter%backend_mode) == "persistent_worker") then
      call ensure_workdir(trim(adapter%persistent_worker_root))
      call write_text_file(trim(adapter%persistent_worker_root) // "/persistent_init.marker", "persistent_backend_init")
    end if
    call write_text_file(trim(work_root) // "/backend_mode_seen.txt", trim(adapter%backend_mode))
    adapter%call_count = 0
    adapter%initialized = .true.
  end subroutine oracle_adapter_init

  subroutine oracle_adapter_shutdown(adapter)
    type(oracle_adapter_state), intent(inout) :: adapter

    if (.not. adapter%initialized) return
    if (trim(adapter%backend_mode) /= "persistent_worker") return
    if (.not. adapter%persistent_worker_started) return
    call remove_file(trim(adapter%persistent_done_path))
    call write_text_file(trim(adapter%persistent_shutdown_path), "shutdown")
    call wait_for_file(trim(adapter%persistent_done_path), adapter%persistent_worker_start_timeout_s * 1000)
    adapter%persistent_worker_started = .false.
  end subroutine oracle_adapter_shutdown

  pure function oracle_adapter_worker_state(adapter) result(state)
    type(oracle_adapter_state), intent(in) :: adapter
    character(len=16) :: state

    if (trim(adapter%backend_mode) == "persistent_worker") then
      if (adapter%persistent_worker_request_count <= 1) then
        state = "startup"
      else
        state = "reused"
      end if
    else
      state = "external"
    end if
  end function oracle_adapter_worker_state

  subroutine oracle_adapter_eval(adapter, positions_A, charge_in, solver_guess_in, bead_id, charge_policy, result, call_slot)
    type(oracle_adapter_state), intent(inout) :: adapter
    real(dp), intent(in) :: positions_A(:, :)
    real(dp), intent(in) :: charge_in(:)
    real(dp), intent(in) :: solver_guess_in(:)
    integer, intent(in), optional :: bead_id
    character(len=*), intent(in), optional :: charge_policy
    type(oracle_adapter_result), intent(inout) :: result
    integer, intent(in), optional :: call_slot
    character(len=shell_path_len) :: work_dir
    integer :: bead_value, status, natom
    real(dp) :: energy_kcal, max_force, charge_sum
    character(len=512) :: summary_line

    if (.not. adapter%initialized) error stop "oracle_adapter_eval: adapter not initialized"
    if (size(positions_A, 1) /= adapter%template%ndim) error stop "oracle_adapter_eval: position rank mismatch"
    if (size(positions_A, 2) /= adapter%template%nparticles) error stop "oracle_adapter_eval: particle count mismatch"
    if (size(charge_in) /= adapter%template%nparticles) error stop "oracle_adapter_eval: charge size mismatch"
    if (size(solver_guess_in) /= adapter%template%nparticles) error stop "oracle_adapter_eval: solver guess size mismatch"
    if (any(.not. ieee_is_finite(positions_A))) error stop "oracle_adapter_eval: non-finite position"
    if (any(.not. ieee_is_finite(charge_in))) error stop "oracle_adapter_eval: non-finite charge"
    if (any(.not. ieee_is_finite(solver_guess_in))) error stop "oracle_adapter_eval: non-finite solver guess"

    bead_value = 1
    if (present(bead_id)) bead_value = bead_id
    if (bead_value < 1) error stop "oracle_adapter_eval: bead_id must be positive"

    if (present(call_slot)) then
      if (call_slot < 1) error stop "oracle_adapter_eval: call_slot must be positive"
      work_dir = trim(adapter%work_root) // "/adapter_call_" // pad_int(call_slot, 4)
    else
      adapter%call_count = adapter%call_count + 1
      work_dir = trim(adapter%work_root) // "/adapter_call_" // pad_int(adapter%call_count, 4)
    end if

    if (trim(adapter%backend_mode) == "persistent_worker") then
      call persistent_worker_eval(adapter, positions_A, charge_in, solver_guess_in, bead_value, &
          energy_kcal, max_force, charge_sum, natom, status, summary_line, result)
    else
      call oracle_adapter_smoke(adapter, positions_A, charge_in, solver_guess_in, bead_value, &
          work_dir, energy_kcal, max_force, charge_sum, natom, status, summary_line, result)
    end if

    result%force%energy_eV = energy_kcal * kcal2eV
    result%force%status = status
    result%force%diagnostic = trim(summary_line)
    result%diagnostic = trim(summary_line)
    result%legacy_summary = trim(summary_line)
    if (present(charge_policy)) then
      result%diagnostic = trim(result%diagnostic) // " policy=" // trim(charge_policy)
    end if
    if (status /= 0) then
      result%diagnostic = "legacy oracle status=" // trim(adjustl(itoa(status)))
    end if
  end subroutine oracle_adapter_eval

  subroutine oracle_adapter_smoke(adapter, positions_A, charge_in, solver_guess_in, bead_id, work_dir, &
      energy_kcal, max_force, charge_sum, natom, status, summary_line, result)
    type(oracle_adapter_state), intent(in) :: adapter
    real(dp), intent(in) :: positions_A(:, :)
    real(dp), intent(in) :: charge_in(:)
    real(dp), intent(in) :: solver_guess_in(:)
    integer, intent(in) :: bead_id
    character(len=*), intent(in) :: work_dir
    real(dp), intent(out) :: energy_kcal, max_force, charge_sum
    integer, intent(out) :: natom, status
    character(len=*), intent(out) :: summary_line
    type(oracle_adapter_result), intent(inout) :: result
    character(len=shell_path_len) :: wrapper_script
    integer :: exitstat
    logical :: has_geo_path

    call ensure_workdir(trim(work_dir))
    call write_legacy_inputs(adapter, charge_in, solver_guess_in, trim(work_dir))

    has_geo_path = len_trim(adapter%geo_path) > 0
    if (has_geo_path) then
      call write_geo_with_positions(adapter, adapter%geo_path, trim(work_dir) // "/geo", positions_A)
    else
      call generate_phase4_geo(trim(work_dir), adapter%make_geo_script, positions_A)
    end if

    wrapper_script = trim(work_dir) // "/run_oracle_wrapper.sh"
    call write_oracle_wrapper(trim(work_dir), trim(wrapper_script), adapter%timeout_runner, &
        adapter%timeout_s, adapter%smoke_binary, adapter%repeat_binary, adapter%repeatability_check)
    call run_shell("sh " // shell_quote(trim(wrapper_script)), exitstat)
    if (exitstat /= 0) then
      error stop "oracle_adapter_smoke: legacy smoke command failed"
    end if

    call parse_single_point_summary(trim(work_dir) // "/phase4_single_point_oracle.out", &
        energy_kcal, max_force, charge_sum, status, natom)
    call parse_single_point_forces(trim(work_dir) // "/single_point_forces.out", result)

    write(summary_line,'(a,1x,f20.12,1x,es12.4,1x,es12.4,1x,i0,1x,i0,1x,i0)') &
        "oracle_adapter_eval", energy_kcal, max_force, charge_sum, status, natom, bead_id
  end subroutine oracle_adapter_smoke

  subroutine ensure_persistent_worker(adapter)
    type(oracle_adapter_state), intent(inout) :: adapter
    integer :: timeout_ms

    if (trim(adapter%backend_mode) /= "persistent_worker") return
    if (adapter%persistent_worker_started) return

    call ensure_workdir(trim(adapter%persistent_worker_root))
    call write_text_file(trim(adapter%persistent_worker_root) // "/ensure_called.txt", "ensure_persistent_worker")
    call start_persistent_worker(adapter)
    timeout_ms = adapter%persistent_worker_start_timeout_s * 1000
    call wait_for_file(trim(adapter%persistent_ready_path), timeout_ms)
    adapter%persistent_worker_started = .true.
  end subroutine ensure_persistent_worker

  subroutine start_persistent_worker(adapter)
    type(oracle_adapter_state), intent(inout) :: adapter
    character(len=shell_path_len) :: launch_script
    integer :: unit, exitstat, ios

    launch_script = trim(adapter%persistent_worker_root) // "/persistent_launch.sh"
    open(newunit=unit, file=trim(launch_script), status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "start_persistent_worker: unable to write launch script"
    write(unit,'(a)') "#!/usr/bin/env sh"
    write(unit,'(a)') "set -eu"
    write(unit,'(a)') "cd '" // trim(adapter%persistent_worker_root) // "'"
    write(unit,'(a)') "nohup '" // trim(adapter%persistent_worker_binary) // "' '" // trim(adapter%control_path) // "' '" // &
        trim(adapter%ffield_path) // "' '" // trim(adapter%iopt_path) // "' '" // trim(adapter%geo_path) // "' '" // &
        trim(adapter%timing_log_path) // "' '" // trim(adapter%persistent_worker_root) // "' '" // &
        trim(adapter%smoke_binary) // "' '" // trim(adapter%make_geo_script) // "' > '" // &
        trim(adapter%persistent_log_path) // "' 2>&1 &"
    write(unit,'(a)') "echo $! > '" // trim(adapter%persistent_pid_path) // "'"
    close(unit)
    call execute_command_line("sh '" // trim(launch_script) // "'", wait=.true., exitstat=exitstat)
    if (exitstat /= 0) error stop "start_persistent_worker: unable to launch worker"
  end subroutine start_persistent_worker

  subroutine write_persistent_request(path, request_id, bead_id, stage_label, positions_A, charge_in, &
      solver_guess_in)
    character(len=*), intent(in) :: path, stage_label
    integer, intent(in) :: request_id, bead_id
    real(dp), intent(in) :: positions_A(:, :)
    real(dp), intent(in) :: charge_in(:), solver_guess_in(:)
    integer :: unit, i

    if (size(positions_A, 2) /= size(charge_in)) error stop "write_persistent_request: size mismatch"
    open(newunit=unit, file=trim(path), status="replace", action="write")
    write(unit,'(i0)') request_id
    write(unit,'(i0,1x,i0,1x,i0)') 0, bead_id, size(charge_in)
    write(unit,'(a)') trim(stage_label)
    do i = 1, size(charge_in)
      write(unit,'(i8,5(1x,es24.15))') i, positions_A(1, i), positions_A(2, i), positions_A(3, i), charge_in(i), &
          solver_guess_in(i)
    end do
    close(unit)
  end subroutine write_persistent_request

  subroutine wait_for_file(path, timeout_ms)
    character(len=*), intent(in) :: path
    integer, intent(in) :: timeout_ms
    integer :: elapsed_ms

    elapsed_ms = 0
    do
      if (file_exists(trim(path))) return
      if (elapsed_ms >= timeout_ms) error stop "wait_for_file: timed out"
      call wait_milliseconds(25)
      elapsed_ms = elapsed_ms + 25
    end do
  end subroutine wait_for_file

  logical function file_exists(path) result(exists)
    character(len=*), intent(in) :: path
    inquire(file=trim(path), exist=exists)
  end function file_exists

  subroutine read_integer_file(path, value, ios_out)
    character(len=*), intent(in) :: path
    integer, intent(out) :: value
    integer, intent(out) :: ios_out
    integer :: unit

    value = -1
    open(newunit=unit, file=trim(path), status="old", action="read", iostat=ios_out)
    if (ios_out /= 0) return
    read(unit, *, iostat=ios_out) value
    close(unit)
  end subroutine read_integer_file

  subroutine remove_file(path)
    character(len=*), intent(in) :: path
    integer :: exitstat

    call execute_command_line("rm -f " // shell_quote(trim(path)), wait=.true., exitstat=exitstat)
  end subroutine remove_file

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

  subroutine persistent_worker_eval(adapter, positions_A, charge_in, solver_guess_in, bead_id, &
      energy_kcal, max_force, charge_sum, natom, status, summary_line, result)
    type(oracle_adapter_state), intent(inout) :: adapter
    real(dp), intent(in) :: positions_A(:, :)
    real(dp), intent(in) :: charge_in(:)
    real(dp), intent(in) :: solver_guess_in(:)
    integer, intent(in) :: bead_id
    real(dp), intent(out) :: energy_kcal, max_force, charge_sum
    integer, intent(out) :: natom, status
    character(len=*), intent(out) :: summary_line
    type(oracle_adapter_result), intent(inout) :: result
    integer :: request_id, response_id, response_ios, start_count, end_count, clock_rate
    character(len=64) :: stage_label, worker_state
    character(len=512) :: request_path, response_path, persistent_work_dir, summary_path, forces_path
    character(len=512) :: completion_path

    call ensure_persistent_worker(adapter)

    request_id = adapter%persistent_worker_request_count + 1
    adapter%persistent_worker_request_count = request_id
    stage_label = "persistent_worker"
    request_path = trim(adapter%persistent_request_path)
    response_path = trim(adapter%persistent_response_path)
    persistent_work_dir = trim(adapter%persistent_worker_root) // "/adapter_call_" // trim(pad_int(request_id, 4))
    summary_path = trim(persistent_work_dir) // "/phase4_single_point_oracle.out"
    forces_path = trim(persistent_work_dir) // "/single_point_forces.out"
    completion_path = trim(persistent_work_dir) // "/persistent_request.done"
    call remove_file(request_path)
    call remove_file(response_path)
    call remove_file(summary_path)
    call remove_file(forces_path)
    call remove_file(completion_path)
    call write_persistent_request(request_path, request_id, bead_id, stage_label, positions_A, charge_in, &
        solver_guess_in)

    call system_clock(start_count, clock_rate)
    if (clock_rate <= 0) error stop "persistent_worker_eval: system clock unavailable"
    call wait_for_file(trim(response_path), adapter%persistent_worker_start_timeout_s * 1000)
    call read_integer_file(trim(response_path), response_id, response_ios)
    if (response_ios /= 0 .or. response_id /= request_id) error stop "persistent_worker_eval: response id mismatch"
    call wait_for_file(trim(summary_path), adapter%persistent_worker_start_timeout_s * 1000)
    call wait_for_file(trim(forces_path), adapter%persistent_worker_start_timeout_s * 1000)
    call system_clock(end_count)

    call parse_single_point_summary(trim(summary_path), &
        energy_kcal, max_force, charge_sum, status, natom)
    call parse_single_point_forces(trim(forces_path), result)

    result%force%energy_eV = energy_kcal * kcal2eV
    result%force%status = status
    result%force%diagnostic = "persistent worker"
    result%diagnostic = "persistent worker"
    result%legacy_summary = "persistent worker"
    if (status /= 0) then
      result%diagnostic = "legacy oracle status=" // trim(adjustl(itoa(status)))
    end if
    if (adapter%persistent_worker_request_count == 1) then
      worker_state = "startup"
    else
      worker_state = "reused"
    end if
    summary_line = "persistent_worker " // trim(worker_state) // " request=" // trim(itoa(request_id))
  end subroutine persistent_worker_eval

  subroutine parse_single_point_summary(path, energy_kcal, max_force, charge_sum, status, natom)
    character(len=*), intent(in) :: path
    real(dp), intent(out) :: energy_kcal, max_force, charge_sum
    integer, intent(out) :: status, natom
    integer :: unit, ios
    character(len=512) :: line

    open(newunit=unit, file=path, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "parse_single_point_summary: unable to open summary"
    read(unit,'(A)',iostat=ios) line
    if (ios /= 0) error stop "parse_single_point_summary: missing header"
    read(unit,*,iostat=ios) energy_kcal, max_force, charge_sum, status, natom
    if (ios /= 0) error stop "parse_single_point_summary: unable to parse data line"
    close(unit)
  end subroutine parse_single_point_summary

  subroutine parse_single_point_forces(path, result)
    character(len=*), intent(in) :: path
    type(oracle_adapter_result), intent(inout) :: result
    integer :: unit, ios, idx, natom, i
    character(len=512) :: line
    real(dp) :: fx, fy, fz, charge, solver

    open(newunit=unit, file=path, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "parse_single_point_forces: unable to open forces file"
    read(unit,'(A)',iostat=ios) line
    if (ios /= 0) error stop "parse_single_point_forces: missing header"
    natom = 0
    do
      read(unit,*,iostat=ios) idx, fx, fy, fz, charge, solver
      if (ios /= 0) exit
      natom = natom + 1
    end do
    rewind(unit)
    read(unit,'(A)',iostat=ios) line
    if (ios /= 0) error stop "parse_single_point_forces: header rewind failed"
    call allocate_force_result(result%force, 3, natom, .true.)
    if (allocated(result%solver_guess)) deallocate(result%solver_guess)
    allocate(result%solver_guess(natom))
    result%solver_guess = 0.0_dp
    do i = 1, natom
      read(unit,*,iostat=ios) idx, fx, fy, fz, charge, solver
      if (ios /= 0) error stop "parse_single_point_forces: malformed force line"
      if (idx /= i) error stop "parse_single_point_forces: atom index mismatch"
      result%force%forces_eV_per_A(:, i) = [fx, fy, fz]
      result%force%charges(i) = charge
      result%solver_guess(i) = solver
    end do
    close(unit)
  end subroutine parse_single_point_forces

  subroutine write_legacy_inputs(adapter, charge_in, solver_guess_in, work_dir)
    type(oracle_adapter_state), intent(in) :: adapter
    real(dp), intent(in) :: charge_in(:)
    real(dp), intent(in) :: solver_guess_in(:)
    character(len=*), intent(in) :: work_dir

    call write_control_with_cell(adapter, trim(work_dir) // "/control")
    call copy_fixture(adapter%ffield_path, trim(work_dir) // "/ffield")
    call copy_fixture(adapter%iopt_path, trim(work_dir) // "/fort.20")
    call write_text_file(trim(work_dir) // "/fort.35", "0")
    call write_charge_trace(trim(work_dir) // "/charge_memory.in", charge_in, solver_guess_in)
  end subroutine write_legacy_inputs

  subroutine write_control_with_cell(adapter, target)
    type(oracle_adapter_state), intent(in) :: adapter
    character(len=*), intent(in) :: target
    integer, parameter :: line_len = 8192
    integer :: src_unit, dst_unit, ios
    character(len=line_len) :: line
    logical :: replaced

    open(newunit=src_unit, file=adapter%control_path, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "write_control_with_cell: unable to open source control"
    open(newunit=dst_unit, file=target, status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "write_control_with_cell: unable to open target control"

    do
      read(src_unit,'(A)',iostat=ios) line
      if (ios /= 0) exit
      replaced = .false.
      if (adapter%template%cell%defined) then
        call maybe_write_control_axis(dst_unit, line, "axis1", adapter%template%cell%crystx(1), replaced)
        if (.not. replaced) call maybe_write_control_axis(dst_unit, line, "axis2", adapter%template%cell%crystx(2), replaced)
        if (.not. replaced) call maybe_write_control_axis(dst_unit, line, "axis3", adapter%template%cell%crystx(3), replaced)
      end if
      if (.not. replaced) write(dst_unit,'(a)') trim(line)
    end do

    close(src_unit)
    close(dst_unit)
  end subroutine write_control_with_cell

  subroutine maybe_write_control_axis(unit, line, axis_label, axis_value, replaced)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: line, axis_label
    real(dp), intent(in) :: axis_value
    logical, intent(inout) :: replaced
    character(len=64) :: label
    real(dp) :: old_value
    integer :: ios

    if (replaced) return
    label = ""
    old_value = 0.0_dp
    read(line, *, iostat=ios) old_value, label
    if (ios == 0 .and. trim(label) == axis_label) then
      write(unit,'(f12.5,1x,a)') axis_value, trim(axis_label)
      replaced = .true.
    end if
  end subroutine maybe_write_control_axis

  subroutine write_oracle_wrapper(work_dir, script_path, timeout_runner, timeout_s, smoke_binary, repeat_binary, &
      repeatability_check)
    character(len=*), intent(in) :: work_dir, script_path, timeout_runner, smoke_binary, repeat_binary
    integer, intent(in) :: timeout_s
    logical, intent(in) :: repeatability_check
    integer :: unit, ios

    open(newunit=unit, file=script_path, status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "write_oracle_wrapper: unable to open wrapper script"
    write(unit,'(a)') "#!/usr/bin/env sh"
    write(unit,'(a)') "set -eu"
    write(unit,'(a)') "cd '" // trim(work_dir) // "'"
    write(unit,'(a)') "python3 '" // trim(timeout_runner) // "' " // trim(itoa(timeout_s)) // " '" // &
        "run_single_point.log' '" // trim(smoke_binary) // "'"
    if (repeatability_check) then
      write(unit,'(a)') "python3 '" // trim(timeout_runner) // "' " // trim(itoa(timeout_s)) // " '" // &
          "run_repeatability.log' '" // trim(repeat_binary) // "' 3"
    else
      write(unit,'(a)') "printf '%s\n' 'repeatability_check skipped' > run_repeatability.log"
    end if
    close(unit)
  end subroutine write_oracle_wrapper

  subroutine generate_phase4_geo(work_dir, make_geo_script, positions_A)
    character(len=*), intent(in) :: work_dir, make_geo_script
    real(dp), intent(in) :: positions_A(:, :)
    character(len=64) :: el_x, el_y, el_z
    character(len=shell_path_len) :: command
    integer :: exitstat

    if (size(positions_A, 1) /= 3) error stop "generate_phase4_geo: positions must be 3D"
    write(el_x,'(es24.15)') positions_A(1, 1)
    write(el_y,'(es24.15)') positions_A(2, 1)
    write(el_z,'(es24.15)') positions_A(3, 1)
    command = "cd " // shell_quote(work_dir) // " && PHASE4_EL_X=" // trim(adjustl(el_x)) // &
        " PHASE4_EL_Y=" // trim(adjustl(el_y)) // " PHASE4_EL_Z=" // trim(adjustl(el_z)) // &
        " python3 " // shell_quote(make_geo_script)
    call run_shell(command, exitstat)
    if (exitstat /= 0) error stop "generate_phase4_geo: geometry generation failed"
  end subroutine generate_phase4_geo

  subroutine write_geo_with_positions(adapter, source, target, positions_A)
    type(oracle_adapter_state), intent(in) :: adapter
    character(len=*), intent(in) :: source, target
    real(dp), intent(in) :: positions_A(:, :)
    integer, parameter :: line_len = 8192
    integer :: src_unit, dst_unit, ios, atom_counter
    character(len=line_len) :: line
    character(len=8) :: record, res3, res5
    character(len=5) :: atom_name, atom_type
    character(len=1) :: res1
    integer :: atom_id_read, dummy_int1, dummy_int2
    real(dp) :: x, y, z, dummy_real
    logical :: wrote_crystx

    if (size(positions_A, 1) /= 3) error stop "write_geo_with_positions: positions must be 3D"
    open(newunit=src_unit, file=source, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "write_geo_with_positions: unable to open source geo"
    open(newunit=dst_unit, file=target, status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "write_geo_with_positions: unable to open target geo"

    atom_counter = 0
    wrote_crystx = .false.
    do
      read(src_unit,'(A)',iostat=ios) line
      if (ios /= 0) exit
      if (starts_with(trim(line), "CRYSTX")) then
        if (adapter%template%cell%defined) then
          call write_crystx_line(dst_unit, adapter%template%cell%crystx)
        else
          write(dst_unit,'(a)') trim(line)
        end if
        wrote_crystx = .true.
      else if (starts_with(trim(line), "HETATM")) then
        if (adapter%template%cell%defined .and. .not. wrote_crystx) then
          call write_crystx_line(dst_unit, adapter%template%cell%crystx)
          wrote_crystx = .true.
        end if
        record = ""
        atom_name = ""
        atom_type = ""
        res3 = ""
        res1 = ""
        res5 = ""
        x = 0.0_dp
        y = 0.0_dp
        z = 0.0_dp
        dummy_int1 = 0
        dummy_int2 = 0
        dummy_real = 0.0_dp
        read(line, '(a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)', iostat=ios) &
            record, atom_id_read, atom_name, res3, res1, atom_type, x, y, z, res5, dummy_int1, dummy_int2, dummy_real
        if (ios /= 0) error stop "write_geo_with_positions: malformed HETATM record"
        atom_counter = atom_counter + 1
        if (atom_counter > size(positions_A, 2)) error stop "write_geo_with_positions: too many HETATM records"
        write(dst_unit,'(a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)') &
            record, atom_id_read, atom_name, res3, res1, atom_type, positions_A(:, atom_counter), &
            res5, dummy_int1, dummy_int2, dummy_real
      else
        write(dst_unit,'(a)') trim(line)
      end if
    end do
    close(src_unit)
    close(dst_unit)
    if (atom_counter /= size(positions_A, 2)) error stop "write_geo_with_positions: HETATM count mismatch"
  end subroutine write_geo_with_positions

  subroutine write_crystx_line(unit, crystx)
    integer, intent(in) :: unit
    real(dp), intent(in) :: crystx(6)

    write(unit,'(a6,6(1x,f12.5))') "CRYSTX", crystx
  end subroutine write_crystx_line

  subroutine write_charge_trace(path, charge_in, solver_guess_in)
    character(len=*), intent(in) :: path
    real(dp), intent(in) :: charge_in(:)
    real(dp), intent(in) :: solver_guess_in(:)
    integer :: unit, i

    open(newunit=unit, file=path, status="replace", action="write")
    do i = 1, size(charge_in)
      write(unit,'(i0,1x,es24.15,1x,es24.15)') i, charge_in(i), solver_guess_in(i)
    end do
    close(unit)
  end subroutine write_charge_trace

  subroutine ensure_workdir(path)
    character(len=*), intent(in) :: path
    character(len=shell_path_len) :: cmd
    integer :: exitstat

    cmd = "mkdir -p " // shell_quote(path)
    call run_shell(cmd, exitstat)
    if (exitstat /= 0) error stop "ensure_workdir: unable to create working directory"
  end subroutine ensure_workdir

  subroutine copy_fixture(source, target)
    character(len=*), intent(in) :: source, target
    integer, parameter :: line_len = 8192
    integer :: src_unit, dst_unit, ios
    character(len=line_len) :: line

    open(newunit=src_unit, file=source, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "copy_fixture: unable to open source file"
    open(newunit=dst_unit, file=target, status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "copy_fixture: unable to open target file"
    do
      read(src_unit,'(A)',iostat=ios) line
      if (ios /= 0) exit
      write(dst_unit,'(a)') trim(line)
    end do
    close(src_unit)
    close(dst_unit)
  end subroutine copy_fixture

  pure logical function starts_with(text, prefix) result(match)
    character(len=*), intent(in) :: text, prefix
    integer :: n

    n = len_trim(prefix)
    if (len_trim(text) < n) then
      match = .false.
    else
      match = text(1:n) == prefix(1:n)
    end if
  end function starts_with

  subroutine write_text_file(path, text)
    character(len=*), intent(in) :: path, text
    integer :: unit

    open(newunit=unit, file=path, status="replace", action="write")
    write(unit,'(a)') trim(text)
    close(unit)
  end subroutine write_text_file

  subroutine run_shell(command, exitstat)
    character(len=*), intent(in) :: command
    integer, intent(out) :: exitstat
    integer :: cmdstat
    character(len=256) :: errmsg

    call execute_command_line(trim(command), wait=.true., exitstat=exitstat, cmdstat=cmdstat, cmdmsg=errmsg)
    if (cmdstat /= 0) error stop "run_shell: " // trim(errmsg)
  end subroutine run_shell

  pure function shell_quote(text) result(quoted)
    character(len=*), intent(in) :: text
    character(len=:), allocatable :: quoted
    integer :: n

    n = len_trim(text)
    allocate(character(len=n + 2) :: quoted)
    quoted = "'" // trim(text) // "'"
  end function shell_quote

  pure function pad_int(value, width) result(text)
    integer, intent(in) :: value, width
    character(len=32) :: text
    integer :: nzeros

    write(text,'(i0)') value
    nzeros = max(0, width - len_trim(text))
    if (nzeros > 0) text = repeat("0", nzeros) // trim(text)
  end function pad_int

  pure function itoa(value) result(text)
    integer, intent(in) :: value
    character(len=32) :: text

    write(text,'(i0)') value
  end function itoa

  subroutine assert_nonempty(text, label)
    character(len=*), intent(in) :: text, label

    if (len_trim(text) == 0) error stop trim(label) // " must not be empty"
  end subroutine assert_nonempty

end module ereaxff_oracle_adapter
