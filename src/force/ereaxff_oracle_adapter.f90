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
    logical :: repeatability_check = .true.
  end type oracle_adapter_state

  type, public :: oracle_adapter_result
    type(force_result) :: force
    real(dp), allocatable :: solver_guess(:)
    character(len=512) :: legacy_summary = ""
    character(len=512) :: diagnostic = ""
  end type oracle_adapter_result

  public :: oracle_adapter_init
  public :: oracle_adapter_eval
  public :: oracle_adapter_smoke

contains

  subroutine oracle_adapter_init(adapter, template, v1_root, v2_root, work_root, control_path, ffield_path, iopt_path, &
      smoke_binary, repeat_binary, make_geo_script, timeout_runner, geo_path, repeatability_check)
    type(oracle_adapter_state), intent(inout) :: adapter
    type(system_state), intent(in) :: template
    character(len=*), intent(in) :: v1_root, v2_root, work_root
    character(len=*), intent(in), optional :: control_path, ffield_path, iopt_path
    character(len=*), intent(in), optional :: smoke_binary, repeat_binary, make_geo_script, timeout_runner, geo_path
    logical, intent(in), optional :: repeatability_check

    call assert_nonempty(v1_root, "oracle_adapter_init: v1_root")
    call assert_nonempty(v2_root, "oracle_adapter_init: v2_root")
    call assert_nonempty(work_root, "oracle_adapter_init: work_root")
    if (template%nparticles < 1) error stop "oracle_adapter_init: template state is empty"

    adapter%template = template
    adapter%v1_root = trim(v1_root)
    adapter%v2_root = trim(v2_root)
    adapter%work_root = trim(work_root)
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
    adapter%call_count = 0
    adapter%initialized = .true.
  end subroutine oracle_adapter_init

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

    call oracle_adapter_smoke(adapter, positions_A, charge_in, solver_guess_in, bead_value, &
        work_dir, energy_kcal, max_force, charge_sum, natom, status, summary_line, result)

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
        adapter%smoke_binary, adapter%repeat_binary, adapter%repeatability_check)
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

  subroutine write_oracle_wrapper(work_dir, script_path, timeout_runner, smoke_binary, repeat_binary, repeatability_check)
    character(len=*), intent(in) :: work_dir, script_path, timeout_runner, smoke_binary, repeat_binary
    logical, intent(in) :: repeatability_check
    integer :: unit, ios

    open(newunit=unit, file=script_path, status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "write_oracle_wrapper: unable to open wrapper script"
    write(unit,'(a)') "#!/usr/bin/env sh"
    write(unit,'(a)') "set -eu"
    write(unit,'(a)') "cd '" // trim(work_dir) // "'"
    write(unit,'(a)') "python3 '" // trim(timeout_runner) // "' 90 '" // &
        "run_single_point.log' '" // trim(smoke_binary) // "'"
    if (repeatability_check) then
      write(unit,'(a)') "python3 '" // trim(timeout_runner) // "' 90 '" // &
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
