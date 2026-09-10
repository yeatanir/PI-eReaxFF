module ereaxff_fast_jax_adapter
  use ereaxff_kinds, only: dp
  use ereaxff_system_state, only: system_state
  implicit none
  private

  integer, parameter :: shell_path_len = 16384

  type, public :: fast_jax_adapter_state
    logical :: initialized = .false.
    logical :: worker_started = .false.
    integer :: request_count = 0
    type(system_state) :: template
    character(len=:), allocatable :: work_root
    character(len=:), allocatable :: command
    character(len=:), allocatable :: ffield_path
    character(len=:), allocatable :: geo_path
    character(len=:), allocatable :: kernel
    real(dp), allocatable :: fixed_charges(:)
  end type fast_jax_adapter_state

  type, public :: fast_jax_batch_result
    real(dp), allocatable :: bead_energies_eV(:)
    real(dp), allocatable :: forces_eV_per_A(:, :, :)
    real(dp), allocatable :: charges(:, :)
    real(dp), allocatable :: solver_guess(:, :)
    integer, allocatable :: bead_status(:)
  end type fast_jax_batch_result

  public :: fast_jax_adapter_init
  public :: fast_jax_adapter_eval_batch
  public :: fast_jax_adapter_shutdown

contains

  subroutine fast_jax_adapter_init(adapter, template, work_root, command, ffield_path, geo_path, kernel, allow_unvalidated)
    type(fast_jax_adapter_state), intent(inout) :: adapter
    type(system_state), intent(in) :: template
    character(len=*), intent(in) :: work_root, command, ffield_path, geo_path, kernel
    logical, intent(in) :: allow_unvalidated

    if (.not. allow_unvalidated) then
      error stop "fast_jax_adapter_init: fast_jax_oracle requires fast_backend_allow_unvalidated=true"
    end if
    if (template%nparticles < 1) error stop "fast_jax_adapter_init: template state is empty"
    if (.not. template%cell%defined) error stop "fast_jax_adapter_init: fast backend requires a defined cell"
    if (abs(template%cell%crystx(4) - 90.0_dp) > 1.0e-8_dp .or. &
        abs(template%cell%crystx(5) - 90.0_dp) > 1.0e-8_dp .or. &
        abs(template%cell%crystx(6) - 90.0_dp) > 1.0e-8_dp) then
      error stop "fast_jax_adapter_init: fast backend currently supports only orthorhombic cells"
    end if
    call assert_electron_only_jamil_system(template)
    if (len_trim(work_root) == 0) error stop "fast_jax_adapter_init: empty work_root"
    if (len_trim(ffield_path) == 0) error stop "fast_jax_adapter_init: empty ffield_path"
    if (len_trim(geo_path) == 0) error stop "fast_jax_adapter_init: empty geo_path"

    adapter%template = template
    adapter%work_root = trim(work_root)
    if (len_trim(command) > 0) then
      adapter%command = trim(command)
    else
      adapter%command = ".venv-fast/bin/python projects/jamil_replication/fast_backend/jax_live_worker.py"
    end if
    adapter%ffield_path = trim(ffield_path)
    adapter%geo_path = trim(geo_path)
    if (len_trim(kernel) > 0) then
      adapter%kernel = trim(kernel)
    else
      adapter%kernel = "reaxff_corrected_bond_nonbonded_valence"
    end if
    call ensure_workdir(trim(adapter%work_root))
    call load_molcharge_constraints(trim(adapter%geo_path), template%nparticles, adapter%fixed_charges)
    adapter%request_count = 0
    adapter%worker_started = .false.
    adapter%initialized = .true.
    call ensure_worker(adapter)
  end subroutine fast_jax_adapter_init

  subroutine fast_jax_adapter_eval_batch(adapter, positions_A, charge_in_by_bead, solver_guess_by_bead, result)
    type(fast_jax_adapter_state), intent(inout) :: adapter
    real(dp), intent(in) :: positions_A(:, :, :)
    real(dp), intent(in) :: charge_in_by_bead(:, :)
    real(dp), intent(in) :: solver_guess_by_bead(:, :)
    type(fast_jax_batch_result), intent(inout) :: result
    integer :: request_id, natoms, nbeads

    if (.not. adapter%initialized) error stop "fast_jax_adapter_eval_batch: adapter not initialized"
    natoms = size(positions_A, 2)
    nbeads = size(positions_A, 3)
    if (size(positions_A, 1) /= 3) error stop "fast_jax_adapter_eval_batch: only 3D systems are supported"
    if (natoms /= adapter%template%nparticles) error stop "fast_jax_adapter_eval_batch: atom count mismatch"
    if (size(charge_in_by_bead, 1) /= natoms .or. size(charge_in_by_bead, 2) /= nbeads) then
      error stop "fast_jax_adapter_eval_batch: charge shape mismatch"
    end if
    if (size(solver_guess_by_bead, 1) /= natoms .or. size(solver_guess_by_bead, 2) /= nbeads) then
      error stop "fast_jax_adapter_eval_batch: solver shape mismatch"
    end if
    call ensure_worker(adapter)
    adapter%request_count = adapter%request_count + 1
    request_id = adapter%request_count
    call cleanup_request_files(adapter)
    call write_batch_request(adapter, request_id, positions_A, charge_in_by_bead, solver_guess_by_bead)
    call write_text_file(trim(adapter%work_root) // "/request.ready", "ready")
    call wait_for_file(trim(adapter%work_root) // "/response.ready", 120000)
    call read_batch_response(trim(adapter%work_root) // "/response.out", request_id, natoms, nbeads, result)
  end subroutine fast_jax_adapter_eval_batch

  subroutine fast_jax_adapter_shutdown(adapter)
    type(fast_jax_adapter_state), intent(inout) :: adapter

    if (.not. adapter%initialized) return
    if (.not. adapter%worker_started) return
    call remove_file(trim(adapter%work_root) // "/worker.done")
    call write_text_file(trim(adapter%work_root) // "/shutdown.flag", "shutdown")
    call wait_for_file(trim(adapter%work_root) // "/worker.done", 60000)
    adapter%worker_started = .false.
  end subroutine fast_jax_adapter_shutdown

  subroutine ensure_worker(adapter)
    type(fast_jax_adapter_state), intent(inout) :: adapter

    if (adapter%worker_started) return
    call cleanup_worker_files(adapter)
    call start_worker(adapter)
    call wait_for_file(trim(adapter%work_root) // "/worker.ready", 120000)
    adapter%worker_started = .true.
  end subroutine ensure_worker

  subroutine start_worker(adapter)
    type(fast_jax_adapter_state), intent(in) :: adapter
    character(len=shell_path_len) :: launch_script
    integer :: unit, ios, exitstat

    launch_script = trim(adapter%work_root) // "/fast_jax_launch.sh"
    open(newunit=unit, file=trim(launch_script), status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "start_worker: unable to write launch script"
    write(unit,'(a)') "#!/usr/bin/env sh"
    write(unit,'(a)') "set -eu"
    write(unit,'(a)') "nohup " // trim(adapter%command) // &
        " --worker-root " // shell_quote(trim(adapter%work_root)) // &
        " --ffield " // shell_quote(trim(adapter%ffield_path)) // &
        " --kernel " // shell_quote(trim(adapter%kernel)) // &
        " > " // shell_quote(trim(adapter%work_root) // "/fast_jax_worker.log") // " 2>&1 &"
    write(unit,'(a)') "echo $! > " // shell_quote(trim(adapter%work_root) // "/fast_jax_worker.pid")
    close(unit)
    call execute_command_line("sh " // shell_quote(trim(launch_script)), wait=.true., exitstat=exitstat)
    if (exitstat /= 0) error stop "start_worker: unable to launch fast JAX worker"
  end subroutine start_worker

  subroutine write_batch_request(adapter, request_id, positions_A, charge_in_by_bead, solver_guess_by_bead)
    type(fast_jax_adapter_state), intent(in) :: adapter
    integer, intent(in) :: request_id
    real(dp), intent(in) :: positions_A(:, :, :)
    real(dp), intent(in) :: charge_in_by_bead(:, :)
    real(dp), intent(in) :: solver_guess_by_bead(:, :)
    character(len=1024) :: tmp_path, final_path
    integer :: unit, ios, atom, bead

    final_path = trim(adapter%work_root) // "/request.in"
    tmp_path = trim(final_path) // ".tmp"
    open(newunit=unit, file=trim(tmp_path), status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "write_batch_request: unable to write request"
    write(unit,'(a,1x,i0)') "request_id", request_id
    write(unit,'(a,1x,i0)') "natoms", size(positions_A, 2)
    write(unit,'(a,1x,i0)') "nbeads", size(positions_A, 3)
    write(unit,'(a,6(1x,es24.15))') "crystx", adapter%template%cell%crystx
    write(unit,'(a)') "atoms"
    do atom = 1, size(positions_A, 2)
      write(unit,'(i0,1x,a,1x,es24.15)') atom, trim(adapter%template%element(atom)), adapter%fixed_charges(atom)
    end do
    do bead = 1, size(positions_A, 3)
      write(unit,'(a,1x,i0)') "bead", bead
      do atom = 1, size(positions_A, 2)
        write(unit,'(i0,5(1x,es24.15))') atom, positions_A(1, atom, bead), positions_A(2, atom, bead), &
            positions_A(3, atom, bead), charge_in_by_bead(atom, bead), solver_guess_by_bead(atom, bead)
      end do
    end do
    close(unit)
    call publish_file_atomic(trim(tmp_path), trim(final_path))
  end subroutine write_batch_request

  subroutine read_batch_response(path, expected_request_id, natoms, nbeads, result)
    character(len=*), intent(in) :: path
    integer, intent(in) :: expected_request_id, natoms, nbeads
    type(fast_jax_batch_result), intent(inout) :: result
    integer :: unit, ios, request_id, status, response_nbeads, response_natoms, bead, bead_index, bead_status, atom, atom_index
    real(dp) :: energy
    character(len=128) :: label
    character(len=512) :: diagnostic

    if (allocated(result%bead_energies_eV)) deallocate(result%bead_energies_eV)
    if (allocated(result%forces_eV_per_A)) deallocate(result%forces_eV_per_A)
    if (allocated(result%charges)) deallocate(result%charges)
    if (allocated(result%solver_guess)) deallocate(result%solver_guess)
    if (allocated(result%bead_status)) deallocate(result%bead_status)
    allocate(result%bead_energies_eV(nbeads))
    allocate(result%forces_eV_per_A(3, natoms, nbeads))
    allocate(result%charges(natoms, nbeads))
    allocate(result%solver_guess(natoms, nbeads))
    allocate(result%bead_status(nbeads))
    result%bead_energies_eV = 0.0_dp
    result%forces_eV_per_A = 0.0_dp
    result%charges = 0.0_dp
    result%solver_guess = 0.0_dp
    result%bead_status = 0

    open(newunit=unit, file=trim(path), status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "read_batch_response: unable to open response"
    read(unit, *, iostat=ios) label, request_id
    if (ios /= 0 .or. trim(label) /= "request_id" .or. request_id /= expected_request_id) then
      error stop "read_batch_response: request id mismatch"
    end if
    read(unit, *, iostat=ios) label, status
    if (ios /= 0 .or. trim(label) /= "status") error stop "read_batch_response: missing status"
    read(unit, '(a)', iostat=ios) diagnostic
    if (status /= 0) error stop "read_batch_response: fast JAX worker failed; see fast_jax_work/last_error.txt"
    read(unit, *, iostat=ios) label, response_nbeads
    if (ios /= 0 .or. trim(label) /= "nbeads" .or. response_nbeads /= nbeads) then
      error stop "read_batch_response: bead count mismatch"
    end if
    read(unit, *, iostat=ios) label, response_natoms
    if (ios /= 0 .or. trim(label) /= "natoms" .or. response_natoms /= natoms) then
      error stop "read_batch_response: atom count mismatch"
    end if
    do bead = 1, nbeads
      read(unit, *, iostat=ios) label, bead_index, bead_status, energy
      if (ios /= 0 .or. trim(label) /= "BEAD" .or. bead_index /= bead) error stop "read_batch_response: bead header mismatch"
      result%bead_status(bead) = bead_status
      result%bead_energies_eV(bead) = energy
      read(unit, *, iostat=ios) label
      if (ios /= 0 .or. trim(label) /= "CHARGES") error stop "read_batch_response: missing charges"
      do atom = 1, natoms
        read(unit, *, iostat=ios) atom_index, result%charges(atom, bead)
        if (ios /= 0 .or. atom_index /= atom) error stop "read_batch_response: charge atom mismatch"
      end do
      read(unit, *, iostat=ios) label
      if (ios /= 0 .or. trim(label) /= "SOLVER") error stop "read_batch_response: missing solver"
      do atom = 1, natoms
        read(unit, *, iostat=ios) atom_index, result%solver_guess(atom, bead)
        if (ios /= 0 .or. atom_index /= atom) error stop "read_batch_response: solver atom mismatch"
      end do
      read(unit, *, iostat=ios) label
      if (ios /= 0 .or. trim(label) /= "FORCES") error stop "read_batch_response: missing forces"
      do atom = 1, natoms
        read(unit, *, iostat=ios) atom_index, result%forces_eV_per_A(1, atom, bead), &
            result%forces_eV_per_A(2, atom, bead), result%forces_eV_per_A(3, atom, bead)
        if (ios /= 0 .or. atom_index /= atom) error stop "read_batch_response: force atom mismatch"
      end do
      read(unit, *, iostat=ios) label
      if (ios /= 0 .or. trim(label) /= "END_BEAD") error stop "read_batch_response: missing END_BEAD"
    end do
    close(unit)
  end subroutine read_batch_response

  subroutine load_molcharge_constraints(path, natoms, fixed_charges)
    character(len=*), intent(in) :: path
    integer, intent(in) :: natoms
    real(dp), allocatable, intent(out) :: fixed_charges(:)
    integer :: unit, ios, first_atom, last_atom
    real(dp) :: total_charge
    character(len=8192) :: line
    character(len=64) :: label
    logical :: found

    allocate(fixed_charges(natoms))
    fixed_charges = 0.0_dp
    found = .false.
    open(newunit=unit, file=trim(path), status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "load_molcharge_constraints: unable to open geo"
    do
      read(unit,'(a)',iostat=ios) line
      if (ios /= 0) exit
      if (.not. starts_with(adjustl(line), "MOLCHARGE")) cycle
      read(line, *, iostat=ios) label, first_atom, last_atom, total_charge
      if (ios /= 0) error stop "load_molcharge_constraints: malformed MOLCHARGE line"
      if (first_atom < 1 .or. last_atom > natoms .or. first_atom > last_atom) then
        error stop "load_molcharge_constraints: MOLCHARGE range out of bounds"
      end if
      fixed_charges(first_atom:last_atom) = total_charge / real(last_atom - first_atom + 1, dp)
      found = .true.
    end do
    close(unit)
    if (.not. found) error stop "load_molcharge_constraints: no MOLCHARGE constraints found"
  end subroutine load_molcharge_constraints

  subroutine assert_electron_only_jamil_system(template)
    type(system_state), intent(in) :: template
    integer :: atom

    do atom = 1, template%nparticles
      select case (trim(template%element(atom)))
      case ("C", "El")
        cycle
      case ("Ho")
        error stop "assert_electron_only_jamil_system: Ho is not supported; use electron-only eReaxFF"
      case default
        error stop "assert_electron_only_jamil_system: fast backend currently supports only C and El"
      end select
    end do
  end subroutine assert_electron_only_jamil_system

  subroutine cleanup_worker_files(adapter)
    type(fast_jax_adapter_state), intent(in) :: adapter

    call remove_file(trim(adapter%work_root) // "/worker.ready")
    call remove_file(trim(adapter%work_root) // "/worker.done")
    call remove_file(trim(adapter%work_root) // "/shutdown.flag")
    call cleanup_request_files(adapter)
  end subroutine cleanup_worker_files

  subroutine cleanup_request_files(adapter)
    type(fast_jax_adapter_state), intent(in) :: adapter

    call remove_file(trim(adapter%work_root) // "/request.ready")
    call remove_file(trim(adapter%work_root) // "/response.ready")
    call remove_file(trim(adapter%work_root) // "/request.in")
    call remove_file(trim(adapter%work_root) // "/response.out")
  end subroutine cleanup_request_files

  subroutine ensure_workdir(path)
    character(len=*), intent(in) :: path
    integer :: exitstat

    call execute_command_line("mkdir -p " // shell_quote(trim(path)), wait=.true., exitstat=exitstat)
    if (exitstat /= 0) error stop "ensure_workdir: unable to create fast JAX working directory"
  end subroutine ensure_workdir

  subroutine write_text_file(path, text)
    character(len=*), intent(in) :: path, text
    integer :: unit

    open(newunit=unit, file=trim(path), status="replace", action="write")
    write(unit,'(a)') trim(text)
    close(unit)
  end subroutine write_text_file

  subroutine remove_file(path)
    character(len=*), intent(in) :: path
    integer :: exitstat

    call execute_command_line("rm -f " // shell_quote(trim(path)), wait=.true., exitstat=exitstat)
  end subroutine remove_file

  subroutine publish_file_atomic(tmp_path, final_path)
    character(len=*), intent(in) :: tmp_path, final_path
    integer :: exitstat

    call execute_command_line("mv -f " // shell_quote(trim(tmp_path)) // " " // shell_quote(trim(final_path)), &
        wait=.true., exitstat=exitstat)
    if (exitstat /= 0) error stop "publish_file_atomic: unable to publish fast JAX request"
  end subroutine publish_file_atomic

  subroutine wait_for_file(path, timeout_ms)
    character(len=*), intent(in) :: path
    integer, intent(in) :: timeout_ms
    integer :: elapsed_ms

    elapsed_ms = 0
    do
      if (file_exists(trim(path))) return
      if (elapsed_ms >= timeout_ms) error stop "wait_for_file: timed out waiting for fast JAX worker"
      call wait_milliseconds(25)
      elapsed_ms = elapsed_ms + 25
    end do
  end subroutine wait_for_file

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

  logical function file_exists(path) result(exists)
    character(len=*), intent(in) :: path

    inquire(file=trim(path), exist=exists)
  end function file_exists

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

  pure function shell_quote(text) result(quoted)
    character(len=*), intent(in) :: text
    character(len=:), allocatable :: quoted

    allocate(character(len=len_trim(text) + 2) :: quoted)
    quoted = '"' // trim(text) // '"'
  end function shell_quote

end module ereaxff_fast_jax_adapter
