program persistent_oracle_worker
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 300)
  integer, parameter :: natmax = 2600
  integer, parameter :: nbondmax = 400000
  integer, parameter :: nbucketmax = 256
  real(dp), parameter :: kcal2eV = 0.043364_dp

  integer :: request_id, step_index, bead_id, natom, status, ios, nbond, bo_natom, nbuckets
  integer :: bond_atoms(2, nbondmax)
  integer :: bucket_left(nbucketmax), bucket_right(nbucketmax), bucket_self(nbucketmax), bucket_count(nbucketmax)
  integer :: start_count, end_count, clock_rate
  real(dp) :: x(3, natmax), force(3, natmax), charge_in(natmax), solver_in(natmax)
  real(dp) :: charge_out(natmax), solver_out(natmax), energy_kcal, max_force, charge_sum
  real(dp) :: energy_terms(14), energy_terms_total
  real(dp) :: abo(natmax), bo(nbondmax), bopi(nbondmax), bopi2(nbondmax)
  real(dp) :: bucket_vdw(nbucketmax), bucket_coulomb(nbucketmax)
  character(len=512) :: control_path, ffield_path, iopt_path, geo_path, timing_log_path, worker_root, smoke_binary
  character(len=512) :: make_geo_script
  character(len=128) :: stage_label
  character(len=512) :: request_file, response_file, ready_file, shutdown_file, done_file
  character(len=512) :: work_dir, summary_path, forces_path, completion_path
  logical :: first_request

  interface
    subroutine ereaxff_force_oracle_init()
    end subroutine ereaxff_force_oracle_init
    subroutine ereaxff_force_eval(x_all, natom_eval, charge_in, acks2_solver_guess_in, energy_kcal, forces_work, &
        charge_out, acks2_solver_guess_out, status)
      integer :: natom_eval, status
      double precision :: x_all(3,*), charge_in(*)
      double precision :: acks2_solver_guess_in(*), energy_kcal
      double precision :: forces_work(3,*), charge_out(*)
      double precision :: acks2_solver_guess_out(*)
    end subroutine ereaxff_force_eval
    subroutine ereaxff_oracle_get_energy_terms(terms_kcal, total_kcal)
      double precision :: terms_kcal(14), total_kcal
    end subroutine ereaxff_oracle_get_energy_terms
    subroutine ereaxff_oracle_get_bond_orders(natom_eval, nbond_eval, abo_out, bond_atoms_out, bo_out, bopi_out, bopi2_out)
      integer :: natom_eval, nbond_eval
      integer :: bond_atoms_out(2,*)
      double precision :: abo_out(*), bo_out(*), bopi_out(*), bopi2_out(*)
    end subroutine ereaxff_oracle_get_bond_orders
    subroutine ereaxff_oracle_get_nonbonded_buckets(max_rows, nrows, left_type_out, right_type_out, self_out, count_out, vdw_out, coulomb_out)
      integer :: max_rows, nrows
      integer :: left_type_out(*), right_type_out(*), self_out(*), count_out(*)
      double precision :: vdw_out(*), coulomb_out(*)
    end subroutine ereaxff_oracle_get_nonbonded_buckets
  end interface

  call get_command_argument(1, control_path)
  call get_command_argument(2, ffield_path)
  call get_command_argument(3, iopt_path)
  call get_command_argument(4, geo_path)
  call get_command_argument(5, timing_log_path)
  call get_command_argument(6, worker_root)
  call get_command_argument(7, smoke_binary)
  call get_command_argument(8, make_geo_script)
  if (len_trim(timing_log_path) == 0) timing_log_path = "persistent_oracle_timing.log"
  if (len_trim(worker_root) == 0) worker_root = "persistent_worker"
  if (len_trim(smoke_binary) == 0) smoke_binary = "phase4_single_point_oracle"
  if (len_trim(make_geo_script) == 0) make_geo_script = "make_geo.py"

  request_file = trim(worker_root) // "/persistent_request.in"
  response_file = trim(worker_root) // "/persistent_response.ready"
  ready_file = trim(worker_root) // "/persistent_worker.ready"
  shutdown_file = trim(worker_root) // "/persistent_shutdown.flag"
  done_file = trim(worker_root) // "/persistent_worker.done"
  first_request = .true.

  call stage_static_inputs(trim(worker_root), trim(control_path), trim(ffield_path), trim(iopt_path), trim(geo_path))
  call ereaxff_force_oracle_init()
  call write_text_file(trim(ready_file), "ready")

  do
    if (file_exists(trim(shutdown_file))) exit
    if (file_exists(trim(request_file))) then
      call process_request(trim(worker_root), trim(request_file), trim(response_file), trim(timing_log_path), trim(smoke_binary), &
      trim(make_geo_script), first_request, request_id, step_index, bead_id, stage_label, natom, energy_kcal, force, charge_out, solver_out, status, charge_sum, &
          max_force, energy_terms, energy_terms_total, bo_natom, nbond, abo, bond_atoms, bo, bopi, bopi2, &
          nbuckets, bucket_left, bucket_right, bucket_self, bucket_count, bucket_vdw, bucket_coulomb, work_dir)
      first_request = .false.
    else
      call wait_milliseconds(10)
    end if
  end do

  call write_text_file(trim(done_file), "done")

contains

  subroutine process_request(worker_root, request_path, response_path, timing_log, smoke_binary, make_geo_script, first_request_inout, request_id_out, step_out, &
      bead_out, stage_out, natom_out, energy_out, force_out, charge_out, solver_out, status_out, charge_sum_out, max_force_out, &
      energy_terms_out, energy_terms_total_out, bo_natom_out, nbond_out, abo_out, bond_atoms_out, bo_out, bopi_out, bopi2_out, &
      nbuckets_out, bucket_left_out, bucket_right_out, bucket_self_out, bucket_count_out, bucket_vdw_out, bucket_coulomb_out, work_dir_out)
    character(len=*), intent(in) :: worker_root, request_path, response_path, timing_log, smoke_binary, make_geo_script
    logical, intent(inout) :: first_request_inout
    integer, intent(out) :: request_id_out, step_out, bead_out, natom_out, status_out
    integer, intent(out) :: bo_natom_out, nbond_out, bond_atoms_out(:, :)
    integer, intent(out) :: nbuckets_out, bucket_left_out(:), bucket_right_out(:), bucket_self_out(:), bucket_count_out(:)
    character(len=*), intent(out) :: stage_out, work_dir_out
    real(dp), intent(out) :: energy_out, force_out(:, :), charge_out(:), solver_out(:), charge_sum_out, max_force_out
    real(dp), intent(out) :: energy_terms_out(14), energy_terms_total_out
    real(dp), intent(out) :: abo_out(:), bo_out(:), bopi_out(:), bopi2_out(:)
    real(dp), intent(out) :: bucket_vdw_out(:), bucket_coulomb_out(:)
    integer :: unit, i, atom_index, start_count_local, end_count_local, rate
    character(len=512) :: line
    character(len=64) :: worker_state

    open(newunit=unit, file=trim(request_path), status="old", action="read", iostat=ios)
    if (ios /= 0) stop "persistent_oracle_worker: unable to open request"
    read(unit, *, iostat=ios) request_id_out
    if (ios /= 0) stop "persistent_oracle_worker: malformed request id"
    read(unit, *, iostat=ios) step_out, bead_out, natom_out
    if (ios /= 0) stop "persistent_oracle_worker: malformed request header"
    read(unit, '(A)', iostat=ios) line
    if (ios /= 0) stop "persistent_oracle_worker: missing stage label"
    stage_out = adjustl(trim(line))
    if (natom_out < 1 .or. natom_out > natmax) stop "persistent_oracle_worker: natom out of range"
    x = 0.0_dp
    charge_in = 0.0_dp
    solver_in = 0.0_dp
    do i = 1, natom_out
      read(unit, *, iostat=ios) atom_index, x(1,i), x(2,i), x(3,i), charge_in(i), solver_in(i)
      if (ios /= 0) stop "persistent_oracle_worker: malformed atom line"
    end do
    close(unit)

    work_dir_out = trim(worker_root) // "/adapter_call_" // trim(pad_int(request_id_out, 4))
    call ensure_workdir(trim(work_dir_out))
    summary_path = trim(work_dir_out) // "/phase4_single_point_oracle.out"
    forces_path = trim(work_dir_out) // "/single_point_forces.out"
    completion_path = trim(work_dir_out) // "/persistent_request.done"
    call write_request_inputs(trim(worker_root), trim(work_dir_out), trim(make_geo_script), x, natom_out, charge_in, solver_in)

    call system_clock(start_count_local, rate)
    if (rate <= 0) stop "persistent_oracle_worker: system clock unavailable"
    call ereaxff_force_eval(x, natom_out, charge_in, solver_in, energy_out, force_out, charge_out, solver_out, status_out)
    call ereaxff_oracle_get_energy_terms(energy_terms_out, energy_terms_total_out)
    call ereaxff_oracle_get_bond_orders(bo_natom_out, nbond_out, abo_out, bond_atoms_out, bo_out, bopi_out, bopi2_out)
    call ereaxff_oracle_get_nonbonded_buckets(nbucketmax, nbuckets_out, bucket_left_out, bucket_right_out, bucket_self_out, &
        bucket_count_out, bucket_vdw_out, bucket_coulomb_out)
    if (bo_natom_out /= natom_out) stop "persistent_oracle_worker: bond-order atom count mismatch"
    if (nbond_out > nbondmax) stop "persistent_oracle_worker: too many bond-order rows"
    call system_clock(end_count_local)

    call summarize_force_result(natom_out, force_out, charge_out, max_force_out, charge_sum_out)
    call write_summary(trim(summary_path), energy_out, max_force_out, charge_sum_out, status_out, natom_out)
    call write_force_file(trim(forces_path), natom_out, force_out, charge_out, solver_out)
    call write_energy_terms(trim(work_dir_out) // "/single_point_energy_terms.out", energy_terms_out, energy_terms_total_out)
    call write_bond_orders(trim(work_dir_out) // "/single_point_bond_orders.out", natom_out, nbond_out, abo_out, bond_atoms_out, bo_out, bopi_out, bopi2_out)
    call write_nonbonded_buckets(trim(work_dir_out) // "/single_point_nonbonded_buckets.out", nbuckets_out, bucket_left_out, bucket_right_out, &
        bucket_self_out, bucket_count_out, bucket_vdw_out, bucket_coulomb_out)
    if (first_request_inout) then
      worker_state = "startup"
    else
      worker_state = "reused"
    end if
    call append_timing(trim(timing_log), request_id_out, step_out, bead_out, trim(stage_out), trim(worker_state), &
        end_count_local - start_count_local, status_out, natom_out)
    call write_text_file(trim(completion_path), "done")
    call write_text_file_atomic(trim(response_path), trim(itoa(request_id_out)))
    call remove_file(trim(request_path))
  end subroutine process_request

  subroutine stage_static_inputs(worker_root, control_src, ffield_src, iopt_src, geo_src)
    character(len=*), intent(in) :: worker_root, control_src, ffield_src, iopt_src, geo_src
    real(dp) :: crystx(6)
    logical :: has_crystx

    call ensure_workdir(trim(worker_root))
    call copy_file(ffield_src, trim(worker_root) // "/ffield")
    call copy_file(iopt_src, trim(worker_root) // "/fort.20")
    call write_text_file(trim(worker_root) // "/fort.35", "0")
    if (len_trim(geo_src) > 0) then
      call copy_file(geo_src, trim(worker_root) // "/geo")
      call extract_crystx(trim(worker_root) // "/geo", crystx, has_crystx)
    else
      has_crystx = .false.
    end if
    call write_control_with_cell(control_src, trim(worker_root) // "/control", crystx, has_crystx)
  end subroutine stage_static_inputs

  subroutine write_request_inputs(worker_root, work_dir, make_geo_script, positions_A, natom, charge_in, solver_in)
    character(len=*), intent(in) :: worker_root, work_dir, make_geo_script
    real(dp), intent(in) :: positions_A(:, :)
    integer, intent(in) :: natom
    real(dp), intent(in) :: charge_in(:), solver_in(:)
    real(dp) :: crystx(6)
    logical :: has_crystx

    call copy_file(trim(worker_root) // "/ffield", trim(work_dir) // "/ffield")
    call copy_file(trim(worker_root) // "/fort.20", trim(work_dir) // "/fort.20")
    call write_text_file(trim(work_dir) // "/fort.35", "0")
    if (file_exists(trim(worker_root) // "/geo")) then
      call write_geo_with_positions(trim(worker_root) // "/geo", trim(work_dir) // "/geo", positions_A, natom)
    else
      call generate_phase4_geo(trim(work_dir), trim(make_geo_script), positions_A)
    end if
    call extract_crystx(trim(work_dir) // "/geo", crystx, has_crystx)
    call write_control_with_cell(trim(worker_root) // "/control", trim(work_dir) // "/control", crystx, has_crystx)
    call write_charge_trace(trim(work_dir) // "/charge_memory.in", natom, charge_in, solver_in)
  end subroutine write_request_inputs

  subroutine write_geo_with_positions(source, target, positions_A, natom)
    character(len=*), intent(in) :: source, target
    real(dp), intent(in) :: positions_A(:, :)
    integer, intent(in) :: natom
    integer, parameter :: line_len = 8192
    integer :: src_unit, dst_unit, ios, atom_counter
    character(len=line_len) :: line
    character(len=8) :: record, res3, res5
    character(len=5) :: atom_name, atom_type
    character(len=1) :: res1
    integer :: atom_id_read, dummy_int1, dummy_int2
    real(dp) :: x, y, z, dummy_real

    open(newunit=src_unit, file=trim(source), status="old", action="read", iostat=ios)
    if (ios /= 0) stop "persistent_oracle_worker: unable to open source geo"
    open(newunit=dst_unit, file=trim(target), status="replace", action="write", iostat=ios)
    if (ios /= 0) stop "persistent_oracle_worker: unable to open target geo"
    atom_counter = 0
    do
      read(src_unit,'(A)',iostat=ios) line
      if (ios /= 0) exit
      if (starts_with(trim(line), "HETATM")) then
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
        if (ios /= 0) stop "persistent_oracle_worker: malformed HETATM record"
        atom_counter = atom_counter + 1
        if (atom_counter > natom) stop "persistent_oracle_worker: too many atoms"
        write(dst_unit,'(a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)') &
            record, atom_id_read, atom_name, res3, res1, atom_type, positions_A(:, atom_counter), &
            res5, dummy_int1, dummy_int2, dummy_real
      else
        write(dst_unit,'(a)') trim(line)
      end if
    end do
    close(src_unit)
    close(dst_unit)
  end subroutine write_geo_with_positions

  subroutine generate_phase4_geo(work_dir, make_geo_script, positions_A)
    character(len=*), intent(in) :: work_dir, make_geo_script
    real(dp), intent(in) :: positions_A(:, :)
    character(len=64) :: el_x, el_y, el_z
    character(len=512) :: command
    integer :: exitstat

    if (size(positions_A, 1) /= 3) stop "persistent_oracle_worker: positions must be 3D"
    write(el_x,'(es24.15)') positions_A(1, 1)
    write(el_y,'(es24.15)') positions_A(2, 1)
    write(el_z,'(es24.15)') positions_A(3, 1)
    command = "cd " // shell_quote(trim(work_dir)) // " && PHASE4_EL_X=" // trim(adjustl(el_x)) // &
        " PHASE4_EL_Y=" // trim(adjustl(el_y)) // " PHASE4_EL_Z=" // trim(adjustl(el_z)) // &
        " python3 " // shell_quote(trim(make_geo_script))
    call run_shell(trim(command), exitstat)
    if (exitstat /= 0) stop "persistent_oracle_worker: geometry generation failed"
  end subroutine generate_phase4_geo

  subroutine extract_crystx(path, crystx, has_crystx)
    character(len=*), intent(in) :: path
    real(dp), intent(out) :: crystx(6)
    logical, intent(out) :: has_crystx
    integer, parameter :: line_len = 8192
    integer :: unit, ios
    character(len=line_len) :: line
    character(len=16) :: label
    real(dp) :: values(6)

    crystx = 0.0_dp
    has_crystx = .false.
    open(newunit=unit, file=trim(path), status="old", action="read", iostat=ios)
    if (ios /= 0) return
    do
      read(unit,'(A)',iostat=ios) line
      if (ios /= 0) exit
      label = ""
      values = 0.0_dp
      read(line, *, iostat=ios) label, values
      if (ios == 0 .and. trim(label) == "CRYSTX") then
        crystx = values
        has_crystx = .true.
        exit
      end if
    end do
    close(unit)
  end subroutine extract_crystx

  subroutine write_control_with_cell(source, target, crystx, has_crystx)
    character(len=*), intent(in) :: source, target
    real(dp), intent(in) :: crystx(6)
    logical, intent(in) :: has_crystx
    integer, parameter :: line_len = 8192
    integer :: src_unit, dst_unit, ios
    character(len=line_len) :: line
    logical :: replaced

    open(newunit=src_unit, file=trim(source), status="old", action="read", iostat=ios)
    if (ios /= 0) stop "persistent_oracle_worker: unable to open source control"
    open(newunit=dst_unit, file=trim(target), status="replace", action="write", iostat=ios)
    if (ios /= 0) stop "persistent_oracle_worker: unable to open target control"

    do
      read(src_unit,'(A)',iostat=ios) line
      if (ios /= 0) exit
      replaced = .false.
      if (has_crystx) then
        call maybe_write_control_axis(dst_unit, line, "axis1", crystx(1), replaced)
        if (.not. replaced) call maybe_write_control_axis(dst_unit, line, "axis2", crystx(2), replaced)
        if (.not. replaced) call maybe_write_control_axis(dst_unit, line, "axis3", crystx(3), replaced)
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

  subroutine write_summary(path, energy_kcal, max_force, charge_sum, status, natom)
    character(len=*), intent(in) :: path
    real(dp), intent(in) :: energy_kcal, max_force, charge_sum
    integer, intent(in) :: status, natom
    integer :: unit

    open(newunit=unit, file=trim(path), status="replace", action="write")
    write(unit,'(a)') '# energy_kcal max_force charge_sum status natom'
    write(unit,'(es24.15,1x,es24.15,1x,es24.15,1x,i8,1x,i8)') energy_kcal, max_force, charge_sum, status, natom
    close(unit)
  end subroutine write_summary

  subroutine write_force_file(path, natom, force_in, charge_in, solver_in)
    character(len=*), intent(in) :: path
    integer, intent(in) :: natom
    real(dp), intent(in) :: force_in(:, :), charge_in(:), solver_in(:)
    integer :: unit, j

    open(newunit=unit, file=trim(path), status="replace", action="write")
    write(unit,'(a)') '# atom Fx Fy Fz charge solver_guess'
    do j = 1, natom
      write(unit,'(i8,5es24.15)') j, force_in(1, j), force_in(2, j), force_in(3, j), charge_in(j), solver_in(j)
    end do
    close(unit)
  end subroutine write_force_file

  subroutine write_energy_terms(path, terms_kcal, total_kcal)
    character(len=*), intent(in) :: path
    real(dp), intent(in) :: terms_kcal(14), total_kcal
    character(len=32), parameter :: names(14) = [ &
        character(len=32) :: &
        'bond', 'overcoordination', 'lone_pair', 'molecule', &
        'valence_penalty', 'angle_conjugation', 'hydrogen_bond', &
        'torsion', 'torsion_conjugation', 'vdw', 'coulomb', &
        'charge', 'external_field', 'electron' &
    ]
    integer :: unit, i

    open(newunit=unit, file=trim(path), status="replace", action="write")
    write(unit,'(a)') '# term energy_kcal'
    do i = 1, size(names)
      write(unit,'(a,1x,es24.15)') trim(names(i)), terms_kcal(i)
    end do
    write(unit,'(a,1x,es24.15)') 'total', total_kcal
    close(unit)
  end subroutine write_energy_terms

  subroutine write_bond_orders(path, natom_in, nbond_in, abo_in, bond_atoms_in, bo_in, bopi_in, bopi2_in)
    character(len=*), intent(in) :: path
    integer, intent(in) :: natom_in, nbond_in, bond_atoms_in(:, :)
    real(dp), intent(in) :: abo_in(:), bo_in(:), bopi_in(:), bopi2_in(:)
    integer :: unit, i

    open(newunit=unit, file=trim(path), status="replace", action="write")
    write(unit,'(a)') '# kind index atom_i atom_j value bopi bopi2'
    do i = 1, natom_in
      write(unit,'(a,1x,i8,1x,es24.15)') 'ATOM', i, abo_in(i)
    end do
    do i = 1, nbond_in
      write(unit,'(a,1x,i8,2(1x,i8),3es24.15)') 'BOND', i, bond_atoms_in(1,i), bond_atoms_in(2,i), &
          bo_in(i), bopi_in(i), bopi2_in(i)
    end do
    close(unit)
  end subroutine write_bond_orders

  subroutine write_nonbonded_buckets(path, nrows, left_type, right_type, self_image, bucket_count, vdw_kcal, coulomb_kcal)
    character(len=*), intent(in) :: path
    integer, intent(in) :: nrows, left_type(:), right_type(:), self_image(:), bucket_count(:)
    real(dp), intent(in) :: vdw_kcal(:), coulomb_kcal(:)
    integer :: unit, i

    open(newunit=unit, file=trim(path), status="replace", action="write")
    write(unit,'(a)') '# left_type right_type self_image count vdw_kcal coulomb_kcal'
    do i = 1, nrows
      write(unit,'(4i8,2es24.15)') left_type(i), right_type(i), self_image(i), bucket_count(i), vdw_kcal(i), coulomb_kcal(i)
    end do
    close(unit)
  end subroutine write_nonbonded_buckets

  subroutine write_charge_trace(path, natom, charge_in, solver_in)
    character(len=*), intent(in) :: path
    integer, intent(in) :: natom
    real(dp), intent(in) :: charge_in(:), solver_in(:)
    integer :: unit, i

    open(newunit=unit, file=trim(path), status="replace", action="write")
    do i = 1, natom
      write(unit,'(i0,1x,es24.15,1x,es24.15)') i, charge_in(i), solver_in(i)
    end do
    close(unit)
  end subroutine write_charge_trace

  subroutine append_timing(path, request_id_in, step_in, bead_in, stage_in, worker_state_in, elapsed_count, status_in, natom_in)
    character(len=*), intent(in) :: path, stage_in, worker_state_in
    integer, intent(in) :: request_id_in, step_in, bead_in, elapsed_count, status_in, natom_in
    integer :: unit, rate, ios
    real(dp) :: elapsed_s
    logical :: exists

    call system_clock(count_rate=rate)
    if (rate <= 0) stop "persistent_oracle_worker: clock rate unavailable"
    elapsed_s = real(elapsed_count, dp) / real(rate, dp)
    inquire(file=trim(path), exist=exists)
    open(newunit=unit, file=trim(path), status="unknown", action="write", position="append", iostat=ios)
    if (ios /= 0) stop "persistent_oracle_worker: unable to append timing log"
    if (.not. exists) then
      write(unit,'(a)') '# kind request step bead stage worker_state elapsed_s status natom'
    end if
    write(unit,'(a,1x,i0,1x,i0,1x,i0,1x,a,1x,a,1x,es24.15,1x,i0,1x,i0)') &
        'REQ', request_id_in, step_in, bead_in, trim(stage_in), trim(worker_state_in), elapsed_s, status_in, natom_in
    close(unit)
  end subroutine append_timing

  subroutine ensure_workdir(path)
    character(len=*), intent(in) :: path
    character(len=1024) :: command
    integer :: exitstat

    command = "mkdir -p " // shell_quote(trim(path))
    call execute_command_line(trim(command), wait=.true., exitstat=exitstat)
    if (exitstat /= 0) stop "persistent_oracle_worker: unable to create workdir"
  end subroutine ensure_workdir

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

  subroutine copy_file(source, target)
    character(len=*), intent(in) :: source, target
    integer :: src_unit, dst_unit, ios
    character(len=8192) :: line

    open(newunit=src_unit, file=trim(source), status="old", action="read", iostat=ios)
    if (ios /= 0) stop "persistent_oracle_worker: unable to open source file"
    open(newunit=dst_unit, file=trim(target), status="replace", action="write", iostat=ios)
    if (ios /= 0) stop "persistent_oracle_worker: unable to open target file"
    do
      read(src_unit,'(A)',iostat=ios) line
      if (ios /= 0) exit
      write(dst_unit,'(a)') trim(line)
    end do
    close(src_unit)
    close(dst_unit)
  end subroutine copy_file

  subroutine write_text_file(path, text)
    character(len=*), intent(in) :: path, text
    integer :: unit

    open(newunit=unit, file=trim(path), status="replace", action="write")
    write(unit,'(a)') trim(text)
    close(unit)
  end subroutine write_text_file

  subroutine write_text_file_atomic(path, text)
    character(len=*), intent(in) :: path, text
    character(len=1024) :: tmp_path

    tmp_path = trim(path) // ".tmp"
    call remove_file(trim(tmp_path))
    call write_text_file(trim(tmp_path), trim(text))
    call publish_file_atomic(trim(tmp_path), trim(path))
  end subroutine write_text_file_atomic

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
    if (exitstat /= 0) stop "persistent_oracle_worker: unable to publish file"
  end subroutine publish_file_atomic

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

  pure function shell_quote(text) result(quoted)
    character(len=*), intent(in) :: text
    character(len=:), allocatable :: quoted
    integer :: n
    n = len_trim(text)
    allocate(character(len=n + 2) :: quoted)
    quoted = "'" // trim(text) // "'"
  end function shell_quote

  subroutine parse_single_point_summary(path, energy_kcal, max_force, charge_sum, status, natom)
    character(len=*), intent(in) :: path
    real(dp), intent(out) :: energy_kcal, max_force, charge_sum
    integer, intent(out) :: status, natom
    integer :: unit, ios
    character(len=512) :: line

    open(newunit=unit, file=trim(path), status="old", action="read", iostat=ios)
    if (ios /= 0) stop "persistent_oracle_worker: unable to open summary"
    read(unit,'(A)',iostat=ios) line
    if (ios /= 0) stop "persistent_oracle_worker: missing summary header"
    read(unit,*,iostat=ios) energy_kcal, max_force, charge_sum, status, natom
    if (ios /= 0) stop "persistent_oracle_worker: malformed summary line"
    close(unit)
  end subroutine parse_single_point_summary

  subroutine parse_single_point_forces(path, force_out, charge_out, solver_out)
    character(len=*), intent(in) :: path
    real(dp), intent(out) :: force_out(:, :), charge_out(:), solver_out(:)
    integer :: unit, ios, idx, natom, i
    character(len=512) :: line
    real(dp) :: fx, fy, fz, charge, solver

    open(newunit=unit, file=trim(path), status="old", action="read", iostat=ios)
    if (ios /= 0) stop "persistent_oracle_worker: unable to open forces"
    read(unit,'(A)',iostat=ios) line
    if (ios /= 0) stop "persistent_oracle_worker: missing forces header"
    natom = 0
    do
      read(unit,*,iostat=ios) idx, fx, fy, fz, charge, solver
      if (ios /= 0) exit
      natom = natom + 1
    end do
    rewind(unit)
    read(unit,'(A)',iostat=ios) line
    if (ios /= 0) stop "persistent_oracle_worker: header rewind failed"
    do i = 1, natom
      read(unit,*,iostat=ios) idx, fx, fy, fz, charge, solver
      if (ios /= 0) stop "persistent_oracle_worker: malformed force row"
      force_out(:, i) = [fx, fy, fz]
      charge_out(i) = charge
      solver_out(i) = solver
    end do
    close(unit)
  end subroutine parse_single_point_forces

  subroutine summarize_force_result(natom, force_in, charge_in, max_force, charge_sum)
    integer, intent(in) :: natom
    real(dp), intent(in) :: force_in(:, :), charge_in(:)
    real(dp), intent(out) :: max_force, charge_sum
    integer :: i

    max_force = 0.0_dp
    charge_sum = 0.0_dp
    do i = 1, natom
      max_force = max(max_force, abs(force_in(1, i)))
      max_force = max(max_force, abs(force_in(2, i)))
      max_force = max(max_force, abs(force_in(3, i)))
      charge_sum = charge_sum + charge_in(i)
    end do
  end subroutine summarize_force_result

  subroutine run_shell(command, exitstat)
    character(len=*), intent(in) :: command
    integer, intent(out) :: exitstat
    integer :: cmdstat
    character(len=256) :: errmsg

    call execute_command_line(trim(command), wait=.true., exitstat=exitstat, cmdstat=cmdstat, cmdmsg=errmsg)
    if (cmdstat /= 0) stop "persistent_oracle_worker: " // trim(errmsg)
  end subroutine run_shell

end program persistent_oracle_worker
