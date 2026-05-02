program ereaxff_mixed_baoab_main
  use ereaxff_kinds, only: dp
  use ereaxff_system_state, only: system_state, init_system_state, set_particle
  use ereaxff_particle_groups, only: particle_group_spec, particle_group_assignment, &
      assign_particle_groups, count_particles_with_dynamics, dynamics_classical, dynamics_pimd
  use ereaxff_input_config, only: input_config, read_input_config
  use ereaxff_ring_polymer, only: ring_polymer_state, init_ring_polymer
  use ereaxff_normal_modes, only: bead_to_nm, nm_to_bead
  use ereaxff_charge_coupling_policy, only: charge_memory_state, init_charge_memory, &
      charge_policy_bead_local_acks2, charge_policy_is_valid, max_charge_difference_across_beads
  use ereaxff_oracle_adapter, only: oracle_adapter_state, oracle_adapter_init
  use ereaxff_bead_force_manager, only: oracle_bead_summary, evaluate_oracle_beads
  use ereaxff_mixed_baoab_integrator, only: mixed_baoab_electron_step
  use ereaxff_mixed_baoab_restart, only: mixed_baoab_restart_state, write_mixed_baoab_restart, &
      read_mixed_baoab_restart
  use ereaxff_mixed_restart, only: write_charge_by_bead_summary, write_charge_centroid_summary, &
      write_electron_localization_summary
  use ereaxff_reaxff_compat_output, only: write_reaxff_compat_snapshot
  use ereaxff_output_manager, only: join_path
  use ereaxff_spring_energy, only: compute_free_ring_energy_nm, compute_spring_energy_nm
  use ereaxff_diagnostics, only: nvt_conserved_like
  implicit none

  integer, parameter :: ndim = 3
  integer, parameter :: nparticles = 61
  real(dp), parameter :: electron_mass_amu = 0.00054858_dp
  real(dp), parameter :: electron_shift_x = 0.08_dp
  real(dp), parameter :: electron_shift_y = -0.03_dp
  type(system_state) :: state
  type(input_config) :: cfg
  type(particle_group_spec), allocatable :: groups(:)
  type(particle_group_assignment) :: assignment
  type(ring_polymer_state) :: rp
  type(oracle_adapter_state) :: adapter
  type(charge_memory_state) :: charge_memory
  type(oracle_bead_summary) :: summary
  type(mixed_baoab_restart_state) :: restart_state
  real(dp), allocatable :: positions_A(:, :, :)
  real(dp), allocatable :: q_nm(:, :)
  real(dp), allocatable :: v_nm(:, :)
  real(dp), allocatable :: solver_guess(:, :)
  real(dp) :: masses(3), kinetic_eV, spring_eV, physical_potential_eV, conserved_like
  real(dp) :: centroid_A(3), bead_spread_A, max_force, charge_diff
  character(len=512) :: arg1, arg2, beads_arg, steps_arg, policy_arg, restart_in, restart_out
  character(len=512) :: config_path, output_dir, adapter_work_root
  character(len=512) :: v1_root, v2_root
  character(len=1100) :: mkdir_command
  integer :: beads, steps, ios, step, start_step, energy_unit, spread_unit, heat_unit, xyz_unit
  integer :: pimd_index, pimd_count, output_stride
  real(dp) :: temperature_K, dt_fs, pimd_mass_amu
  logical :: config_mode, exists

  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  inquire(file=trim(arg1), exist=exists)
  config_mode = exists .and. (index(trim(arg1), ".yaml") > 0 .or. index(trim(arg1), ".yml") > 0)
  if (config_mode) then
    config_path = trim(arg1)
    call read_input_config(trim(config_path), cfg)
    output_dir = trim(arg2)
    if (len_trim(output_dir) == 0) output_dir = trim(cfg%simulation%output_dir)
    if (len_trim(output_dir) == 0) output_dir = "ereaxff_mixed_baoab_output"
    beads = cfg%ring_polymer%default_beads
    steps = cfg%simulation%steps
    policy_arg = cfg%charge_coupling%mode
    dt_fs = cfg%simulation%dt_fs
    temperature_K = select_mixed_temperature(cfg)
    call get_command_argument(3, restart_in)
    call get_command_argument(4, restart_out)
  else
    output_dir = trim(arg1)
    if (len_trim(output_dir) == 0) output_dir = "ereaxff_mixed_baoab_output"
    call get_command_argument(2, beads_arg)
    beads = 4
    if (len_trim(beads_arg) > 0) then
      read(beads_arg, *, iostat=ios) beads
      if (ios /= 0 .or. beads < 1) error stop "ereaxff_mixed_baoab_main: invalid beads argument"
    end if
    call get_command_argument(3, steps_arg)
    steps = 8
    if (len_trim(steps_arg) > 0) then
      read(steps_arg, *, iostat=ios) steps
      if (ios /= 0 .or. steps < 0) error stop "ereaxff_mixed_baoab_main: invalid steps argument"
    end if
    call get_command_argument(4, policy_arg)
    call get_command_argument(5, restart_in)
    call get_command_argument(6, restart_out)
    if (len_trim(policy_arg) == 0) policy_arg = charge_policy_bead_local_acks2
    temperature_K = 300.0_dp
    dt_fs = 0.1_dp
  end if
  if (.not. charge_policy_is_valid(policy_arg)) error stop "ereaxff_mixed_baoab_main: invalid charge policy"
  call get_environment_variable("EREAXFF_V1_ROOT", v1_root)
  call get_environment_variable("EREAXFF_V2_ROOT", v2_root)
  if (len_trim(v1_root) == 0) v1_root = "/Users/yeatanir/PhD/PIMD/ereaxff_pimd/V1"
  if (len_trim(v2_root) == 0) v2_root = "/Users/yeatanir/PhD/PIMD/ereaxff_pimd/V2"
  if (config_mode) then
    if (len_trim(cfg%force_provider%v1_path) > 0) v1_root = trim(cfg%force_provider%v1_path)
    if (len_trim(cfg%force_provider%v2_path) > 0) v2_root = trim(cfg%force_provider%v2_path)
  end if

  mkdir_command = 'mkdir -p "' // trim(output_dir) // '"'
  call execute_command_line(trim(mkdir_command))

  write(*,'(a)') "PI-eReaxFF mixed eReaxFF BAOAB driver"
  write(*,'(a)') "build/version: development"
  write(*,'(a)') "ring-polymer convention: standard_ipi"
  write(*,'(a,i0)') "P: ", beads
  write(*,'(a,es12.4)') "T_K: ", temperature_K
  write(*,'(a)') "thermostat mode: none"
  write(*,'(a)') "integrator mode: mixed_baoab"
  write(*,'(a)') "charge coupling mode: " // trim(policy_arg)
  write(*,'(a)') "force provider: eReaxFF adapter"
  if (config_mode) then
    write(*,'(a)') "input boundary: config-driven legacy V1/V2 fixtures"
    write(*,'(a)') "input file: " // trim(config_path)
  else
    write(*,'(a)') "input boundary: legacy V1/V2 fixtures copied by adapter"
  end if
  write(*,'(a)') "output directory: " // trim(output_dir)

  call init_system_state(state, ndim, nparticles)
  call build_phase4_fixture(state)

  start_step = 0
  if (len_trim(restart_in) > 0 .and. trim(restart_in) /= "none") then
    call read_mixed_baoab_restart(trim(restart_in), restart_state)
    beads = restart_state%beads
    policy_arg = restart_state%charge_policy
    if (.not. charge_policy_is_valid(policy_arg)) error stop "ereaxff_mixed_baoab_main: restart charge policy invalid"
    start_step = restart_state%step
  end if

  if (config_mode) then
    if (.not. allocated(cfg%groups)) error stop "ereaxff_mixed_baoab_main: config groups missing"
    groups = cfg%groups
  else
    call build_default_groups(groups, beads, temperature_K)
  end if
  call assign_particle_groups(state, groups, assignment)
  pimd_count = count_particles_with_dynamics(assignment, dynamics_pimd)
  if (pimd_count /= 1) error stop "ereaxff_mixed_baoab_main: exactly one PIMD particle is supported"
  pimd_index = first_particle_with_dynamics(assignment, dynamics_pimd)
  beads = assignment%beads(pimd_index)

  call init_ring_polymer(rp, beads, temperature_K, "standard_ipi")
  allocate(positions_A(ndim, nparticles, beads))
  allocate(q_nm(ndim, beads))
  allocate(v_nm(ndim, beads))
  allocate(solver_guess(nparticles, beads))

  if (len_trim(restart_in) > 0 .and. trim(restart_in) /= "none") then
    positions_A = restart_state%positions_A
    q_nm = restart_state%q_nm
    v_nm = restart_state%v_nm
    solver_guess = restart_state%solver_guess_by_bead
  else
    call build_positions(state, positions_A, beads)
    call perturb_pimd_beads(positions_A, pimd_index)
    q_nm = 0.0_dp
    v_nm = 0.0_dp
    call bead_positions_to_nm(positions_A, pimd_index, rp, q_nm)
    solver_guess = 0.0_dp
  end if

  call init_charge_memory(charge_memory, nparticles, beads, policy_arg)
  if (config_mode) then
    adapter_work_root = join_path(output_dir, "adapter_work")
    if (len_trim(cfg%force_provider%work_root) > 0) adapter_work_root = trim(cfg%force_provider%work_root)
    call oracle_adapter_init(adapter, state, trim(v1_root), trim(v2_root), trim(adapter_work_root), &
        control_path=join_path(v1_root, trim(cfg%force_provider%control)), &
        ffield_path=join_path(v1_root, trim(cfg%force_provider%ffield)), &
        iopt_path=join_path(v1_root, trim(cfg%force_provider%iopt)), &
        smoke_binary=trim(cfg%force_provider%smoke_binary), repeat_binary=trim(cfg%force_provider%repeat_binary), &
        make_geo_script=trim(cfg%force_provider%make_geo_script), timeout_runner=trim(cfg%force_provider%timeout_runner))
  else
    call oracle_adapter_init(adapter, state, trim(v1_root), trim(v2_root), join_path(output_dir, "adapter_work"))
  end if
  pimd_mass_amu = state%mass_amu(pimd_index)
  masses = [pimd_mass_amu, pimd_mass_amu, pimd_mass_amu]

  call evaluate_oracle_beads(adapter, positions_A, charge_memory, solver_guess, policy_arg, summary)
  output_stride = 1
  if (config_mode) output_stride = max(1, cfg%simulation%output_stride)

  open(newunit=energy_unit, file=join_path(output_dir, "pimd_energy.out"), status="replace", action="write")
  open(newunit=spread_unit, file=join_path(output_dir, "bead_spread.out"), status="replace", action="write")
  open(newunit=heat_unit, file=join_path(output_dir, "thermostat_heat.out"), status="replace", action="write")
  open(newunit=xyz_unit, file=join_path(output_dir, "centroid.xyz"), status="replace", action="write")

  write(energy_unit,'(a)') "# step kinetic_eV physical_potential_eV spring_eV conserved_like_eV charge_diff max_force"
  write(spread_unit,'(a)') "# step centroid_x centroid_y centroid_z bead_rms_spread_A"
  write(heat_unit,'(a)') "# step thermostat_heat_eV"

  if (mod(start_step, output_stride) == 0 .or. start_step == steps) call write_step_outputs(start_step)
  step = start_step
  do while (step < steps)
    call mixed_baoab_electron_step(positions_A, q_nm, v_nm, pimd_index, pimd_mass_amu, rp, dt_fs, &
        adapter, charge_memory, solver_guess, policy_arg, summary)
    step = step + 1
    if (mod(step, output_stride) == 0 .or. step == steps) call write_step_outputs(step)
  end do

  close(energy_unit)
  close(spread_unit)
  close(heat_unit)
  close(xyz_unit)

  call write_final_restart_and_compat()

  write(*,'(a,i0)') "classical_count = ", count_particles_with_dynamics(assignment, dynamics_classical)
  write(*,'(a,i0)') "pimd_count = ", count_particles_with_dynamics(assignment, dynamics_pimd)

contains

  subroutine build_phase4_fixture(state_out)
    type(system_state), intent(inout) :: state_out
    integer :: ix, iy, unit_index, idx
    real(dp), parameter :: a_cc = 1.4210_dp
    real(dp), parameter :: z_c = 10.0_dp
    real(dp), parameter :: el_x = 7.38300_dp
    real(dp), parameter :: el_y = 5.68433_dp
    real(dp), parameter :: el_z = 10.001_dp
    real(dp) :: ax, ay, x, y
    real(dp), dimension(4, 2) :: unit

    unit = reshape([0.0_dp, 0.0_dp, 0.5_dp, 1.0_dp/6.0_dp, 0.5_dp, 3.0_dp/6.0_dp, &
        0.0_dp, 4.0_dp/6.0_dp], shape(unit))
    ax = sqrt(3.0_dp) * a_cc
    ay = 3.0_dp * a_cc

    idx = 1
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
    call set_particle(state_out, 61, 61, "El", electron_mass_amu, [el_x, el_y, el_z])
  end subroutine build_phase4_fixture

  subroutine build_positions(state_in, positions_out, beads_in)
    type(system_state), intent(in) :: state_in
    real(dp), intent(out) :: positions_out(:, :, :)
    integer, intent(in) :: beads_in
    integer :: bead

    do bead = 1, beads_in
      positions_out(:, :, bead) = state_in%position_A(:, :)
    end do
  end subroutine build_positions

  subroutine perturb_pimd_beads(positions_inout, pimd_index)
    real(dp), intent(inout) :: positions_inout(:, :, :)
    integer, intent(in) :: pimd_index
    integer :: bead
    real(dp) :: scale

    do bead = 1, size(positions_inout, 3)
      if (size(positions_inout, 3) == 1) cycle
      scale = real(bead - 1, dp) / real(size(positions_inout, 3) - 1, dp)
      positions_inout(1, pimd_index, bead) = positions_inout(1, pimd_index, bead) + electron_shift_x * scale
      positions_inout(2, pimd_index, bead) = positions_inout(2, pimd_index, bead) + electron_shift_y * scale
    end do
  end subroutine perturb_pimd_beads

  subroutine build_default_groups(groups_out, beads_in, temperature_in)
    type(particle_group_spec), allocatable, intent(out) :: groups_out(:)
    integer, intent(in) :: beads_in
    real(dp), intent(in) :: temperature_in

    allocate(groups_out(2))
    groups_out(1)%name = "graphene_carbons"
    groups_out(1)%selector = "element == C"
    groups_out(1)%dynamics = dynamics_classical
    groups_out(1)%mobile = .false.
    groups_out(1)%beads = 1
    groups_out(1)%thermostat = "none"
    groups_out(1)%temperature_K = 0.0_dp

    groups_out(2)%name = "excess_electron"
    groups_out(2)%selector = "element == El and id == 61"
    groups_out(2)%dynamics = dynamics_pimd
    groups_out(2)%mobile = .true.
    groups_out(2)%beads = beads_in
    groups_out(2)%thermostat = "none"
    groups_out(2)%temperature_K = temperature_in
  end subroutine build_default_groups

  integer function first_particle_with_dynamics(assignment_in, dynamics) result(index_value)
    type(particle_group_assignment), intent(in) :: assignment_in
    character(len=*), intent(in) :: dynamics
    integer :: i

    index_value = 0
    do i = 1, size(assignment_in%dynamics)
      if (trim(assignment_in%dynamics(i)) == trim(dynamics)) then
        index_value = i
        return
      end if
    end do
    error stop "first_particle_with_dynamics: particle not found"
  end function first_particle_with_dynamics

  real(dp) function select_mixed_temperature(cfg_in) result(temperature_out)
    type(input_config), intent(in) :: cfg_in
    integer :: i

    temperature_out = 300.0_dp
    do i = 1, size(cfg_in%groups)
      if (trim(cfg_in%groups(i)%dynamics) == dynamics_pimd) then
        temperature_out = cfg_in%groups(i)%temperature_K
        return
      end if
    end do
  end function select_mixed_temperature

  subroutine bead_positions_to_nm(positions_in, electron_index, rp_in, q_nm_out)
    real(dp), intent(in) :: positions_in(:, :, :)
    integer, intent(in) :: electron_index
    type(ring_polymer_state), intent(in) :: rp_in
    real(dp), intent(out) :: q_nm_out(:, :)
    real(dp) :: q_bead(size(q_nm_out, 2))
    integer :: coord

    do coord = 1, size(q_nm_out, 1)
      q_bead = positions_in(coord, electron_index, :)
      call bead_to_nm(rp_in%cmat, q_bead, q_nm_out(coord, :))
    end do
  end subroutine bead_positions_to_nm

  subroutine write_step_outputs(step_in)
    integer, intent(in) :: step_in
    real(dp) :: q_bead_local(3, beads)
    real(dp) :: free_ring_eV
    integer :: coord

    do coord = 1, 3
      call nm_to_bead(rp%cmat, q_nm(coord, :), q_bead_local(coord, :))
    end do
    centroid_A = sum(q_bead_local, dim=2) / real(beads, dp)
    bead_spread_A = sqrt(sum((q_bead_local(1, :) - centroid_A(1))**2 + &
        (q_bead_local(2, :) - centroid_A(2))**2 + (q_bead_local(3, :) - centroid_A(3))**2) / real(beads, dp))

    free_ring_eV = compute_free_ring_energy_nm(q_nm, v_nm, masses, rp%omega_k)
    spring_eV = compute_spring_energy_nm(q_nm, masses, rp%omega_k)
    kinetic_eV = free_ring_eV - spring_eV
    physical_potential_eV = summary%pimd_potential_eV
    max_force = maxval(abs(summary%scaled_forces_eV_per_A(:, pimd_index, :)))
    charge_diff = max_charge_difference_across_beads(charge_memory)
    conserved_like = nvt_conserved_like(kinetic_eV, physical_potential_eV, spring_eV, 0.0_dp)

    write(energy_unit,'(i8,6(1x,es22.14))') step_in, kinetic_eV, physical_potential_eV, spring_eV, &
        conserved_like, charge_diff, max_force
    write(spread_unit,'(i8,4(1x,es22.14))') step_in, centroid_A(1), centroid_A(2), centroid_A(3), bead_spread_A
    write(heat_unit,'(i8,1x,es22.14)') step_in, 0.0_dp
    call write_centroid_xyz_frame(xyz_unit, step_in, centroid_A)
    call write_charge_by_bead_summary(join_path(output_dir, "charge_by_bead.out"), charge_memory%charges_by_bead)
    call write_charge_centroid_summary(join_path(output_dir, "charge_centroid.out"), &
        charge_memory%charges_by_bead, solver_guess)
    call write_electron_localization_summary(join_path(output_dir, "electron_localization.out"), step_in, &
        charge_diff, 0.0_dp, max_force, physical_potential_eV)
  end subroutine write_step_outputs

  subroutine write_final_restart_and_compat()
    type(mixed_baoab_restart_state) :: restart_out_state
    real(dp) :: centroid_positions(ndim, nparticles)
    real(dp) :: centroid_charges(nparticles)
    real(dp) :: final_charge_sum, final_max_force

    restart_out_state%step = step
    restart_out_state%ndim = ndim
    restart_out_state%nparticles = nparticles
    restart_out_state%beads = beads
    restart_out_state%charge_policy = trim(policy_arg)
    allocate(restart_out_state%positions_A(ndim, nparticles, beads))
    allocate(restart_out_state%q_nm(ndim, beads))
    allocate(restart_out_state%v_nm(ndim, beads))
    allocate(restart_out_state%charges_by_bead(nparticles, beads))
    allocate(restart_out_state%solver_guess_by_bead(nparticles, beads))
    restart_out_state%positions_A = positions_A
    restart_out_state%q_nm = q_nm
    restart_out_state%v_nm = v_nm
    restart_out_state%charges_by_bead = charge_memory%charges_by_bead
    restart_out_state%solver_guess_by_bead = solver_guess

    if (len_trim(restart_out) == 0) restart_out = join_path(output_dir, "restart.dat")
    if (trim(restart_out) /= "none") then
      call write_mixed_baoab_restart(trim(restart_out), restart_out_state)
    end if

    centroid_positions = sum(positions_A, dim=3) / real(beads, dp)
    centroid_charges = sum(charge_memory%charges_by_bead, dim=2) / real(beads, dp)
    final_charge_sum = sum(centroid_charges)
    final_max_force = maxval(abs(summary%scaled_forces_eV_per_A(:, pimd_index, :)))
    call write_reaxff_compat_snapshot(output_dir, step, centroid_positions, centroid_charges, &
        summary%pimd_potential_eV, final_max_force, final_charge_sum, beads)
  end subroutine write_final_restart_and_compat

  subroutine write_centroid_xyz_frame(unit, step_in, centroid)
    integer, intent(in) :: unit, step_in
    real(dp), intent(in) :: centroid(:)

    write(unit,'(i0)') 1
    write(unit,'(a,i0)') "step=", step_in
    write(unit,'(a,3(1x,es22.14))') "El", centroid(1), centroid(2), centroid(3)
  end subroutine write_centroid_xyz_frame

end program ereaxff_mixed_baoab_main
