program ereaxff_mixed_baoab_main
  use ereaxff_kinds, only: dp, i8
  use ereaxff_system_state, only: system_state, init_system_state, set_particle, set_cell_from_crystx, &
      set_orthorhombic_cell, cell_is_defined, cell_lengths_A, wrap_position_into_cell
  use ereaxff_legacy_geo_loader, only: load_legacy_geo_file
  use ereaxff_particle_groups, only: particle_group_spec, particle_group_assignment, &
      assign_particle_groups, count_particles_with_dynamics, dynamics_classical, dynamics_pimd
  use ereaxff_input_config, only: input_config, read_input_config
  use ereaxff_rng, only: rng_state, rng_seed, rng_normal
  use ereaxff_ring_polymer, only: ring_polymer_state, init_ring_polymer
  use ereaxff_normal_modes, only: bead_to_nm, nm_to_bead
  use ereaxff_charge_coupling_policy, only: charge_memory_state, init_charge_memory, &
      charge_policy_bead_local_acks2, charge_policy_is_valid, max_charge_difference_across_beads
  use ereaxff_oracle_adapter, only: oracle_adapter_state, oracle_adapter_init
  use ereaxff_bead_force_manager, only: oracle_bead_summary, evaluate_oracle_beads
  use ereaxff_mixed_baoab_restart, only: mixed_baoab_restart_state, write_mixed_baoab_restart, &
      read_mixed_baoab_restart
  use ereaxff_mixed_restart, only: write_charge_by_bead_summary, write_charge_centroid_summary, &
      write_electron_localization_summary
  use ereaxff_reaxff_compat_output, only: write_reaxff_compat_snapshot
  use ereaxff_output_manager, only: join_path
  use ereaxff_spring_energy, only: compute_free_ring_energy_nm, compute_spring_energy_nm
  use ereaxff_diagnostics, only: nvt_conserved_like
  use ereaxff_constants, only: boltzmann_eV_per_K
  use ereaxff_units, only: amu_angstrom2_per_fs2_to_eV
  use ereaxff_thermostat_base, only: kinetic_energy_1d_eV, maxwell_velocity_sigma
  use ereaxff_mixed_baoab_integrator, only: mixed_baoab_quantum_step
  use ereaxff_thermostat_pile_l, only: pile_l_thermostat_state, pile_l_init
  implicit none

  integer, parameter :: ndim = 3
  real(dp), parameter :: electron_mass_amu = 0.00054858_dp
  real(dp), parameter :: electron_shift_x = 0.30_dp
  real(dp), parameter :: electron_shift_y = -0.12_dp
  real(dp), parameter :: electron_shift_z = 0.08_dp
  type(system_state) :: state
  type(input_config) :: cfg
  type(particle_group_spec), allocatable :: groups(:)
  type(particle_group_assignment) :: assignment
  type(ring_polymer_state), allocatable :: rp(:)
  type(oracle_adapter_state) :: adapter
  type(charge_memory_state) :: charge_memory
  type(oracle_bead_summary) :: summary
  type(mixed_baoab_restart_state) :: restart_state
  real(dp), allocatable :: positions_A(:, :, :)
  real(dp), allocatable :: q_nm(:, :, :)
  real(dp), allocatable :: v_nm(:, :, :)
  real(dp), allocatable :: classical_velocities(:, :)
  real(dp), allocatable :: solver_guess(:, :)
  real(dp), allocatable :: quantum_masses_amu(:)
  real(dp), allocatable :: quantum_heat_eV(:)
  integer, allocatable :: quantum_bead_counts(:)
  integer, allocatable :: quantum_indices(:)
  type(pile_l_thermostat_state), allocatable :: quantum_thermostats(:)
  real(dp) :: masses(3), kinetic_eV, spring_eV, physical_potential_eV, conserved_like
  real(dp) :: classical_heat_eV, thermostat_heat_eV
  real(dp) :: centroid_A(3), bead_spread_A, max_force, charge_diff
  character(len=512) :: arg1, arg2, beads_arg, steps_arg, policy_arg, restart_in, restart_out
  character(len=512) :: config_path, output_dir, adapter_work_root
  character(len=512) :: v1_root, v2_root
  character(len=512) :: geometry_path
  character(len=1100) :: mkdir_command
  integer :: beads, steps, ios, step, start_step, energy_unit, spread_unit, heat_unit, xyz_unit
  integer :: carbon_temp_unit, electron_temp_unit, electron_mode_temp_unit
  integer :: nparticles
  integer :: pimd_count, output_stride, nquantum, q, primary_quantum, primary_beads
  real(dp) :: temperature_K, dt_fs, pimd_mass_amu
  type(rng_state) :: rng
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
  geometry_path = resolve_geometry_path(cfg)
  if (config_mode) then
    call rng_seed(rng, int(cfg%simulation%seed, kind=8))
  else
    call rng_seed(rng, 1_8)
  end if

  mkdir_command = 'mkdir -p "' // trim(output_dir) // '"'
  call execute_command_line(trim(mkdir_command))

  write(*,'(a)') "PI-eReaxFF mixed eReaxFF BAOAB driver"
  write(*,'(a)') "build/version: development"
  write(*,'(a)') "ring-polymer convention: standard_ipi"
  write(*,'(a,i0)') "P: ", beads
  write(*,'(a,es12.4)') "T_K: ", temperature_K
  write(*,'(a)') "thermostat mode: mixed_group"
  write(*,'(a)') "integrator mode: mixed_baoab"
  write(*,'(a)') "charge coupling mode: " // trim(policy_arg)
  write(*,'(a)') "force provider: eReaxFF adapter"
  write(*,'(a)') "geometry template: " // trim(cfg%geometry%template)
  if (config_mode) then
    write(*,'(a)') "input boundary: config-driven legacy V1/V2 fixtures"
    write(*,'(a)') "input file: " // trim(config_path)
    if (len_trim(geometry_path) > 0) write(*,'(a)') "geometry source: " // trim(geometry_path)
  else
    write(*,'(a)') "input boundary: legacy V1/V2 fixtures copied by adapter"
  end if
  write(*,'(a)') "output directory: " // trim(output_dir)

  if (trim(cfg%geometry%template) == "user_geo_structure") then
    call load_legacy_geo_file(trim(geometry_path), state)
  else if (config_mode) then
    call init_system_state(state, ndim, 61)
    call build_phase4_fixture(state, cfg%geometry%electron_position_A)
  else
    call init_system_state(state, ndim, 61)
    call build_phase4_fixture(state, [7.38300_dp, 5.68433_dp, 10.001_dp])
  end if
  call apply_configured_cell(state, cfg%geometry%crystx, cfg%geometry%periodic)
  if (cell_is_defined(state)) then
    write(*,'(a,l1,1x,a,l1,1x,a,6(1x,f12.5))') "cell_defined=", cell_is_defined(state), &
        "periodic=", state%cell%periodic, "crystx:", state%cell%crystx
    write(*,'(a,3(1x,f12.5))') "cell_lengths_A:", cell_lengths_A(state)
  else
    write(*,'(a)') "cell_defined=F"
  end if
  nparticles = state%nparticles

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
  do q = 1, state%nparticles
    if (assignment%mass_amu(q) > 0.0_dp) state%mass_amu(q) = assignment%mass_amu(q)
  end do
  pimd_count = count_particles_with_dynamics(assignment, dynamics_pimd)
  if (pimd_count < 1) error stop "ereaxff_mixed_baoab_main: at least one PIMD particle is required"
  quantum_indices = collect_particles_with_dynamics(assignment, dynamics_pimd)
  nquantum = size(quantum_indices)
  if (trim(cfg%geometry%template) == "phase4_fixed_graphene" .and. nquantum > 1) then
    error stop "ereaxff_mixed_baoab_main: phase4_fixed_graphene supports only one PIMD particle; use user_geo_structure"
  end if
  if (len_trim(restart_in) > 0 .and. trim(restart_in) /= "none") then
    if (restart_state%nquantum > 0 .and. restart_state%nquantum /= nquantum) then
      error stop "ereaxff_mixed_baoab_main: restart quantum count mismatch"
    end if
  end if
  allocate(quantum_bead_counts(nquantum))
  do q = 1, nquantum
    quantum_bead_counts(q) = assignment%beads(quantum_indices(q))
  end do
  if (len_trim(restart_in) > 0 .and. trim(restart_in) /= "none") then
    if (allocated(restart_state%quantum_beads)) then
      if (any(restart_state%quantum_beads /= quantum_bead_counts)) then
        error stop "ereaxff_mixed_baoab_main: restart quantum bead-count mismatch"
      end if
    end if
  end if
  beads = maxval(quantum_bead_counts)
  if (len_trim(restart_in) > 0 .and. trim(restart_in) /= "none") then
    if (restart_state%beads /= beads) error stop "ereaxff_mixed_baoab_main: restart max bead-count mismatch"
  end if
  call print_group_thermostats(groups)

  allocate(rp(nquantum))
  do q = 1, nquantum
    call init_ring_polymer(rp(q), quantum_bead_counts(q), &
        get_pimd_group_temperature(groups, assignment, quantum_indices(q)), "standard_ipi")
  end do
  allocate(positions_A(ndim, nparticles, beads))
  allocate(q_nm(ndim, nquantum, beads))
  allocate(v_nm(ndim, nquantum, beads))
  allocate(classical_velocities(ndim, nparticles))
  allocate(solver_guess(nparticles, beads))
  allocate(quantum_masses_amu(nquantum))
  allocate(quantum_heat_eV(nquantum))
  allocate(quantum_thermostats(nquantum))
  classical_velocities = 0.0_dp
  classical_heat_eV = 0.0_dp
  quantum_heat_eV = 0.0_dp

  if (len_trim(restart_in) > 0 .and. trim(restart_in) /= "none") then
    positions_A = restart_state%positions_A
    q_nm = restart_state%q_nm
    v_nm = restart_state%v_nm
    solver_guess = restart_state%solver_guess_by_bead
    if (allocated(restart_state%classical_velocities)) classical_velocities = restart_state%classical_velocities
    classical_heat_eV = restart_state%classical_heat_eV
    if (allocated(restart_state%quantum_heat_eV)) quantum_heat_eV = restart_state%quantum_heat_eV
    if (restart_state%rng_state_i8 /= 0_i8) rng%state = restart_state%rng_state_i8
    if (allocated(restart_state%quantum_indices)) quantum_indices = restart_state%quantum_indices
    if (allocated(restart_state%quantum_masses_amu)) quantum_masses_amu = restart_state%quantum_masses_amu
  else
    if (config_mode) then
      call build_positions(state, positions_A, beads)
      if (trim(cfg%geometry%template) == "phase4_fixed_graphene") then
        do q = 1, nquantum
          call set_quantum_position(positions_A, quantum_indices(q), cfg%geometry%electron_position_A + &
              quantum_initial_offset(real(q, dp)))
        end do
      end if
    else
      call build_positions(state, positions_A, beads)
      do q = 1, nquantum
        call perturb_pimd_beads(positions_A, quantum_indices(q), real(q - 1, dp))
      end do
    end if
    q_nm = 0.0_dp
    v_nm = 0.0_dp
    do q = 1, nquantum
      call bead_positions_to_nm(positions_A, quantum_indices(q), rp(q), q_nm(:, q, 1:quantum_bead_counts(q)))
    end do
    call initialize_group_velocities(assignment, groups, quantum_indices, quantum_bead_counts, state, &
        classical_velocities, v_nm, rng)
    solver_guess = 0.0_dp
  end if
  if (.not. allocated(quantum_indices)) quantum_indices = collect_particles_with_dynamics(assignment, dynamics_pimd)
  if (.not. allocated(quantum_masses_amu)) then
    do q = 1, nquantum
      quantum_masses_amu(q) = state%mass_amu(quantum_indices(q))
    end do
  end if
  write(*,'(a,3(1x,es12.4))') "electron_start_A:", state%position_A(:, quantum_indices(1))

  call init_charge_memory(charge_memory, nparticles, beads, policy_arg)
  if (config_mode) then
    adapter_work_root = join_path(output_dir, "adapter_work")
    if (len_trim(cfg%force_provider%work_root) > 0) adapter_work_root = trim(cfg%force_provider%work_root)
    call oracle_adapter_init(adapter, state, trim(v1_root), trim(v2_root), trim(adapter_work_root), &
        control_path=join_path(v1_root, trim(cfg%force_provider%control)), &
        ffield_path=join_path(v1_root, trim(cfg%force_provider%ffield)), &
        iopt_path=join_path(v1_root, trim(cfg%force_provider%iopt)), &
        smoke_binary=trim(cfg%force_provider%smoke_binary), repeat_binary=trim(cfg%force_provider%repeat_binary), &
        make_geo_script=trim(cfg%force_provider%make_geo_script), timeout_runner=trim(cfg%force_provider%timeout_runner), &
        geo_path=trim(geometry_path))
  else
    call oracle_adapter_init(adapter, state, trim(v1_root), trim(v2_root), join_path(output_dir, "adapter_work"))
  end if
  primary_quantum = quantum_indices(1)
  primary_beads = quantum_bead_counts(1)
  do q = 1, nquantum
    quantum_masses_amu(q) = state%mass_amu(quantum_indices(q))
    call init_quantum_thermostat(quantum_thermostats(q), groups(assignment%group_index(quantum_indices(q)))%thermostat, &
        get_pimd_group_temperature(groups, assignment, quantum_indices(q)), &
        get_pimd_group_centroid_gamma(groups, assignment, quantum_indices(q)), dt_fs, rp(q)%omega_k)
    if (len_trim(restart_in) > 0 .and. trim(restart_in) /= "none") quantum_thermostats(q)%heat_eV = quantum_heat_eV(q)
  end do
  pimd_mass_amu = quantum_masses_amu(1)
  masses = [pimd_mass_amu, pimd_mass_amu, pimd_mass_amu]

  call evaluate_oracle_beads(adapter, positions_A, charge_memory, solver_guess, policy_arg, summary)
  output_stride = 1
  if (config_mode) output_stride = max(1, cfg%simulation%output_stride)

  open(newunit=energy_unit, file=join_path(output_dir, "pimd_energy.out"), status="replace", action="write")
  open(newunit=spread_unit, file=join_path(output_dir, "bead_spread.out"), status="replace", action="write")
  open(newunit=heat_unit, file=join_path(output_dir, "thermostat_heat.out"), status="replace", action="write")
  open(newunit=xyz_unit, file=join_path(output_dir, "centroid.xyz"), status="replace", action="write")
  open(newunit=carbon_temp_unit, file=join_path(output_dir, "carbon_temperature.out"), status="replace", action="write")
  open(newunit=electron_temp_unit, file=join_path(output_dir, "electron_temperature.out"), status="replace", action="write")
  open(newunit=electron_mode_temp_unit, &
      file=join_path(output_dir, "electron_mode_temperature.out"), status="replace", action="write")

  write(energy_unit,'(a)') "# step kinetic_eV physical_potential_eV spring_eV conserved_like_eV charge_diff max_force"
  write(spread_unit,'(a)') "# step centroid_x centroid_y centroid_z bead_rms_spread_A"
  write(heat_unit,'(a)') "# step thermostat_heat_eV classical_heat_eV quantum_heat_eV"
  write(carbon_temp_unit,'(a)') "# step carbon_temperature_K carbon_ke_eV"
  write(electron_temp_unit,'(a)') "# step electron_temperature_K electron_ke_eV"
  write(electron_mode_temp_unit,'(a)') "# step mode mode_temperature_K mode_ke_eV"

  if (mod(start_step, output_stride) == 0 .or. start_step == steps) call write_step_outputs(start_step)
  step = start_step
  do while (step < steps)
    call mixed_group_baoab_step(positions_A, classical_velocities, q_nm, v_nm, assignment, groups, quantum_indices, &
        quantum_bead_counts, state, rp, dt_fs, adapter, charge_memory, solver_guess, policy_arg, summary, rng, &
        quantum_thermostats, quantum_heat_eV, classical_heat_eV)
    step = step + 1
    if (mod(step, output_stride) == 0 .or. step == steps) call write_step_outputs(step)
  end do

  close(energy_unit)
  close(spread_unit)
  close(heat_unit)
  close(xyz_unit)
  close(carbon_temp_unit)
  close(electron_temp_unit)
  close(electron_mode_temp_unit)

  call write_final_restart_and_compat()

  write(*,'(a,i0)') "classical_count = ", count_particles_with_dynamics(assignment, dynamics_classical)
  write(*,'(a,i0)') "pimd_count = ", count_particles_with_dynamics(assignment, dynamics_pimd)
  write(*,'(a,es12.4)') "classical_thermostat_heat_eV = ", classical_heat_eV
  write(*,'(a,es12.4)') "quantum_thermostat_heat_eV = ", sum(quantum_heat_eV)

contains

  subroutine build_phase4_fixture(state_out, electron_position_A)
    type(system_state), intent(inout) :: state_out
    real(dp), intent(in) :: electron_position_A(3)
    integer :: ix, iy, unit_index, idx
    real(dp), parameter :: a_cc = 1.4210_dp
    real(dp), parameter :: z_c = 10.0_dp
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
    call set_particle(state_out, 61, 61, "El", electron_mass_amu, electron_position_A)
    call set_orthorhombic_cell(state_out, [12.5_dp, 12.5_dp, 20.0_dp], .false.)
  end subroutine build_phase4_fixture

  subroutine apply_configured_cell(state_inout, crystx, periodic)
    type(system_state), intent(inout) :: state_inout
    real(dp), intent(in) :: crystx(6)
    logical, intent(in) :: periodic

    if (any(abs(crystx) > 0.0_dp)) call set_cell_from_crystx(state_inout, crystx, periodic)
  end subroutine apply_configured_cell

  pure function resolve_geometry_path(cfg_in) result(path_out)
    type(input_config), intent(in) :: cfg_in
    character(len=512) :: path_out

    path_out = ""
    if (len_trim(cfg_in%geometry%geo_path) > 0) then
      path_out = trim(cfg_in%geometry%geo_path)
    else if (len_trim(cfg_in%geometry%structure_path) > 0) then
      path_out = trim(cfg_in%geometry%structure_path)
    end if
  end function resolve_geometry_path

  subroutine build_positions(state_in, positions_out, beads_in)
    type(system_state), intent(in) :: state_in
    real(dp), intent(out) :: positions_out(:, :, :)
    integer, intent(in) :: beads_in
    integer :: bead

    do bead = 1, beads_in
      positions_out(:, :, bead) = state_in%position_A(:, :)
    end do
  end subroutine build_positions

  subroutine set_quantum_position(positions_inout, pimd_index, electron_position_A)
    real(dp), intent(inout) :: positions_inout(:, :, :)
    integer, intent(in) :: pimd_index
    real(dp), intent(in) :: electron_position_A(3)
    integer :: bead

    do bead = 1, size(positions_inout, 3)
      positions_inout(:, pimd_index, bead) = electron_position_A
    end do
  end subroutine set_quantum_position

  subroutine perturb_pimd_beads(positions_inout, pimd_index, seed_offset)
    real(dp), intent(inout) :: positions_inout(:, :, :)
    integer, intent(in) :: pimd_index
    real(dp), intent(in) :: seed_offset
    integer :: bead
    real(dp) :: scale

    do bead = 1, size(positions_inout, 3)
      if (size(positions_inout, 3) == 1) cycle
      scale = real(bead - 1, dp) / real(size(positions_inout, 3) - 1, dp)
      positions_inout(1, pimd_index, bead) = positions_inout(1, pimd_index, bead) + (electron_shift_x + seed_offset) * scale
      positions_inout(2, pimd_index, bead) = positions_inout(2, pimd_index, bead) + (electron_shift_y - seed_offset) * scale
      positions_inout(3, pimd_index, bead) = positions_inout(3, pimd_index, bead) + electron_shift_z * scale
    end do
  end subroutine perturb_pimd_beads

  pure function quantum_initial_offset(q_index) result(offset_A)
    real(dp), intent(in) :: q_index
    real(dp) :: offset_A(3)

    offset_A = [0.08_dp * (q_index - 1.0_dp), -0.03_dp * (q_index - 1.0_dp), 0.02_dp * (q_index - 1.0_dp)]
  end function quantum_initial_offset

  subroutine build_default_groups(groups_out, beads_in, temperature_in)
    type(particle_group_spec), allocatable, intent(out) :: groups_out(:)
    integer, intent(in) :: beads_in
    real(dp), intent(in) :: temperature_in

    allocate(groups_out(2))
    groups_out(1)%name = "graphene_carbons"
    groups_out(1)%selector = "element == C"
    groups_out(1)%dynamics = dynamics_classical
    groups_out(1)%mobile = .true.
    groups_out(1)%beads = 1
    groups_out(1)%thermostat = "langevin"
    groups_out(1)%temperature_K = temperature_in
    groups_out(1)%gamma_fs_inv = 0.001_dp

    groups_out(2)%name = "excess_electron"
    groups_out(2)%selector = "element == El and id == 61"
    groups_out(2)%dynamics = dynamics_pimd
    groups_out(2)%mobile = .true.
    groups_out(2)%beads = beads_in
    groups_out(2)%thermostat = "pile_l"
    groups_out(2)%temperature_K = temperature_in
    groups_out(2)%centroid_gamma_fs_inv = 0.001_dp
    groups_out(2)%mass_amu = electron_mass_amu
  end subroutine build_default_groups

  function collect_particles_with_dynamics(assignment_in, dynamics) result(indices)
    type(particle_group_assignment), intent(in) :: assignment_in
    character(len=*), intent(in) :: dynamics
    integer, allocatable :: indices(:)
    integer :: i, count

    count = 0
    do i = 1, size(assignment_in%dynamics)
      if (trim(assignment_in%dynamics(i)) == trim(dynamics)) count = count + 1
    end do
    if (count < 1) then
      allocate(indices(0))
      return
    end if
    allocate(indices(count))
    count = 0
    do i = 1, size(assignment_in%dynamics)
      if (trim(assignment_in%dynamics(i)) == trim(dynamics)) then
        count = count + 1
        indices(count) = i
      end if
    end do
  end function collect_particles_with_dynamics

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
      q_bead = positions_in(coord, electron_index, 1:size(q_nm_out, 2))
      call bead_to_nm(rp_in%cmat, q_bead, q_nm_out(coord, :))
    end do
  end subroutine bead_positions_to_nm

  subroutine initialize_group_velocities(assignment_in, groups_in, quantum_indices, quantum_bead_counts, state_in, &
      classical_velocities_out, v_nm_out, rng_inout)
    type(particle_group_assignment), intent(in) :: assignment_in
    type(particle_group_spec), intent(in) :: groups_in(:)
    integer, intent(in) :: quantum_indices(:), quantum_bead_counts(:)
    type(system_state), intent(in) :: state_in
    real(dp), intent(inout) :: classical_velocities_out(:, :)
    real(dp), intent(inout) :: v_nm_out(:, :, :)
    type(rng_state), intent(inout) :: rng_inout
    integer :: i, q, coord, mode, group_index
    real(dp) :: sigma, kBT_eV

    classical_velocities_out = 0.0_dp
    do i = 1, state_in%nparticles
      if (any(quantum_indices == i)) cycle
      if (.not. assignment_in%mobile(i)) cycle
      if (trim(assignment_in%dynamics(i)) /= dynamics_classical) cycle
      group_index = assignment_in%group_index(i)
      kBT_eV = boltzmann_eV_per_K * groups_in(group_index)%temperature_K
      sigma = maxwell_velocity_sigma(kBT_eV, state_in%mass_amu(i))
      do coord = 1, size(classical_velocities_out, 1)
        classical_velocities_out(coord, i) = sigma * rng_normal(rng_inout)
      end do
    end do

    v_nm_out = 0.0_dp
    do q = 1, size(quantum_indices)
      group_index = assignment_in%group_index(quantum_indices(q))
      kBT_eV = boltzmann_eV_per_K * groups_in(group_index)%temperature_K
      sigma = maxwell_velocity_sigma(kBT_eV, state_in%mass_amu(quantum_indices(q)))
      do mode = 1, quantum_bead_counts(q)
        do coord = 1, size(v_nm_out, 1)
          v_nm_out(coord, q, mode) = sigma * rng_normal(rng_inout)
        end do
      end do
    end do
  end subroutine initialize_group_velocities

  subroutine write_step_outputs(step_in)
    integer, intent(in) :: step_in
    real(dp) :: q_bead_local(3, primary_beads)
    real(dp) :: free_ring_eV
    integer :: coord

    do coord = 1, 3
      call nm_to_bead(rp(1)%cmat, q_nm(coord, 1, 1:primary_beads), q_bead_local(coord, :))
    end do
    centroid_A = sum(q_bead_local, dim=2) / real(primary_beads, dp)
    bead_spread_A = sqrt(sum((q_bead_local(1, :) - centroid_A(1))**2 + &
        (q_bead_local(2, :) - centroid_A(2))**2 + (q_bead_local(3, :) - centroid_A(3))**2) / real(primary_beads, dp))

    free_ring_eV = compute_free_ring_energy_nm(q_nm(:, 1, 1:primary_beads), v_nm(:, 1, 1:primary_beads), &
        masses, rp(1)%omega_k)
    spring_eV = compute_spring_energy_nm(q_nm(:, 1, 1:primary_beads), masses, rp(1)%omega_k)
    kinetic_eV = free_ring_eV - spring_eV
    physical_potential_eV = summary%pimd_potential_eV
    max_force = maxval(abs(summary%scaled_forces_eV_per_A(:, primary_quantum, :)))
    charge_diff = max_charge_difference_across_beads(charge_memory)
    call write_temperature_diagnostics(step_in, carbon_temp_unit, electron_temp_unit, electron_mode_temp_unit)
    thermostat_heat_eV = classical_heat_eV + sum(quantum_heat_eV)
    conserved_like = nvt_conserved_like(kinetic_eV, physical_potential_eV, spring_eV, thermostat_heat_eV)

    write(energy_unit,'(i8,6(1x,es22.14))') step_in, kinetic_eV, physical_potential_eV, spring_eV, &
        conserved_like, charge_diff, max_force
    write(spread_unit,'(i8,4(1x,es22.14))') step_in, centroid_A(1), centroid_A(2), centroid_A(3), bead_spread_A
    write(heat_unit,'(i8,3(1x,es22.14))') step_in, thermostat_heat_eV, classical_heat_eV, sum(quantum_heat_eV)
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
    restart_out_state%nquantum = nquantum
    restart_out_state%beads = beads
    restart_out_state%charge_policy = trim(policy_arg)
    restart_out_state%classical_heat_eV = classical_heat_eV
    restart_out_state%rng_state_i8 = rng%state
    allocate(restart_out_state%positions_A(ndim, nparticles, beads))
    allocate(restart_out_state%classical_velocities(ndim, nparticles))
    allocate(restart_out_state%quantum_indices(nquantum))
    allocate(restart_out_state%quantum_beads(nquantum))
    allocate(restart_out_state%quantum_masses_amu(nquantum))
    allocate(restart_out_state%quantum_heat_eV(nquantum))
    allocate(restart_out_state%q_nm(ndim, nquantum, beads))
    allocate(restart_out_state%v_nm(ndim, nquantum, beads))
    allocate(restart_out_state%charges_by_bead(nparticles, beads))
    allocate(restart_out_state%solver_guess_by_bead(nparticles, beads))
    restart_out_state%positions_A = positions_A
    restart_out_state%classical_velocities = classical_velocities
    restart_out_state%quantum_indices = quantum_indices
    restart_out_state%quantum_beads = quantum_bead_counts
    restart_out_state%quantum_masses_amu = quantum_masses_amu
    restart_out_state%quantum_heat_eV = quantum_heat_eV
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
    final_max_force = maxval(abs(summary%scaled_forces_eV_per_A(:, primary_quantum, :)))
    call write_reaxff_compat_snapshot(output_dir, step, centroid_positions, centroid_charges, &
        summary%pimd_potential_eV, final_max_force, final_charge_sum, beads, state)
  end subroutine write_final_restart_and_compat

  subroutine write_centroid_xyz_frame(unit, step_in, centroid)
    integer, intent(in) :: unit, step_in
    real(dp), intent(in) :: centroid(:)

    write(unit,'(i0)') 1
    write(unit,'(a,i0)') "step=", step_in
    write(unit,'(a,3(1x,es22.14))') "El", centroid(1), centroid(2), centroid(3)
  end subroutine write_centroid_xyz_frame

  subroutine write_temperature_diagnostics(step_in, carbon_unit, electron_unit, mode_unit)
    integer, intent(in) :: step_in, carbon_unit, electron_unit, mode_unit
    real(dp) :: carbon_ke_eV, carbon_temp_K
    real(dp) :: electron_ke_eV, electron_temp_K, mode_ke_eV, mode_temp_K
    integer :: i, mode, carbon_count

    carbon_ke_eV = 0.0_dp
    carbon_count = 0
    do i = 1, nparticles
      if (any(quantum_indices == i)) cycle
      if (trim(assignment%dynamics(i)) /= dynamics_classical) cycle
      if (.not. assignment%mobile(i)) cycle
      carbon_count = carbon_count + 1
      carbon_ke_eV = carbon_ke_eV + kinetic_energy_1d_eV(classical_velocities(:, i), &
          [state%mass_amu(i), state%mass_amu(i), state%mass_amu(i)])
    end do
    if (carbon_count > 0) then
      carbon_temp_K = 2.0_dp * carbon_ke_eV / (3.0_dp * real(carbon_count, dp) * boltzmann_eV_per_K)
    else
      carbon_temp_K = 0.0_dp
    end if

    electron_ke_eV = 0.0_dp
    do mode = 1, primary_beads
      electron_ke_eV = electron_ke_eV + kinetic_energy_1d_eV(v_nm(:, 1, mode), masses)
    end do
    electron_temp_K = 2.0_dp * electron_ke_eV / (3.0_dp * real(primary_beads, dp) * boltzmann_eV_per_K)

    write(carbon_unit,'(i8,2(1x,es22.14))') step_in, carbon_temp_K, carbon_ke_eV
    write(electron_unit,'(i8,2(1x,es22.14))') step_in, electron_temp_K, electron_ke_eV
    do mode = 1, primary_beads
      mode_ke_eV = kinetic_energy_1d_eV(v_nm(:, 1, mode), masses)
      mode_temp_K = 2.0_dp * mode_ke_eV / (3.0_dp * boltzmann_eV_per_K)
      write(mode_unit,'(i8,1x,i8,2(1x,es22.14))') step_in, mode, mode_temp_K, mode_ke_eV
    end do
  end subroutine write_temperature_diagnostics

  subroutine print_group_thermostats(groups_in)
    type(particle_group_spec), intent(in) :: groups_in(:)
    integer :: i

    do i = 1, size(groups_in)
      write(*,'(a,1x,a,1x,a,1x,a,1x,a,1x,es12.4,1x,a,1x,es12.4)') &
          "group thermostat:", trim(groups_in(i)%name), &
          "dynamics=" // trim(groups_in(i)%dynamics), &
          "thermostat=" // trim(groups_in(i)%thermostat), &
          "T_K=", groups_in(i)%temperature_K, &
          "gamma_fs_inv=", merge(groups_in(i)%gamma_fs_inv, groups_in(i)%centroid_gamma_fs_inv, &
          trim(groups_in(i)%dynamics) == dynamics_classical)
    end do
  end subroutine print_group_thermostats

  subroutine init_quantum_thermostat(state_inout, thermostat_name, temperature_K, centroid_gamma_fs_inv, dt_fs, omega_k)
    type(pile_l_thermostat_state), intent(inout) :: state_inout
    character(len=*), intent(in) :: thermostat_name
    real(dp), intent(in) :: temperature_K, centroid_gamma_fs_inv, dt_fs, omega_k(:)

    if (trim(thermostat_name) == "none") then
      state_inout%temperature_K = temperature_K
      state_inout%centroid_gamma_fs_inv = 0.0_dp
      state_inout%lambda = 1.0_dp
      state_inout%dt_fs = dt_fs
      state_inout%heat_eV = 0.0_dp
      if (allocated(state_inout%gamma_mode_fs_inv)) deallocate(state_inout%gamma_mode_fs_inv)
      if (allocated(state_inout%c1)) deallocate(state_inout%c1)
      return
    end if
    if (trim(thermostat_name) == "pile_g") error stop "init_quantum_thermostat: PILE-G requested but not implemented"
    if (trim(thermostat_name) /= "pile_l") error stop "init_quantum_thermostat: unsupported thermostat for PIMD group"
    if (centroid_gamma_fs_inv < 0.0_dp) error stop "init_quantum_thermostat: centroid gamma must be non-negative"
    call pile_l_init(state_inout, temperature_K, centroid_gamma_fs_inv, 1.0_dp, dt_fs, omega_k)
  end subroutine init_quantum_thermostat

  real(dp) function get_pimd_group_temperature(groups_in, assignment_in, particle_index) result(temperature_out)
    type(particle_group_spec), intent(in) :: groups_in(:)
    type(particle_group_assignment), intent(in) :: assignment_in
    integer, intent(in) :: particle_index
    integer :: group_index

    group_index = assignment_in%group_index(particle_index)
    temperature_out = groups_in(group_index)%temperature_K
  end function get_pimd_group_temperature

  real(dp) function get_pimd_group_centroid_gamma(groups_in, assignment_in, particle_index) result(gamma_out)
    type(particle_group_spec), intent(in) :: groups_in(:)
    type(particle_group_assignment), intent(in) :: assignment_in
    integer, intent(in) :: particle_index
    integer :: group_index

    group_index = assignment_in%group_index(particle_index)
    gamma_out = groups_in(group_index)%centroid_gamma_fs_inv
  end function get_pimd_group_centroid_gamma

  subroutine mixed_group_baoab_step(positions_A, classical_velocities, q_nm, v_nm, assignment, groups_in, quantum_indices, &
      quantum_bead_counts, state_in, rp_in, dt_fs, adapter, charge_memory, solver_guess_by_bead, charge_policy, summary, rng, &
      quantum_thermostats, quantum_heat_eV, classical_heat_eV)
    real(dp), intent(inout) :: positions_A(:, :, :)
    real(dp), intent(inout) :: classical_velocities(:, :)
    real(dp), intent(inout) :: q_nm(:, :, :)
    real(dp), intent(inout) :: v_nm(:, :, :)
    type(particle_group_assignment), intent(in) :: assignment
    type(particle_group_spec), intent(in) :: groups_in(:)
    integer, intent(in) :: quantum_indices(:)
    integer, intent(in) :: quantum_bead_counts(:)
    type(system_state), intent(in) :: state_in
    type(ring_polymer_state), intent(in) :: rp_in(:)
    real(dp), intent(in) :: dt_fs
    type(oracle_adapter_state), intent(inout) :: adapter
    type(charge_memory_state), intent(inout) :: charge_memory
    real(dp), intent(inout) :: solver_guess_by_bead(:, :)
    character(len=*), intent(in) :: charge_policy
    type(oracle_bead_summary), intent(inout) :: summary
    type(rng_state), intent(inout) :: rng
    type(pile_l_thermostat_state), intent(inout) :: quantum_thermostats(:)
    real(dp), intent(inout) :: quantum_heat_eV(:)
    real(dp), intent(inout) :: classical_heat_eV
    integer :: i, g, q
    real(dp) :: force_particle(3)
    real(dp) :: mass_amu
    real(dp), allocatable :: quantum_masses(:)

    if (size(quantum_indices) /= size(quantum_thermostats)) error stop "mixed_group_baoab_step: thermostat count mismatch"
    if (size(quantum_indices) /= size(quantum_heat_eV)) error stop "mixed_group_baoab_step: heat count mismatch"
    if (size(quantum_indices) /= size(quantum_bead_counts)) error stop "mixed_group_baoab_step: bead-count mismatch"
    allocate(quantum_masses(size(quantum_indices)))
    do q = 1, size(quantum_indices)
      quantum_masses(q) = state_in%mass_amu(quantum_indices(q))
    end do

    call replicate_classical_positions(positions_A, assignment, quantum_indices)
    call evaluate_oracle_beads(adapter, positions_A, charge_memory, solver_guess_by_bead, charge_policy, summary)

    do i = 1, state_in%nparticles
      if (any(quantum_indices == i)) cycle
      if (.not. assignment%mobile(i)) cycle
      g = assignment%group_index(i)
      if (trim(assignment%thermostat(i)) == "pile_g") error stop "mixed_group_baoab_step: pile_g not implemented"
      force_particle = sum(summary%scaled_forces_eV_per_A(:, i, :), dim=2)
      mass_amu = state_in%mass_amu(i)
      classical_velocities(:, i) = classical_velocities(:, i) + 0.5_dp * dt_fs * force_particle / &
          (mass_amu * amu_angstrom2_per_fs2_to_eV)
      positions_A(:, i, 1) = positions_A(:, i, 1) + 0.5_dp * dt_fs * classical_velocities(:, i)
      call wrap_position_into_cell(state_in, positions_A(:, i, 1))
    end do

    call apply_classical_group_thermostats(groups_in, assignment, quantum_indices, state_in, classical_velocities, &
        dt_fs, rng, classical_heat_eV)

    do i = 1, state_in%nparticles
      if (any(quantum_indices == i)) cycle
      if (.not. assignment%mobile(i)) cycle
      force_particle = sum(summary%scaled_forces_eV_per_A(:, i, :), dim=2)
      mass_amu = state_in%mass_amu(i)
      positions_A(:, i, 1) = positions_A(:, i, 1) + 0.5_dp * dt_fs * classical_velocities(:, i)
      call wrap_position_into_cell(state_in, positions_A(:, i, 1))
      classical_velocities(:, i) = classical_velocities(:, i) + 0.5_dp * dt_fs * force_particle / &
          (mass_amu * amu_angstrom2_per_fs2_to_eV)
      positions_A(:, i, :) = spread(positions_A(:, i, 1), dim=2, ncopies=size(positions_A, 3))
    end do

    call mixed_baoab_quantum_step(positions_A, q_nm, v_nm, quantum_indices, quantum_masses, rp_in, &
        quantum_bead_counts, dt_fs, adapter, charge_memory, solver_guess_by_bead, charge_policy, summary, &
        quantum_thermostats, rng, state_in)
    do q = 1, size(quantum_indices)
      quantum_heat_eV(q) = quantum_thermostats(q)%heat_eV
    end do

    call replicate_classical_positions(positions_A, assignment, quantum_indices)
    call evaluate_oracle_beads(adapter, positions_A, charge_memory, solver_guess_by_bead, charge_policy, summary)

    do i = 1, state_in%nparticles
      if (any(quantum_indices == i)) cycle
      if (.not. assignment%mobile(i)) cycle
      g = assignment%group_index(i)
      force_particle = sum(summary%scaled_forces_eV_per_A(:, i, :), dim=2)
      mass_amu = state_in%mass_amu(i)
      classical_velocities(:, i) = classical_velocities(:, i) + 0.5_dp * dt_fs * force_particle / &
          (mass_amu * amu_angstrom2_per_fs2_to_eV)
    end do
    deallocate(quantum_masses)
  end subroutine mixed_group_baoab_step

  subroutine replicate_classical_positions(positions_A, assignment_in, quantum_indices)
    real(dp), intent(inout) :: positions_A(:, :, :)
    type(particle_group_assignment), intent(in) :: assignment_in
    integer, intent(in) :: quantum_indices(:)
    integer :: i, bead

    do i = 1, size(assignment_in%group_index)
      if (any(quantum_indices == i)) cycle
      if (trim(assignment_in%dynamics(i)) /= dynamics_classical) cycle
      do bead = 1, size(positions_A, 3)
        positions_A(:, i, bead) = positions_A(:, i, 1)
      end do
    end do
  end subroutine replicate_classical_positions

  subroutine apply_classical_thermostat(group, velocity_vec, mass_amu, dt_fs, rng_inout, heat_eV)
    type(particle_group_spec), intent(in) :: group
    real(dp), intent(inout) :: velocity_vec(:)
    real(dp), intent(in) :: mass_amu, dt_fs
    type(rng_state), intent(inout) :: rng_inout
    real(dp), intent(inout) :: heat_eV
    real(dp) :: v_before(size(velocity_vec)), ke_before, ke_after, c1, sigma
    integer :: coord

    if (trim(group%thermostat) == "none") return
    if (trim(group%thermostat) /= "langevin") error stop "apply_classical_thermostat: unsupported thermostat"
    if (group%temperature_K < 0.0_dp) error stop "apply_classical_thermostat: temperature must be non-negative"
    if (group%gamma_fs_inv < 0.0_dp) error stop "apply_classical_thermostat: gamma must be non-negative"
    c1 = exp(-group%gamma_fs_inv * dt_fs)
    sigma = maxwell_velocity_sigma(boltzmann_eV_per_K * group%temperature_K, mass_amu) * &
        sqrt(max(0.0_dp, 1.0_dp - c1**2))
    v_before = velocity_vec
    ke_before = kinetic_energy_1d_eV(v_before, [mass_amu, mass_amu, mass_amu])
    do coord = 1, size(velocity_vec)
      velocity_vec(coord) = c1 * velocity_vec(coord) + sigma * rng_normal(rng_inout)
    end do
    ke_after = kinetic_energy_1d_eV(velocity_vec, [mass_amu, mass_amu, mass_amu])
    heat_eV = heat_eV + (ke_before - ke_after)
  end subroutine apply_classical_thermostat

  subroutine apply_classical_group_thermostats(groups_in, assignment_in, quantum_indices, state_in, velocities, &
      dt_fs, rng_inout, heat_eV)
    type(particle_group_spec), intent(in) :: groups_in(:)
    type(particle_group_assignment), intent(in) :: assignment_in
    integer, intent(in) :: quantum_indices(:)
    type(system_state), intent(in) :: state_in
    real(dp), intent(inout) :: velocities(:, :)
    real(dp), intent(in) :: dt_fs
    type(rng_state), intent(inout) :: rng_inout
    real(dp), intent(inout) :: heat_eV
    integer :: g, i, count

    do g = 1, size(groups_in)
      select case (trim(groups_in(g)%thermostat))
      case ("none", "pile_l")
        cycle
      case ("langevin")
        do i = 1, state_in%nparticles
          if (assignment_in%group_index(i) /= g) cycle
          if (any(quantum_indices == i)) cycle
          if (.not. assignment_in%mobile(i)) cycle
          call apply_classical_thermostat(groups_in(g), velocities(:, i), state_in%mass_amu(i), dt_fs, rng_inout, heat_eV)
        end do
      case ("berendsen")
        count = count_mobile_classical_group_members(assignment_in, quantum_indices, g)
        if (count > 0) call apply_berendsen_group(groups_in(g), assignment_in, quantum_indices, state_in, velocities, &
            g, count, dt_fs, heat_eV)
      case ("pile_g")
        error stop "apply_classical_group_thermostats: pile_g not implemented"
      case default
        error stop "apply_classical_group_thermostats: unsupported thermostat"
      end select
    end do
  end subroutine apply_classical_group_thermostats

  integer function count_mobile_classical_group_members(assignment_in, quantum_indices, group_index) result(count)
    type(particle_group_assignment), intent(in) :: assignment_in
    integer, intent(in) :: quantum_indices(:), group_index
    integer :: i

    count = 0
    do i = 1, size(assignment_in%group_index)
      if (assignment_in%group_index(i) /= group_index) cycle
      if (any(quantum_indices == i)) cycle
      if (.not. assignment_in%mobile(i)) cycle
      if (trim(assignment_in%dynamics(i)) /= dynamics_classical) cycle
      count = count + 1
    end do
  end function count_mobile_classical_group_members

  subroutine apply_berendsen_group(group, assignment_in, quantum_indices, state_in, velocities, group_index, count, &
      dt_fs, heat_eV)
    type(particle_group_spec), intent(in) :: group
    type(particle_group_assignment), intent(in) :: assignment_in
    integer, intent(in) :: quantum_indices(:), group_index, count
    type(system_state), intent(in) :: state_in
    real(dp), intent(inout) :: velocities(:, :)
    real(dp), intent(in) :: dt_fs
    real(dp), intent(inout) :: heat_eV
    real(dp) :: ke_before, ke_after, temp_now_K, lambda
    integer :: i

    if (group%temperature_K < 0.0_dp) error stop "apply_berendsen_group: temperature must be non-negative"
    if (group%tau_fs <= 0.0_dp) error stop "apply_berendsen_group: tau must be positive"

    ke_before = classical_group_ke_eV(assignment_in, quantum_indices, state_in, velocities, group_index)
    if (ke_before <= 0.0_dp) return
    temp_now_K = 2.0_dp * ke_before / (3.0_dp * real(count, dp) * boltzmann_eV_per_K)
    if (temp_now_K <= 0.0_dp) return
    lambda = sqrt(max(0.0_dp, 1.0_dp + dt_fs / group%tau_fs * (group%temperature_K / temp_now_K - 1.0_dp)))
    do i = 1, state_in%nparticles
      if (assignment_in%group_index(i) /= group_index) cycle
      if (any(quantum_indices == i)) cycle
      if (.not. assignment_in%mobile(i)) cycle
      if (trim(assignment_in%dynamics(i)) /= dynamics_classical) cycle
      velocities(:, i) = lambda * velocities(:, i)
    end do
    ke_after = classical_group_ke_eV(assignment_in, quantum_indices, state_in, velocities, group_index)
    heat_eV = heat_eV + (ke_before - ke_after)
  end subroutine apply_berendsen_group

  real(dp) function classical_group_ke_eV(assignment_in, quantum_indices, state_in, velocities, group_index) result(ke_eV)
    type(particle_group_assignment), intent(in) :: assignment_in
    integer, intent(in) :: quantum_indices(:), group_index
    type(system_state), intent(in) :: state_in
    real(dp), intent(in) :: velocities(:, :)
    integer :: i

    ke_eV = 0.0_dp
    do i = 1, state_in%nparticles
      if (assignment_in%group_index(i) /= group_index) cycle
      if (any(quantum_indices == i)) cycle
      if (.not. assignment_in%mobile(i)) cycle
      if (trim(assignment_in%dynamics(i)) /= dynamics_classical) cycle
      ke_eV = ke_eV + kinetic_energy_1d_eV(velocities(:, i), &
          [state_in%mass_amu(i), state_in%mass_amu(i), state_in%mass_amu(i)])
    end do
  end function classical_group_ke_eV

end program ereaxff_mixed_baoab_main
