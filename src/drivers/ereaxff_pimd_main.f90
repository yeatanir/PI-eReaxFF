program ereaxff_pimd_main
  use ereaxff_kinds, only: dp
  use ereaxff_system_state, only: system_state, init_system_state, set_particle
  use ereaxff_particle_groups, only: particle_group_spec, particle_group_assignment, &
      assign_particle_groups, count_particles_with_dynamics, dynamics_classical, dynamics_pimd
  use ereaxff_charge_coupling_policy, only: charge_memory_state, init_charge_memory, &
      charge_policy_bead_local_acks2, charge_policy_is_valid, &
      max_charge_difference_across_beads
  use ereaxff_oracle_adapter, only: oracle_adapter_state, oracle_adapter_init, oracle_adapter_shutdown
  use ereaxff_bead_force_manager, only: oracle_bead_summary, evaluate_oracle_beads
  use ereaxff_output_manager, only: join_path
  use ereaxff_mixed_restart, only: mixed_restart_state, write_mixed_restart, read_mixed_restart, &
      write_charge_by_bead_summary, write_charge_centroid_summary, write_electron_localization_summary
  implicit none

  integer, parameter :: ndim = 3
  integer, parameter :: nparticles = 61
  integer, parameter :: default_beads = 4
  character(len=*), parameter :: default_policy = charge_policy_bead_local_acks2
  real(dp), parameter :: electron_shift_x = 0.08_dp
  real(dp), parameter :: electron_shift_y = -0.03_dp

  type(system_state) :: state
  type(particle_group_spec) :: groups(2)
  type(particle_group_assignment) :: assignment
  type(oracle_adapter_state) :: adapter
  type(charge_memory_state) :: charge_memory
  type(oracle_bead_summary) :: summary
  type(mixed_restart_state) :: restart_state
  real(dp), allocatable :: positions(:, :, :)
  real(dp), allocatable :: solver_guess(:, :)
  real(dp) :: local_charge_diff, max_force, centroid_charge_diff, energy_eV
  character(len=512) :: output_dir, beads_arg, policy_arg, restart_in, restart_out
  character(len=512) :: v1_root, v2_root
  character(len=1100) :: mkdir_command
  integer :: beads, ios

  call get_command_argument(1, output_dir)
  if (len_trim(output_dir) == 0) output_dir = "ereaxff_mixed_output"
  call get_command_argument(2, beads_arg)
  beads = default_beads
  if (len_trim(beads_arg) > 0) then
    read(beads_arg, *, iostat=ios) beads
    if (ios /= 0 .or. beads < 1) error stop "ereaxff_pimd_main: invalid beads argument"
  end if
  call get_command_argument(3, policy_arg)
  if (len_trim(policy_arg) == 0) policy_arg = default_policy
  if (.not. charge_policy_is_valid(policy_arg)) error stop "ereaxff_pimd_main: invalid charge policy"
  call get_command_argument(4, restart_in)
  call get_command_argument(5, restart_out)
  call get_environment_variable("EREAXFF_V1_ROOT", v1_root)
  call get_environment_variable("EREAXFF_V2_ROOT", v2_root)
  if (len_trim(v1_root) == 0) v1_root = "/Users/yeatanir/PhD/PIMD/ereaxff_pimd/V1"
  if (len_trim(v2_root) == 0) v2_root = "/Users/yeatanir/PhD/PIMD/ereaxff_pimd/V2"

  mkdir_command = 'mkdir -p "' // trim(output_dir) // '"'
  call execute_command_line(trim(mkdir_command))

  write(*,'(a)') "PI-eReaxFF eReaxFF PIMD smoke driver"
  write(*,'(a)') "build/version: development"
  write(*,'(a)') "ring-polymer convention: standard_ipi"
  write(*,'(a,i0)') "P: ", beads
  write(*,'(a)') "thermostat mode: none"
  write(*,'(a)') "charge coupling mode: " // trim(policy_arg)
  write(*,'(a)') "force provider: eReaxFF adapter"
  write(*,'(a)') "output directory: " // trim(output_dir)

  if (len_trim(restart_in) > 0 .and. trim(restart_in) /= "none") then
    call read_mixed_restart(trim(restart_in), restart_state)
    beads = restart_state%beads
    policy_arg = restart_state%charge_policy
    call build_phase4_fixture(state)
    call assign_mixed_groups(state, groups, assignment)
    allocate(positions(ndim, nparticles, beads))
    allocate(solver_guess(nparticles, beads))
    positions = restart_state%positions_A
    solver_guess = restart_state%solver_guess_by_bead
  else
    call build_phase4_fixture(state)
    call assign_mixed_groups(state, groups, assignment)
    allocate(positions(ndim, nparticles, beads))
    allocate(solver_guess(nparticles, beads))
    call build_positions(state, positions)
    call perturb_electron_beads(positions)
    solver_guess = 0.0_dp
  end if

  call oracle_adapter_init(adapter, state, trim(v1_root), trim(v2_root), join_path(output_dir, "adapter_work"), &
      backend_mode="external_process")
  call init_charge_memory(charge_memory, nparticles, beads, policy_arg)
  call evaluate_oracle_beads(adapter, positions, charge_memory, solver_guess, policy_arg, summary, 0, "initial_force")

  local_charge_diff = max_charge_difference_across_beads(charge_memory)
  centroid_charge_diff = 0.0_dp
  max_force = maxval(abs(summary%scaled_forces_eV_per_A))
  energy_eV = summary%pimd_potential_eV

  call write_charge_by_bead_summary(join_path(output_dir, "charge_by_bead.out"), charge_memory%charges_by_bead)
  call write_charge_centroid_summary(join_path(output_dir, "charge_centroid.out"), &
      charge_memory%charges_by_bead, solver_guess)
  call write_electron_localization_summary(join_path(output_dir, "electron_localization.out"), 0, &
      local_charge_diff, centroid_charge_diff, max_force, energy_eV)

  call write_mixed_energy_summary(join_path(output_dir, "pimd_energy.out"), energy_eV, local_charge_diff, &
      centroid_charge_diff, max_force, beads)

  restart_state%step = 0
  restart_state%ndim = ndim
  restart_state%nparticles = nparticles
  restart_state%beads = beads
  restart_state%charge_policy = trim(policy_arg)
  if (allocated(restart_state%positions_A)) deallocate(restart_state%positions_A)
  if (allocated(restart_state%charges_by_bead)) deallocate(restart_state%charges_by_bead)
  if (allocated(restart_state%solver_guess_by_bead)) deallocate(restart_state%solver_guess_by_bead)
  allocate(restart_state%positions_A(ndim, nparticles, beads))
  allocate(restart_state%charges_by_bead(nparticles, beads))
  allocate(restart_state%solver_guess_by_bead(nparticles, beads))
  restart_state%positions_A = positions
  restart_state%charges_by_bead = charge_memory%charges_by_bead
  restart_state%solver_guess_by_bead = solver_guess

  if (len_trim(restart_out) == 0) restart_out = join_path(output_dir, "restart.dat")
  if (trim(restart_out) /= "none") then
    call write_mixed_restart(trim(restart_out), restart_state)
  end if

  write(*,'(a,es24.15)') "local_charge_diff = ", local_charge_diff
  write(*,'(a,es24.15)') "energy_eV = ", energy_eV
  write(*,'(a,i0)') "classical_count = ", count_particles_with_dynamics(assignment, dynamics_classical)
  write(*,'(a,i0)') "pimd_count = ", count_particles_with_dynamics(assignment, dynamics_pimd)

  call oracle_adapter_shutdown(adapter)

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
  end subroutine build_phase4_fixture

  subroutine build_positions(state_in, positions_out)
    type(system_state), intent(in) :: state_in
    real(dp), intent(out) :: positions_out(:, :, :)
    integer :: bead

    do bead = 1, size(positions_out, 3)
      positions_out(:, :, bead) = state_in%position_A(:, :)
    end do
  end subroutine build_positions

  subroutine perturb_electron_beads(positions_inout)
    real(dp), intent(inout) :: positions_inout(:, :, :)
    integer :: bead
    real(dp) :: scale

    do bead = 1, size(positions_inout, 3)
      if (size(positions_inout, 3) == 1) cycle
      scale = real(bead - 1, dp) / real(size(positions_inout, 3) - 1, dp)
      positions_inout(1, 1, bead) = positions_inout(1, 1, bead) + electron_shift_x * scale
      positions_inout(2, 1, bead) = positions_inout(2, 1, bead) + electron_shift_y * scale
    end do
  end subroutine perturb_electron_beads

  subroutine assign_mixed_groups(state_in, groups_out, assignment_out)
    type(system_state), intent(in) :: state_in
    type(particle_group_spec), intent(out) :: groups_out(:)
    type(particle_group_assignment), intent(inout) :: assignment_out

    groups_out(1)%name = "graphene_carbons"
    groups_out(1)%selector = "element == C"
    groups_out(1)%dynamics = dynamics_classical
    groups_out(1)%mobile = .true.
    groups_out(1)%beads = 1
    groups_out(1)%thermostat = "none"
    groups_out(1)%temperature_K = 0.0_dp

    groups_out(2)%name = "excess_electron"
    groups_out(2)%selector = "element == El and id == 1"
    groups_out(2)%dynamics = dynamics_pimd
    groups_out(2)%mobile = .true.
    groups_out(2)%beads = 16
    groups_out(2)%thermostat = "none"
    groups_out(2)%temperature_K = 0.0_dp

    call assign_particle_groups(state_in, groups_out, assignment_out)
  end subroutine assign_mixed_groups

  subroutine write_mixed_energy_summary(path, energy_eV, local_charge_diff, centroid_charge_diff, max_force, beads)
    character(len=*), intent(in) :: path
    real(dp), intent(in) :: energy_eV, local_charge_diff, centroid_charge_diff, max_force
    integer, intent(in) :: beads
    integer :: unit

    open(newunit=unit, file=trim(path), status="replace", action="write")
    write(unit,'(a)') "# beads energy_eV local_charge_diff centroid_charge_diff max_force_eV_per_A"
    write(unit,'(i8,4(1x,es24.16))') beads, energy_eV, local_charge_diff, centroid_charge_diff, max_force
    close(unit)
  end subroutine write_mixed_energy_summary

end program ereaxff_pimd_main
