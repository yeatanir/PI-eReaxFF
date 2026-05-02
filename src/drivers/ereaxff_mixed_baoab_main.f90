program ereaxff_mixed_baoab_main
  use ereaxff_kinds, only: dp
  use ereaxff_system_state, only: system_state, init_system_state, set_particle
  use ereaxff_particle_groups, only: particle_group_spec, particle_group_assignment, &
      assign_particle_groups, count_particles_with_dynamics, dynamics_classical, dynamics_pimd
  use ereaxff_ring_polymer, only: ring_polymer_state, init_ring_polymer
  use ereaxff_normal_modes, only: bead_to_nm, nm_to_bead
  use ereaxff_charge_coupling_policy, only: charge_memory_state, init_charge_memory, &
      charge_policy_bead_local_acks2, charge_policy_is_valid, max_charge_difference_across_beads
  use ereaxff_oracle_adapter, only: oracle_adapter_state, oracle_adapter_init
  use ereaxff_bead_force_manager, only: oracle_bead_summary, evaluate_oracle_beads
  use ereaxff_mixed_baoab_integrator, only: mixed_baoab_electron_step
  use ereaxff_mixed_restart, only: write_charge_by_bead_summary, write_charge_centroid_summary, &
      write_electron_localization_summary
  use ereaxff_output_manager, only: join_path
  use ereaxff_spring_energy, only: compute_free_ring_energy_nm, compute_spring_energy_nm
  use ereaxff_diagnostics, only: nvt_conserved_like
  implicit none

  integer, parameter :: ndim = 3
  integer, parameter :: nparticles = 61
  integer, parameter :: default_beads = 4
  integer, parameter :: default_steps = 8
  real(dp), parameter :: temperature_K = 300.0_dp
  real(dp), parameter :: electron_mass_amu = 0.00054858_dp
  real(dp), parameter :: dt_fs = 0.1_dp
  real(dp), parameter :: electron_shift_x = 0.08_dp
  real(dp), parameter :: electron_shift_y = -0.03_dp
  type(system_state) :: state
  type(particle_group_spec) :: groups(2)
  type(particle_group_assignment) :: assignment
  type(ring_polymer_state) :: rp
  type(oracle_adapter_state) :: adapter
  type(charge_memory_state) :: charge_memory
  type(oracle_bead_summary) :: summary
  real(dp), allocatable :: positions_A(:, :, :)
  real(dp), allocatable :: q_nm(:, :)
  real(dp), allocatable :: v_nm(:, :)
  real(dp), allocatable :: solver_guess(:, :)
  real(dp) :: masses(3), kinetic_eV, spring_eV, physical_potential_eV, conserved_like
  real(dp) :: centroid_A(3), bead_spread_A, max_force, charge_diff
  character(len=512) :: output_dir, beads_arg, steps_arg, policy_arg
  character(len=512) :: v1_root, v2_root
  character(len=1100) :: mkdir_command
  integer :: beads, steps, ios, step, energy_unit, spread_unit, heat_unit, xyz_unit

  call get_command_argument(1, output_dir)
  if (len_trim(output_dir) == 0) output_dir = "ereaxff_mixed_baoab_output"
  call get_command_argument(2, beads_arg)
  beads = default_beads
  if (len_trim(beads_arg) > 0) then
    read(beads_arg, *, iostat=ios) beads
    if (ios /= 0 .or. beads < 1) error stop "ereaxff_mixed_baoab_main: invalid beads argument"
  end if
  call get_command_argument(3, steps_arg)
  steps = default_steps
  if (len_trim(steps_arg) > 0) then
    read(steps_arg, *, iostat=ios) steps
    if (ios /= 0 .or. steps < 0) error stop "ereaxff_mixed_baoab_main: invalid steps argument"
  end if
  call get_command_argument(4, policy_arg)
  if (len_trim(policy_arg) == 0) policy_arg = charge_policy_bead_local_acks2
  if (.not. charge_policy_is_valid(policy_arg)) error stop "ereaxff_mixed_baoab_main: invalid charge policy"
  call get_environment_variable("EREAXFF_V1_ROOT", v1_root)
  call get_environment_variable("EREAXFF_V2_ROOT", v2_root)
  if (len_trim(v1_root) == 0) v1_root = "/Users/yeatanir/PhD/PIMD/ereaxff_pimd/V1"
  if (len_trim(v2_root) == 0) v2_root = "/Users/yeatanir/PhD/PIMD/ereaxff_pimd/V2"

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
  write(*,'(a)') "input boundary: legacy V1/V2 fixtures copied by adapter"
  write(*,'(a)') "output directory: " // trim(output_dir)

  call init_system_state(state, ndim, nparticles)
  call build_phase4_fixture(state)
  call assign_mixed_groups(state, groups, assignment)
  call init_ring_polymer(rp, beads, temperature_K, "standard_ipi")

  allocate(positions_A(ndim, nparticles, beads))
  allocate(q_nm(ndim, beads))
  allocate(v_nm(ndim, beads))
  allocate(solver_guess(nparticles, beads))

  call build_positions(state, positions_A, beads)
  call perturb_electron_beads(positions_A)
  q_nm = 0.0_dp
  v_nm = 0.0_dp
  call bead_positions_to_nm(positions_A, 1, rp, q_nm)
  call init_charge_memory(charge_memory, nparticles, beads, policy_arg)
  solver_guess = 0.0_dp

  call oracle_adapter_init(adapter, state, trim(v1_root), trim(v2_root), join_path(output_dir, "adapter_work"))
  masses = [electron_mass_amu, electron_mass_amu, electron_mass_amu]

  call evaluate_oracle_beads(adapter, positions_A, charge_memory, solver_guess, policy_arg, summary)

  open(newunit=energy_unit, file=join_path(output_dir, "pimd_energy.out"), status="replace", action="write")
  open(newunit=spread_unit, file=join_path(output_dir, "bead_spread.out"), status="replace", action="write")
  open(newunit=heat_unit, file=join_path(output_dir, "thermostat_heat.out"), status="replace", action="write")
  open(newunit=xyz_unit, file=join_path(output_dir, "centroid.xyz"), status="replace", action="write")

  write(energy_unit,'(a)') "# step kinetic_eV physical_potential_eV spring_eV conserved_like_eV charge_diff max_force"
  write(spread_unit,'(a)') "# step centroid_x centroid_y centroid_z bead_rms_spread_A"
  write(heat_unit,'(a)') "# step thermostat_heat_eV"

  call write_step_outputs(0)
  do step = 1, steps
    call mixed_baoab_electron_step(positions_A, q_nm, v_nm, 1, electron_mass_amu, rp, dt_fs, &
        adapter, charge_memory, solver_guess, policy_arg, summary)
    call write_step_outputs(step)
  end do

  close(energy_unit)
  close(spread_unit)
  close(heat_unit)
  close(xyz_unit)

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

    call set_particle(state_out, 1, 1, "El", electron_mass_amu, [el_x, el_y, el_z])
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

  subroutine build_positions(state_in, positions_out, beads_in)
    type(system_state), intent(in) :: state_in
    real(dp), intent(out) :: positions_out(:, :, :)
    integer, intent(in) :: beads_in
    integer :: bead

    do bead = 1, beads_in
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
    groups_out(1)%mobile = .false.
    groups_out(1)%beads = 1
    groups_out(1)%thermostat = "none"
    groups_out(1)%temperature_K = 0.0_dp

    groups_out(2)%name = "excess_electron"
    groups_out(2)%selector = "element == El and id == 1"
    groups_out(2)%dynamics = dynamics_pimd
    groups_out(2)%mobile = .true.
    groups_out(2)%beads = beads
    groups_out(2)%thermostat = "none"
    groups_out(2)%temperature_K = temperature_K

    call assign_particle_groups(state_in, groups_out, assignment_out)
  end subroutine assign_mixed_groups

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
    max_force = maxval(abs(summary%scaled_forces_eV_per_A(:, 1, :)))
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

  subroutine write_centroid_xyz_frame(unit, step_in, centroid)
    integer, intent(in) :: unit, step_in
    real(dp), intent(in) :: centroid(:)

    write(unit,'(i0)') 1
    write(unit,'(a,i0)') "step=", step_in
    write(unit,'(a,3(1x,es22.14))') "El", centroid(1), centroid(2), centroid(3)
  end subroutine write_centroid_xyz_frame

end program ereaxff_mixed_baoab_main
