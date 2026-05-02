program test_mixed_oracle_manager_matrix
  use ereaxff_kinds, only: dp
  use ereaxff_system_state, only: system_state, init_system_state, set_particle
  use ereaxff_particle_groups, only: particle_group_spec, particle_group_assignment, &
      assign_particle_groups, count_particles_with_dynamics, dynamics_classical, dynamics_pimd
  use ereaxff_charge_coupling_policy, only: charge_memory_state, init_charge_memory, &
      max_charge_difference_across_beads, charge_policy_bead_local_acks2, &
      charge_policy_centroid_acks2
  use ereaxff_oracle_adapter, only: oracle_adapter_state, oracle_adapter_init
  use ereaxff_bead_force_manager, only: oracle_bead_summary, evaluate_oracle_beads
  implicit none

  integer, parameter :: ndim = 3
  integer, parameter :: nparticles = 61
  integer, parameter :: ncases = 3
  integer, parameter :: p_values(ncases) = [1, 2, 4]
  type(system_state) :: state
  type(particle_group_spec) :: groups(2)
  type(particle_group_assignment) :: assignment
  type(oracle_adapter_state) :: adapter
  type(charge_memory_state) :: local_memory, centroid_memory
  type(oracle_bead_summary) :: local_summary, centroid_summary
  real(dp), allocatable :: positions(:, :, :)
  real(dp), allocatable :: local_solver(:, :), centroid_solver(:, :)
  character(len=512) :: v1_root, v2_root, work_root
  integer :: icase, p
  real(dp) :: local_diff, centroid_diff, energy_diff

  call get_environment_variable("EREAXFF_V1_ROOT", v1_root)
  call get_environment_variable("EREAXFF_V2_ROOT", v2_root)
  call get_environment_variable("EREAXFF_ADAPTER_WORK_ROOT", work_root)
  if (len_trim(v1_root) == 0) v1_root = "/Users/yeatanir/PhD/PIMD/ereaxff_pimd/V1"
  if (len_trim(v2_root) == 0) v2_root = "/Users/yeatanir/PhD/PIMD/ereaxff_pimd/V2"
  if (len_trim(work_root) == 0) work_root = "build/test_adapter_work"

  call build_phase4_template(state)
  call assign_mixed_groups(state, groups, assignment)
  call assert_int("classical count", count_particles_with_dynamics(assignment, dynamics_classical), 60)
  call assert_int("pimd count", count_particles_with_dynamics(assignment, dynamics_pimd), 1)
  call assert_int("electron beads", assignment%beads(1), 16)

  call oracle_adapter_init(adapter, state, trim(v1_root), trim(v2_root), trim(work_root))

  do icase = 1, ncases
    p = p_values(icase)
    allocate(positions(ndim, nparticles, p))
    allocate(local_solver(nparticles, p))
    allocate(centroid_solver(nparticles, p))

    call build_positions(state, positions)
    call perturb_electron_beads(positions, p)

    call init_charge_memory(local_memory, nparticles, p, charge_policy_bead_local_acks2)
    call init_charge_memory(centroid_memory, nparticles, p, charge_policy_centroid_acks2)
    local_solver = 0.0_dp
    centroid_solver = 0.0_dp

    call evaluate_oracle_beads(adapter, positions, local_memory, local_solver, &
        charge_policy_bead_local_acks2, local_summary)
    call evaluate_oracle_beads(adapter, positions, centroid_memory, centroid_solver, &
        charge_policy_centroid_acks2, centroid_summary)

    local_diff = max_charge_difference_across_beads(local_memory)
    centroid_diff = max_charge_difference_across_beads(centroid_memory)
    energy_diff = abs(local_summary%pimd_potential_eV - centroid_summary%pimd_potential_eV)

    call assert_true("case status zero", all(local_summary%bead_status == 0) .and. &
        all(centroid_summary%bead_status == 0))
    call assert_true("classical subset replicated", classical_subset_replicated(positions))

    if (p == 1) then
      call assert_close("P=1 charge diff", local_diff, 0.0_dp, 1.0e-12_dp)
      call assert_close("P=1 centroid diff", centroid_diff, 0.0_dp, 1.0e-12_dp)
      call assert_close("P=1 energy diff", energy_diff, 0.0_dp, 1.0e-12_dp)
    else
      call assert_true("P>1 local charge diff positive", local_diff > 0.0_dp)
      call assert_close("P>1 centroid charge diff", centroid_diff, 0.0_dp, 1.0e-12_dp)
      call assert_true("P>1 energy diff positive", energy_diff > 1.0e-12_dp)
    end if

    deallocate(positions, local_solver, centroid_solver)
  end do

  write(*,'(a)') "test_15_mixed_oracle_manager_matrix: PASS"
  write(*,'(a,i0)') "  classical_count = ", count_particles_with_dynamics(assignment, dynamics_classical)
  write(*,'(a,i0)') "  pimd_count      = ", count_particles_with_dynamics(assignment, dynamics_pimd)
  write(*,'(a,3(i0,1x))') "  p_values        = ", p_values

contains

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
  end subroutine build_phase4_template

  subroutine build_positions(state_in, positions_out)
    type(system_state), intent(in) :: state_in
    real(dp), intent(out) :: positions_out(:, :, :)
    integer :: bead

    do bead = 1, size(positions_out, 3)
      positions_out(:, :, bead) = state_in%position_A(:, :)
    end do
  end subroutine build_positions

  subroutine perturb_electron_beads(positions_inout, p)
    real(dp), intent(inout) :: positions_inout(:, :, :)
    integer, intent(in) :: p
    integer :: bead
    real(dp) :: scale

    do bead = 1, p
      if (p == 1) cycle
      scale = real(bead - 1, dp) / real(p - 1, dp)
      positions_inout(1, 1, bead) = positions_inout(1, 1, bead) + 0.08_dp * scale
      positions_inout(2, 1, bead) = positions_inout(2, 1, bead) - 0.03_dp * scale
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
    groups_out(1)%thermostat = "langevin"
    groups_out(1)%temperature_K = 300.0_dp

    groups_out(2)%name = "excess_electron"
    groups_out(2)%selector = "element == El and id == 1"
    groups_out(2)%dynamics = dynamics_pimd
    groups_out(2)%mobile = .true.
    groups_out(2)%beads = 16
    groups_out(2)%thermostat = "pile_l"
    groups_out(2)%temperature_K = 300.0_dp

    call assign_particle_groups(state_in, groups_out, assignment_out)
  end subroutine assign_mixed_groups

  subroutine assert_true(label, condition)
    character(len=*), intent(in) :: label
    logical, intent(in) :: condition

    if (.not. condition) then
      write(*,'(a)') "FAIL: " // trim(label)
      error stop "mixed oracle-manager matrix regression failed"
    end if
  end subroutine assert_true

  subroutine assert_int(label, value, expected)
    character(len=*), intent(in) :: label
    integer, intent(in) :: value, expected

    if (value /= expected) then
      write(*,'(a)') "FAIL: " // trim(label)
      write(*,'(a,i0)') "  value    = ", value
      write(*,'(a,i0)') "  expected = ", expected
      error stop "mixed oracle-manager matrix regression failed"
    end if
  end subroutine assert_int

  subroutine assert_close(label, value, expected, threshold)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: value, expected, threshold

    if (abs(value - expected) > threshold) then
      write(*,'(a)') "FAIL: " // trim(label)
      write(*,'(a,es24.16)') "  value    = ", value
      write(*,'(a,es24.16)') "  expected = ", expected
      write(*,'(a,es24.16)') "  tol      = ", threshold
      error stop "mixed oracle-manager matrix regression failed"
    end if
  end subroutine assert_close

  logical function classical_subset_replicated(positions_in) result(is_replicated)
    real(dp), intent(in) :: positions_in(:, :, :)
    integer :: bead

    is_replicated = .true.
    do bead = 2, size(positions_in, 3)
      if (maxval(abs(positions_in(:, 2:nparticles, bead) - positions_in(:, 2:nparticles, 1))) > 1.0e-12_dp) then
        is_replicated = .false.
        return
      end if
    end do
  end function classical_subset_replicated

end program test_mixed_oracle_manager_matrix
