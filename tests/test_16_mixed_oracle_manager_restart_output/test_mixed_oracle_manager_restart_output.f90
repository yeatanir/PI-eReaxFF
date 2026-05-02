program test_mixed_oracle_manager_restart_output
  use ereaxff_kinds, only: dp
  use ereaxff_system_state, only: system_state, init_system_state, set_particle
  use ereaxff_particle_groups, only: particle_group_spec, particle_group_assignment, &
      assign_particle_groups, count_particles_with_dynamics, dynamics_classical, dynamics_pimd
  use ereaxff_charge_coupling_policy, only: charge_memory_state, init_charge_memory, &
      charge_policy_bead_local_acks2
  use ereaxff_oracle_adapter, only: oracle_adapter_state, oracle_adapter_init
  use ereaxff_bead_force_manager, only: oracle_bead_summary, evaluate_oracle_beads
  use ereaxff_mixed_restart, only: mixed_restart_state, write_mixed_restart, read_mixed_restart, &
      write_charge_by_bead_summary, write_charge_centroid_summary
  implicit none

  integer, parameter :: ndim = 3
  integer, parameter :: nparticles = 61
  integer, parameter :: p = 4
  type(system_state) :: state
  type(particle_group_spec) :: groups(2)
  type(particle_group_assignment) :: assignment
  type(oracle_adapter_state) :: adapter
  type(charge_memory_state) :: charge_memory
  type(oracle_bead_summary) :: summary
  type(mixed_restart_state) :: restart_out, restart_in
  real(dp) :: positions(ndim, nparticles, p)
  real(dp) :: solver_guess(nparticles, p)
  character(len=512) :: v1_root, v2_root, work_root, out_dir, restart_path, bead_path, centroid_path
  integer :: bead
  real(dp) :: local_charge_diff

  call get_environment_variable("EREAXFF_V1_ROOT", v1_root)
  call get_environment_variable("EREAXFF_V2_ROOT", v2_root)
  call get_environment_variable("EREAXFF_ADAPTER_WORK_ROOT", work_root)
  if (len_trim(v1_root) == 0) v1_root = "/Users/yeatanir/PhD/PIMD/ereaxff_pimd/V1"
  if (len_trim(v2_root) == 0) v2_root = "/Users/yeatanir/PhD/PIMD/ereaxff_pimd/V2"
  if (len_trim(work_root) == 0) work_root = "build/test_adapter_work"
  out_dir = "build/test_mixed_restart_output"
  restart_path = trim(out_dir) // "/mixed_restart.dat"
  bead_path = trim(out_dir) // "/charge_by_bead.out"
  centroid_path = trim(out_dir) // "/charge_centroid.out"

  call execute_command_line("mkdir -p '" // trim(out_dir) // "'")

  call build_phase4_template(state, positions(:, :, 1))
  do bead = 2, p
    positions(:, :, bead) = positions(:, :, 1)
  end do
  do bead = 1, p
    if (p > 1) then
      positions(1, 1, bead) = positions(1, 1, bead) + 0.05_dp * real(bead - 1, dp) / real(p - 1, dp)
      positions(2, 1, bead) = positions(2, 1, bead) - 0.02_dp * real(bead - 1, dp) / real(p - 1, dp)
    end if
  end do

  call assign_mixed_groups(state, groups, assignment)
  call assert_int("classical count", count_particles_with_dynamics(assignment, dynamics_classical), 60)
  call assert_int("pimd count", count_particles_with_dynamics(assignment, dynamics_pimd), 1)

  call oracle_adapter_init(adapter, state, trim(v1_root), trim(v2_root), trim(work_root))
  call init_charge_memory(charge_memory, nparticles, p, charge_policy_bead_local_acks2)
  solver_guess = 0.0_dp

  call evaluate_oracle_beads(adapter, positions, charge_memory, solver_guess, &
      charge_policy_bead_local_acks2, summary)

  restart_out%step = 777
  restart_out%ndim = ndim
  restart_out%nparticles = nparticles
  restart_out%beads = p
  restart_out%charge_policy = charge_policy_bead_local_acks2
  allocate(restart_out%positions_A(ndim, nparticles, p))
  allocate(restart_out%charges_by_bead(nparticles, p))
  allocate(restart_out%solver_guess_by_bead(nparticles, p))
  restart_out%positions_A = positions
  restart_out%charges_by_bead = charge_memory%charges_by_bead
  restart_out%solver_guess_by_bead = solver_guess

  call write_mixed_restart(restart_path, restart_out)
  call write_charge_by_bead_summary(bead_path, charge_memory%charges_by_bead)
  call write_charge_centroid_summary(centroid_path, charge_memory%charges_by_bead, solver_guess)

  call assert_file_exists(restart_path)
  call assert_file_exists(bead_path)
  call assert_file_exists(centroid_path)

  call read_mixed_restart(restart_path, restart_in)
  call assert_int("restart step", restart_in%step, restart_out%step)
  call assert_int("restart ndim", restart_in%ndim, restart_out%ndim)
  call assert_int("restart nparticles", restart_in%nparticles, restart_out%nparticles)
  call assert_int("restart beads", restart_in%beads, restart_out%beads)
  call assert_true("restart charge policy", trim(restart_in%charge_policy) == trim(restart_out%charge_policy))
  call assert_close("restart position roundtrip", maxval(abs(restart_in%positions_A - restart_out%positions_A)), &
      0.0_dp, 1.0e-12_dp)
  call assert_close("restart charge roundtrip", maxval(abs(restart_in%charges_by_bead - restart_out%charges_by_bead)), &
      0.0_dp, 1.0e-12_dp)
  call assert_close("restart solver roundtrip", maxval(abs(restart_in%solver_guess_by_bead - restart_out%solver_guess_by_bead)), &
      0.0_dp, 1.0e-12_dp)
  call assert_true("mixed output local split present", maxval(abs(charge_memory%charges_by_bead(:, 1) - &
      charge_memory%charges_by_bead(:, p))) > 0.0_dp)

  local_charge_diff = maxval(abs(charge_memory%charges_by_bead(:, 1) - charge_memory%charges_by_bead(:, p)))
  write(*,'(a)') "test_16_mixed_oracle_manager_restart_output: PASS"
  write(*,'(a,i0)') "  restart_step = ", restart_in%step
  write(*,'(a,es24.15)') "  local_charge_diff = ", local_charge_diff

contains

  subroutine build_phase4_template(state_out, positions_out)
    type(system_state), intent(inout) :: state_out
    real(dp), intent(out) :: positions_out(ndim, nparticles)
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
    positions_out(:, 1) = [el_x, el_y, el_z]

    idx = 2
    do ix = 0, 4
      do iy = 0, 2
        do unit_index = 1, 4
          x = (real(ix, dp) + unit(unit_index, 1)) * ax
          y = (real(iy, dp) + unit(unit_index, 2)) * ay
          call set_particle(state_out, idx, idx, "C", 12.011_dp, [x, y, z_c])
          positions_out(:, idx) = [x, y, z_c]
          idx = idx + 1
        end do
      end do
    end do
  end subroutine build_phase4_template

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

  subroutine assert_file_exists(path)
    character(len=*), intent(in) :: path
    logical :: exists

    inquire(file=trim(path), exist=exists)
    if (.not. exists) error stop "expected output file missing"
  end subroutine assert_file_exists

  subroutine assert_true(label, condition)
    character(len=*), intent(in) :: label
    logical, intent(in) :: condition

    if (.not. condition) then
      write(*,'(a)') "FAIL: " // trim(label)
      error stop "mixed restart/output regression failed"
    end if
  end subroutine assert_true

  subroutine assert_int(label, value, expected)
    character(len=*), intent(in) :: label
    integer, intent(in) :: value, expected

    if (value /= expected) then
      write(*,'(a)') "FAIL: " // trim(label)
      write(*,'(a,i0)') "  value    = ", value
      write(*,'(a,i0)') "  expected = ", expected
      error stop "mixed restart/output regression failed"
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
      error stop "mixed restart/output regression failed"
    end if
  end subroutine assert_close

end program test_mixed_oracle_manager_restart_output
