program test_multi_quantum_restart
  use ereaxff_kinds, only: dp, i8
  use ereaxff_system_state, only: system_state, init_system_state, set_particle
  use ereaxff_particle_groups, only: particle_group_spec, particle_group_assignment, &
      assign_particle_groups, count_particles_with_dynamics, dynamics_classical, dynamics_pimd
  use ereaxff_mixed_baoab_restart, only: mixed_baoab_restart_state, write_mixed_baoab_restart, &
      read_mixed_baoab_restart
  implicit none

  type(system_state) :: state
  type(particle_group_spec) :: groups(2)
  type(particle_group_assignment) :: assignment
  type(mixed_baoab_restart_state) :: restart_a, restart_b
  character(len=512) :: path
  integer :: i

  call get_command_argument(1, path)
  if (len_trim(path) == 0) error stop "test_multi_quantum_restart: output path required"

  call init_system_state(state, 3, 3)
  call set_particle(state, 1, 1, "C", 12.011_dp, [0.0_dp, 0.0_dp, 0.0_dp])
  call set_particle(state, 2, 2, "El", 1.0_dp, [1.0_dp, 0.0_dp, 0.0_dp])
  call set_particle(state, 3, 3, "El", 1.0_dp, [2.0_dp, 0.0_dp, 0.0_dp])

  groups(1)%name = "carbon"
  groups(1)%selector = "element == C"
  groups(1)%dynamics = dynamics_classical
  groups(1)%mobile = .true.
  groups(1)%beads = 1
  groups(1)%thermostat = "langevin"
  groups(1)%temperature_K = 300.0_dp
  groups(1)%mass_amu = 12.011_dp

  groups(2)%name = "electrons"
  groups(2)%selector = "element == El"
  groups(2)%dynamics = dynamics_pimd
  groups(2)%mobile = .true.
  groups(2)%beads = 4
  groups(2)%thermostat = "pile_l"
  groups(2)%temperature_K = 300.0_dp
  groups(2)%centroid_gamma_fs_inv = 0.001_dp
  groups(2)%mass_amu = 1.0_dp

  call assign_particle_groups(state, groups, assignment)
  call assert_int("classical_count", count_particles_with_dynamics(assignment, dynamics_classical), 1)
  call assert_int("pimd_count", count_particles_with_dynamics(assignment, dynamics_pimd), 2)

  restart_a%step = 12
  restart_a%ndim = 3
  restart_a%nparticles = 3
  restart_a%nquantum = 2
  restart_a%beads = 4
  restart_a%charge_policy = "bead_local_acks2"
  restart_a%classical_heat_eV = 0.25_dp
  restart_a%rng_state_i8 = int(42, kind=i8)
  allocate(restart_a%positions_A(3, 3, 4))
  allocate(restart_a%classical_velocities(3, 3))
  allocate(restart_a%quantum_indices(2))
  allocate(restart_a%quantum_masses_amu(2))
  allocate(restart_a%quantum_heat_eV(2))
  allocate(restart_a%q_nm(3, 2, 4))
  allocate(restart_a%v_nm(3, 2, 4))
  allocate(restart_a%charges_by_bead(3, 4))
  allocate(restart_a%solver_guess_by_bead(3, 4))
  restart_a%quantum_indices = [2, 3]
  restart_a%quantum_masses_amu = [1.0_dp, 1.0_dp]
  restart_a%quantum_heat_eV = [0.1_dp, 0.2_dp]
  restart_a%positions_A = 0.0_dp
  restart_a%classical_velocities = 0.0_dp
  restart_a%q_nm = 0.0_dp
  restart_a%v_nm = 0.0_dp
  restart_a%charges_by_bead = 0.0_dp
  restart_a%solver_guess_by_bead = 0.0_dp
  do i = 1, 4
    restart_a%q_nm(:, 1, i) = [0.1_dp * real(i, dp), 0.0_dp, 0.0_dp]
    restart_a%q_nm(:, 2, i) = [0.2_dp * real(i, dp), 0.1_dp, 0.0_dp]
    restart_a%v_nm(:, 1, i) = [0.01_dp * real(i, dp), 0.0_dp, 0.0_dp]
    restart_a%v_nm(:, 2, i) = [0.02_dp * real(i, dp), 0.01_dp, 0.0_dp]
    restart_a%positions_A(:, 1, i) = [0.0_dp, 0.0_dp, 0.0_dp]
    restart_a%positions_A(:, 2, i) = [1.0_dp, 0.0_dp, 0.0_dp]
    restart_a%positions_A(:, 3, i) = [2.0_dp, 0.0_dp, 0.0_dp]
    restart_a%charges_by_bead(:, i) = [0.0_dp, -0.5_dp, -0.5_dp]
    restart_a%solver_guess_by_bead(:, i) = [0.0_dp, -1.0_dp, -1.0_dp]
  end do

  call write_mixed_baoab_restart(trim(path), restart_a)
  call read_mixed_baoab_restart(trim(path), restart_b)

  call assert_int("restart_nquantum", restart_b%nquantum, 2)
  call assert_int("restart_quantum_index_1", restart_b%quantum_indices(1), 2)
  call assert_int("restart_quantum_index_2", restart_b%quantum_indices(2), 3)
  call assert_close("restart_heat_1", restart_b%quantum_heat_eV(1), 0.1_dp, 1.0e-12_dp)
  call assert_close("restart_heat_2", restart_b%quantum_heat_eV(2), 0.2_dp, 1.0e-12_dp)
  call assert_close("restart_qnm", restart_b%q_nm(1, 2, 4), restart_a%q_nm(1, 2, 4), 1.0e-12_dp)
  call assert_close("restart_vnm", restart_b%v_nm(2, 1, 3), restart_a%v_nm(2, 1, 3), 1.0e-12_dp)
  call assert_close("restart_charge", restart_b%charges_by_bead(3, 2), restart_a%charges_by_bead(3, 2), 1.0e-12_dp)

  write(*,'(a)') "test_27_multi_quantum_restart: PASS"
  write(*,'(a,i0)') "  pimd_count = ", count_particles_with_dynamics(assignment, dynamics_pimd)
  write(*,'(a,i0)') "  nquantum   = ", restart_b%nquantum

contains

  subroutine assert_int(label, got, want)
    character(len=*), intent(in) :: label
    integer, intent(in) :: got, want

    if (got /= want) error stop "assert_int failed: " // trim(label)
  end subroutine assert_int

  subroutine assert_close(label, got, want, tol)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: got, want, tol

    if (abs(got - want) > tol) error stop "assert_close failed: " // trim(label)
  end subroutine assert_close

end program test_multi_quantum_restart
