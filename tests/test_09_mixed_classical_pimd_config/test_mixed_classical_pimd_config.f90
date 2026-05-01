program test_mixed_classical_pimd_config
  use ereaxff_kinds, only: dp
  use ereaxff_system_state, only: system_state, init_system_state, set_particle
  use ereaxff_selectors, only: selector_matches, selector_is_supported
  use ereaxff_particle_groups, only: particle_group_spec, particle_group_assignment, &
      assign_particle_groups, count_particles_with_dynamics, dynamics_classical, &
      dynamics_pimd, dynamics_frozen
  implicit none

  type(system_state) :: state
  type(particle_group_spec) :: groups(3)
  type(particle_group_assignment) :: assignment

  call init_system_state(state, 3, 6)
  call set_particle(state, 1, 1, "C", 12.011_dp, [0.0_dp, 0.0_dp, 0.0_dp])
  call set_particle(state, 2, 2, "C", 12.011_dp, [1.4_dp, 0.0_dp, 0.0_dp])
  call set_particle(state, 3, 3, "C", 12.011_dp, [0.0_dp, 1.4_dp, 0.0_dp])
  call set_particle(state, 4, 61, "El", 0.00054858_dp, [0.7_dp, 0.7_dp, 2.0_dp])
  call set_particle(state, 5, 62, "El", 0.00054858_dp, [2.0_dp, 2.0_dp, 2.0_dp])
  call set_particle(state, 6, 63, "El", 0.00054858_dp, [3.0_dp, 3.0_dp, 2.0_dp])

  groups(1)%name = "graphene_carbons"
  groups(1)%selector = "element == C"
  groups(1)%dynamics = dynamics_classical
  groups(1)%mobile = .true.
  groups(1)%beads = 1
  groups(1)%thermostat = "langevin"
  groups(1)%temperature_K = 300.0_dp

  groups(2)%name = "excess_electron"
  groups(2)%selector = "element == El and id == 61"
  groups(2)%dynamics = dynamics_pimd
  groups(2)%mobile = .true.
  groups(2)%beads = 16
  groups(2)%thermostat = "pile_l"
  groups(2)%temperature_K = 300.0_dp

  groups(3)%name = "background_electrons"
  groups(3)%selector = "element == El and id != 61"
  groups(3)%dynamics = dynamics_frozen
  groups(3)%mobile = .true.
  groups(3)%beads = 1
  groups(3)%thermostat = "none"
  groups(3)%temperature_K = 0.0_dp

  call assign_particle_groups(state, groups, assignment)

  call assert_int("classical count", count_particles_with_dynamics(assignment, dynamics_classical), 3)
  call assert_int("pimd count", count_particles_with_dynamics(assignment, dynamics_pimd), 1)
  call assert_int("frozen count", count_particles_with_dynamics(assignment, dynamics_frozen), 2)

  call assert_int("carbon group", assignment%group_index(1), 1)
  call assert_int("excess electron group", assignment%group_index(4), 2)
  call assert_int("background electron group", assignment%group_index(5), 3)
  call assert_int("PIMD bead count", assignment%beads(4), 16)
  call assert_true("frozen group forces mobile false", .not. assignment%mobile(5))
  call assert_true("PIMD thermostat selected", trim(assignment%thermostat(4)) == "pile_l")
  call assert_true("list selector matches id", selector_matches(state, 4, "element == El and id in [61,64]"))
  call assert_true("list selector rejects other id", .not. selector_matches(state, 5, "element == El and id in [61,64]"))
  call assert_true("invalid selector rejected", .not. selector_is_supported("element = C"))

  write(*,'(a)') "test_09_mixed_classical_pimd_config: PASS"
  write(*,'(a,i0)') "  classical_count = ", count_particles_with_dynamics(assignment, dynamics_classical)
  write(*,'(a,i0)') "  pimd_count      = ", count_particles_with_dynamics(assignment, dynamics_pimd)
  write(*,'(a,i0)') "  frozen_count    = ", count_particles_with_dynamics(assignment, dynamics_frozen)
  write(*,'(a,i0)') "  excess_beads    = ", assignment%beads(4)

contains

  subroutine assert_true(label, condition)
    character(len=*), intent(in) :: label
    logical, intent(in) :: condition

    if (.not. condition) then
      write(*,'(a)') "FAIL: " // trim(label)
      error stop "mixed config regression failed"
    end if
  end subroutine assert_true

  subroutine assert_int(label, value, expected)
    character(len=*), intent(in) :: label
    integer, intent(in) :: value, expected

    if (value /= expected) then
      write(*,'(a)') "FAIL: " // trim(label)
      write(*,'(a,i0)') "  value    = ", value
      write(*,'(a,i0)') "  expected = ", expected
      error stop "mixed config regression failed"
    end if
  end subroutine assert_int

end program test_mixed_classical_pimd_config

