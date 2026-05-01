program test_invalid_group_overlap
  use ereaxff_kinds, only: dp
  use ereaxff_system_state, only: system_state, init_system_state, set_particle
  use ereaxff_particle_groups, only: particle_group_spec, particle_group_assignment, &
      assign_particle_groups, dynamics_classical
  implicit none

  type(system_state) :: state
  type(particle_group_spec) :: groups(2)
  type(particle_group_assignment) :: assignment

  call init_system_state(state, 3, 1)
  call set_particle(state, 1, 1, "C", 12.011_dp, [0.0_dp, 0.0_dp, 0.0_dp])

  groups(1)%name = "all_particles"
  groups(1)%selector = "all"
  groups(1)%dynamics = dynamics_classical
  groups(1)%beads = 1

  groups(2)%name = "carbons"
  groups(2)%selector = "element == C"
  groups(2)%dynamics = dynamics_classical
  groups(2)%beads = 1

  call assign_particle_groups(state, groups, assignment)
end program test_invalid_group_overlap

