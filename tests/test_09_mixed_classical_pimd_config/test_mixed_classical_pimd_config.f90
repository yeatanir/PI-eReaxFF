program test_mixed_classical_pimd_config
  use ereaxff_kinds, only: dp
  use ereaxff_system_state, only: system_state, init_system_state, set_particle
  use ereaxff_selectors, only: selector_matches, selector_is_supported
  use ereaxff_particle_groups, only: particle_group_spec, particle_group_assignment, &
      assign_particle_groups, count_particles_with_dynamics, dynamics_classical, &
      dynamics_pimd, dynamics_frozen
  use ereaxff_input_config, only: input_config, read_input_config
  implicit none

  type(system_state) :: state
  type(particle_group_spec) :: groups(4)
  type(particle_group_assignment) :: assignment
  type(input_config) :: cfg
  character(len=*), parameter :: parsed_config_path = "build/test_09_multi_pimd_groups.yaml"

  call init_system_state(state, 3, 7)
  call set_particle(state, 1, 1, "C", 12.011_dp, [0.0_dp, 0.0_dp, 0.0_dp])
  call set_particle(state, 2, 2, "C", 12.011_dp, [1.4_dp, 0.0_dp, 0.0_dp])
  call set_particle(state, 3, 3, "C", 12.011_dp, [0.0_dp, 1.4_dp, 0.0_dp])
  call set_particle(state, 4, 61, "El", 0.00054858_dp, [0.7_dp, 0.7_dp, 2.0_dp])
  call set_particle(state, 5, 62, "El", 0.00054858_dp, [2.0_dp, 2.0_dp, 2.0_dp])
  call set_particle(state, 6, 63, "El", 0.00054858_dp, [3.0_dp, 3.0_dp, 2.0_dp])
  call set_particle(state, 7, 64, "El", 0.00054858_dp, [4.0_dp, 4.0_dp, 2.0_dp])

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
  groups(2)%mass_amu = 0.00054858_dp

  groups(3)%name = "second_quantum_electron"
  groups(3)%selector = "element == El and id == 64"
  groups(3)%dynamics = dynamics_pimd
  groups(3)%mobile = .true.
  groups(3)%beads = 8
  groups(3)%thermostat = "pile_l"
  groups(3)%temperature_K = 300.0_dp
  groups(3)%mass_amu = 0.00109716_dp

  groups(4)%name = "background_electrons"
  groups(4)%selector = "element == El and id in [62,63]"
  groups(4)%dynamics = dynamics_frozen
  groups(4)%mobile = .true.
  groups(4)%beads = 1
  groups(4)%thermostat = "none"
  groups(4)%temperature_K = 0.0_dp

  call assign_particle_groups(state, groups, assignment)

  call assert_int("classical count", count_particles_with_dynamics(assignment, dynamics_classical), 3)
  call assert_int("pimd count", count_particles_with_dynamics(assignment, dynamics_pimd), 2)
  call assert_int("frozen count", count_particles_with_dynamics(assignment, dynamics_frozen), 2)

  call assert_int("carbon group", assignment%group_index(1), 1)
  call assert_int("excess electron group", assignment%group_index(4), 2)
  call assert_int("background electron group", assignment%group_index(5), 4)
  call assert_int("second PIMD electron group", assignment%group_index(7), 3)
  call assert_int("first PIMD bead count", assignment%beads(4), 16)
  call assert_int("second PIMD bead count", assignment%beads(7), 8)
  call assert_close("first PIMD explicit mass", assignment%mass_amu(4), 0.00054858_dp, 1.0e-12_dp)
  call assert_close("second PIMD explicit mass", assignment%mass_amu(7), 0.00109716_dp, 1.0e-12_dp)
  call assert_true("frozen group forces mobile false", .not. assignment%mobile(5))
  call assert_true("PIMD thermostat selected", trim(assignment%thermostat(4)) == "pile_l")
  call assert_true("list selector matches id", selector_matches(state, 4, "element == El and id in [61,64]"))
  call assert_true("list selector rejects other id", .not. selector_matches(state, 5, "element == El and id in [61,64]"))
  call assert_true("invalid selector rejected", .not. selector_is_supported("element = C"))

  call write_multi_pimd_config(parsed_config_path)
  call read_input_config(parsed_config_path, cfg)
  call assert_int("parsed group count", size(cfg%groups), 4)
  call assert_int("parsed first PIMD beads", cfg%groups(2)%beads, 16)
  call assert_int("parsed second PIMD beads", cfg%groups(3)%beads, 8)
  call assert_close("parsed first PIMD mass", cfg%groups(2)%mass_amu, 0.00054858_dp, 1.0e-12_dp)
  call assert_close("parsed second PIMD mass", cfg%groups(3)%mass_amu, 0.00109716_dp, 1.0e-12_dp)
  call assert_true("parsed production oracle mode", .not. cfg%force_provider%repeatability_check)

  write(*,'(a)') "test_09_mixed_classical_pimd_config: PASS"
  write(*,'(a,i0)') "  classical_count = ", count_particles_with_dynamics(assignment, dynamics_classical)
  write(*,'(a,i0)') "  pimd_count      = ", count_particles_with_dynamics(assignment, dynamics_pimd)
  write(*,'(a,i0)') "  frozen_count    = ", count_particles_with_dynamics(assignment, dynamics_frozen)
  write(*,'(a,i0)') "  excess_beads    = ", assignment%beads(4)
  write(*,'(a,i0)') "  shell_beads     = ", assignment%beads(7)

contains

  subroutine write_multi_pimd_config(path)
    character(len=*), intent(in) :: path
    integer :: unit

    open(newunit=unit, file=trim(path), status="replace", action="write")
    write(unit,'(a)') "simulation:"
    write(unit,'(a)') "  dt_fs: 0.1"
    write(unit,'(a)') "  steps: 0"
    write(unit,'(a)') "force_provider:"
    write(unit,'(a)') "  type: ereaxff_oracle"
    write(unit,'(a)') "  repeatability_check: false"
    write(unit,'(a)') "ring_polymer:"
    write(unit,'(a)') "  default_beads: 16"
    write(unit,'(a)') "charge_coupling:"
    write(unit,'(a)') "  mode: bead_local_acks2"
    write(unit,'(a)') "geometry:"
    write(unit,'(a)') "  template: phase4_fixed_graphene"
    write(unit,'(a)') "groups:"
    write(unit,'(a)') "  - name: carbons"
    write(unit,'(a)') "    selector: ""element == C"""
    write(unit,'(a)') "    dynamics: classical"
    write(unit,'(a)') "    beads: 1"
    write(unit,'(a)') "    thermostat: langevin"
    write(unit,'(a)') "    temperature_K: 300.0"
    write(unit,'(a)') "  - name: excess"
    write(unit,'(a)') "    selector: ""element == El and id == 61"""
    write(unit,'(a)') "    dynamics: pimd"
    write(unit,'(a)') "    beads: 16"
    write(unit,'(a)') "    thermostat: pile_l"
    write(unit,'(a)') "    temperature_K: 300.0"
    write(unit,'(a)') "    mass_amu: 0.00054858"
    write(unit,'(a)') "  - name: selected_pair"
    write(unit,'(a)') "    selector: ""element == El and id == 64"""
    write(unit,'(a)') "    dynamics: pimd"
    write(unit,'(a)') "    beads: 8"
    write(unit,'(a)') "    thermostat: pile_l"
    write(unit,'(a)') "    temperature_K: 300.0"
    write(unit,'(a)') "    mass_amu: 0.00109716"
    write(unit,'(a)') "  - name: background"
    write(unit,'(a)') "    selector: ""element == El and id in [62,63]"""
    write(unit,'(a)') "    dynamics: frozen"
    write(unit,'(a)') "    beads: 1"
    write(unit,'(a)') "    thermostat: none"
    write(unit,'(a)') "output:"
    write(unit,'(a)') "  diagnostics: true"
    close(unit)
  end subroutine write_multi_pimd_config

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

  subroutine assert_close(label, value, expected, tol)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: value, expected, tol

    if (abs(value - expected) > tol) then
      write(*,'(a)') "FAIL: " // trim(label)
      write(*,'(a,es22.14)') "  value    = ", value
      write(*,'(a,es22.14)') "  expected = ", expected
      error stop "mixed config regression failed"
    end if
  end subroutine assert_close

end program test_mixed_classical_pimd_config
