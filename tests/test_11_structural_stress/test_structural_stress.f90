program test_structural_stress
  use ereaxff_kinds, only: dp
  use ereaxff_normal_modes, only: build_cmat, bead_to_nm, nm_to_bead
  use ereaxff_ring_polymer, only: ring_polymer_state, init_ring_polymer
  use ereaxff_harmonic_force_provider, only: harmonic_force_provider_state, harmonic_force_provider_init
  use ereaxff_bead_force_manager, only: bead_force_summary, evaluate_harmonic_beads
  use ereaxff_charge_coupling_policy, only: charge_memory_state, init_charge_memory, set_bead_charges, &
      max_charge_difference_across_beads, charge_policy_bead_local_acks2
  use ereaxff_system_state, only: system_state, init_system_state, set_particle
  use ereaxff_particle_groups, only: particle_group_spec, particle_group_assignment, assign_particle_groups, &
      count_particles_with_dynamics, dynamics_classical, dynamics_pimd, dynamics_frozen
  implicit none

  call stress_normal_modes()
  call stress_repeated_ring_initialization()
  call stress_large_force_manager()
  call stress_large_charge_memory()
  call stress_large_group_assignment()

  write(*,'(a)') "test_11_structural_stress: PASS"

contains

  subroutine stress_normal_modes()
    integer, parameter :: ncases = 9
    integer, parameter :: p_values(ncases) = [1, 2, 3, 4, 5, 16, 31, 64, 128]
    integer :: icase, p, i
    real(dp), allocatable :: cmat(:, :), ident(:, :), x(:), x_nm(:), x_round(:)
    real(dp) :: max_orth, max_round

    do icase = 1, ncases
      p = p_values(icase)
      allocate(cmat(p, p), ident(p, p), x(p), x_nm(p), x_round(p))
      call build_cmat(p, cmat)
      ident = matmul(cmat, transpose(cmat))
      do i = 1, p
        ident(i, i) = ident(i, i) - 1.0_dp
        x(i) = sin(0.19_dp * real(i, dp)) + cos(0.07_dp * real(i * i, dp))
      end do
      max_orth = maxval(abs(ident))
      call bead_to_nm(cmat, x, x_nm)
      call nm_to_bead(cmat, x_nm, x_round)
      max_round = maxval(abs(x_round - x))
      call assert_less("stress normal-mode orthogonality", max_orth, 5.0e-12_dp)
      call assert_less("stress normal-mode roundtrip", max_round, 5.0e-12_dp)
      deallocate(cmat, ident, x, x_nm, x_round)
    end do
  end subroutine stress_normal_modes

  subroutine stress_repeated_ring_initialization()
    type(ring_polymer_state) :: rp
    integer, parameter :: p_values(6) = [1, 7, 16, 32, 64, 8]
    integer :: i

    call init_ring_polymer(rp, p_values(1), 300.0_dp, "standard_ipi")
    call assert_int("repeated ring init bead count", rp%beads, p_values(1))
    call assert_int("repeated ring init omega size", size(rp%omega_k), p_values(1))
    call assert_int("repeated ring init cmat rows", size(rp%cmat, 1), p_values(1))

    do i = 2, size(p_values)
      call init_ring_polymer(rp, p_values(i), 300.0_dp, "standard_ipi")
      call assert_int("repeated ring init bead count", rp%beads, p_values(i))
      call assert_int("repeated ring init omega size", size(rp%omega_k), p_values(i))
      call assert_int("repeated ring init cmat rows", size(rp%cmat, 1), p_values(i))
    end do
  end subroutine stress_repeated_ring_initialization

  subroutine stress_large_force_manager()
    integer, parameter :: ndim = 3
    integer, parameter :: nparticles = 1024
    integer, parameter :: beads = 32
    real(dp), parameter :: k_eV_per_A2 = 0.75_dp
    type(harmonic_force_provider_state) :: provider
    type(bead_force_summary) :: summary
    real(dp), allocatable :: positions(:, :, :)
    real(dp) :: expected_energy
    integer :: i, bead

    allocate(positions(ndim, nparticles, beads))
    do bead = 1, beads
      do i = 1, nparticles
        positions(1, i, bead) = 0.001_dp * real(i, dp)
        positions(2, i, bead) = -0.002_dp * real(bead, dp)
        positions(3, i, bead) = 0.0005_dp * real(i + bead, dp)
      end do
    end do

    call harmonic_force_provider_init(provider, k_eV_per_A2)
    call evaluate_harmonic_beads(provider, positions, summary)

    expected_energy = 0.0_dp
    do bead = 1, beads
      expected_energy = expected_energy + 0.5_dp * k_eV_per_A2 * sum(positions(:, :, bead)**2)
    end do
    expected_energy = expected_energy / real(beads, dp)

    call assert_close("large force manager PIMD potential", summary%pimd_potential_eV, expected_energy, 1.0e-10_dp)
    call assert_close("large force manager force spot", summary%scaled_forces_eV_per_A(1, 100, 7), &
        -k_eV_per_A2 * positions(1, 100, 7) / real(beads, dp), 1.0e-14_dp)
  end subroutine stress_large_force_manager

  subroutine stress_large_charge_memory()
    integer, parameter :: nparticles = 2048
    integer, parameter :: beads = 32
    type(charge_memory_state) :: memory
    real(dp), allocatable :: charges(:)
    integer :: i, bead

    allocate(charges(nparticles))
    call init_charge_memory(memory, nparticles, beads, charge_policy_bead_local_acks2)
    do bead = 1, beads
      do i = 1, nparticles
        charges(i) = 1.0e-4_dp * real(mod(i, 17), dp) + 1.0e-3_dp * real(bead, dp)
      end do
      call set_bead_charges(memory, bead, charges)
    end do

    call assert_close("large charge memory max diff", max_charge_difference_across_beads(memory), &
        1.0e-3_dp * real(beads - 1, dp), 1.0e-14_dp)
  end subroutine stress_large_charge_memory

  subroutine stress_large_group_assignment()
    integer, parameter :: nparticles = 10000
    type(system_state) :: state
    type(particle_group_spec) :: groups(3)
    type(particle_group_assignment) :: assignment
    integer :: i
    character(len=8) :: element

    call init_system_state(state, 3, nparticles)
    do i = 1, nparticles
      if (i <= 6000) then
        element = "C"
        call set_particle(state, i, i, element, 12.011_dp, [0.01_dp * real(i, dp), 0.0_dp, 0.0_dp])
      else
        element = "El"
        call set_particle(state, i, i, element, 0.00054858_dp, [0.0_dp, 0.01_dp * real(i, dp), 0.0_dp])
      end if
    end do

    groups(1)%name = "carbons"
    groups(1)%selector = "element == C"
    groups(1)%dynamics = dynamics_classical
    groups(1)%mobile = .true.
    groups(1)%thermostat = "langevin"
    groups(1)%temperature_K = 300.0_dp

    groups(2)%name = "selected_electron"
    groups(2)%selector = "element == El and id == 9999"
    groups(2)%dynamics = dynamics_pimd
    groups(2)%mobile = .true.
    groups(2)%beads = 16
    groups(2)%thermostat = "pile_l"
    groups(2)%temperature_K = 300.0_dp

    groups(3)%name = "background_electrons"
    groups(3)%selector = "element == El and id != 9999"
    groups(3)%dynamics = dynamics_frozen
    groups(3)%mobile = .true.
    groups(3)%thermostat = "none"

    call assign_particle_groups(state, groups, assignment)
    call assert_int("large group classical count", count_particles_with_dynamics(assignment, dynamics_classical), 6000)
    call assert_int("large group PIMD count", count_particles_with_dynamics(assignment, dynamics_pimd), 1)
    call assert_int("large group frozen count", count_particles_with_dynamics(assignment, dynamics_frozen), 3999)
    call assert_int("large group selected beads", assignment%beads(9999), 16)
  end subroutine stress_large_group_assignment

  subroutine assert_less(label, value, threshold)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: value, threshold

    if (value > threshold) then
      write(*,'(a)') "FAIL: " // trim(label)
      write(*,'(a,es24.16)') "  value = ", value
      write(*,'(a,es24.16)') "  max   = ", threshold
      error stop "structural stress regression failed"
    end if
  end subroutine assert_less

  subroutine assert_close(label, value, expected, threshold)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: value, expected, threshold

    if (abs(value - expected) > threshold) then
      write(*,'(a)') "FAIL: " // trim(label)
      write(*,'(a,es24.16)') "  value    = ", value
      write(*,'(a,es24.16)') "  expected = ", expected
      write(*,'(a,es24.16)') "  tol      = ", threshold
      error stop "structural stress regression failed"
    end if
  end subroutine assert_close

  subroutine assert_int(label, value, expected)
    character(len=*), intent(in) :: label
    integer, intent(in) :: value, expected

    if (value /= expected) then
      write(*,'(a)') "FAIL: " // trim(label)
      write(*,'(a,i0)') "  value    = ", value
      write(*,'(a,i0)') "  expected = ", expected
      error stop "structural stress regression failed"
    end if
  end subroutine assert_int

end program test_structural_stress
