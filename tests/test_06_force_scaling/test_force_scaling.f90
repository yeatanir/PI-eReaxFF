program test_force_scaling
  use ereaxff_kinds, only: dp
  use ereaxff_force_provider_base, only: force_result, batched_force_result, allocate_batched_force_result
  use ereaxff_harmonic_force_provider, only: harmonic_force_provider_state, &
      harmonic_force_provider_init, harmonic_force_provider_evaluate
  use ereaxff_dummy_force_provider, only: dummy_force_provider_state, &
      dummy_force_provider_init, dummy_force_provider_evaluate
  use ereaxff_bead_force_manager, only: bead_force_summary, evaluate_harmonic_beads, &
      evaluate_dummy_beads, average_forces_over_beads
  implicit none

  integer, parameter :: ndim = 3
  integer, parameter :: nparticles = 2
  integer, parameter :: p = 4
  real(dp), parameter :: k_eV_per_A2 = 2.5_dp
  type(harmonic_force_provider_state) :: harmonic
  type(dummy_force_provider_state) :: dummy
  type(force_result) :: result
  type(batched_force_result) :: batched_result
  type(bead_force_summary) :: harmonic_summary, dummy_summary
  real(dp) :: positions(ndim, nparticles)
  real(dp) :: bead_positions(ndim, nparticles, p)
  real(dp) :: mean_force(ndim, nparticles)
  real(dp) :: expected_energy, expected_pimd_energy
  integer :: bead

  positions(:, 1) = [1.0_dp, -2.0_dp, 0.5_dp]
  positions(:, 2) = [-0.25_dp, 0.75_dp, 1.25_dp]

  call harmonic_force_provider_init(harmonic, k_eV_per_A2)
  call harmonic_force_provider_evaluate(harmonic, positions, result)

  expected_energy = 0.5_dp * k_eV_per_A2 * sum(positions**2)
  call assert_close("harmonic provider energy", result%energy_eV, expected_energy, 1.0e-14_dp)
  call assert_array_close("harmonic provider force", result%forces_eV_per_A, -k_eV_per_A2 * positions, 1.0e-14_dp)

  call dummy_force_provider_init(dummy)
  call dummy_force_provider_evaluate(dummy, positions, result)
  call assert_close("dummy provider energy", result%energy_eV, 0.0_dp, 0.0_dp)
  call assert_array_close("dummy provider force", result%forces_eV_per_A, 0.0_dp * positions, 0.0_dp)

  do bead = 1, p
    bead_positions(:, :, bead) = positions * real(bead, dp)
  end do

  call evaluate_harmonic_beads(harmonic, bead_positions, harmonic_summary)
  expected_pimd_energy = 0.0_dp
  do bead = 1, p
    expected_pimd_energy = expected_pimd_energy + 0.5_dp * k_eV_per_A2 * sum(bead_positions(:, :, bead)**2)
    call assert_array_close("bead force scaling", harmonic_summary%scaled_forces_eV_per_A(:, :, bead), &
        -k_eV_per_A2 * bead_positions(:, :, bead) / real(p, dp), 1.0e-14_dp)
  end do
  expected_pimd_energy = expected_pimd_energy / real(p, dp)
  call assert_close("PIMD bead-averaged potential", harmonic_summary%pimd_potential_eV, expected_pimd_energy, 1.0e-14_dp)

  call average_forces_over_beads(harmonic_summary%scaled_forces_eV_per_A, mean_force)
  call assert_array_close("sum scaled bead forces for classical projection", mean_force, &
      -k_eV_per_A2 * sum(bead_positions, dim=3) / real(p, dp), 1.0e-14_dp)

  call evaluate_dummy_beads(dummy, bead_positions, dummy_summary)
  call assert_close("dummy bead PIMD potential", dummy_summary%pimd_potential_eV, 0.0_dp, 0.0_dp)
  call assert_array_close("dummy bead forces", dummy_summary%scaled_forces_eV_per_A, &
      0.0_dp * dummy_summary%scaled_forces_eV_per_A, 0.0_dp)

  call allocate_batched_force_result(batched_result, ndim, nparticles, p, with_charges=.true.)
  call assert_int("batched result bead count", batched_result%nbeads, p)
  call assert_int("batched energy shape", size(batched_result%energies_eV), p)
  call assert_int("batched force bead shape", size(batched_result%forces_eV_per_A, 3), p)
  call assert_int("batched charge bead shape", size(batched_result%charges, 2), p)

  write(*,'(a)') "test_06_force_scaling: PASS"
  write(*,'(a,es12.4)') "  harmonic_energy_eV       = ", expected_energy
  write(*,'(a,es12.4)') "  pimd_average_energy_eV   = ", harmonic_summary%pimd_potential_eV
  write(*,'(a,es12.4)') "  max_scaled_force_abs     = ", maxval(abs(harmonic_summary%scaled_forces_eV_per_A))

contains

  subroutine assert_close(label, value, expected, threshold)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: value, expected, threshold

    if (abs(value - expected) > threshold) then
      write(*,'(a)') "FAIL: " // trim(label)
      write(*,'(a,es24.16)') "  value    = ", value
      write(*,'(a,es24.16)') "  expected = ", expected
      write(*,'(a,es24.16)') "  tol      = ", threshold
      error stop "force scaling regression failed"
    end if
  end subroutine assert_close

  subroutine assert_int(label, value, expected)
    character(len=*), intent(in) :: label
    integer, intent(in) :: value, expected

    if (value /= expected) then
      write(*,'(a)') "FAIL: " // trim(label)
      write(*,'(a,i0)') "  value    = ", value
      write(*,'(a,i0)') "  expected = ", expected
      error stop "force integer regression failed"
    end if
  end subroutine assert_int

  subroutine assert_array_close(label, value, expected, threshold)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: value(..)
    real(dp), intent(in) :: expected(..)
    real(dp), intent(in) :: threshold
    real(dp) :: max_error

    select rank (value)
    rank (2)
      select rank (expected)
      rank (2)
        if (any(shape(value) /= shape(expected))) error stop "assert_array_close: rank-2 shape mismatch"
        max_error = maxval(abs(value - expected))
      rank default
        error stop "assert_array_close: expected rank mismatch"
      end select
    rank (3)
      select rank (expected)
      rank (3)
        if (any(shape(value) /= shape(expected))) error stop "assert_array_close: rank-3 shape mismatch"
        max_error = maxval(abs(value - expected))
      rank default
        error stop "assert_array_close: expected rank mismatch"
      end select
    rank default
      error stop "assert_array_close: unsupported rank"
    end select

    if (max_error > threshold) then
      write(*,'(a)') "FAIL: " // trim(label)
      write(*,'(a,es24.16)') "  max_error = ", max_error
      write(*,'(a,es24.16)') "  tol       = ", threshold
      error stop "force array regression failed"
    end if
  end subroutine assert_array_close

end program test_force_scaling
