program test_harmonic_quantum_variance
  use ereaxff_kinds, only: dp, i8
  use ereaxff_constants, only: pi, boltzmann_eV_per_K, hbar_eV_fs
  use ereaxff_rng, only: rng_state, rng_seed
  use ereaxff_ring_polymer, only: ring_polymer_state, init_ring_polymer
  use ereaxff_normal_modes, only: nm_to_bead
  use ereaxff_harmonic_quantum_variance, only: harmonic_quantum_variance_A2, &
      harmonic_finite_p_bead_variance_A2
  use ereaxff_thermostat_pile_l, only: pile_l_thermostat_state, pile_l_init
  use ereaxff_baoab_integrator, only: baoab_harmonic_pimd_step_1d
  implicit none

  integer, parameter :: np = 6
  integer, parameter :: p_values(np) = [1, 2, 4, 8, 16, 32]
  integer, parameter :: p_sample = 16
  integer, parameter :: nsteps = 180000
  real(dp), parameter :: temperature_K = 300.0_dp
  real(dp), parameter :: mass_amu = 1.0_dp
  real(dp), parameter :: k_eV_per_A2 = 1.0_dp
  real(dp), parameter :: dt_fs = 0.2_dp
  real(dp), parameter :: centroid_gamma_fs_inv = 0.02_dp
  type(ring_polymer_state) :: rp
  type(pile_l_thermostat_state) :: pile_l
  type(rng_state) :: rng
  real(dp) :: quantum_var, finite_var(np), errors(np), legacy_var
  real(dp) :: q_nm(p_sample), v_nm(p_sample), q_bead(p_sample)
  real(dp) :: sample_sum, sampled_var, expected_sample_var, rel_sample_error
  integer :: i, p, step, sample_count

  quantum_var = harmonic_quantum_variance_A2(k_eV_per_A2, mass_amu, temperature_K)

  do i = 1, np
    p = p_values(i)
    call init_ring_polymer(rp, p, temperature_K, "standard_ipi")
    finite_var(i) = harmonic_finite_p_bead_variance_A2(k_eV_per_A2, mass_amu, temperature_K, rp%omega_k)
    errors(i) = abs(finite_var(i) - quantum_var) / quantum_var
  end do

  do i = 2, np
    call assert_true("finite-P harmonic variance converges toward quantum limit", errors(i) < errors(i - 1))
  end do

  call init_ring_polymer(rp, p_sample, temperature_K, "standard_ipi")
  legacy_var = legacy_sqrt_prefactor_variance(p_sample)
  call assert_true("standard omega_P is closer to quantum variance than legacy sqrt(P)", &
      abs(finite_var(5) - quantum_var) < abs(legacy_var - quantum_var))

  call rng_seed(rng, 20260430_i8)
  call pile_l_init(pile_l, temperature_K * real(p_sample, dp), centroid_gamma_fs_inv, 1.0_dp, dt_fs, rp%omega_k)
  q_nm = 0.0_dp
  v_nm = 0.0_dp
  sample_sum = 0.0_dp
  sample_count = 0

  do step = 1, nsteps
    call baoab_harmonic_pimd_step_1d(q_nm, v_nm, rp%cmat, rp%omega_k, mass_amu, &
        k_eV_per_A2, dt_fs, pile_l, rng, physical_force_scale=1.0_dp)
    if (step > nsteps / 5 .and. mod(step, 5) == 0) then
      call nm_to_bead(rp%cmat, q_nm, q_bead)
      sample_sum = sample_sum + sum(q_bead**2) / real(p_sample, dp)
      sample_count = sample_count + 1
    end if
  end do

  sampled_var = sample_sum / real(sample_count, dp)
  expected_sample_var = finite_var(5)
  rel_sample_error = abs(sampled_var - expected_sample_var) / expected_sample_var
  call assert_true("sampled BAOAB variance is finite", sampled_var > 0.0_dp .and. sampled_var < huge(sampled_var))
  call assert_less("sampled BAOAB harmonic variance relative error", rel_sample_error, 0.15_dp)
  call assert_true("BAOAB path used thermostat heat accounting", abs(pile_l%heat_eV) > 1.0e-8_dp)

  write(*,'(a)') "test_03_harmonic_quantum_variance: PASS"
  write(*,'(a,es12.4)') "  quantum_variance_A2        = ", quantum_var
  write(*,'(a,es12.4)') "  finite_P16_variance_A2     = ", finite_var(5)
  write(*,'(a,es12.4)') "  legacy_sqrt_P16_variance   = ", legacy_var
  write(*,'(a,es12.4)') "  sampled_BAOAB_variance_A2  = ", sampled_var
  write(*,'(a,es12.4)') "  sampled_relative_error     = ", rel_sample_error

contains

  real(dp) function legacy_sqrt_prefactor_variance(p) result(var_A2)
    integer, intent(in) :: p
    real(dp) :: omega_k(p), omega_p, stiffness
    integer :: mode

    omega_p = sqrt(real(p, dp)) * boltzmann_eV_per_K * temperature_K / hbar_eV_fs
    omega_k(1) = 0.0_dp
    do mode = 2, p
      omega_k(mode) = 2.0_dp * omega_p * sin(pi * real(mode - 1, dp) / real(p, dp))
    end do
    var_A2 = 0.0_dp
    do mode = 1, p
      stiffness = mass_amu * 103.642696671608_dp * omega_k(mode)**2 + k_eV_per_A2
      var_A2 = var_A2 + real(p, dp) * boltzmann_eV_per_K * temperature_K / stiffness
    end do
    var_A2 = var_A2 / real(p, dp)
  end function legacy_sqrt_prefactor_variance

  subroutine assert_true(label, condition)
    character(len=*), intent(in) :: label
    logical, intent(in) :: condition

    if (.not. condition) then
      write(*,'(a)') "FAIL: " // trim(label)
      error stop "harmonic quantum variance regression failed"
    end if
  end subroutine assert_true

  subroutine assert_less(label, value, threshold)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: value, threshold

    if (value > threshold) then
      write(*,'(a)') "FAIL: " // trim(label)
      write(*,'(a,es24.16)') "  value = ", value
      write(*,'(a,es24.16)') "  max   = ", threshold
      error stop "harmonic quantum variance regression failed"
    end if
  end subroutine assert_less

end program test_harmonic_quantum_variance
