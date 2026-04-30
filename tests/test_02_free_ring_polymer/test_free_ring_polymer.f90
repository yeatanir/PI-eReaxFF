program test_free_ring_polymer
  use ereaxff_kinds, only: dp
  use ereaxff_ring_polymer, only: ring_polymer_state, init_ring_polymer
  use ereaxff_free_ring_propagator, only: exact_free_propagation
  use ereaxff_spring_energy, only: compute_free_ring_energy_nm, compute_spring_energy_nm
  implicit none

  integer, parameter :: ndim = 3
  integer, parameter :: p = 16
  integer, parameter :: nsteps = 20000
  real(dp), parameter :: temperature_K = 300.0_dp
  real(dp), parameter :: dt_fs = 0.05_dp
  real(dp), parameter :: mass_amu(ndim) = [1.0_dp, 2.0_dp, 12.0_dp]
  type(ring_polymer_state) :: rp
  real(dp) :: q_nm(ndim, p), v_nm(ndim, p)
  real(dp) :: e0, e1, drift, spring_e
  integer :: coord, mode, step

  call init_ring_polymer(rp, p, temperature_K, "standard_ipi")

  do mode = 1, p
    do coord = 1, ndim
      q_nm(coord, mode) = 0.03_dp * sin(0.17_dp * real(coord * mode, dp)) + &
          0.01_dp * cos(0.11_dp * real(coord + 2 * mode, dp))
      v_nm(coord, mode) = 0.002_dp * cos(0.13_dp * real(coord * mode, dp)) - &
          0.001_dp * sin(0.07_dp * real(3 * coord + mode, dp))
    end do
  end do

  e0 = compute_free_ring_energy_nm(q_nm, v_nm, mass_amu, rp%omega_k)
  spring_e = compute_spring_energy_nm(q_nm, mass_amu, rp%omega_k)

  do step = 1, nsteps
    call exact_free_propagation(q_nm, v_nm, dt_fs, rp%omega_k)
  end do

  e1 = compute_free_ring_energy_nm(q_nm, v_nm, mass_amu, rp%omega_k)
  drift = abs(e1 - e0)
  call assert_less("free ring-polymer energy drift", drift, 1.0e-10_dp)
  call assert_less("centroid frequency is zero", abs(rp%omega_k(1)), 0.0_dp)

  write(*,'(a)') "test_02_free_ring_polymer: PASS"
  write(*,'(a,es12.4)') "  initial_free_energy_eV = ", e0
  write(*,'(a,es12.4)') "  initial_spring_eV      = ", spring_e
  write(*,'(a,es12.4)') "  abs_energy_drift_eV    = ", drift

contains

  subroutine assert_less(label, value, threshold)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: value, threshold

    if (value > threshold) then
      write(*,'(a)') "FAIL: " // trim(label)
      write(*,'(a,es24.16)') "  value = ", value
      write(*,'(a,es24.16)') "  max   = ", threshold
      error stop "free ring-polymer regression failed"
    end if
  end subroutine assert_less

end program test_free_ring_polymer

