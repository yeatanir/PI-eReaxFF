module ereaxff_harmonic_quantum_variance
  use ereaxff_kinds, only: dp
  use ereaxff_constants, only: boltzmann_eV_per_K, hbar_eV_fs
  use ereaxff_units, only: amu_angstrom2_per_fs2_to_eV
  implicit none
  private

  public :: harmonic_frequency_fs_inv
  public :: harmonic_quantum_variance_A2
  public :: harmonic_finite_p_bead_variance_A2

contains

  pure real(dp) function harmonic_frequency_fs_inv(k_eV_per_A2, mass_amu) result(omega_fs_inv)
    real(dp), intent(in) :: k_eV_per_A2, mass_amu

    if (k_eV_per_A2 <= 0.0_dp) error stop "harmonic_frequency_fs_inv: k must be positive"
    if (mass_amu <= 0.0_dp) error stop "harmonic_frequency_fs_inv: mass must be positive"

    omega_fs_inv = sqrt(k_eV_per_A2 / (mass_amu * amu_angstrom2_per_fs2_to_eV))
  end function harmonic_frequency_fs_inv

  pure real(dp) function harmonic_quantum_variance_A2(k_eV_per_A2, mass_amu, temperature_K) result(var_A2)
    real(dp), intent(in) :: k_eV_per_A2, mass_amu, temperature_K
    real(dp) :: omega, x

    if (temperature_K <= 0.0_dp) error stop "harmonic_quantum_variance_A2: temperature must be positive"

    omega = harmonic_frequency_fs_inv(k_eV_per_A2, mass_amu)
    x = hbar_eV_fs * omega / (2.0_dp * boltzmann_eV_per_K * temperature_K)
    var_A2 = hbar_eV_fs / (2.0_dp * mass_amu * omega * amu_angstrom2_per_fs2_to_eV) * coth(x)
  end function harmonic_quantum_variance_A2

  pure real(dp) function harmonic_finite_p_bead_variance_A2(k_eV_per_A2, mass_amu, temperature_K, omega_k) result(var_A2)
    real(dp), intent(in) :: k_eV_per_A2, mass_amu, temperature_K
    real(dp), intent(in) :: omega_k(:)
    integer :: mode, p
    real(dp) :: stiffness_eV_per_A2

    if (size(omega_k) < 1) error stop "harmonic_finite_p_bead_variance_A2: omega_k must not be empty"
    if (temperature_K <= 0.0_dp) error stop "harmonic_finite_p_bead_variance_A2: temperature must be positive"

    p = size(omega_k)
    var_A2 = 0.0_dp
    do mode = 1, p
      stiffness_eV_per_A2 = mass_amu * amu_angstrom2_per_fs2_to_eV * omega_k(mode)**2 + &
          k_eV_per_A2
      var_A2 = var_A2 + real(p, dp) * boltzmann_eV_per_K * temperature_K / stiffness_eV_per_A2
    end do
    var_A2 = var_A2 / real(p, dp)
  end function harmonic_finite_p_bead_variance_A2

  pure real(dp) function coth(x) result(value)
    real(dp), intent(in) :: x

    value = cosh(x) / sinh(x)
  end function coth

end module ereaxff_harmonic_quantum_variance
