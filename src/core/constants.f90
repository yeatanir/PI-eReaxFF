module ereaxff_constants
  use ereaxff_kinds, only: dp
  implicit none
  private

  real(dp), parameter, public :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter, public :: twopi = 2.0_dp * pi

  ! Units used by the initial PIMD modules: eV, Angstrom, fs, amu, K.
  real(dp), parameter, public :: boltzmann_eV_per_K = 8.617333262145179e-5_dp
  real(dp), parameter, public :: hbar_eV_fs = 0.6582119569509067_dp
  real(dp), parameter, public :: kcalmol_to_eV = 0.0433641153087705_dp
  real(dp), parameter, public :: eV_to_kcalmol = 1.0_dp / kcalmol_to_eV
  real(dp), parameter, public :: angstrom_per_fs_to_m_per_s = 1.0e5_dp
end module ereaxff_constants

