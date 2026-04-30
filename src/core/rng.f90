module ereaxff_rng
  use ereaxff_kinds, only: dp, i8
  implicit none
  private

  type, public :: rng_state
    integer(i8) :: state = 1_i8
  end type rng_state

  public :: rng_seed
  public :: rng_uniform
  public :: rng_normal

contains

  subroutine rng_seed(rng, seed)
    type(rng_state), intent(inout) :: rng
    integer(i8), intent(in) :: seed

    rng%state = modulo(seed, 2147483647_i8)
    if (rng%state <= 0_i8) rng%state = rng%state + 2147483646_i8
  end subroutine rng_seed

  real(dp) function rng_uniform(rng) result(u)
    type(rng_state), intent(inout) :: rng

    rng%state = modulo(16807_i8 * rng%state, 2147483647_i8)
    u = real(rng%state, dp) / 2147483647.0_dp
  end function rng_uniform

  real(dp) function rng_normal(rng) result(z)
    type(rng_state), intent(inout) :: rng
    real(dp) :: u1, u2

    u1 = max(rng_uniform(rng), tiny(1.0_dp))
    u2 = rng_uniform(rng)
    z = sqrt(-2.0_dp * log(u1)) * cos(6.283185307179586476925286766559_dp * u2)
  end function rng_normal

end module ereaxff_rng

