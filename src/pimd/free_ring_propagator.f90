module ereaxff_free_ring_propagator
  use ereaxff_kinds, only: dp
  implicit none
  private

  public :: exact_free_propagation

contains

  subroutine exact_free_propagation(q_nm, v_nm, dt_fs, omega_k)
    real(dp), intent(inout) :: q_nm(:, :)
    real(dp), intent(inout) :: v_nm(:, :)
    real(dp), intent(in) :: dt_fs
    real(dp), intent(in) :: omega_k(:)
    integer :: coord, mode
    real(dp) :: omega, coswt, sinwt, q_old, v_old

    if (any(shape(v_nm) /= shape(q_nm))) error stop "exact_free_propagation: velocity shape mismatch"
    if (size(q_nm, 2) /= size(omega_k)) error stop "exact_free_propagation: omega_k size mismatch"

    do mode = 1, size(q_nm, 2)
      omega = omega_k(mode)
      if (abs(omega) <= tiny(1.0_dp)) then
        q_nm(:, mode) = q_nm(:, mode) + dt_fs * v_nm(:, mode)
      else
        coswt = cos(omega * dt_fs)
        sinwt = sin(omega * dt_fs)
        do coord = 1, size(q_nm, 1)
          q_old = q_nm(coord, mode)
          v_old = v_nm(coord, mode)
          q_nm(coord, mode) = q_old * coswt + (v_old / omega) * sinwt
          v_nm(coord, mode) = v_old * coswt - q_old * omega * sinwt
        end do
      end if
    end do
  end subroutine exact_free_propagation

end module ereaxff_free_ring_propagator

