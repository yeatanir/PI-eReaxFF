module ereaxff_pimd_output
  use ereaxff_kinds, only: dp
  implicit none
  private

  public :: write_centroid_xyz_frame
  public :: compute_centroid_1d
  public :: compute_bead_rms_spread_1d

contains

  subroutine write_centroid_xyz_frame(unit, step, centroid_A)
    integer, intent(in) :: unit, step
    real(dp), intent(in) :: centroid_A

    write(unit,'(i0)') 1
    write(unit,'(a,i0)') "step=", step
    write(unit,'(a,3(1x,es22.14))') "H", centroid_A, 0.0_dp, 0.0_dp
  end subroutine write_centroid_xyz_frame

  pure real(dp) function compute_centroid_1d(q_bead) result(centroid_A)
    real(dp), intent(in) :: q_bead(:)

    if (size(q_bead) < 1) error stop "compute_centroid_1d: empty bead vector"
    centroid_A = sum(q_bead) / real(size(q_bead), dp)
  end function compute_centroid_1d

  pure real(dp) function compute_bead_rms_spread_1d(q_bead) result(spread_A)
    real(dp), intent(in) :: q_bead(:)
    real(dp) :: centroid_A

    centroid_A = compute_centroid_1d(q_bead)
    spread_A = sqrt(sum((q_bead - centroid_A)**2) / real(size(q_bead), dp))
  end function compute_bead_rms_spread_1d

end module ereaxff_pimd_output

