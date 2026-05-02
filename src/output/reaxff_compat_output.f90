module ereaxff_reaxff_compat_output
  use ereaxff_kinds, only: dp
  implicit none
  private

  public :: write_reaxff_compat_snapshot

contains

  subroutine write_reaxff_compat_snapshot(output_dir, step, centroid_positions_A, centroid_charges, energy_eV, &
      max_force_eV_per_A, charge_sum, beads)
    character(len=*), intent(in) :: output_dir
    integer, intent(in) :: step, beads
    real(dp), intent(in) :: centroid_positions_A(:, :)
    real(dp), intent(in) :: centroid_charges(:)
    real(dp), intent(in) :: energy_eV, max_force_eV_per_A, charge_sum
    character(len=:), allocatable :: fort7_path, fort8_path, molfra_path, xmol_path
    integer :: unit, i

    fort7_path = trim(output_dir) // "/fort.7"
    fort8_path = trim(output_dir) // "/fort.8"
    molfra_path = trim(output_dir) // "/molfra.out"
    xmol_path = trim(output_dir) // "/xmolout"

    open(newunit=unit, file=fort7_path, status="replace", action="write")
    write(unit,'(a)') "# step energy_eV max_force_eV_per_A charge_sum beads"
    write(unit,'(i8,3(1x,es24.16),1x,i8)') step, energy_eV, max_force_eV_per_A, charge_sum, beads
    close(unit)

    open(newunit=unit, file=fort8_path, status="replace", action="write")
    write(unit,'(a)') "# atom x_A y_A z_A charge"
    do i = 1, size(centroid_charges)
      write(unit,'(i8,4(1x,es24.16))') i, centroid_positions_A(:, i), centroid_charges(i)
    end do
    close(unit)

    open(newunit=unit, file=molfra_path, status="replace", action="write")
    write(unit,'(a)') "# centroid-projected molecular summary"
    write(unit,'(a,1x,i0)') "step", step
    write(unit,'(a,1x,i0)') "beads", beads
    write(unit,'(a,1x,es24.16)') "energy_eV", energy_eV
    write(unit,'(a,1x,es24.16)') "max_force_eV_per_A", max_force_eV_per_A
    write(unit,'(a,1x,es24.16)') "charge_sum", charge_sum
    close(unit)

    open(newunit=unit, file=xmol_path, status="replace", action="write")
    write(unit,'(i0)') size(centroid_charges)
    write(unit,'(a,i0)') "step=", step
    do i = 1, size(centroid_charges)
      write(unit,'(a,3(1x,es22.14),1x,es22.14)') "El", centroid_positions_A(:, i), centroid_charges(i)
    end do
    close(unit)
  end subroutine write_reaxff_compat_snapshot

end module ereaxff_reaxff_compat_output
