module ereaxff_mixed_restart
  use ereaxff_kinds, only: dp
  implicit none
  private

  type, public :: mixed_restart_state
    integer :: step = 0
    integer :: ndim = 0
    integer :: nparticles = 0
    integer :: beads = 0
    character(len=32) :: charge_policy = ""
    real(dp), allocatable :: positions_A(:, :, :)
    real(dp), allocatable :: charges_by_bead(:, :)
    real(dp), allocatable :: solver_guess_by_bead(:, :)
  end type mixed_restart_state

  public :: write_mixed_restart
  public :: read_mixed_restart
  public :: write_charge_by_bead_summary
  public :: write_charge_centroid_summary
  public :: write_electron_localization_summary

contains

  subroutine write_mixed_restart(path, state)
    character(len=*), intent(in) :: path
    type(mixed_restart_state), intent(in) :: state
    integer :: unit, bead, i

    if (.not. allocated(state%positions_A)) error stop "write_mixed_restart: positions not allocated"
    if (.not. allocated(state%charges_by_bead)) error stop "write_mixed_restart: charges not allocated"
    if (.not. allocated(state%solver_guess_by_bead)) error stop "write_mixed_restart: solver memory not allocated"

    open(newunit=unit, file=trim(path), status="replace", action="write")
    write(unit,'(a)') "# PI-eReaxFF mixed restart v1"
    write(unit,'(a,1x,i0)') "step", state%step
    write(unit,'(a,1x,i0)') "ndim", state%ndim
    write(unit,'(a,1x,i0)') "nparticles", state%nparticles
    write(unit,'(a,1x,i0)') "beads", state%beads
    write(unit,'(a,1x,a)') "charge_policy", trim(state%charge_policy)
    write(unit,'(a)') "positions_A"
    do bead = 1, state%beads
      write(unit,'(a,1x,i0)') "bead", bead
      do i = 1, state%nparticles
        write(unit,'(i8,3(1x,es24.16))') i, state%positions_A(:, i, bead)
      end do
    end do
    write(unit,'(a)') "charges_by_bead"
    do bead = 1, state%beads
      write(unit,'(a,1x,i0)') "bead", bead
      do i = 1, state%nparticles
        write(unit,'(i8,2(1x,es24.16))') i, state%charges_by_bead(i, bead), state%solver_guess_by_bead(i, bead)
      end do
    end do
    close(unit)
  end subroutine write_mixed_restart

  subroutine read_mixed_restart(path, state)
    character(len=*), intent(in) :: path
    type(mixed_restart_state), intent(inout) :: state
    integer :: unit, i, bead, index_value
    character(len=64) :: key
    character(len=32) :: charge_policy

    open(newunit=unit, file=trim(path), status="old", action="read")
    read(unit,*)
    read(unit,*) key, state%step
    if (trim(key) /= "step") error stop "read_mixed_restart: expected step"
    read(unit,*) key, state%ndim
    if (trim(key) /= "ndim") error stop "read_mixed_restart: expected ndim"
    read(unit,*) key, state%nparticles
    if (trim(key) /= "nparticles") error stop "read_mixed_restart: expected nparticles"
    read(unit,*) key, state%beads
    if (trim(key) /= "beads") error stop "read_mixed_restart: expected beads"
    read(unit,*) key, charge_policy
    if (trim(key) /= "charge_policy") error stop "read_mixed_restart: expected charge_policy"
    state%charge_policy = trim(charge_policy)
    if (state%ndim < 1 .or. state%nparticles < 1 .or. state%beads < 1) error stop "read_mixed_restart: invalid shape"

    call allocate_state(state)
    read(unit,*) key
    if (trim(key) /= "positions_A") error stop "read_mixed_restart: expected positions_A"
    do bead = 1, state%beads
      read(unit,*) key, index_value
      if (trim(key) /= "bead" .or. index_value /= bead) error stop "read_mixed_restart: expected bead section"
      do i = 1, state%nparticles
        read(unit,*) index_value, state%positions_A(:, i, bead)
        if (index_value /= i) error stop "read_mixed_restart: position index mismatch"
      end do
    end do

    read(unit,*) key
    if (trim(key) /= "charges_by_bead") error stop "read_mixed_restart: expected charges_by_bead"
    do bead = 1, state%beads
      read(unit,*) key, index_value
      if (trim(key) /= "bead" .or. index_value /= bead) error stop "read_mixed_restart: expected charge bead section"
      do i = 1, state%nparticles
        read(unit,*) index_value, state%charges_by_bead(i, bead), state%solver_guess_by_bead(i, bead)
        if (index_value /= i) error stop "read_mixed_restart: charge index mismatch"
      end do
    end do
    close(unit)
  end subroutine read_mixed_restart

  subroutine write_charge_by_bead_summary(path, charges_by_bead)
    character(len=*), intent(in) :: path
    real(dp), intent(in) :: charges_by_bead(:, :)
    integer :: unit, bead

    open(newunit=unit, file=trim(path), status="replace", action="write")
    write(unit,'(a)') "# bead total_charge max_abs_charge"
    do bead = 1, size(charges_by_bead, 2)
      write(unit,'(i8,2(1x,es24.16))') bead, sum(charges_by_bead(:, bead)), maxval(abs(charges_by_bead(:, bead)))
    end do
    close(unit)
  end subroutine write_charge_by_bead_summary

  subroutine write_charge_centroid_summary(path, charges_by_bead, solver_guess_by_bead)
    character(len=*), intent(in) :: path
    real(dp), intent(in) :: charges_by_bead(:, :)
    real(dp), intent(in) :: solver_guess_by_bead(:, :)
    integer :: unit, i
    real(dp), allocatable :: centroid_charges(:), centroid_solver(:)

    allocate(centroid_charges(size(charges_by_bead, 1)))
    allocate(centroid_solver(size(solver_guess_by_bead, 1)))
    centroid_charges = sum(charges_by_bead, dim=2) / real(size(charges_by_bead, 2), dp)
    centroid_solver = sum(solver_guess_by_bead, dim=2) / real(size(solver_guess_by_bead, 2), dp)

    open(newunit=unit, file=trim(path), status="replace", action="write")
    write(unit,'(a)') "# atom centroid_charge centroid_solver_guess"
    do i = 1, size(centroid_charges)
      write(unit,'(i8,2(1x,es24.16))') i, centroid_charges(i), centroid_solver(i)
    end do
    close(unit)
  end subroutine write_charge_centroid_summary

  subroutine write_electron_localization_summary(path, step, local_charge_diff, centroid_charge_diff, max_force, energy_eV)
    character(len=*), intent(in) :: path
    integer, intent(in) :: step
    real(dp), intent(in) :: local_charge_diff, centroid_charge_diff, max_force, energy_eV
    integer :: unit

    open(newunit=unit, file=trim(path), status="replace", action="write")
    write(unit,'(a)') "# step local_charge_diff centroid_charge_diff max_force_eV_per_A energy_eV"
    write(unit,'(i8,4(1x,es24.16))') step, local_charge_diff, centroid_charge_diff, max_force, energy_eV
    close(unit)
  end subroutine write_electron_localization_summary

  subroutine allocate_state(state)
    type(mixed_restart_state), intent(inout) :: state

    if (allocated(state%positions_A)) deallocate(state%positions_A)
    if (allocated(state%charges_by_bead)) deallocate(state%charges_by_bead)
    if (allocated(state%solver_guess_by_bead)) deallocate(state%solver_guess_by_bead)
    allocate(state%positions_A(state%ndim, state%nparticles, state%beads))
    allocate(state%charges_by_bead(state%nparticles, state%beads))
    allocate(state%solver_guess_by_bead(state%nparticles, state%beads))
    state%positions_A = 0.0_dp
    state%charges_by_bead = 0.0_dp
    state%solver_guess_by_bead = 0.0_dp
  end subroutine allocate_state

end module ereaxff_mixed_restart
