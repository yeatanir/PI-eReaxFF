module ereaxff_restart
  use ereaxff_kinds, only: dp, i8
  implicit none
  private

  type, public :: harmonic_restart_state
    integer :: step = 0
    integer :: beads = 0
    integer(i8) :: rng_state_i8 = 1_i8
    real(dp) :: thermostat_heat_eV = 0.0_dp
    real(dp), allocatable :: q_nm(:)
    real(dp), allocatable :: v_nm(:)
  end type harmonic_restart_state

  public :: write_harmonic_restart
  public :: read_harmonic_restart

contains

  subroutine write_harmonic_restart(path, step, q_nm, v_nm, thermostat_heat_eV, rng_state_i8)
    character(len=*), intent(in) :: path
    integer, intent(in) :: step
    real(dp), intent(in) :: q_nm(:), v_nm(:)
    real(dp), intent(in) :: thermostat_heat_eV
    integer(i8), intent(in) :: rng_state_i8
    integer :: unit, i

    if (size(q_nm) /= size(v_nm)) error stop "write_harmonic_restart: q/v size mismatch"

    open(newunit=unit, file=trim(path), status="replace", action="write")
    write(unit,'(a)') "# PI-eReaxFF harmonic restart v1"
    write(unit,'(a,1x,i0)') "step", step
    write(unit,'(a,1x,i0)') "beads", size(q_nm)
    write(unit,'(a,1x,i0)') "rng_state", rng_state_i8
    write(unit,'(a,1x,es24.16)') "thermostat_heat_eV", thermostat_heat_eV
    write(unit,'(a)') "q_nm"
    do i = 1, size(q_nm)
      write(unit,'(i8,1x,es24.16)') i, q_nm(i)
    end do
    write(unit,'(a)') "v_nm"
    do i = 1, size(v_nm)
      write(unit,'(i8,1x,es24.16)') i, v_nm(i)
    end do
    close(unit)
  end subroutine write_harmonic_restart

  subroutine read_harmonic_restart(path, restart)
    character(len=*), intent(in) :: path
    type(harmonic_restart_state), intent(inout) :: restart
    integer :: unit, i, index_value
    character(len=64) :: key

    open(newunit=unit, file=trim(path), status="old", action="read")
    read(unit,*)
    read(unit,*) key, restart%step
    if (trim(key) /= "step") error stop "read_harmonic_restart: expected step"
    read(unit,*) key, restart%beads
    if (trim(key) /= "beads") error stop "read_harmonic_restart: expected beads"
    read(unit,*) key, restart%rng_state_i8
    if (trim(key) /= "rng_state") error stop "read_harmonic_restart: expected rng_state"
    read(unit,*) key, restart%thermostat_heat_eV
    if (trim(key) /= "thermostat_heat_eV") error stop "read_harmonic_restart: expected thermostat heat"
    if (restart%beads < 1) error stop "read_harmonic_restart: invalid bead count"

    if (allocated(restart%q_nm)) deallocate(restart%q_nm)
    if (allocated(restart%v_nm)) deallocate(restart%v_nm)
    allocate(restart%q_nm(restart%beads))
    allocate(restart%v_nm(restart%beads))

    read(unit,*) key
    if (trim(key) /= "q_nm") error stop "read_harmonic_restart: expected q_nm"
    do i = 1, restart%beads
      read(unit,*) index_value, restart%q_nm(i)
      if (index_value /= i) error stop "read_harmonic_restart: q_nm index mismatch"
    end do

    read(unit,*) key
    if (trim(key) /= "v_nm") error stop "read_harmonic_restart: expected v_nm"
    do i = 1, restart%beads
      read(unit,*) index_value, restart%v_nm(i)
      if (index_value /= i) error stop "read_harmonic_restart: v_nm index mismatch"
    end do
    close(unit)
  end subroutine read_harmonic_restart

end module ereaxff_restart
