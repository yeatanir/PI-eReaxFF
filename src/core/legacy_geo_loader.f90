module ereaxff_legacy_geo_loader
  use ereaxff_kinds, only: dp
  use ereaxff_system_state, only: system_state, init_system_state, set_particle, set_cell_from_crystx
  implicit none
  private

  real(dp), parameter :: electron_mass_amu = 0.00054858_dp
  real(dp), parameter :: carbon_mass_amu = 12.011_dp

  public :: load_legacy_geo_file

contains

  subroutine load_legacy_geo_file(path, state)
    character(len=*), intent(in) :: path
    type(system_state), intent(inout) :: state
    integer :: unit, ios, natom, idx
    character(len=1024) :: line
    character(len=8) :: record, res3, res5
    character(len=5) :: atom_name, atom_type
    character(len=1) :: res1
    real(dp) :: x, y, z
    integer :: atom_id_read, dummy_int1, dummy_int2
    real(dp) :: dummy_real
    real(dp) :: crystx(6)
    logical :: found_crystx
    character(len=8) :: tag

    natom = 0
    crystx = 0.0_dp
    found_crystx = .false.
    open(newunit=unit, file=trim(path), status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "load_legacy_geo_file: unable to open geo file"
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (starts_with(trim(line), "CRYSTX")) then
        read(line, *, iostat=ios) tag, crystx
        if (ios == 0 .and. trim(tag) == "CRYSTX") found_crystx = .true.
      end if
      if (starts_with(trim(line), "HETATM")) natom = natom + 1
    end do
    close(unit)
    if (natom < 1) error stop "load_legacy_geo_file: no HETATM records found"

    call init_system_state(state, 3, natom)
    open(newunit=unit, file=trim(path), status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "load_legacy_geo_file: unable to reopen geo file"

    idx = 0
    do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (.not. starts_with(trim(line), "HETATM")) cycle
      record = ""
      atom_name = ""
      atom_type = ""
      res3 = ""
      res1 = ""
      res5 = ""
      x = 0.0_dp
      y = 0.0_dp
      z = 0.0_dp
      read(line, '(a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)', iostat=ios) &
          record, atom_id_read, atom_name, res3, res1, atom_type, x, y, z, res5, dummy_int1, dummy_int2, dummy_real
      if (ios /= 0) error stop "load_legacy_geo_file: malformed HETATM record"
      idx = idx + 1
      call set_particle(state, idx, idx, trim(atom_type), infer_mass(trim(atom_type), trim(atom_name)), [x, y, z])
    end do
    close(unit)

    if (idx /= natom) error stop "load_legacy_geo_file: particle count mismatch"
    if (found_crystx) call set_cell_from_crystx(state, crystx, .true.)
  end subroutine load_legacy_geo_file

  pure logical function starts_with(text, prefix) result(match)
    character(len=*), intent(in) :: text, prefix
    integer :: n

    n = len_trim(prefix)
    if (len_trim(text) < n) then
      match = .false.
    else
      match = text(1:n) == prefix(1:n)
    end if
  end function starts_with

  real(dp) function infer_mass(atom_type, atom_name) result(mass_amu)
    character(len=*), intent(in) :: atom_type, atom_name
    character(len=32) :: key

    key = adjustl(trim(atom_type))
    if (len_trim(key) == 0) key = adjustl(trim(atom_name))
    select case (trim(key))
    case ("C")
      mass_amu = carbon_mass_amu
    case ("El")
      mass_amu = electron_mass_amu
    case default
      error stop "infer_mass: unsupported element in legacy geo"
    end select
  end function infer_mass

end module ereaxff_legacy_geo_loader
