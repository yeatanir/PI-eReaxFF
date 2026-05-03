module ereaxff_mixed_baoab_restart
  use ereaxff_kinds, only: dp, i8
  implicit none
  private

  type, public :: mixed_baoab_restart_state
    integer :: step = 0
    integer :: ndim = 0
    integer :: nparticles = 0
    integer :: nquantum = 0
    integer :: beads = 0
    character(len=32) :: charge_policy = ""
    real(dp) :: classical_heat_eV = 0.0_dp
    real(dp), allocatable :: quantum_heat_eV(:)
    integer(i8) :: rng_state_i8 = 0_i8
    real(dp), allocatable :: positions_A(:, :, :)
    real(dp), allocatable :: classical_velocities(:, :)
    integer, allocatable :: quantum_indices(:)
    real(dp), allocatable :: quantum_masses_amu(:)
    real(dp), allocatable :: q_nm(:, :, :)
    real(dp), allocatable :: v_nm(:, :, :)
    real(dp), allocatable :: charges_by_bead(:, :)
    real(dp), allocatable :: solver_guess_by_bead(:, :)
  end type mixed_baoab_restart_state

  public :: write_mixed_baoab_restart
  public :: read_mixed_baoab_restart

contains

  subroutine write_mixed_baoab_restart(path, state)
    character(len=*), intent(in) :: path
    type(mixed_baoab_restart_state), intent(in) :: state
    integer :: unit, bead, i, coord, q

    if (.not. allocated(state%positions_A)) error stop "write_mixed_baoab_restart: positions not allocated"
    if (.not. allocated(state%quantum_indices)) error stop "write_mixed_baoab_restart: quantum indices not allocated"
    if (.not. allocated(state%quantum_masses_amu)) error stop "write_mixed_baoab_restart: quantum masses not allocated"
    if (.not. allocated(state%quantum_heat_eV)) error stop "write_mixed_baoab_restart: quantum heat not allocated"
    if (.not. allocated(state%q_nm)) error stop "write_mixed_baoab_restart: q_nm not allocated"
    if (.not. allocated(state%v_nm)) error stop "write_mixed_baoab_restart: v_nm not allocated"
    if (.not. allocated(state%charges_by_bead)) error stop "write_mixed_baoab_restart: charges not allocated"
    if (.not. allocated(state%solver_guess_by_bead)) error stop "write_mixed_baoab_restart: solver memory not allocated"
    if (size(state%quantum_indices) /= state%nquantum) error stop "write_mixed_baoab_restart: nquantum mismatch"
    if (size(state%quantum_masses_amu) /= state%nquantum) error stop "write_mixed_baoab_restart: mass count mismatch"
    if (size(state%quantum_heat_eV) /= state%nquantum) error stop "write_mixed_baoab_restart: heat count mismatch"

    open(newunit=unit, file=trim(path), status="replace", action="write")
    write(unit,'(a)') "# PI-eReaxFF mixed BAOAB restart v2"
    write(unit,'(a,1x,i0)') "step", state%step
    write(unit,'(a,1x,i0)') "ndim", state%ndim
    write(unit,'(a,1x,i0)') "nparticles", state%nparticles
    write(unit,'(a,1x,i0)') "nquantum", state%nquantum
    write(unit,'(a,1x,i0)') "beads", state%beads
    write(unit,'(a,1x,a)') "charge_policy", trim(state%charge_policy)
    write(unit,'(a,1x,es24.16)') "classical_heat_eV", state%classical_heat_eV
    write(unit,'(a,1x,i0)') "rng_state_i8", state%rng_state_i8
    write(unit,'(a)') "quantum_indices"
    do q = 1, state%nquantum
      write(unit,'(i8,1x,i0,1x,es24.16,1x,es24.16)') q, state%quantum_indices(q), &
          state%quantum_masses_amu(q), state%quantum_heat_eV(q)
    end do
    write(unit,'(a)') "positions_A"
    do bead = 1, state%beads
      write(unit,'(a,1x,i0)') "bead", bead
      do i = 1, state%nparticles
        write(unit,'(i8,3(1x,es24.16))') i, state%positions_A(:, i, bead)
      end do
    end do
    write(unit,'(a)') "classical_velocities"
    do i = 1, state%nparticles
      write(unit,'(i8,3(1x,es24.16))') i, state%classical_velocities(:, i)
    end do
    write(unit,'(a)') "q_nm"
    do q = 1, state%nquantum
      write(unit,'(a,1x,i0)') "quantum", q
      write(unit,'(a,1x,i0)') "particle", state%quantum_indices(q)
      do coord = 1, state%ndim
        write(unit,'(a,1x,i0)') "coord", coord
        do bead = 1, state%beads
          write(unit,'(i8,1x,es24.16)') bead, state%q_nm(coord, q, bead)
        end do
      end do
    end do
    write(unit,'(a)') "v_nm"
    do q = 1, state%nquantum
      write(unit,'(a,1x,i0)') "quantum", q
      write(unit,'(a,1x,i0)') "particle", state%quantum_indices(q)
      do coord = 1, state%ndim
        write(unit,'(a,1x,i0)') "coord", coord
        do bead = 1, state%beads
          write(unit,'(i8,1x,es24.16)') bead, state%v_nm(coord, q, bead)
        end do
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
  end subroutine write_mixed_baoab_restart

  subroutine read_mixed_baoab_restart(path, state)
    character(len=*), intent(in) :: path
    type(mixed_baoab_restart_state), intent(inout) :: state
    integer :: unit, i, bead, coord, index_value, q, quantum_index
    character(len=64) :: key
    character(len=32) :: charge_policy

    open(newunit=unit, file=trim(path), status="old", action="read")
    read(unit,*)
    read(unit,*) key, state%step
    if (trim(key) /= "step") error stop "read_mixed_baoab_restart: expected step"
    read(unit,*) key, state%ndim
    if (trim(key) /= "ndim") error stop "read_mixed_baoab_restart: expected ndim"
    read(unit,*) key, state%nparticles
    if (trim(key) /= "nparticles") error stop "read_mixed_baoab_restart: expected nparticles"
    read(unit,*) key, state%nquantum
    if (trim(key) /= "nquantum") error stop "read_mixed_baoab_restart: expected nquantum"
    read(unit,*) key, state%beads
    if (trim(key) /= "beads") error stop "read_mixed_baoab_restart: expected beads"
    read(unit,*) key, charge_policy
    if (trim(key) /= "charge_policy") error stop "read_mixed_baoab_restart: expected charge_policy"
    state%charge_policy = trim(charge_policy)
    read(unit,*) key, state%classical_heat_eV
    if (trim(key) /= "classical_heat_eV") error stop "read_mixed_baoab_restart: expected classical_heat_eV"
    read(unit,*) key, state%rng_state_i8
    if (trim(key) /= "rng_state_i8") error stop "read_mixed_baoab_restart: expected rng_state_i8"
    call allocate_state(state)

    read(unit,*) key
    if (trim(key) /= "quantum_indices") error stop "read_mixed_baoab_restart: expected quantum_indices"
    do q = 1, state%nquantum
      read(unit,*) index_value, state%quantum_indices(q), state%quantum_masses_amu(q), state%quantum_heat_eV(q)
      if (index_value /= q) error stop "read_mixed_baoab_restart: quantum index mismatch"
    end do

    read(unit,*) key
    if (trim(key) /= "positions_A") error stop "read_mixed_baoab_restart: expected positions_A"
    do bead = 1, state%beads
      read(unit,*) key, index_value
      if (trim(key) /= "bead" .or. index_value /= bead) error stop "read_mixed_baoab_restart: expected bead section"
      do i = 1, state%nparticles
        read(unit,*) index_value, state%positions_A(:, i, bead)
        if (index_value /= i) error stop "read_mixed_baoab_restart: position index mismatch"
      end do
    end do

    read(unit,*) key
    if (trim(key) /= "classical_velocities") error stop "read_mixed_baoab_restart: expected classical_velocities"
    do i = 1, state%nparticles
      read(unit,*) index_value, state%classical_velocities(:, i)
      if (index_value /= i) error stop "read_mixed_baoab_restart: velocity index mismatch"
    end do

    read(unit,*) key
    if (trim(key) /= "q_nm") error stop "read_mixed_baoab_restart: expected q_nm"
    do q = 1, state%nquantum
      read(unit,*) key, index_value
      if (trim(key) /= "quantum" .or. index_value /= q) error stop "read_mixed_baoab_restart: expected q_nm quantum"
      read(unit,*) key, quantum_index
      if (trim(key) /= "particle" .or. quantum_index /= state%quantum_indices(q)) &
          error stop "read_mixed_baoab_restart: q_nm particle mismatch"
      do coord = 1, state%ndim
        read(unit,*) key, index_value
        if (trim(key) /= "coord" .or. index_value /= coord) error stop "read_mixed_baoab_restart: expected q_nm coord"
        do bead = 1, state%beads
          read(unit,*) index_value, state%q_nm(coord, q, bead)
          if (index_value /= bead) error stop "read_mixed_baoab_restart: q_nm index mismatch"
        end do
      end do
    end do

    read(unit,*) key
    if (trim(key) /= "v_nm") error stop "read_mixed_baoab_restart: expected v_nm"
    do q = 1, state%nquantum
      read(unit,*) key, index_value
      if (trim(key) /= "quantum" .or. index_value /= q) error stop "read_mixed_baoab_restart: expected v_nm quantum"
      read(unit,*) key, quantum_index
      if (trim(key) /= "particle" .or. quantum_index /= state%quantum_indices(q)) &
          error stop "read_mixed_baoab_restart: v_nm particle mismatch"
      do coord = 1, state%ndim
        read(unit,*) key, index_value
        if (trim(key) /= "coord" .or. index_value /= coord) error stop "read_mixed_baoab_restart: expected v_nm coord"
        do bead = 1, state%beads
          read(unit,*) index_value, state%v_nm(coord, q, bead)
          if (index_value /= bead) error stop "read_mixed_baoab_restart: v_nm index mismatch"
        end do
      end do
    end do

    read(unit,*) key
    if (trim(key) /= "charges_by_bead") error stop "read_mixed_baoab_restart: expected charges_by_bead"
    do bead = 1, state%beads
      read(unit,*) key, index_value
      if (trim(key) /= "bead" .or. index_value /= bead) error stop "read_mixed_baoab_restart: expected charge bead section"
      do i = 1, state%nparticles
        read(unit,*) index_value, state%charges_by_bead(i, bead), state%solver_guess_by_bead(i, bead)
        if (index_value /= i) error stop "read_mixed_baoab_restart: charge index mismatch"
      end do
    end do
    close(unit)
  end subroutine read_mixed_baoab_restart

  subroutine allocate_state(state)
    type(mixed_baoab_restart_state), intent(inout) :: state

    if (allocated(state%positions_A)) deallocate(state%positions_A)
    if (allocated(state%classical_velocities)) deallocate(state%classical_velocities)
    if (allocated(state%quantum_indices)) deallocate(state%quantum_indices)
    if (allocated(state%quantum_masses_amu)) deallocate(state%quantum_masses_amu)
    if (allocated(state%quantum_heat_eV)) deallocate(state%quantum_heat_eV)
    if (allocated(state%q_nm)) deallocate(state%q_nm)
    if (allocated(state%v_nm)) deallocate(state%v_nm)
    if (allocated(state%charges_by_bead)) deallocate(state%charges_by_bead)
    if (allocated(state%solver_guess_by_bead)) deallocate(state%solver_guess_by_bead)
    allocate(state%positions_A(state%ndim, state%nparticles, state%beads))
    allocate(state%classical_velocities(state%ndim, state%nparticles))
    allocate(state%quantum_indices(state%nquantum))
    allocate(state%quantum_masses_amu(state%nquantum))
    allocate(state%quantum_heat_eV(state%nquantum))
    allocate(state%q_nm(state%ndim, state%nquantum, state%beads))
    allocate(state%v_nm(state%ndim, state%nquantum, state%beads))
    allocate(state%charges_by_bead(state%nparticles, state%beads))
    allocate(state%solver_guess_by_bead(state%nparticles, state%beads))
    state%positions_A = 0.0_dp
    state%classical_velocities = 0.0_dp
    state%quantum_indices = 0
    state%quantum_masses_amu = 0.0_dp
    state%quantum_heat_eV = 0.0_dp
    state%q_nm = 0.0_dp
    state%v_nm = 0.0_dp
    state%charges_by_bead = 0.0_dp
    state%solver_guess_by_bead = 0.0_dp
  end subroutine allocate_state

end module ereaxff_mixed_baoab_restart
