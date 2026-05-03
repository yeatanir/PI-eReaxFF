module ereaxff_reaxff_compat_output
  use ereaxff_kinds, only: dp
  use ereaxff_system_state, only: system_state, cell_is_defined, minimum_image_displacement
  implicit none
  private

  public :: write_reaxff_compat_snapshot

  real(dp), parameter :: bond_order_cutoff = 0.3000_dp
  real(dp), parameter :: raw_bond_scale = 1.18_dp
  real(dp), parameter :: corrected_bond_scale = 1.12_dp

contains

  subroutine write_reaxff_compat_snapshot(output_dir, step, centroid_positions_A, centroid_charges, energy_eV, &
      max_force_eV_per_A, charge_sum, beads, state)
    character(len=*), intent(in) :: output_dir
    integer, intent(in) :: step, beads
    real(dp), intent(in) :: centroid_positions_A(:, :)
    real(dp), intent(in) :: centroid_charges(:)
    real(dp), intent(in) :: energy_eV, max_force_eV_per_A, charge_sum
    type(system_state), intent(in) :: state
    character(len=:), allocatable :: fort7_path, fort8_path, molfra_path, xmol_path
    integer, allocatable :: raw_components(:), corrected_components(:)
    integer :: raw_ncomponents, corrected_ncomponents

    fort7_path = trim(output_dir) // "/fort.7"
    fort8_path = trim(output_dir) // "/fort.8"
    molfra_path = trim(output_dir) // "/molfra.out"
    xmol_path = trim(output_dir) // "/xmolout"

    call build_components(state, centroid_positions_A, raw_bond_scale, raw_components, raw_ncomponents)
    call build_components(state, centroid_positions_A, corrected_bond_scale, corrected_components, &
        corrected_ncomponents)

    call write_connection_table(fort7_path, "fort.7 centroid-projected connection table (raw)", step, state, &
        centroid_positions_A, centroid_charges, raw_components, raw_bond_scale)
    call write_connection_table(fort8_path, "fort.8 centroid-projected connection table (corrected)", step, state, &
        centroid_positions_A, centroid_charges, corrected_components, corrected_bond_scale)
    call write_molecule_list(molfra_path, step, state, centroid_positions_A, corrected_components, &
        corrected_ncomponents)
    call write_xmolout(xmol_path, step, state, centroid_positions_A, centroid_charges)
  end subroutine write_reaxff_compat_snapshot

  subroutine write_connection_table(path, title, step, state, positions_A, charges, components, bond_scale)
    character(len=*), intent(in) :: path, title
    integer, intent(in) :: step
    type(system_state), intent(in) :: state
    real(dp), intent(in) :: positions_A(:, :)
    real(dp), intent(in) :: charges(:)
    integer, intent(in) :: components(:)
    real(dp), intent(in) :: bond_scale
    integer :: unit, i, j, nn, molecule_id, bond_count
    real(dp) :: total_bo, bond_order
    integer, allocatable :: neighbors(:)
    real(dp), allocatable :: bond_orders(:)

    bond_count = count_bonds(state, positions_A, bond_scale)
    open(newunit=unit, file=trim(path), status="replace", action="write")
    write(unit,'(a)') "# " // trim(title)
    if (cell_is_defined(state)) then
      write(unit,'(a,1x,l1,1x,a,1x,l1,1x,a,6(1x,f12.5))') "# cell_defined", cell_is_defined(state), &
          "# periodic", state%cell%periodic, "# crystx", state%cell%crystx
    else
      write(unit,'(a)') "# cell_defined F"
    end if
    write(unit,'(i6,1x,a,1x,i0,1x,a,1x,i0)') state%nparticles, "Iteration:", step, "#Bonds:", bond_count
    write(unit,'(a)') "# idx element molecule_id x_A y_A z_A charge total_bo n_neighbors neighbors..."

    do i = 1, state%nparticles
      allocate(neighbors(state%nparticles - 1))
      allocate(bond_orders(state%nparticles - 1))
      nn = 0
      total_bo = 0.0_dp
      molecule_id = components(i)
      do j = 1, state%nparticles
        if (i == j) cycle
        bond_order = bond_order_proxy(state%element(i), state%element(j), &
            pair_distance(state, positions_A(:, i), positions_A(:, j)), bond_scale)
        if (bond_order >= bond_order_cutoff) then
          nn = nn + 1
          neighbors(nn) = j
          bond_orders(nn) = bond_order
          total_bo = total_bo + bond_order
        end if
      end do
      write(unit,'(i8,1x,a,1x,i8,1x,3(1x,es16.8),1x,es16.8,1x,es16.8,1x,i8)', advance="no") &
          i, trim(adjustl(state%element(i))), molecule_id, positions_A(:, i), charges(i), total_bo, nn
      do j = 1, nn
        write(unit,'(1x,i8,1x,es16.8)', advance="no") neighbors(j), bond_orders(j)
      end do
      write(unit,*)
      deallocate(neighbors)
      deallocate(bond_orders)
    end do

    write(unit,'(a,1x,es16.8,1x,a,1x,es16.8,1x,a,1x,es16.8)') "# energy_eV", 0.0_dp, &
        "max_force_eV_per_A", 0.0_dp, "charge_sum", sum(charges)
    close(unit)
  end subroutine write_connection_table

  subroutine write_molecule_list(path, step, state, positions_A, components, ncomponents)
    character(len=*), intent(in) :: path
    integer, intent(in) :: step, ncomponents
    type(system_state), intent(in) :: state
    real(dp), intent(in) :: positions_A(:, :)
    integer, intent(in) :: components(:)
    integer :: unit, comp, i, natoms, unique_count, pos
    integer, allocatable :: members(:), frequencies(:)
    character(len=128), allocatable :: formulas(:)
    real(dp), allocatable :: masses(:)

    allocate(formulas(ncomponents))
    allocate(masses(ncomponents))
    allocate(frequencies(ncomponents))
    formulas = ""
    masses = 0.0_dp
    frequencies = 0

    do comp = 1, ncomponents
      natoms = count(components == comp)
      allocate(members(natoms))
      pos = 0
      do i = 1, state%nparticles
        if (components(i) == comp) then
          pos = pos + 1
          members(pos) = i
        end if
      end do
      formulas(comp) = fragment_formula(state, members(:pos))
      masses(comp) = fragment_mass(state, members(:pos))
      deallocate(members)
    end do

    unique_count = 0
    do comp = 1, ncomponents
      if (frequencies(comp) /= 0) cycle
      unique_count = unique_count + 1
      frequencies(comp) = 1
      do i = comp + 1, ncomponents
        if (frequencies(i) == 0 .and. trim(formulas(i)) == trim(formulas(comp))) then
          frequencies(comp) = frequencies(comp) + 1
          frequencies(i) = -1
        end if
      end do
    end do

    open(newunit=unit, file=trim(path), status="replace", action="write")
    write(unit,'(a)') "# molecule list per frame"
    if (cell_is_defined(state)) then
      write(unit,'(a,1x,l1,1x,a,1x,l1,1x,a,6(1x,f12.5))') "# cell_defined", cell_is_defined(state), &
          "# periodic", state%cell%periodic, "# crystx", state%cell%crystx
    else
      write(unit,'(a)') "# cell_defined F"
    end if
    write(unit,'(a,f8.4)') "Bond order cutoff:", bond_order_cutoff
    write(unit,'(a)') "Iteration Freq. Molecular formula               Molecular mass"
    do comp = 1, ncomponents
      if (frequencies(comp) <= 0) cycle
      write(unit,'(i8,3x,i2,1x,a2,1x,a32,1x,f14.4)') step, frequencies(comp), "x", trim(formulas(comp)), masses(comp)
    end do
    write(unit,'(a,i12)') " Total number of molecules:", ncomponents
    write(unit,'(a,i12)') " Total number of atoms:", state%nparticles
    write(unit,'(a,f20.12)') " Total system mass:", sum(state%mass_amu)
    write(unit,'(a)') "# molecule list"
    do comp = 1, ncomponents
      if (frequencies(comp) <= 0) cycle
      call write_component_line(unit, step, comp, state, formulas, masses, components)
    end do
    close(unit)
    deallocate(formulas)
    deallocate(masses)
    deallocate(frequencies)
  end subroutine write_molecule_list

  subroutine write_component_line(unit, step, comp, state, formulas, masses, components)
    integer, intent(in) :: unit, step, comp
    type(system_state), intent(in) :: state
    character(len=*), intent(in) :: formulas(:)
    real(dp), intent(in) :: masses(:)
    integer, intent(in) :: components(:)
    integer :: i, natoms, seen_count

    natoms = count(components == comp)
    if (natoms < 1) return
    write(unit,'(a,1x,i0,1x,a,1x,i0,1x,a,1x,a,1x,f12.4,1x,a)', advance="no") &
        "molecule", comp, "step", step, "label=", trim(formulas(comp)), masses(comp), "atoms="
    seen_count = 0
    do i = 1, state%nparticles
      if (components(i) /= comp) cycle
      seen_count = seen_count + 1
      write(unit,'(i0)', advance="no") i
      if (seen_count < natoms) write(unit,'(a)', advance="no") ","
    end do
    write(unit,*)
  end subroutine write_component_line

  subroutine write_xmolout(path, step, state, positions_A, charges)
    character(len=*), intent(in) :: path
    integer, intent(in) :: step
    type(system_state), intent(in) :: state
    real(dp), intent(in) :: positions_A(:, :)
    real(dp), intent(in) :: charges(:)
    integer :: unit, i

    open(newunit=unit, file=trim(path), status="replace", action="write")
    write(unit,'(i0)') size(charges)
    write(unit,'(a,i0)') "step=", step
    do i = 1, size(charges)
      write(unit,'(a,3(1x,es22.14),1x,es22.14)') trim(state%element(i)), positions_A(:, i), charges(i)
    end do
    close(unit)
  end subroutine write_xmolout

  subroutine build_components(state, positions_A, bond_scale, components, ncomponents)
    type(system_state), intent(in) :: state
    real(dp), intent(in) :: positions_A(:, :)
    real(dp), intent(in) :: bond_scale
    integer, allocatable, intent(out) :: components(:)
    integer, intent(out) :: ncomponents
    logical, allocatable :: visited(:)
    integer, allocatable :: stack(:)
    integer :: atom, top, current, other

    allocate(components(state%nparticles))
    allocate(visited(state%nparticles))
    allocate(stack(state%nparticles))
    components = 0
    visited = .false.
    ncomponents = 0

    do atom = 1, state%nparticles
      if (visited(atom)) cycle
      ncomponents = ncomponents + 1
      top = 1
      stack(top) = atom
      visited(atom) = .true.
      components(atom) = ncomponents
      do while (top > 0)
        current = stack(top)
        top = top - 1
        do other = 1, state%nparticles
          if (visited(other)) cycle
          if (bond_order_proxy(state%element(current), state%element(other), &
              pair_distance(state, positions_A(:, current), positions_A(:, other)), bond_scale) &
              >= bond_order_cutoff) then
            top = top + 1
            stack(top) = other
            visited(other) = .true.
            components(other) = ncomponents
          end if
        end do
      end do
    end do

    deallocate(visited)
    deallocate(stack)
  end subroutine build_components

  pure real(dp) function bond_order_proxy(element_i, element_j, distance, bond_scale) result(bond_order)
    character(len=*), intent(in) :: element_i, element_j
    real(dp), intent(in) :: distance, bond_scale
    real(dp) :: cutoff, scale

    if (trim(element_i) == "El" .or. trim(element_j) == "El") then
      bond_order = 0.0_dp
      return
    end if

    cutoff = bond_scale * (covalent_radius_A(element_i) + covalent_radius_A(element_j))
    if (distance > cutoff) then
      bond_order = 0.0_dp
    else
      scale = max(0.2_dp, covalent_radius_A(element_i) + covalent_radius_A(element_j))
      bond_order = exp(-3.0_dp * max(0.0_dp, distance / scale - 1.0_dp))
    end if
  end function bond_order_proxy

  pure real(dp) function pair_distance(state, a, b) result(distance)
    type(system_state), intent(in) :: state
    real(dp), intent(in) :: a(:), b(:)
    real(dp) :: delta(size(a))

    call minimum_image_displacement(state, a, b, delta)
    distance = sqrt(sum(delta**2))
  end function pair_distance

  pure real(dp) function covalent_radius_A(element) result(radius)
    character(len=*), intent(in) :: element

    select case (trim(element))
    case ("H")
      radius = 0.31_dp
    case ("C")
      radius = 0.76_dp
    case ("N")
      radius = 0.71_dp
    case ("O")
      radius = 0.66_dp
    case ("El")
      radius = 0.0_dp
    case default
      radius = 0.75_dp
    end select
  end function covalent_radius_A

  pure real(dp) function fragment_mass(state, members) result(mass)
    type(system_state), intent(in) :: state
    integer, intent(in) :: members(:)
    integer :: i

    mass = 0.0_dp
    do i = 1, size(members)
      mass = mass + state%mass_amu(members(i))
    end do
  end function fragment_mass

  character(len=128) function fragment_formula(state, members) result(formula)
    type(system_state), intent(in) :: state
    integer, intent(in) :: members(:)
    integer :: i, count_c, count_el

    count_c = 0
    count_el = 0
    do i = 1, size(members)
      select case (trim(state%element(members(i))))
      case ("C")
        count_c = count_c + 1
      case ("El")
        count_el = count_el + 1
      case default
        count_el = count_el + 0
      end select
    end do

    formula = ""
    if (count_c > 0) formula = trim(formula) // "C" // trim(int_to_string(count_c))
    if (count_el > 0) formula = trim(formula) // "El" // trim(int_to_string(count_el))
    if (len_trim(formula) == 0) formula = "fragment"
  end function fragment_formula

  pure integer function count_bonds(state, positions_A, bond_scale) result(nbonds)
    type(system_state), intent(in) :: state
    real(dp), intent(in) :: positions_A(:, :)
    real(dp), intent(in) :: bond_scale
    integer :: i, j

    nbonds = 0
    do i = 1, state%nparticles - 1
      do j = i + 1, state%nparticles
        if (bond_order_proxy(state%element(i), state%element(j), &
            pair_distance(state, positions_A(:, i), positions_A(:, j)), bond_scale) >= bond_order_cutoff) then
          nbonds = nbonds + 1
        end if
      end do
    end do
  end function count_bonds

  character(len=32) function int_to_string(i) result(text)
    integer, intent(in) :: i
    write(text,'(i0)') i
  end function int_to_string

end module ereaxff_reaxff_compat_output
