module ereaxff_selectors
  use ereaxff_system_state, only: system_state
  implicit none
  private

  public :: selector_matches
  public :: selector_is_supported

contains

  logical function selector_matches(state, particle_index, selector) result(matches)
    type(system_state), intent(in) :: state
    integer, intent(in) :: particle_index
    character(len=*), intent(in) :: selector
    character(len=:), allocatable :: compact
    integer :: and_pos
    logical :: lhs_matches, rhs_matches

    if (particle_index < 1 .or. particle_index > state%nparticles) then
      error stop "selector_matches: particle index out of range"
    end if

    compact = remove_spaces(selector)
    if (compact == "all") then
      matches = .true.
    else
      and_pos = index(compact, "and")
      if (and_pos > 0) then
        lhs_matches = match_clause(state, particle_index, compact(:and_pos - 1))
        rhs_matches = match_clause(state, particle_index, compact(and_pos + 3:))
        matches = lhs_matches .and. rhs_matches
      else
        matches = match_clause(state, particle_index, compact)
      end if
    end if
  end function selector_matches

  logical function selector_is_supported(selector) result(is_supported)
    character(len=*), intent(in) :: selector
    character(len=:), allocatable :: compact
    integer :: and_pos

    compact = remove_spaces(selector)
    if (compact == "all") then
      is_supported = .true.
      return
    end if

    and_pos = index(compact, "and")
    if (and_pos > 0) then
      is_supported = clause_is_supported(compact(:and_pos - 1)) .and. &
          clause_is_supported(compact(and_pos + 3:))
    else
      is_supported = clause_is_supported(compact)
    end if
  end function selector_is_supported

  logical function match_clause(state, particle_index, clause) result(matches)
    type(system_state), intent(in) :: state
    integer, intent(in) :: particle_index
    character(len=*), intent(in) :: clause
    character(len=:), allocatable :: rhs
    integer :: value

    if (index(clause, "element==") == 1) then
      rhs = clause(10:)
      matches = trim(state%element(particle_index)) == rhs
    else if (index(clause, "id==") == 1) then
      read(clause(5:), *) value
      matches = state%id(particle_index) == value
    else if (index(clause, "id!=") == 1) then
      read(clause(5:), *) value
      matches = state%id(particle_index) /= value
    else if (index(clause, "idin[") == 1 .and. clause(len_trim(clause):len_trim(clause)) == "]") then
      matches = id_in_list(state%id(particle_index), clause(6:len_trim(clause) - 1))
    else
      error stop "selector_matches: unsupported selector clause"
    end if
  end function match_clause

  logical function clause_is_supported(clause) result(is_supported)
    character(len=*), intent(in) :: clause
    integer :: value, ios

    is_supported = .false.
    if (index(clause, "element==") == 1 .and. len_trim(clause) > 9) then
      is_supported = .true.
    else if (index(clause, "id==") == 1) then
      read(clause(5:), *, iostat=ios) value
      is_supported = ios == 0
    else if (index(clause, "id!=") == 1) then
      read(clause(5:), *, iostat=ios) value
      is_supported = ios == 0
    else if (index(clause, "idin[") == 1 .and. clause(len_trim(clause):len_trim(clause)) == "]") then
      is_supported = len_trim(clause) > 6
    end if
  end function clause_is_supported

  logical function id_in_list(id, list_text) result(found)
    integer, intent(in) :: id
    character(len=*), intent(in) :: list_text
    integer :: start, comma, value, ios
    character(len=:), allocatable :: token

    found = .false.
    start = 1
    do
      comma = index(list_text(start:), ",")
      if (comma == 0) then
        token = list_text(start:)
      else
        token = list_text(start:start + comma - 2)
      end if
      read(token, *, iostat=ios) value
      if (ios /= 0) error stop "id_in_list: invalid id list"
      if (id == value) found = .true.
      if (comma == 0) exit
      start = start + comma
    end do
  end function id_in_list

  pure function remove_spaces(text) result(compact)
    character(len=*), intent(in) :: text
    character(len=:), allocatable :: compact
    integer :: i

    compact = ""
    do i = 1, len_trim(text)
      if (text(i:i) /= " " .and. text(i:i) /= char(9)) compact = compact // text(i:i)
    end do
  end function remove_spaces

end module ereaxff_selectors
