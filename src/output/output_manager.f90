module ereaxff_output_manager
  implicit none
  private

  public :: join_path

contains

  pure function join_path(directory, filename) result(path)
    character(len=*), intent(in) :: directory, filename
    character(len=:), allocatable :: path

    if (len_trim(directory) == 0 .or. trim(directory) == ".") then
      path = trim(filename)
    else if (directory(len_trim(directory):len_trim(directory)) == "/") then
      path = trim(directory) // trim(filename)
    else
      path = trim(directory) // "/" // trim(filename)
    end if
  end function join_path

end module ereaxff_output_manager

