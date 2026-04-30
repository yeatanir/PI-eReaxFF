module ereaxff_errors
  implicit none
  private

  public :: require
  public :: fail

contains

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop trim(message)
  end subroutine require

  subroutine fail(message)
    character(len=*), intent(in) :: message

    error stop trim(message)
  end subroutine fail

end module ereaxff_errors

