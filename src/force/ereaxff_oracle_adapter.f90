module ereaxff_oracle_adapter
  implicit none
  private

  public :: ereaxff_oracle_not_implemented

contains

  subroutine ereaxff_oracle_not_implemented()
    error stop "eReaxFF oracle adapter is not implemented yet; validate harmonic/dummy force paths first."
  end subroutine ereaxff_oracle_not_implemented

end module ereaxff_oracle_adapter

