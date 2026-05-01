module ereaxff_reaxff_compat_output
  implicit none
  private

  public :: reaxff_compat_not_implemented

contains

  subroutine reaxff_compat_not_implemented()
    error stop "ReaxFF compatibility outputs are not implemented yet; current harmonic driver writes PIMD-native outputs only."
  end subroutine reaxff_compat_not_implemented

end module ereaxff_reaxff_compat_output

