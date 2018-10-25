module mod_ikapmue

contains

  !-------------------------------------------------------------------------------
  !> Summary: Get combined (kappa,mue) index from kappa and mue
  !> Author: 
  !> Category: KKRhost, dirac, special-functions
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> INDEXING OF MATRIX-ELEMENTS:                                    *
  !>                                                                 *
  !> I = 2*L*(J+1/2) + J + MUE + 1
  !-------------------------------------------------------------------------------
  integer function ikapmue(kappa, muem05)

    implicit none

    ! Dummy arguments
    integer :: kappa, muem05

    ! Local variables
    integer :: iabs
    integer :: jp05, l

    jp05 = iabs(kappa)

    if (kappa<0) then
      l = -kappa - 1
    else
      l = kappa
    end if

    ikapmue = 2*l*jp05 + jp05 + muem05 + 1

  end function ikapmue

end module mod_ikapmue
