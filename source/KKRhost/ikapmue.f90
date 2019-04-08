!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

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
