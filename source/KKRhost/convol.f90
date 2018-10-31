!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_convol

  private
  public :: convol

contains

  !-------------------------------------------------------------------------------
  !> Summary: Convolutes potentials with shape functions
  !> Author: 
  !> Category: KKRhost, shape-functions 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Calculate convolution of potential with shapefunction.
  !-------------------------------------------------------------------------------
  subroutine convol(imt1, irc1, icell, imaxsh, ilm_map, ifunm, lmpot, gsh, thetas, thesme, z, rfpi, r, vons, vspsmo, lmsp)

    use :: mod_datatypes, only: dp
    use :: global_variables, only: irmd, irid, nfund, natypd, ngshd, lmpotd
    implicit none

    ! .. Local Scalars ..
    real (kind=dp) :: rfpi, z
    integer :: icell, imaxsh, imt1, irc1, lmpot

    ! .. Local Arrays ..
    real (kind=dp) :: gsh(*), r(*), thetas(irid, nfund, *), vons(irmd, *), thesme(irid, nfund, *), vspsmo(irmd)
    integer :: ifunm(natypd, *), ilm_map(ngshd, 3), lmsp(natypd, *)

    real (kind=dp) :: zzor
    integer :: i, ifun, ir, irh, lm, lm1, lm2, lm3

    real (kind=dp) :: vstore(irid, lmpotd), vstsme(irid, lmpotd)


    do lm = 1, lmpot
      do ir = 1, irc1 - imt1
        vstore(ir, lm) = 0.0e0_dp
        vstsme(ir, lm) = 0.0e0_dp
      end do
    end do
    ! COPY THE PART INSIDE THE MT SPHERE
    do ir = imt1 + 1, irc1
      zzor = 2.0e0_dp*z/r(ir)*rfpi
      vons(ir, 1) = vons(ir, 1) - zzor
    end do

    do i = 1, imaxsh
      lm1 = ilm_map(i, 1)
      lm2 = ilm_map(i, 2)
      lm3 = ilm_map(i, 3)
      if (lmsp(icell,lm3)>0) then
        ifun = ifunm(icell, lm3)
        do ir = imt1 + 1, irc1
          irh = ir - imt1
          vstore(irh, lm1) = vstore(irh, lm1) + gsh(i)*vons(ir, lm2)*thetas(irh, ifun, icell)
          vstsme(irh, lm1) = vstsme(irh, lm1) + gsh(i)*vons(ir, lm2)*thesme(irh, ifun, icell)
        end do
      end if
    end do

    do ir = imt1 + 1, irc1
      irh = ir - imt1
      zzor = 2.0e0_dp*z/r(ir)*rfpi
      vons(ir, 1) = vstore(irh, 1) + zzor
      vspsmo(ir) = (vstsme(irh,1)+zzor)/rfpi
    end do

    ! ************************************************************************
    do ir = 1, imt1
      vspsmo(ir) = vons(ir, 1)/rfpi
    end do
    ! ************************************************************************

    do lm = 2, lmpot
      do ir = imt1 + 1, irc1
        irh = ir - imt1
        vons(ir, lm) = vstore(irh, lm)
      end do
    end do

    return

  end subroutine convol

end module mod_convol
