!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_brysh2

contains

  !-------------------------------------------------------------------------------
  !> Summary: Broyden mixing tool
  !> Author: S. Bluegel
  !> Date: 1987
  !> Category: KKRhost, mixing
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Maps the density or potential back from one single vector into
  !> the proper bins of each single mt-cell. The magnetization
  !> density is also added in.
  !> S. Bluegel , KFA , 1987
  !-------------------------------------------------------------------------------
  subroutine brysh2(y, x, xsme, ins, irmin, irc, natps, natyp, nspin, imap, lmpot, lsmear)

    use :: mod_datatypes, only: dp
    use :: global_variables, only: irmd, lmpotd
    implicit none
    ! ..
    ! .. Local Scalars ..
    integer :: imap, ins, lmpot, natps, natyp, nspin, lsmear
    ! ..
    real (kind=dp) :: x(irmd, lmpotd, *), y(*), xsme(irmd, *)
    integer :: irc(*), irmin(*)

    integer :: ia, ip, ir, irc1, irmin1, is, lm

    imap = 0

    do is = 1, nspin
      do ia = natps, natyp

        ip = nspin*(ia-1) + is
        irc1 = irc(ia)
        do ir = 1, irc1
          imap = imap + 1
          x(ir, 1, ip) = y(imap)
        end do

        ! ************************************************************************
        ! Next for SMEARed spherical potential
        if (lsmear>0) then
          do ir = 1, irc1
            imap = imap + 1
            xsme(ir, ip) = y(imap)
          end do
        end if
        ! *********************************************************************

        if (ins>0 .and. lmpot>1) then
          irmin1 = irmin(ia)
          do lm = 2, lmpot
            do ir = irmin1, irc1
              imap = imap + 1
              x(ir, lm, ip) = y(imap)
            end do
          end do
        end if

      end do
    end do

  end subroutine brysh2

end module mod_brysh2
