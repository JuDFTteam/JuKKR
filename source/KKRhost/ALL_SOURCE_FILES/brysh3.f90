module mod_brysh3

contains

  !-------------------------------------------------------------------------------
  !> Summary: Broyden mixing tool
  !> Author: S. Bluegel
  !> Date: 1987
  !> Category: KKRhost, mixing
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Shifts the density or potential of all mt-cell into one single
  !> vector and projects out the coulomb part only.
  !>
  !> S. Bluegel , KFA , 1987
  !-------------------------------------------------------------------------------
  subroutine brysh3(y, x, z, xsme, ins, irmin, irc, natps, natyp, nspin, imap, lmpot, lsmear)

    use :: mod_datatypes, only: dp
    use :: global_variables
    implicit none
    ! ..
    ! .. Local Scalars ..
    integer :: imap, ins, lmpot, natps, natyp, nspin, lsmear
    ! ..

    real (kind=dp) :: x(irmd, *), y(*), z(irmind:irmd, lmpotd, *), xsme(irmd, *)
    integer :: irc(*), irmin(*)

    integer :: ia, ip, ir, irc1, irmin1, is, lm

    imap = 0
    do is = 1, nspin
      do ia = natps, natyp

        ip = nspin*(ia-1) + is
        irc1 = irc(ia)
        do ir = 1, irc1
          imap = imap + 1
          y(imap) = x(ir, ip)
        end do

        ! ************************************************************************
        ! SMEARed spherical potential
        if (lsmear>0) then
          do ir = 1, irc1
            imap = imap + 1
            y(imap) = xsme(ir, ip)
          end do
        end if
        ! *********************************************************************

        if (ins>0 .and. lmpot>1) then
          irmin1 = irmin(ia)
          do lm = 2, lmpot
            do ir = irmin1, irc1
              imap = imap + 1
              y(imap) = z(ir, lm, ip)
            end do
          end do
        end if

      end do
    end do

  end subroutine brysh3

end module mod_brysh3
