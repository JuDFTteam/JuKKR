module mod_brysh2

contains

! ************************************************************************
subroutine brysh2(y, x, xsme, ins, irmin, irc, natps, natyp, nspin, imap, &
  lmpot, lsmear)
  ! *********************************************************************
  ! maps the density or potential back from one single vector into
  ! the proper bins of each single mt-cell . the magnetization
  ! density is also added in.
  ! s. bluegel , kfa , 1987

  use :: mod_datatypes, only: dp
  use global_variables
  implicit none
  ! ------------------------------------------------------------------------
  ! ..
  ! .. Local Scalars ..
  integer :: imap, ins, lmpot, natps, natyp, nspin, lsmear
  ! ..

  real (kind=dp) :: x(irmd, lmpotd, *), y(*), xsme(irmd, *)
  integer :: irc(*), irmin(*)

  ! Next for SMEARed spherical potential
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
      ! maps the density or potential back from one single vector into
    end do
  end do
  ! the proper bins of each single mt-cell . the magnetization
end subroutine brysh2

end module mod_brysh2
