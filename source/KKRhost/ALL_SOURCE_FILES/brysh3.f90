! ************************************************************************
subroutine brysh3(y, x, z, xsme, ins, irmin, irc, natps, natyp, nspin, imap, &
  lmpot, lsmear)
!*********************************************************************
!     shifts the density or potential of all mt-cell into one single
!     vector and projects out the coulomb part only.

!                                    s. bluegel , kfa , 1987

! ------------------------------------------------------------------------
!.. Parameters ..
  include 'inc.p'
  integer :: lmpotd
  parameter (lmpotd=(lpotd+1)**2)
  integer :: irmind
  parameter (irmind=irmd-irnsd)
!..
!.. Local Scalars ..
  integer :: imap, ins, lmpot, natps, natyp, nspin, lsmear
!     ..

  double precision :: x(irmd, *), y(*), z(irmind:irmd, lmpotd, *), &
    xsme(irmd, *)
  integer :: irc(*), irmin(*)
!     SMEARed spherical potential

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
      if (lsmear>0) then
        do ir = 1, irc1
          imap = imap + 1
          y(imap) = xsme(ir, ip)
        end do
      end if
!*********************************************************************
      if (ins>0 .and. lmpot>1) then
        irmin1 = irmin(ia)
        do lm = 2, lmpot
          do ir = irmin1, irc1
            imap = imap + 1
            y(imap) = z(ir, lm, ip)
          end do
        end do
      end if
!     shifts the density or potential of all mt-cell into one single
    end do
  end do
!     vector and projects out the coulomb part only.
end subroutine
