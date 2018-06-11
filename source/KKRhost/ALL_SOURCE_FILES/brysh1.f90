! ************************************************************************
    Subroutine brysh1(y, x, xsme, ins, irmin, irc, natps, natyp, nspin, imap, &
      lmpot, lsmear)
!*********************************************************************
!     shifts the density or potential of all mt-cell into one single
!     vector and projects out the coulomb part only.
!                                    s. bluegel , kfa , 1987

! ------------------------------------------------------------------------
      Use mod_datatypes, Only: dp
      implicit none
!     .. Parameters ..
      Include 'inc.p'
      Integer :: lmpotd
      Parameter (lmpotd=(lpotd+1)**2)
!..
!.. Local Scalars ..
      Integer :: imap, ins, lmpot, natps, natyp, nspin, lsmear
!..
!     ..
      Real (Kind=dp) :: x(irmd, lmpotd, *), y(*), xsme(irmd, *)
      Integer :: irc(*), irmin(*)

!     Next for SMEARed spherical potential
      Integer :: ia, ip, ir, irc1, irmin1, is, lm


      imap = 0
      Do is = 1, nspin
        Do ia = natps, natyp
          ip = nspin*(ia-1) + is
          irc1 = irc(ia)
          Do ir = 1, irc1
            imap = imap + 1
            y(imap) = x(ir, 1, ip)
          End Do

! ************************************************************************
          If (lsmear>0) Then
            Do ir = 1, irc1
              imap = imap + 1
              y(imap) = xsme(ir, ip)
            End Do
          End If
!*********************************************************************
          If (ins>0 .And. lmpot>1) Then
            irmin1 = irmin(ia)
            Do lm = 2, lmpot
              Do ir = irmin1, irc1
                imap = imap + 1
                y(imap) = x(ir, lm, ip)
              End Do
            End Do
          End If
!     shifts the density or potential of all mt-cell into one single
        End Do
      End Do
!     vector and projects out the coulomb part only.
    End Subroutine
