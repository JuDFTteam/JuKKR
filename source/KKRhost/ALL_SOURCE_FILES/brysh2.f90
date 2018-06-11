! ************************************************************************
    Subroutine brysh2(y, x, xsme, ins, irmin, irc, natps, natyp, nspin, imap, &
      lmpot, lsmear)
!*********************************************************************
!     maps the density or potential back from one single vector into
!     the proper bins of each single mt-cell . the magnetization
!     density is also added in.
!                                    s. bluegel , kfa , 1987

      Use mod_datatypes, Only: dp
      implicit none
! ------------------------------------------------------------------------
!.. Parameters ..
      Include 'inc.p'
!..
!.. Array Arguments ..
      Integer :: lmpotd
      Parameter (lmpotd=(lpotd+1)**2)
!..
!.. Local Scalars ..
      Integer :: imap, ins, lmpot, natps, natyp, nspin, lsmear
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
            x(ir, 1, ip) = y(imap)
          End Do

! ************************************************************************
          If (lsmear>0) Then
            Do ir = 1, irc1
              imap = imap + 1
              xsme(ir, ip) = y(imap)
            End Do
          End If
!*********************************************************************
          If (ins>0 .And. lmpot>1) Then
            irmin1 = irmin(ia)
            Do lm = 2, lmpot
              Do ir = irmin1, irc1
                imap = imap + 1
                x(ir, lm, ip) = y(imap)
              End Do
            End Do
          End If
!     maps the density or potential back from one single vector into
        End Do
      End Do
!     the proper bins of each single mt-cell . the magnetization
    End Subroutine
