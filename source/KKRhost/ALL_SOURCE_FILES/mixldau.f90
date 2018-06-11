    Subroutine mixldau(mmaxd, nspind, natypd, natyp, nspin, lopt, wldauold, &
      wldau)
      Use mod_datatypes, Only: dp
      Implicit None
! Input:
      Integer :: natypd, nspind, mmaxd
      Integer :: lopt(natypd)
      Integer :: natyp, nspin
      Real (Kind=dp) :: wldauold(mmaxd, mmaxd, nspind, natypd)
! Input/Output:
      Real (Kind=dp) :: wldau(mmaxd, mmaxd, nspind, natypd)
! Inside:
      Integer :: iat, is, m1, m2, mmax
      Integer :: ier
      Real (Kind=dp) :: xmix, xmix2, rmserr
      Character (Len=256) :: uio ! NCOLIO=256


      External :: ioinput


! First calculate rms error in interaction matrix

      Do iat = 1, natyp
        rmserr = 0.E0_dp
        If (lopt(iat)>=0) Then
          mmax = 2*lopt(iat) + 1
          Do is = 1, nspin
            Do m2 = 1, mmax
              Do m1 = 1, mmax
                rmserr = rmserr + (wldau(m1,m2,is,iat)-wldauold(m1,m2,is,iat)) &
                  **2
              End Do
            End Do
          End Do
          rmserr = sqrt(rmserr)
          Write (1337, 100) iat, rmserr
100       Format ('LDA+U interaction matrix rms error for atom', I6, ' = ', &
            E10.2)
        End If
      End Do

! Now mix old/new interaction matrices
      ier = 0
      Call ioinput('MIXLDAU         ', uio, 1, 7, ier)
      If (ier/=0) Then
        Write (*, *) 'MIXLDAU not found, setting to 1.'
        Return
      Else
        Read (Unit=uio, Fmt=*) xmix
        Write (1337, *) 'Using MIXLDAU = ', xmix
      End If

      xmix2 = 1.E0_dp - xmix

      Do iat = 1, natyp
        If (lopt(iat)>=0) Then
          Do is = 1, nspin
            Do m2 = 1, mmaxd
              Do m1 = 1, mmaxd
                wldau(m1, m2, is, iat) = xmix*wldau(m1, m2, is, iat) + &
                  xmix2*wldauold(m1, m2, is, iat)
              End Do
            End Do
          End Do
        End If
      End Do

    End Subroutine
