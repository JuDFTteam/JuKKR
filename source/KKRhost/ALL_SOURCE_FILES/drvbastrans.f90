    Subroutine drvbastrans(rc, crel, rrel, srrel, nrrel, irrel, nlmax, nkmmax, &
      nmuemax, nkmpmax, nkmax, linmax)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   *                                                                  *
!   ********************************************************************
      Implicit None

! Dummy arguments
      Integer :: linmax, nkmax, nkmmax, nkmpmax, nlmax, nmuemax
      Complex (Kind=dp) :: crel(nkmmax, nkmmax), rc(nkmmax, nkmmax), &
        rrel(nkmmax, nkmmax), srrel(2, 2, nkmmax)
      Integer :: irrel(2, 2, nkmmax), nrrel(2, nkmmax)

! Local variables
      Real (Kind=dp) :: cgc(nkmpmax, 2)
      Integer :: i, ikm1lin(linmax), ikm2lin(linmax), il, imue, iprint, &
        kaptab(nmuemax), ltab(nmuemax), mmax, nmuetab(nmuemax), &
        nsollm(nlmax, nmuemax)

      If (nkmmax/=2*nlmax**2) Stop ' Check NLMAX,NKMMAX in < DRVBASTRANS > '
      If (nmuemax/=2*nlmax) Stop ' Check NLMAX,NMUEMAX in < DRVBASTRANS > '
      If (nkmpmax/=(nkmmax+2*nlmax)) Stop &
        ' Check NLMAX,NKMMAX,NKMPMAX in < DRVBASTRANS > '
      If (nkmax/=2*nlmax-1) Stop ' Check NLMAX,NKMAX in < DRVBASTRANS > '
      If (linmax/=(2*nlmax*(2*nlmax-1))) Stop &
        ' Check NLMAX,LINMAX in < DRVBASTRANS > '

      iprint = 0

      Do i = 1, nmuemax
        ltab(i) = i/2
        If (2*ltab(i)==i) Then
          kaptab(i) = ltab(i)
        Else
          kaptab(i) = -ltab(i) - 1
        End If
        nmuetab(i) = 2*abs(kaptab(i))
      End Do

      Do il = 1, nlmax
        mmax = 2*il
        Do imue = 1, mmax
          If ((imue==1) .Or. (imue==mmax)) Then
            nsollm(il, imue) = 1
          Else
            nsollm(il, imue) = 2
          End If
        End Do
      End Do

      Call ikmlin(iprint, nsollm, ikm1lin, ikm2lin, nlmax, nmuemax, linmax, &
        nlmax)

      Call calccgc(ltab, kaptab, nmuetab, cgc, nkmax, nmuemax, nkmpmax)

! ---------------------------- now calculate the transformation matrices

      Call strsmat(nlmax-1, cgc, srrel, nrrel, irrel, nkmmax, nkmpmax)

      Call bastrmat(nlmax-1, cgc, rc, crel, rrel, nkmmax, nkmpmax)

      Return
    End Subroutine
