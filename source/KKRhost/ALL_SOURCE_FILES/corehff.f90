    Subroutine corehff(kap1, kap2, mj, s, nsol, bhf, gck, fck, rc, drdic, &
      rnuc, nzero, nrc)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   *   CALCULATE THE RELATIVISTIC HYPERFINEFIELDS FOR THE             *
!   *                  CURRENT  CORE STATE S                           *
!   *                                                                  *
!   *   THE WAVE FUNCTION  {G(K,S),F(K,S)}  IS NORMALIZED TO 1         *
!   *                                                                  *
!   ********************************************************************

      Implicit None


! PARAMETER definitions
      Real (Kind=dp) :: e0, a0, cautog

!CONVERSION FACTOR FOR HYPERFINE FIELDS FROM A.U. TO GAUSS
!                                 ELECTRON CHARGE     IN ESU
!                                 BOHR-RADIUS         IN CM


      Parameter (e0=1.6021892E-19_dp*2.997930E+09_dp, a0=0.52917706E-08_dp, &
        cautog=e0/(a0*a0))

! Dummy arguments
      Integer :: kap1, kap2, nrc, nsol, nzero, s
      Real (Kind=dp) :: mj, rnuc
      Real (Kind=dp) :: bhf(2, 2), drdic(nrc), fck(2, 2, nrc), gck(2, 2, nrc), &
        rc(nrc)

! Local variables
      Real (Kind=dp) :: ame(2, 2), xx(5), yi(nrc), yy(5), zi(nrc)
      Real (Kind=dp) :: dble, dsqrt
      Integer :: i, k1, k2, n
      Real (Kind=dp) :: ylag

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   ANGULAR HYPERFINE MATRIX ELEMENTS   SEE E.G.  E.M.ROSE
!        THE FACTOR  I  HAS BEEN OMITTED

      ame(1, 1) = 4.0E0_dp*kap1*mj/(4.0E0_dp*kap1*kap1-1.0E0_dp)
      If (nsol==2) Then
        ame(2, 2) = 4.0E0_dp*kap2*mj/(4.0E0_dp*kap2*kap2-1.0E0_dp)
        ame(1, 2) = dsqrt(0.25E0_dp-(mj/dble(kap1-kap2))**2)
        ame(2, 1) = ame(1, 2)
      End If
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      Do k1 = 1, nsol
        Do k2 = 1, nsol
          Do n = 1, nzero
            yi(n) = drdic(n)*(gck(k1,s,n)*fck(k2,s,n)+fck(k1,s,n)*gck(k2,s,n))
          End Do
          Call rint4pts(yi, nzero, zi)
          If (rnuc/=0.0E0_dp) Then
            Do i = 1, 5
              xx(i) = rc(nzero-5+i)
              yy(i) = zi(nzero-5+i)
            End Do
            zi(nzero) = ylag(rnuc, xx, yy, 0, 4, 5)
          End If
          xx(1) = 1.0E0_dp
!                      !RC( 1)
          xx(2) = 6.0E0_dp
!                      !RC( 6)
          xx(3) = 11.0E0_dp
!                      !RC(11)
          yy(1) = zi(nzero) - zi(1)
          yy(2) = zi(nzero) - zi(6)
          yy(3) = zi(nzero) - zi(11)
          bhf(k1, k2) = cautog*ame(k1, k2)*ylag(0.0E0_dp, xx, yy, 0, 2, 3)
        End Do
      End Do

    End Subroutine
