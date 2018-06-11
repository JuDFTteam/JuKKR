    Subroutine calccgc(ltab, kaptab, nmuetab, cgc, nkmax, nmuemax, nkmpmax)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   *   CLEBSCH-GORDON-COEFFICIENTS     CGC(IKM,IS)                    *
!   *                                                                  *
!   *   IKM NUMBERS  CGC  FOR INCREASING  K  AND  MUE                  *
!   *   IKM  = L*2*(J+1/2) + J + MUE + 1                               *
!   *   IS= 1/2  SPIN DOWN/UP                                          *
!   *                                                                  *
!   ********************************************************************

      Implicit None

! Dummy arguments
      Integer :: nkmax, nkmpmax, nmuemax
      Real (Kind=dp) :: cgc(nkmpmax, 2)
      Integer :: kaptab(nmuemax), ltab(nmuemax), nmuetab(nmuemax)

! Local variables
      Integer :: ikm, k, kappa, m
      Real (Kind=dp) :: j, l, mue, twolp1

      ikm = 0
      Do k = 1, (nkmax+1)
        l = ltab(k)
        kappa = kaptab(k)
        j = abs(kappa) - 0.5E0_dp
        mue = -j - 1.0E0_dp
        twolp1 = 2.0E0_dp*l + 1.0E0_dp

        If (kappa<0) Then

!     J = L + 1/2
          Do m = 1, nmuetab(k)

            mue = mue + 1.0E0_dp
            ikm = ikm + 1
            cgc(ikm, 1) = sqrt((l-mue+0.5E0_dp)/twolp1)
            cgc(ikm, 2) = sqrt((l+mue+0.5E0_dp)/twolp1)
          End Do
        Else
!     J = L - 1/2
          Do m = 1, nmuetab(k)

            mue = mue + 1.0E0_dp
            ikm = ikm + 1
            cgc(ikm, 1) = sqrt((l+mue+0.5E0_dp)/twolp1)
            cgc(ikm, 2) = -sqrt((l-mue+0.5E0_dp)/twolp1)

          End Do
        End If


      End Do

    End Subroutine
