    Subroutine setfactl(factl, lmax, krel, lmmaxd)
      Use mod_datatypes, Only: dp
      Implicit None
!..
!.. Parameters
      Complex (Kind=dp) :: ci
      Parameter (ci=(0E0_dp,1E0_dp))
!..
!.. Arguments
      Integer :: krel, lmax, lmmaxd
      Complex (Kind=dp) :: factl(lmmaxd, lmmaxd)
!..
!.. Locals
      Integer :: ii1, ii2, l1, l2, lm1, lm2, mm1, mm2, imu1, imu2
      Integer :: kap1(2), kap2(2), nsol1, nsol2
      Real (Kind=dp) :: mu1, mu2, mu1m05, mu2m05
!..
!.. Externals
      External :: cinit
!     ..
      Call cinit(lmmaxd*lmmaxd, factl)

! ----------------------------------------------------------------------
      If (krel==0) Then
        lm1 = 0
        Do l1 = 0, lmax
          Do mm1 = 1, 2*l1 + 1
            lm1 = lm1 + 1
            lm2 = 0
            Do l2 = 0, lmax
              Do mm2 = 1, 2*l2 + 1
                lm2 = lm2 + 1
                factl(lm1, lm2) = (-1)**(l1+l2)
              End Do
            End Do
          End Do
        End Do
! ----------------------------------------------------------------------
      Else ! KREL.EQ.1
! ----------------------------------------------------------------------
        lm1 = 0
! l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1
        Do l1 = 0, lmax
          kap1(1) = l1
          kap1(2) = -l1 - 1
          nsol1 = 2
          If (l1==0) Then
            kap1(1) = kap1(2)
            nsol1 = 1
          End If
          Do ii1 = 1, nsol1
            mu1m05 = abs(kap1(ii1)) - 0.5E0_dp
! m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1
!DO MU1 = -MU1M05,MU1M05,1D0
            Do imu1 = 1, 2*nint(mu1m05) + 1
              mu1 = -mu1m05 + real(imu1-1, kind=dp)
              lm1 = lm1 + 1
              lm2 = 0
! l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2
              Do l2 = 0, lmax
                kap2(1) = l2
                kap2(2) = -l2 - 1
                nsol2 = 2
                If (l2==0) Then
                  kap2(1) = kap2(2)
                  nsol2 = 1
                End If
                Do ii2 = 1, nsol2
                  mu2m05 = abs(kap2(ii2)) - 0.5E0_dp
! m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2
!DO MU2 = -MU2M05,MU2M05,1D0
                  Do imu2 = 1, 2*nint(mu2m05) + 1
                    mu2 = -mu2m05 + real(imu2-1, kind=dp)
                    lm2 = lm2 + 1
                    mm1 = int(mu2-mu1)

                    factl(lm1, lm2) = (-1)**(l1+l2)*ci**mm1

                  End Do
! m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2
                End Do
              End Do
! l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2
            End Do
! m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1
          End Do
        End Do
! l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1
      End If
! ----------------------------------------------------------------------
    End Subroutine
