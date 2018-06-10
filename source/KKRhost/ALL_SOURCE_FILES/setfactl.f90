SUBROUTINE setfactl(factl,lmax,krel,lmmaxd)
      IMPLICIT NONE
!..
!.. Parameters
      DOUBLE COMPLEX CI
      PARAMETER (CI=(0D0,1D0))
!..
!.. Arguments
      INTEGER KREL,LMAX,LMMAXD
      DOUBLE COMPLEX FACTL(LMMAXD,LMMAXD)
!..
!.. Locals
      INTEGER II1,II2,L1,L2,LM1,LM2,MM1,MM2,IMU1,IMU2
      INTEGER KAP1(2),KAP2(2),NSOL1,NSOL2
      DOUBLE PRECISION MU1,MU2,MU1M05,MU2M05
!..
!.. Externals
      EXTERNAL CINIT
!     ..
CALL cinit(lmmaxd*lmmaxd,factl)

! ----------------------------------------------------------------------
IF (krel == 0) THEN
  lm1 = 0
  DO l1 = 0,lmax
    DO mm1 = 1,2*l1+1
      lm1 = lm1 + 1
      lm2 = 0
      DO l2 = 0,lmax
        DO mm2 = 1,2*l2+1
          lm2 = lm2 + 1
          factl(lm1,lm2) = (-1)**(l1+l2)
        END DO
      END DO
    END DO
  END DO
! ----------------------------------------------------------------------
ELSE                      ! KREL.EQ.1
! ----------------------------------------------------------------------
  lm1 = 0
! l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1
  DO l1=0,lmax
    kap1(1) = l1
    kap1(2) = -l1 - 1
    nsol1 = 2
    IF ( l1 == 0 ) THEN
      kap1(1) = kap1(2)
      nsol1 = 1
    END IF
    DO ii1 = 1,nsol1
      mu1m05 = ABS(kap1(ii1))-0.5D0
! m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1
!DO MU1 = -MU1M05,MU1M05,1D0
      DO imu1 = 1,2*nint(mu1m05)+1
        mu1 = -mu1m05 + DBLE(imu1-1)
        lm1 = lm1 + 1
        lm2 = 0
! l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2
        DO l2 = 0,lmax
          kap2(1) = l2
          kap2(2) = -l2 - 1
          nsol2 = 2
          IF ( l2 == 0 ) THEN
            kap2(1) = kap2(2)
            nsol2 = 1
          END IF
          DO ii2 = 1,nsol2
            mu2m05 = ABS(kap2(ii2))-0.5D0
! m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2
!DO MU2 = -MU2M05,MU2M05,1D0
            DO imu2 = 1, 2*nint(mu2m05)+1
              mu2 = -mu2m05 + DBLE(imu2-1)
              lm2 = lm2 + 1
              mm1 = INT(mu2-mu1)
              
              factl(lm1,lm2) = (-1)** (l1+l2) * ci**mm1
              
            END DO
! m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2
          END DO
        END DO
! l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2
      END DO
! m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1
    END DO
  END DO
! l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1
END IF
! ----------------------------------------------------------------------
END SUBROUTINE setfactl
