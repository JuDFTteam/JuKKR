      SUBROUTINE SETFACTL(FACTL,LMAX,KREL,LMMAXD)
      IMPLICIT NONE
C     ..
C     .. Parameters
      DOUBLE COMPLEX CI
      PARAMETER (CI=(0D0,1D0))
C     ..
C     .. Arguments
      INTEGER KREL,LMAX,LMMAXD
      DOUBLE COMPLEX FACTL(LMMAXD,LMMAXD)
C     ..
C     .. Locals
      INTEGER II1,II2,L1,L2,LM1,LM2,MM1,MM2,IMU1,IMU2
      INTEGER KAP1(2),KAP2(2),NSOL1,NSOL2
      DOUBLE PRECISION MU1,MU2,MU1M05,MU2M05
C     ..
C     .. Externals
      EXTERNAL CINIT
C     ..
      CALL CINIT(LMMAXD*LMMAXD,FACTL)
C
C ----------------------------------------------------------------------
      IF (KREL.EQ.0) THEN
         LM1 = 0
         DO L1 = 0,LMAX
            DO MM1 = 1,2*L1+1
               LM1 = LM1 + 1
               LM2 = 0
               DO L2 = 0,LMAX
                  DO MM2 = 1,2*L2+1
                     LM2 = LM2 + 1
                     FACTL(LM1,LM2) = (-1)**(L1+L2)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C ----------------------------------------------------------------------
      ELSE                      ! KREL.EQ.1
C ----------------------------------------------------------------------
         LM1 = 0
C l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1
         DO L1=0,LMAX
            KAP1(1) = L1
            KAP1(2) = -L1 - 1
            NSOL1 = 2
            IF ( L1.EQ.0 ) THEN
               KAP1(1) = KAP1(2)
               NSOL1 = 1
            END IF
            DO II1 = 1,NSOL1
               MU1M05 = ABS(KAP1(II1))-0.5D0
C m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1
               !DO MU1 = -MU1M05,MU1M05,1D0
               DO IMU1 = 1,2*NINT(MU1M05)+1
                  MU1 = -MU1M05 + DBLE(IMU1-1)
                  LM1 = LM1 + 1
                  LM2 = 0
C l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2
                  DO L2 = 0,LMAX
                     KAP2(1) = L2
                     KAP2(2) = -L2 - 1
                     NSOL2 = 2
                     IF ( L2.EQ.0 ) THEN
                        KAP2(1) = KAP2(2)
                        NSOL2 = 1
                     END IF
                     DO II2 = 1,NSOL2
                        MU2M05 = ABS(KAP2(II2))-0.5D0
C m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2
                        !DO MU2 = -MU2M05,MU2M05,1D0
                        DO IMU2 = 1, 2*NINT(MU2M05)+1
                           MU2 = -MU2M05 + DBLE(IMU2-1)
                           LM2 = LM2 + 1
                           MM1 = INT(MU2-MU1)
C     
                           FACTL(LM1,LM2) = (-1)** (L1+L2) * CI**MM1
C
                        END DO
C m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2
                     END DO
                  END DO       
C l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2
               END DO
C m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1
            END DO
         END DO
C l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1
      END IF
C ----------------------------------------------------------------------
      END
