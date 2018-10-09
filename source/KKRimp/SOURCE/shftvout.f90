MODULE MOD_SHFTVOUT

  CONTAINS

  !-------------------------------------------------------------------------------
  !> Summary: 
  !> Author: 
  !> Category: KKRimp, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !-------------------------------------------------------------------------------
  SUBROUTINE SHFTVOUT(VIN, VOUT, SN, LPOT, WG, YRG, LMMAXD, LASSLD, LPOTD, LM3D, NCLEB)
    use mod_ymy, only: YMY
    implicit none
    !-----------------------------------------------------------------------
    !     Transform the 'outer' potential from the old expansion around
    !     the ideal positions to the new expansion around the shifted
    !     positions in case of lattice relaxations.
    !     The 'outer' potential is the madelung potential of the ideal
    !     host minus the intracell potential of the host for the perturbed
    !     cluster.
    !
    !                 gaunt2 has to be called bevor to set up the common
    !                 block assleg
    !
    !-----------------------------------------------------------------------
    !     .. Parameters ..
    !
    !---> attention : ncleb is an empirical factor - it has to be optimized
    !
    !       include 'parameters.file'
    !      INTEGER NTPERD,NTREFD,NATYPD,NATOMD
    !      PARAMETER (NATYPD=19,NTREFD=4,NATOMD=15,NTPERD=NATYPD-NTREFD)
    !      INTEGER LMAXD,LMX,LPOTD
    INTEGER LPOTD
    !      PARAMETER (LMAXD=3,LMX=LMAXD+1,LPOTD=6)
    INTEGER LMMAXD,LM3D !,L3D
    !       PARAMETER (LMMAXD= (LPOTD+1)**2,L3D=2*LPOTD,LM3D= (L3D+1)**2)
    INTEGER LASSLD
    !       PARAMETER (LASSLD=4*LMAXD)
    INTEGER NCLEB
    !       PARAMETER (NCLEB=LM3D*LMMAXD)
    !     ..
    !     .. Scalar Arguments ..
    INTEGER LPOT
    !     ..
    !     .. Array Arguments ..
    REAL*8 SN(3),VIN(*),VOUT(*) !bauer
    !       REAL*8 A(LMMAXD,LMMAXD),SN(3),VIN(*),VOUT(*)
    !     ..
    !     .. Arrays in Common ..
    REAL*8  WG(LASSLD),YRG(LASSLD,0:LASSLD,0:LASSLD)
    !     ..
    !     .. Local Scalars ..
    REAL*8 CLECG,EPI,FACTOR,FPI,PI,R,R1,R2,R3,S,R0
    INTEGER I,IEND,J,L,L1,L2,L3,L3MAX,LM1,LM2,LM3,LMMAX,LX,LY,M,M1,M1A,M1S,M2,M2A,M2S,M3,M3A,M3S,LM
    !     ..
    !     .. Local Arrays ..
    REAL*8 CLEB(NCLEB),DFAC(0:LPOTD,0:LPOTD),Y(LM3D)
    INTEGER ICLEB(NCLEB,3),LOFLM(LM3D)
    !     ..
    !     ..
    !     .. Save statement ..
    SAVE PI
    !     ..
    !     .. External Subroutines ..
    !       EXTERNAL YMY
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC ABS,REAL,SIGN
    !     ..
    PI = 4.D0*DATAN(1.D0)
    !
    !---> determine the l-value for given lm
    
    I = 1
    DO L = 0,2*LPOT !Bauer L3D
      DO M = -L,L
        LOFLM(I) = L
        I = I + 1
      END DO
    END DO
    DO LM=1,LMMAXD
      VOUT(LM) = 0.D0
    END DO
    
    FPI = 4.0D0*PI
    EPI = 8.0D0*PI
    L3MAX = 2*LPOT
    LMMAX = (LPOT+1)**2
    !
    !--->calculate:                  (2*(l+l')+1)!!
    !                dfac(l,l')= ----------------------
    !                            (2*l+1)!! * (2*l'+1)!!
    
    DFAC(0,0) = 1.D0
    DO LX = 1,LPOT
      DFAC(LX,0) = DFAC(LX-1,0)
      DFAC(0,LX) = DFAC(LX,0)
      DO LY = 1,LX
        DFAC(LX,LY) = DFAC(LX,LY-1)*REAL(2* (LX+LY)+1)/REAL(2*LY+1)
        DFAC(LY,LX) = DFAC(LX,LY)
      END DO
    END DO
    !
    !---> set up of the gaunt coefficients with an index field
    !     recognize that they are needed here only for l3=l1+l2
    !
    I = 1
      DO L1 = 0,LPOT
        DO L2 = 0,L1
        L3 = L1 - L2
        DO M1 = -L1,L1
          DO M2 = -L2,L2
            DO M3 = -L3,L3
              M1S = SIGN(1,M1)
              M2S = SIGN(1,M2)
              M3S = SIGN(1,M3)
              
              IF (M1S*M2S*M3S.GE.0) THEN
              
                M1A = ABS(M1)
                M2A = ABS(M2)
                M3A = ABS(M3)
                
                FACTOR = 0.0D0
                
                IF (M1A+M2A.EQ.M3A) FACTOR = FACTOR + REAL(3*M3S+SIGN(1,-M3))/8.0D0
                IF (M1A-M2A.EQ.M3A) FACTOR = FACTOR + REAL(M1S)/4.0D0
                IF (M2A-M1A.EQ.M3A) FACTOR = FACTOR + REAL(M2S)/4.0D0
                
                IF (FACTOR.NE.0.0D0) THEN
                
                  IF (M1S*M2S.NE.1 .OR. M2S*M3S.NE.1 .OR. M1S*M3S.NE.1) FACTOR = -FACTOR
                  
                  S = 0.0D0
                  DO J = 1,LASSLD
                    S = S + WG(J)*YRG(J,L1,M1A)*YRG(J,L2,M2A)*YRG(J,L3,M3A)
                  END DO
                  CLECG = S*FACTOR
                  IF (ABS(CLECG).GT.1.D-10) THEN
                    CLEB(I) = CLECG
                    ICLEB(I,1) = L1* (L1+1) + M1 + 1
                    ICLEB(I,2) = L2* (L2+1) + M2 + 1
                    ICLEB(I,3) = L3* (L3+1) + M3 + 1
                    I = I + 1
                  END IF

                END IF

              END IF

            END DO
          END DO
        END DO
      END DO
    END DO
    IEND = I - 1
    IF (NCLEB.LT.IEND) THEN
      STOP 13

    ELSE
      !        WRITE (6,FMT='(I10)') IEND

      DO LM1 = 1,LMMAX
        ! DO 140 LM2 = 1,LMMAX
        !    A(LM1,LM2) = 0.0D0
        ! END DO
        VOUT(LM1)=0.0D0
      END DO
      

      R1 = SN(1)
      R2 = SN(2)
      R3 = SN(3)
      
      R0 = R1**2+R2**2 +R3**2
      IF (R0.GT.1.D-10) THEN
        CALL YMY(R1,R2,R3,R,Y,L3MAX)
      ELSE
        DO LM=1,LMMAX
          VOUT(LM) = VIN(LM)
        END DO
        RETURN
      END IF
      
      DO I = 1,IEND
        LM1 = ICLEB(I,1)
        LM2 = ICLEB(I,2)
        LM3 = ICLEB(I,3)
        L1 = LOFLM(LM1)
        L2 = LOFLM(LM2)
        L3 = LOFLM(LM3)

        VOUT(LM2) = VOUT(LM2) + FPI*(-1.D0)**L3*DFAC(L2,L3)*CLEB(I)*VIN(LM1)*R**L3*Y(LM3)

      END DO

    END IF
  END SUBROUTINE
END MODULE MOD_SHFTVOUT
