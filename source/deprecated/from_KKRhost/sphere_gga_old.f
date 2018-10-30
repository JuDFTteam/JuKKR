      SUBROUTINE SPHERE_GGA(LMAX,YR,WTYR,RIJ,IJD,LMMAXD,THET,YLM,
     +                      DYLMT1,DYLMT2,DYLMF1,DYLMF2,DYLMTF)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     generate an angular mesh and spherical harmonics at those
c     mesh points. For an angular integration the weights are ge-
c     rated .
c
c     R. Zeller      Feb. 1996
c     Small change for GGA implementation
c     Nikos          Dec. 1996
c-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER IJD,LMAX,LMMAXD
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DX1,DX2,DX3,F0,PI,R,R1,R2,R3
      INTEGER IJ,LM1
C     ..
C     .. External Subroutines ..
      EXTERNAL CYLM02,YMY
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DYLMF1(IJD,LMMAXD),DYLMF2(IJD,LMMAXD),
     +                 DYLMT1(IJD,LMMAXD),DYLMT2(IJD,LMMAXD),
     +                 DYLMTF(IJD,LMMAXD),RIJ(IJD,3),THET(IJD),
     +                 WTYR(IJD,*),YLM(IJD,LMMAXD),YR(IJD,*)
C     ..
C     .. Local Arrays ..
C
      DOUBLE PRECISION COSFI(IJD),COSX(IJD),FAI(IJD),ND(3,3),SINFI(IJD),
     +                 WGHT,Y(1000)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ACOS,ATAN,COS,SIN,SQRT
C     ..
      PI = 4.D0*ATAN(1.D0)
      WRITE (6,*) 'SPHERE for GGA: read LEBEDEV mesh'
      IF (IJD.GT.1000) STOP 'SPHERE'
c
c
      DO 30 IJ = 1,IJD
        CALL LEBEDEV (IJ,R1,R2,R3,WGHT)
c
c      make a small rotation
c
        F0 = 0.08d0
        ND(1,1) = COS(F0)
        ND(1,2) = 0D0
        ND(1,3) = SIN(F0)
        ND(2,1) = 0D0
        ND(2,2) = 1D0
        ND(2,3) = 0D0
        ND(3,1) = -SIN(F0)
        ND(3,2) = 0D0
        ND(3,3) = COS(F0)
C
        DX1 = ND(1,1)*R1 + ND(2,1)*R2 + ND(3,1)*R3
        DX2 = ND(1,2)*R1 + ND(2,2)*R2 + ND(3,2)*R3
        DX3 = ND(1,3)*R1 + ND(2,3)*R2 + ND(3,3)*R3
C
        R1 = DX1
        R2 = DX2
        R3 = DX3
C
        RIJ(IJ,1) = R1
        RIJ(IJ,2) = R2
        RIJ(IJ,3) = R3
C
        CALL YMY(R1,R2,R3,R,Y,LMAX)
        DO 10 LM1 = 1, (LMAX+1)**2
          YR(IJ,LM1) = Y(LM1)
   10   CONTINUE
c
c---> multiply the spherical harmonics with the weights
c
        DO 20 LM1 = 1, (LMAX+1)**2
          WTYR(IJ,LM1) = YR(IJ,LM1)*WGHT*PI*4.D0
   20   CONTINUE
c
c---> produce what is needed for GGA
c
        COSX(IJ) = R3
        IF (ABS(R3).NE.1.d0) THEN
          COSFI(IJ) = R1/SQRT(1.d0-R3*R3)
          SINFI(IJ) = R2/SQRT(1.d0-R3*R3)
          IF (ABS(COSFI(IJ)).GT.1.D0) COSFI(IJ) = COSFI(IJ)/
     +        ABS(COSFI(IJ))
          IF (ABS(SINFI(IJ)).GT.1.D0) SINFI(IJ) = SINFI(IJ)/
     +        ABS(SINFI(IJ))
          FAI(IJ) = ACOS(COSFI(IJ))
        ELSE IF (SINFI(IJ).EQ.0.d0) THEN
          FAI(IJ) = PI/2.d0
        ELSE
          COSFI(IJ) = R1
          SINFI(IJ) = R2
          IF (ABS(COSFI(IJ)).GT.1.D0) COSFI(IJ) = COSFI(IJ)/
     +        ABS(COSFI(IJ))
        END IF
        FAI(IJ) = ACOS(COSFI(IJ))
   30 CONTINUE

      CALL CYLM02(LMAX,COSX,FAI,2*LMAX+1,LMMAXD,THET,YLM,DYLMT1,DYLMT2,
     +            DYLMF1,DYLMF2,DYLMTF)
C
      END
