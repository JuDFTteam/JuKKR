      SUBROUTINE WFMESH(E,EK,CVLIGHT,NSRA,Z,R,S,RS,IRM,IRMD,LMAXD)
C     .. Scalar Arguments ..
      DOUBLE COMPLEX E,EK
      DOUBLE PRECISION CVLIGHT,Z
      INTEGER IRM,IRMD,LMAXD,NSRA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE,SQRT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION R(IRMD),RS(IRMD,0:LMAXD),S(0:LMAXD)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION S1
      INTEGER IR,L
C     ..
      IF (NSRA.EQ.1) EK = SQRT(E)
      IF (NSRA.EQ.2) EK = SQRT(E+E*E/ (CVLIGHT*CVLIGHT))
      DO L = 0,LMAXD

        IF (NSRA.EQ.2) THEN
          S1 = SQRT(DBLE(L*L+L+1)-4.0D0*Z*Z/ (CVLIGHT*CVLIGHT))
          IF (Z.EQ.0.0D0) S1 = DBLE(L)
        ELSE
          S1 = DBLE(L)
        END IF
        S(L) = S1
        RS(1,L) = 0.0D0
        DO IR = 2,IRM
          RS(IR,L) = R(IR)**S1
        END DO
        DO IR = IRM+1, IRMD
           RS(IR,L) = 0.0D0
        END DO

      END DO
      RETURN
      END
