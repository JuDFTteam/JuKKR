      SUBROUTINE POTCUT(IMT1,IRC1,INS,LMPOT,R,VM2Z,VSPSME,VINS,Z1,
     &                  IRMD,IRMIND)
C **********************************************************************
C * set potential equal zero between muffin-tin and outer sphere       *
C **********************************************************************
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION Z1
      INTEGER IRMD,IRMIND
      INTEGER IMT1,INS,IRC1,LMPOT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION R(*),VINS(IRMIND:IRMD,*),VM2Z(*),VSPSME(*)
C     ..
C     .. Local Scalars ..
      INTEGER IR,IST,LM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX
C     ..
      write(1337,*) 'potcut: potential equal 2*Z/R between MT ',
     +           'and outer sphere'
      DO IR = IMT1 + 1,IRC1
         VM2Z(IR) = 2.0D0*Z1/R(IR)
         VSPSME(IR) = 2.0D0*Z1/R(IR)
      END DO
C
      IF (INS.GE.1) THEN
         IST = MAX(IRMIND,IMT1+1)
         DO IR = IST,IRC1
            DO LM = 2,LMPOT
               VINS(IR,LM) = 0.0D0
            END DO
         END DO
      END IF
      END                           ! SUBROUTINE POTCUT
