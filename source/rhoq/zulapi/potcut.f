C ************************************************************************
      SUBROUTINE POTCUT(IMT1,IRC1,INS,LMPOT,R,VM2Z,VINS,Z1)
C ************************************************************************
c
c     set potential equal zero between muffin tin sphere and
c       outer sphere
c
C ************************************************************************
C     .. Parameters ..
      include 'inc.p'
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION Z1
      INTEGER IMT1,INS,IRC1,LMPOT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION R(*),VINS(IRMIND:IRMD,*),VM2Z(*)
C     ..
C     .. Local Scalars ..
      INTEGER IR,IST,LM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX
C     ..
      write(6,*) 'potcut: potential equal 2*Z/R between MT ',
     +           'and outer sphere'
      DO 10 IR = IMT1 + 1,IRC1
        VM2Z(IR) = 2.0D0*Z1/R(IR)
   10 CONTINUE
c
      IF (INS.GE.1) THEN
        IST = MAX(IRMIND,IMT1+1)
        DO 30 IR = IST,IRC1
          DO 20 LM = 2,LMPOT
            VINS(IR,LM) = 0.0D0
   20     CONTINUE
   30   CONTINUE
      END IF

      END                           ! SUBROUTINE POTCUT
