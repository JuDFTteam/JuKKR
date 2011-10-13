      SUBROUTINE RHOVAL0(EZ,WEZ,DRDI,R,IPAN,IRCUT,
     &                  THETAS,DOS0,DOS1)
C
      IMPLICIT NONE
C
C     .. Parameters ..
      INCLUDE 'inc.p'
      INTEGER LMXSPD
      PARAMETER (LMXSPD= (2*LPOTD+1)**2)
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (LMAXD+1)**2)
      INTEGER LMAXD1
      PARAMETER (LMAXD1= LMAXD+1)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
      DOUBLE COMPLEX CONE,CZERO,CI
      PARAMETER ( CONE=(1.D0,0.D0),CZERO=(0.D0,0.D0),CI=(0.D0,1.D0) )
C     ..
C     .. Scalar Arguments ..
      INTEGER IPAN
      DOUBLE COMPLEX EZ,WEZ,DOS0,DOS1
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DRDI(IRMD),
     +                 R(IRMD),
     +                 THETAS(IRID,NFUND)
      INTEGER IRCUT(0:IPAND)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX EK,CIEK,DENL
      DOUBLE PRECISION C0LL
      INTEGER IR,L,L1,IMT1
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX PZ(IRMD,0:LMAXD),QZ(IRMD,0:LMAXD)
      DOUBLE COMPLEX BESSJW(0:LMAXD1),BESSYW(0:LMAXD1),HANKWS(0:LMAXD1)
      DOUBLE COMPLEX CDEN0(IRMD,0:LMAXD1),CDEN1(IRMD,0:LMAXD1)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,DBLE,SQRT
C     ..
C
         EK = SQRT(EZ)
         C0LL = 1.0d0/SQRT(16.0D0*ATAN(1.0D0))
         CIEK=(0.0D0,1.0D0)*EK
C
C initialize arrays ...
C
      DO IR = 1, IRMD
        DO L1 = 0, LMAXD1
          CDEN0(IR,L1) = CZERO
          CDEN1(IR,L1) = CZERO
        ENDDO
      ENDDO
C
C
C=======================================================================
      DO IR = 2,IRMD
        CALL BESHAN(HANKWS,BESSJW,BESSYW,R(IR)*EK,LMAXD1)
        DO L = 0,LMAXD
          PZ(IR,L) = BESSJW(L)*R(IR)
          QZ(IR,L) = (BESSYW(L) - CI*BESSJW(L))*R(IR)
        ENDDO
      ENDDO
      IMT1=IRCUT(1)
      DO L1 = 0,LMAXD1
          CDEN0(1,L1) = (0.0D0,0.0D0)
          CDEN1(1,L1) = (0.0D0,0.0D0)
      END DO
        DO IR = 2,IRMD
          CDEN0(IR,0) = EK*PZ(IR,0)*QZ(IR,0)
          CDEN1(IR,0) = EK*PZ(IR,0)**2*(0.D0,-1.D0)
          CDEN1(IR,LMAXD1) = CIEK*R(IR)**2
          END DO  
      DO L1 = 1,LMAXD
        DO IR = 2,IRMD
          CDEN0(IR,L1) = EK*PZ(IR,L1)*QZ(IR,L1)*(L1+L1+1)
          CDEN1(IR,L1) = EK*PZ(IR,L1)**2*(0.D0,-1.D0)*(L1+L1+1)
        END DO  
      END DO
c
      DO L1 = 0,LMAXD1
        IF (IPAN.GT.1) THEN
          DO IR = IMT1 + 1,IRMD
            CDEN0(IR,L1) = CDEN0(IR,L1)*THETAS(IR-IMT1,1)*C0LL
            CDEN1(IR,L1) = CDEN1(IR,L1)*THETAS(IR-IMT1,1)*C0LL
          END DO
        END IF
        CALL CSIMPK(CDEN0(1,L1),DENL,IPAN,IRCUT,DRDI)
        DOS0 = DOS0 + WEZ*DENL
        CALL CSIMPK(CDEN1(1,L1),DENL,IPAN,IRCUT,DRDI)
        DOS1 = DOS1 + WEZ*DENL
      END DO
C
C
      END
