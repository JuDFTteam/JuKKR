      SUBROUTINE RITES(IFILE,NSPIN,Z,ALAT,RMT,RMTNEW,RWS,
     +                 ITITLE,R,DRDI,VISP,A,B,KXC,IRNS,
     +                 LPOT,VINS,IRC,EFERMI,VBC,ECORE,
     +                 LCORE,NCORE,
C                      new after inc.p removal
     &                 irmd, irnsd)
c ************************************************************************
C      Formatted output of potential.
c      this subroutine stores in 'ifile' the necessary results
c      (potentials e.t.c.) to start self-consistency iterations
c
c       if the sum of absolute values of an lm component of vins (non
c       spher. potential) is less than qbound this
c       component will not be stored .
c
c-----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER irmd
      INTEGER irnsd

      double precision, parameter :: QBOUND = 1d-10

C     .. Scalar Arguments ..
      DOUBLE PRECISION ALAT
      INTEGER IFILE,KXC,LPOT,NSPIN

      DOUBLE PRECISION A
      DOUBLE PRECISION B
      DOUBLE PRECISION DRDI(IRMD)
      DOUBLE PRECISION ECORE(20,2)
      DOUBLE PRECISION EFERMI
      DOUBLE PRECISION R(IRMD)
      DOUBLE PRECISION RMT
      DOUBLE PRECISION RMTNEW
      DOUBLE PRECISION RWS
      DOUBLE PRECISION VBC(2)
      DOUBLE PRECISION VINS(IRMD-IRNSD:IRMD,(LPOT+1)**2,2)
      DOUBLE PRECISION VISP(IRMD,2)
      DOUBLE PRECISION Z
      INTEGER IRC
      INTEGER IRNS
      INTEGER ITITLE(20,NSPIN)
      INTEGER LCORE(20,NSPIN)
      INTEGER NCORE(NSPIN)

      CHARACTER(LEN=24) TXC(4)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION RV,SUM
      INTEGER I,ICORE,INEW,IR,IRMIN,IS,ISAVE,LM,LMNR,
     +        LMPOT,NR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     ..
c -------------------------------------------------------------------
      ISAVE = 1
      INEW  = 1
C
      TXC(1) = ' Morruzi,Janak,Williams '
      TXC(2) = ' von Barth,Hedin        '
      TXC(3) = ' Vosko,Wilk,Nusair      '
      TXC(4) = ' GGA PW91               '
c
      LMPOT = (LPOT+1)* (LPOT+1)

      DO 50 IS = 1,NSPIN

        NR = IRC

        IRMIN = NR - IRNS
c
        WRITE (IFILE,FMT=9000) (ITITLE(I,IS),I=1,7),TXC(KXC+1)
        WRITE (IFILE,FMT=9010) RMT,ALAT,RMTNEW
        WRITE (IFILE,FMT=9020) Z,RWS,EFERMI,VBC(IS)
        WRITE (IFILE,FMT=9030) NR,A,B,NCORE(IS),INEW
C
        IF (NCORE(IS).GE.1) THEN
             WRITE (IFILE,FMT=9040) (LCORE(ICORE,IS),
     +            ECORE(ICORE,IS),ICORE=1,NCORE(IS))
        END IF
c
c
c--->   store the full potential , but the non spherical contribution
c       only from irns1 up to irws1 ;
c       remember that the lm = 1 contribution is multiplied
c       by a factor 1/sqrt(4 pi)
c
        WRITE (IFILE,FMT=9060) NR,IRNS,LMPOT,ISAVE
        WRITE (IFILE,FMT=9070) (VISP(IR,IS),IR=1,NR)
        IF (LPOT.GT.0) THEN
          LMNR = 1
          DO 40 LM = 2,LMPOT
            SUM = 0.0D0
            DO 30 IR = IRMIN,NR
              RV = VINS(IR,LM,IS)*R(IR)
              SUM = SUM + RV*RV*DRDI(IR)
   30       CONTINUE

           IF (SQRT(SUM).GT.QBOUND) THEN
              LMNR = LMNR + 1
              WRITE (IFILE,FMT=9060) LM
              WRITE (IFILE,FMT=9070) (VINS(IR,LM,IS),IR=IRMIN,NR)
           END IF

   40     CONTINUE
c
c--->         write a one to mark the end
c
          IF (LMNR.LT.LMPOT) WRITE (IFILE,FMT=9060) ISAVE
        END IF

   50 CONTINUE



 9000 FORMAT (7a4,6x,'  exc:',a24,3x,a10)
 9010 FORMAT (3f12.8)
 9020 FORMAT (f10.5,/,f10.5,2f15.10)
 9030 FORMAT (i3,/,2d15.8,/,2i2)
 9040 FORMAT (i5,1p,d20.11)
 9060 FORMAT (10i5)
 9070 FORMAT (1p,4d20.13)
      END
