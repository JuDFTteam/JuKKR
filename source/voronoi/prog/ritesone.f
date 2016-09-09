      SUBROUTINE RITESONE(IFILE,IS,Z,ALAT,RMT,RMTNEW,RWS,
     +                 ITITLE,R,DRDI,VM2Z,IRWS,A,B,TXC,KXC,INS,IRNS,
     +                 LPOT,VINS,QBOUND,IRC,KSHAPE,EFERMI,VBC,ECORE,
     +                 LCORE,NCORE,ELEM_NAME,NSPIN)
c ************************************************************************
c      this subroutine stores in 'ifile' the necessary results
c      (potentials e.t.c.) to start self-consistency iterations
c
c      modified for the full potential case - if ins .gt. 0 there
c       is written a different potential card
c       if the sum of absolute values of an lm component of vins (non
c       spher. potential) is less than the given rms error qbound this
c       component will not be stored .
c
c        (see to subroutine start , where most of the arrays are
c         described)
c
c                            modified by b. drittler  aug. 1988
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.geometry'
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
C     ..
C     .. Scalar Arguments ..
      REAL*8           ALAT,QBOUND
      INTEGER IFILE,INS,KSHAPE,KXC,LPOT,NATPS,NATYP,NSPIN
C     ..
C     .. Array Arguments ..
      REAL*8           A,B,DRDI(IRMD),ECORE(20),EFERMI,
     +                 R(IRMD),RMT,RMTNEW,RWS,VBC(2),
     +                 VINS(IRMIND:IRMD,LMPOTD),VM2Z(IRMD),Z
      INTEGER IRC,IRNS,IRWS,ITITLE(20),LCORE(20),NCORE
      CHARACTER*124 TXC(3)
      CHARACTER*4 ELEM_NAME
C     ..
C     .. Local Scalars ..
      REAL*8           A1,B1,RMAX,RMT1,RMTNW1,RV,SIGN,SUM,Z1
      INTEGER I,ICORE,IH,INEW,IP,IR,IRMIN,IRNS1,IS,ISAVE,J,LM,LMNR,LMPOT,
     +        NCORE1,NR
C     ..
C     .. Local Arrays ..
      REAL*8           DRADI(IRMD),ECORE1(20),RA(IRMD),VM2ZA(IRMD)
      INTEGER LCORE1(20)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     ..
c -------------------------------------------------------------------
      ISAVE = 1
      INEW  = 1
c
c
         
          LMPOT = (LPOT+1)* (LPOT+1)

          RMT1 = RMT
          RMTNW1 = RMTNEW
          Z1 = Z
          RMAX = RWS
          IF (KSHAPE.EQ.0) THEN
            NR = IRWS

          ELSE
            NR = IRC
          END IF

          IRNS1 = IRNS
          IRMIN = NR - IRNS1
          A1 = A
          B1 = B
          NCORE1 = NCORE
c
          DO 10 J = 1,NR
            RA(J) = R(J)
            DRADI(J) = DRDI(J)
c
c--->       store only lm=1 component of the potential
c
            VM2ZA(J) = VM2Z(J)
   10     CONTINUE
c
          IF (NCORE1.GE.1) THEN
c
            DO 20 J = 1,NCORE1
              LCORE1(J) = LCORE(J)
              ECORE1(J) = ECORE(J)
   20       CONTINUE
          END IF
c
c
          IF (NSPIN.EQ.1) THEN
             WRITE (IFILE,FMT=8999) ELEM_NAME,TXC(KXC+1)
          ELSE
             IF (IS.EQ.1) THEN
                WRITE (IFILE,FMT=8995) ELEM_NAME,TXC(KXC+1)
             ELSE
                WRITE (IFILE,FMT=8996) ELEM_NAME,TXC(KXC+1)
             END IF
          END IF
c     WRITE (IFILE,FMT=9000) (ITITLE(I),I=1,7),TXC(KXC+1)
          WRITE (IFILE,FMT=9010) RMT1,ALAT,RMTNW1
          WRITE (IFILE,FMT=9020) Z1,RMAX,EFERMI,VBC(IS)
          IF (NR.LE.999) THEN
             WRITE (IFILE,FMT=9030) NR,A1,B1,NCORE1,INEW
          ELSE
             WRITE (IFILE,FMT=9031) NR,A1,B1,NCORE1,INEW
          ENDIF
          IF (NCORE1.GE.1) WRITE (IFILE,FMT=9040) (LCORE1(ICORE),
     +        ECORE1(ICORE),ICORE=1,NCORE1)
c
          
          IF (INS.EQ.0 ) THEN
c
c--->       store only the spherically averaged potential 
c           (in mt or as - case) 
c           this is done always for the host
c
            IF (INEW.EQ.0) THEN
              WRITE (IFILE,FMT=9050)
     +             (RA(IR),DRADI(IR),VM2ZA(IR),IR=1,NR)
            ELSE
              WRITE (IFILE,FMT=9051) (VM2ZA(IR),IR=1,NR)
            END IF

          ELSE
c
c--->     store the full potential , but the non spherical contribution
c         only from irns1 up to irws1 ;
c         remember that the lm = 1 contribution is multiplied
c         by a factor 1/sqrt(4 pi)
c
            WRITE (IFILE,FMT=9060) NR,IRNS1,LMPOT,ISAVE
            WRITE (IFILE,FMT=9070) (VM2ZA(IR),IR=1,NR)
            IF (LPOT.GT.0) THEN
              LMNR = 1
              DO 40 LM = 2,LMPOT
                SUM = 0.0D0
                DO 30 IR = IRMIN,NR
                  RV = VINS(IR,LM)*RA(IR)
                  SUM = SUM + RV*RV*DRADI(IR)
   30           CONTINUE

                IF (SQRT(SUM).GT.QBOUND) THEN
                  LMNR = LMNR + 1
                  WRITE (IFILE,FMT=9060) LM
                  WRITE (IFILE,FMT=9070) (VINS(IR,LM),IR=IRMIN,NR)
                END IF

   40         CONTINUE
c
c--->         write a one to mark the end
c
              IF (LMNR.LT.LMPOT) WRITE (IFILE,FMT=9060) ISAVE
            END IF

          END IF
          !write(6,*) ' Potential finished go on'
   50   CONTINUE
   60 CONTINUE
 8995 format (A4,' POTENTIAL SPIN DOWN',10X,'  exc:',a124)
 8996 format (A4,' POTENTIAL SPIN UP  ',10X,'  exc:',a124)  
 8999 format (A4,' POTENTIAL ',19X,'  exc:',a124)
 9000 FORMAT (7a4,6x,'  exc:',a24,3x,a10)
 9010 FORMAT (3f12.8)
 9020 FORMAT (f10.5,/,f10.5,2f15.10)
 9030 FORMAT (i3,/,2d15.8,/,2i2)
 9031 FORMAT (i4,/,2d15.8,/,2i2)
 9040 FORMAT (i5,1p,d20.11)
 9050 FORMAT (1p,2d15.6,1p,d15.8)
 9051 FORMAT (1p,4d20.12)
 9060 FORMAT (10i5)
 9070 FORMAT (1p,4d20.13)
      END


