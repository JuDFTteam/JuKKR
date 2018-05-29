      SUBROUTINE RITES(IFILE,NATPS,NATYP,NSPIN,Z,ALAT,RMT,RMTNEW,RWS,
     +                 ITITLE,R,DRDI,VM2Z,IRWS,A,B,TXC,KXC,INS,IRNS,
     +                 LPOT,VINS,QBOUND,IRC,KSHAPE,EFERMI,VBC,ECORE,
     +                 LCORE,NCORE,ECOREREL,NKCORE,KAPCORE)
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
      include 'inc.p'
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALAT,QBOUND
      INTEGER IFILE,INS,KSHAPE,KXC,LPOT,NATPS,NATYP,NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*),B(*),DRDI(IRMD,*),ECORE(20,*),EFERMI, !(2), 22.5,2000
     +                 R(IRMD,*),RMT(*),RMTNEW(*),RWS(*),VBC(2),
     +                 VINS(IRMIND:IRMD,LMPOTD,*),VM2Z(IRMD,*),Z(*)
C ===================================================================
      DOUBLE PRECISION ECOREREL(KREL*20+(1-KREL),
     +                          2*NATYPD)  ! relativistic core energies
      INTEGER NKCORE(20,NATYPD),KAPCORE(20,2*NATYPD)
C ===================================================================
      INTEGER IRC(*),IRNS(*),IRWS(*),ITITLE(20,*),LCORE(20,*),NCORE(*)
      CHARACTER*124 TXC(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A1,B1,RMAX,RMT1,RMTNW1,RV,SIGN,SUM,Z1
      INTEGER I,ICORE,IH,INEW,IP,IR,IRMIN,IRNS1,IS,ISAVE,J,LM,LMNR,
     +        LMPOT,NCORE1,NR
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DRADI(IRMD),ECORE1(20),RA(IRMD),VM2ZA(IRMD)
      DOUBLE PRECISION ECORE2(20,2)
      INTEGER LCORE1(20)
      CHARACTER*3 TXTK(4)
      CHARACTER*1 TXTL(0:3)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     ..
      DATA TXTL/'s','p','d','f'/
      DATA TXTK/'1/2','3/2','5/2','7/2'/
C     ..
c -------------------------------------------------------------------
      ISAVE = 1
      INEW  = 1
c
c
      LMPOT = (LPOT+1)* (LPOT+1)
      DO 60 IH = 1,NATYP
        DO 50 IS = 1,NSPIN
          IF (IS.EQ.NSPIN) THEN
            SIGN = 1.0D0

          ELSE
            SIGN = -1.0D0
          END IF
          IP = NSPIN* (IH-1) + IS

          RMT1 = RMT(IH)
          RMTNW1 = RMTNEW(IH)
          Z1 = Z(IH)
          RMAX = RWS(IH)
          IF (KSHAPE.EQ.0) THEN
            NR = IRWS(IH)

          ELSE
            NR = IRC(IH)
          END IF

          IRNS1 = IRNS(IH)
          IRMIN = NR - IRNS1
          A1 = A(IH)
          B1 = B(IH)
          NCORE1 = NCORE(IP)
c
          DO 10 J = 1,NR
            RA(J) = R(J,IH)
            DRADI(J) = DRDI(J,IH)
c
c--->       store only lm=1 component of the potential
c
            VM2ZA(J) = VM2Z(J,IP)
   10     CONTINUE
C
          OPEN(IFILE, FILE='out_potential', FORM='formatted')

          WRITE (IFILE,FMT=9000) (ITITLE(I,IP),I=1,7),TXC(KXC+1)
          WRITE (IFILE,FMT=9010) RMT1,ALAT,RMTNW1
          WRITE (IFILE,FMT=9020) Z1,RMAX,EFERMI,VBC(IS)
          WRITE (IFILE,FMT=9030) NR,A1,B1,NCORE1,INEW
C
          IF (NCORE1.GE.1) THEN
C
            IF (KREL.EQ.0) THEN 
               DO J = 1,NCORE1
                 LCORE1(J) = LCORE(J,IP)
                 ECORE1(J) = ECORE(J,IP)
               END DO
               WRITE (IFILE,FMT=9040) (LCORE1(ICORE),
     +              ECORE1(ICORE),ICORE=1,NCORE1)
            ELSE
              DO J = 1,NCORE1
                 LCORE1(J) = LCORE(J,IP)
                 ECORE2(J,1) = ECOREREL(J,2*IH-1)
                 ECORE2(J,2) = ECOREREL(J,2*IH)
              END DO
C
C --> independent of spin, the \mu-averaged relativistic core energies
C     are written out for \kappa = -l-1,l 
C     format compatible with the non-(scalar) relativistic mode
C     however, the next read in has no meaning for the REL core-solver
C     a detailed output of the core energies is supplied by < CORE >
C
              DO ICORE=1,NCORE1
                 WRITE (IFILE,FMT=9041) LCORE1(ICORE),(
     +                ECORE2(ICORE,I+1),
     +                TXTL(LCORE1(ICORE)),
     +                TXTK(IABS(KAPCORE(ICORE,2*IH-1+I))),
     +                I=0,NKCORE(ICORE,IH)-1)
              END DO
            END IF
C
          END IF
c
          IF (INS.EQ.0 .OR. (IH.LT.NATPS.AND.INS.LE.2)) THEN
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
                  RV = VINS(IR,LM,IP)*RA(IR)
                  SUM = SUM + RV*RV*DRADI(IR)
   30           CONTINUE

                IF (SQRT(SUM).GT.QBOUND) THEN
                  LMNR = LMNR + 1
                  WRITE (IFILE,FMT=9060) LM
                  WRITE (IFILE,FMT=9070) (VINS(IR,LM,IP),IR=IRMIN,NR)
                END IF

   40         CONTINUE
c
c--->         write a one to mark the end
c
              IF (LMNR.LT.LMPOT) WRITE (IFILE,FMT=9060) ISAVE
            END IF

          END IF

   50   CONTINUE
   60 CONTINUE


      close(IFILE)

 9000 FORMAT (7a4,6x,'  exc:',a124,3x,a10)
 9010 FORMAT (3f12.8)
! 9010 FORMAT (3F) !f12.8)
 9020 FORMAT (f10.5,/,f10.5,2f20.15)
 9030 FORMAT (i3,/,2d15.8,/,2i2)
 9040 FORMAT (i5,1p,d20.11)
 9041 FORMAT (i5,2(1p,d20.11,2x,A1,A3))
 9050 FORMAT (1p,2d15.6,1p,d15.8)
 9051 FORMAT (1p,4d20.12)
 9060 FORMAT (10i5)
 9070 FORMAT (1p,4d20.13)
      END
