C ************************************************************************
      SUBROUTINE WFINT(QNS,CDER,DDER,QZEKDR,PZEKDR,VNSPLL,IRMIN,IRC1,
     +                 NSRA)
C ************************************************************************
c      determines the integrands CDER, DDER or ADER, BDER in the
c        integral equations for the non-spherical wavefunctions from
c        the non-spherical contributions of the potential vinsPLL.
c
c      R. Zeller      Aug. 1994
C ************************************************************************
C     .. Parameters ..
      include 'inc.p'
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (LMAXD+1)**2)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
C     ..
C     .. Scalar Arguments ..
      INTEGER IRC1,IRMIN,NSRA
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX CDER(LMMAXD,LMMAXD,IRMIND:IRMD),
     +               DDER(LMMAXD,LMMAXD,IRMIND:IRMD),
     +               PZEKDR(LMMAXD,IRMIND:IRMD,2),
c     +               QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),
     +               QZEKDR(LMMAXD,IRMIND:IRMD,2)
      DOUBLE COMPLEX,Intent(in) ::  QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2)
      DOUBLE PRECISION VNSPLL(LMMAXD,LMMAXD,IRMIND:IRMD)
C     ..
C     .. Local Scalars ..
      INTEGER IR,LM1,LM2
C     ..
C     .. External Subroutines ..
      EXTERNAL DGEMM
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION QNSI(LMMAXD,LMMAXD),QNSR(LMMAXD,LMMAXD),
     +                 VTQNSI(LMMAXD,LMMAXD),VTQNSR(LMMAXD,LMMAXD)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DCMPLX,DIMAG,REAL
C     ..
c      CDER(1:LMMAXD,1:LMMAXD,IRMIND:IRMD)=0d0
c      DDER(1:LMMAXD,1:LMMAXD,IRMIND:IRMD)=0d0

c       WRITE(135, *) "IRMIN",IRMIN
c       WRITE(135, *) "IRC1",IRC1 
c       DO IR=IRMIN,IRC1
c         DO LM1=1,LMMAXD
c           WRITE(135,"((2I5),(8e17.9))") IR,LM1,
c    +          QZEKDR(LM1,IR,1),QZEKDR(LM1,IR,2),
c    +          PZEKDR(LM1,IR,1),PZEKDR(LM1,IR,2)
c           DO LM2=1,LMMAXD
c             WRITE(136,"((3I5),(6e17.9))") IR,LM1,LM2,DDER(LM2,LM1,IR),
c    +                              CDER(LM2,LM1,IR),VNSPLL(LM2,LM1,IR)
c             WRITE(137,"((3I5),(6e17.9))") IR,LM1,LM2,
c    +                              QNS(LM2,LM1,IR,1),QNS(LM2,LM1,IR,2)
c           END DO
c         END DO
c       END DO
c        STOP " " 

      DO 90 IR = IRMIN,IRC1
        DO 20 LM2 = 1,LMMAXD
          DO 10 LM1 = 1,LMMAXD
c            QNSR(LM1,LM2) = REAL(QNS(LM1,LM2,IR,1))
c change for linux
             QNSR(LM1,LM2) = QNS(LM1,LM2,IR,1)
            QNSI(LM1,LM2) = DIMAG(QNS(LM1,LM2,IR,1))
   10     CONTINUE
   20   CONTINUE
        CALL DGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,1.D0,VNSPLL(1,1,IR),
     +             LMMAXD,QNSR,LMMAXD,0.D0,VTQNSR,LMMAXD)
        CALL DGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,1.D0,VNSPLL(1,1,IR),
     +             LMMAXD,QNSI,LMMAXD,0.D0,VTQNSI,LMMAXD)
        DO 40 LM1 = 1,LMMAXD
          DO 30 LM2 = 1,LMMAXD
            CDER(LM1,LM2,IR) = QZEKDR(LM1,IR,1)*
     +                         DCMPLX(VTQNSR(LM1,LM2),VTQNSI(LM1,LM2))
            DDER(LM1,LM2,IR) = PZEKDR(LM1,IR,1)*
     +                         DCMPLX(VTQNSR(LM1,LM2),VTQNSI(LM1,LM2))
   30     CONTINUE
   40   CONTINUE
        IF (NSRA.EQ.2) THEN
          DO 60 LM2 = 1,LMMAXD
            DO 50 LM1 = 1,LMMAXD
c              QNSR(LM1,LM2) = REAL(QNS(LM1,LM2,IR,2))
c change for linux
              QNSR(LM1,LM2) = QNS(LM1,LM2,IR,2)              
              QNSI(LM1,LM2) = DIMAG(QNS(LM1,LM2,IR,2))
   50       CONTINUE
   60     CONTINUE
          CALL DGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,1.D0,VNSPLL(1,1,IR),
     +               LMMAXD,QNSR,LMMAXD,0.D0,VTQNSR,LMMAXD)
          CALL DGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,1.D0,VNSPLL(1,1,IR),
     +               LMMAXD,QNSI,LMMAXD,0.D0,VTQNSI,LMMAXD)
          DO 80 LM2 = 1,LMMAXD
            DO 70 LM1 = 1,LMMAXD
              CDER(LM1,LM2,IR) = CDER(LM1,LM2,IR) +
     +                           QZEKDR(LM1,IR,2)*DCMPLX(VTQNSR(LM1,
     +                           LM2),VTQNSI(LM1,LM2))
              DDER(LM1,LM2,IR) = DDER(LM1,LM2,IR) +
     +                           PZEKDR(LM1,IR,2)*DCMPLX(VTQNSR(LM1,
     +                           LM2),VTQNSI(LM1,LM2))
   70       CONTINUE
   80     CONTINUE
        END IF

   90 CONTINUE


c       WRITE(134, *) "IRMIN",IRMIN
c       WRITE(134, *) "IRC1",IRC1 
c       DO IR=IRMIN,IRC1
c         DO LM1=1,LMMAXD
c           DO LM2=1,LMMAXD
c             WRITE(134,"((3I5),(4e17.9))") IR,LM1,LM2,DDER(LM2,LM1,IR),
c    +                                     CDER(LM2,LM1,IR)
c           END DO
c         END DO
c       END DO

c        STOP " " 

      END


