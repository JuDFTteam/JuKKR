      SUBROUTINE REGNS_VOLTERRA(AR,BR,EFAC,PNS,VNSPLL,ICST,IPAN,IRCUT,
     +              PZLM,QZLM,PZEKDR,QZEKDR,EK,AMAT,BMAT,NSRA,
     +              IRMIND,IRMD,IPAND,LMMAXD)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     determines the regular non spherical wavefunctions , the
c       alpha matrix and the t - matrix in the n-th. born appro-
c       ximation ( n given by input parameter icst )
c
c
c     using the wave functions pz and qz ( regular and irregular
c       solution ) of the spherically averaged potential , the
c       regular wavefunction pns is determined by
c
c           pns(ir,lm1,lm2) = ar(ir,lm1,lm2)*pz(ir,l1)
c                                   + br(ir,lm1,lm2)*qz(ir,l1)
c
c      the matrices ar and br are determined by integral equations
c        containing pns and only the non spherical contributions of
c        the potential , stored in vinspll . these integral equations
c        are  solved iteratively with born approximation up to given n.
c
c     the original way of writing the cr and dr matrices in the equa-
c        tions above caused numerical troubles . therefore here are used
c        rescaled ar and br matrices :
c
c              ~
c              ar(ir,lm1,lm2) = sqrt(e)**(l1-l2)
c                             * ar(ir,lm1,lm2)*((2*l2-1)!!/(2*l1-1)!!)
c
c              ~
c              br(ir,lm1,lm2) = sqrt(e)**(-l1-l2)
c                             * br(ir,lm1,lm2)/((2*l1-1)!!*(2*l2-1)!!)
c
c     for lloyd's formular is only the determinant of the alpha -
c        matrix is needed which is identical with the determinant
c        of the rescaled ar - matrix at the innerst point .
c
c     the non spherical t - matrix is the br matrix at r(irc)
c
c     modified for the use of shape functions
c
c                              (see notes by b.drittler)
c
c                                b.drittler   mar.  1989
c-----------------------------------------------------------------------
c     modified by R. Zeller      Aug. 1994
c-----------------------------------------------------------------------
c     added Volterra equation by M. Ogura      Jan. 2006
c     FRED: true -> use fredholm equation
c           false -> volterra equation
c-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE COMPLEX EK
      INTEGER ICST,IPAN,IPAND,IRMD,IRMIND,LMMAXD,NSRA
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX ADER(LMMAXD,LMMAXD,IRMIND:IRMD),
     +               AMAT(LMMAXD,LMMAXD,IRMIND:IRMD),AR(LMMAXD,LMMAXD),
     +               BDER(LMMAXD,LMMAXD,IRMIND:IRMD),
     +               BMAT(LMMAXD,LMMAXD,IRMIND:IRMD),BR(LMMAXD,LMMAXD),
     +               EFAC(*),PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),
     +               PZEKDR(LMMAXD,IRMIND:IRMD,2),
     +               PZLM(LMMAXD,IRMIND:IRMD,2),
     +               QZEKDR(LMMAXD,IRMIND:IRMD,2),
     +               QZLM(LMMAXD,IRMIND:IRMD,2)
      DOUBLE PRECISION VNSPLL(LMMAXD,LMMAXD,IRMIND:IRMD)
      INTEGER IRCUT(0:IPAND)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX EFAC1,EFAC2
      DOUBLE PRECISION ERR
      INTEGER I,IR,IRC1,J,LM1,LM2,LM3,IP
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX PNS0(LMMAXD,LMMAXD,IRMIND:IRMD,2),
     +               PNS1(LMMAXD,LMMAXD,IRMIND:IRMD)
      INTEGER IPIV(LMMAXD)
C     ..
C     .. External Subroutines ..
      EXTERNAL CSINWD,CSOUT,WFINT,WFINT0,ZGEINV1
C     ..
C     .. Parameters ..
      DOUBLE COMPLEX CONE
      PARAMETER (CONE= (1.0D0,0.0D0))
C     ..
      LOGICAL FRED
      DATA FRED/.false./
c      DATA FRED/.true. /
C     ..
c      WRITE(6 , *) "in REGNS_volterra"
c      write(*,*)ek
      IRC1 = IRCUT(IPAN)
      

c     WRITE(6,*) "IRMD",IRMD
c     WRITE(6,*) "IRC1",IRC1
c     DO IR=IRMIND,IRC1
c       DO LM1=1,LMMAXD
c         WRITE(36,"((2I5),(8e17.9))") IR,LM1,
c    +            QZEKDR(LM1,IR,1),QZEKDR(LM1,IR,2),
c    +            PZEKDR(LM1,IR,1),PZEKDR(LM1,IR,2)
c       END DO
c     END DO

cc      STOP " " 

      DO 1 J = 1,NSRA
        DO 2 IR = IRMIND,IRC1
          DO 3 LM1 = 1,LMMAXD
            DO 4 LM2 = 1,LMMAXD
              IF(LM1.EQ.LM2)THEN
                PNS0(LM1,LM2,IR,J) =  PZLM(LM1,IR,J)
              ELSE
                PNS0(LM1,LM2,IR,J) = (0D0,0D0)
              ENDIF
 4          CONTINUE
 3        CONTINUE
 2      CONTINUE
 1    CONTINUE

c      FRED=.TRUE.

      IF(FRED)THEN
      DO 70 I = 0,ICST
c---> set up integrands for i-th born approximation
c        WRITE(6,*) "I",I
        IF (I.EQ.0) THEN
          CALL WFINT0(ADER,BDER,PZLM,QZEKDR,PZEKDR,VNSPLL,IRMIND,
     +                  IRMD,NSRA)
c         OPEN (unit=30, file="ADER_NOT_ZERO", form="formatted")
c         OPEN (unit=31, file="ADER", form="formatted")
c         OPEN (unit=32, file="BDER", form="formatted")
c         OPEN (unit=33, file="VNSPLL0", form="formatted")
c         DO IR=IRMIND,IRC1
c           DO LM2=1,LMMAXD
c             DO LM1=1,LMMAXD
c               IF( ABS(ADER(LM1,LM2,IR)-ADER(LM2,LM1,IR)) .NE. 0 ) THEN
c                 WRITE(30,"((3I5),(8e17.9))") LM2,LM1,IR,
c    +              ADER(LM1,LM2,IR),
c    +              (ADER(LM1,LM2,IR)-ADER(LM2,LM1,IR)),
c    +              ABS(ADER(LM1,LM2,IR)-ADER(LM2,LM1,IR))
c               END IF
c               IF( ABS(VNSPLL(LM1,LM2,IR)-VNSPLL(LM2,LM1,IR)).NE.0) 
c    +                                                            THEN
c                 WRITE(33,"((3I5),(8e17.9))") LM2,LM1,IR,
c    +                VNSPLL(LM1,LM2,IR),
c    +                (VNSPLL(LM1,LM2,IR)-VNSPLL(LM2,LM1,IR)),
c    +                ABS(VNSPLL(LM1,LM2,IR)-VNSPLL(LM2,LM1,IR))
c               END IF
c               WRITE(31,"((3I5),(8e17.9))") LM2,LM1,IR,
c    +              ADER(LM1,LM2,IR),
c    +              (ADER(LM1,LM2,IR)-ADER(LM2,LM1,IR)),
c    +              ABS(ADER(LM1,LM2,IR)-ADER(LM2,LM1,IR))
c               WRITE(32,"((3I5),(8e17.9))") LM2,LM1,IR,
c    +              BDER(LM1,LM2,IR),
c    +              (BDER(LM1,LM2,IR)-BDER(LM2,LM1,IR)),
c    +              ABS(BDER(LM1,LM2,IR)-BDER(LM2,LM1,IR))
c             END DO
c           END DO
c         END DO
c         CLOSE(30)
c         CLOSE(31)
c         CLOSE(32)
c         CLOSE(33)
c         WRITE(6,*) "after WFINT0_SO"
c           STOP " " 

        ELSE
c         WRITE(6,*) "before wfint"
          CALL WFINT(PNS,ADER,BDER,QZEKDR,PZEKDR,VNSPLL,IRMIND,
     +                 IRMD,NSRA)
c         WRITE(6,*) "after wfint"
c         OPEN (unit=30, file="ADER_NOT_ZERO", form="formatted")
c         OPEN (unit=31, file="ADER", form="formatted")
c         OPEN (unit=32, file="BDER", form="formatted")
c         OPEN (unit=33, file="VNSPLL0", form="formatted")
c         DO IR=IRMIND,IRC1
c           DO LM2=1,LMMAXD
c             DO LM1=1,LMMAXD
c               IF( ABS(ADER(LM1,LM2,IR)-ADER(LM2,LM1,IR)) .NE. 0 ) THEN
c                 WRITE(30,"((3I5),(8e17.9))") LM2,LM1,IR,
c    +              ADER(LM1,LM2,IR),
c    +              (ADER(LM1,LM2,IR)-ADER(LM2,LM1,IR)),
c    +              ABS(ADER(LM1,LM2,IR)-ADER(LM2,LM1,IR))
c               END IF
c               IF( ABS(VNSPLL(LM1,LM2,IR)-VNSPLL(LM2,LM1,IR)).NE.0) 
c    +                                                            THEN
c                 WRITE(33,"((3I5),(8e17.9))") LM2,LM1,IR,
c    +                VNSPLL(LM1,LM2,IR),
c    +                (VNSPLL(LM1,LM2,IR)-VNSPLL(LM2,LM1,IR)),
c    +                ABS(VNSPLL(LM1,LM2,IR)-VNSPLL(LM2,LM1,IR))
c               END IF
c               WRITE(31,"((3I5),(8e17.9))") LM2,LM1,IR,
c    +              ADER(LM1,LM2,IR),
c    +              (ADER(LM1,LM2,IR)-ADER(LM2,LM1,IR)),
c    +              ABS(ADER(LM1,LM2,IR)-ADER(LM2,LM1,IR))
c               WRITE(32,"((3I5),(8e17.9))") LM2,LM1,IR,
c    +              BDER(LM1,LM2,IR),
c    +              (BDER(LM1,LM2,IR)-BDER(LM2,LM1,IR)),
c    +              ABS(BDER(LM1,LM2,IR)-BDER(LM2,LM1,IR))
c             END DO
c           END DO
c         END DO
c         CLOSE(30)
c         CLOSE(31)
c         CLOSE(32)
c         CLOSE(33)
c         WRITE(6,*) "after WFINT_SO"
c          STOP " " 
        END IF
c---> call integration subroutines
        CALL CSINWD(ADER,AMAT,LMMAXD**2,IRMIND,IRMD,IPAN,IRCUT)
        CALL CSOUT(BDER,BMAT,LMMAXD**2,IRMIND,IRMD,IPAN,IRCUT)
        DO 20 IR = IRMIND,IRC1
          DO 10 LM2 = 1,LMMAXD
            AMAT(LM2,LM2,IR) = CONE + AMAT(LM2,LM2,IR)
   10     CONTINUE
   20   CONTINUE
c---> calculate non sph. wft. in i-th born approximation
c       OPEN (unit=30, file="PNS0", form="formatted")
c       OPEN (unit=31, file="AMAT", form="formatted")
c       OPEN (unit=32, file="PZLM", form="formatted")
c       WRITE(30, * ) "I",I
c       WRITE(31, * ) "I",I
        DO 60 J = 1,NSRA
          DO 50 IR = IRMIND,IRC1
            DO 40 LM1 = 1,LMMAXD
c             WRITE(32,"((3I5),(4e17.9))") J,IR,LM1, 
c    +                           PZLM(LM1,IR,J),QZLM(LM1,IR,J)
              DO 30 LM2 = 1,LMMAXD
                PNS(LM1,LM2,IR,J) = (AMAT(LM1,LM2,IR)*PZLM(LM1,IR,J)+
     +                              BMAT(LM1,LM2,IR)*QZLM(LM1,IR,J))
c               WRITE(30,"((4I5),(2e17.9))") J,IR,LM1,LM2, 
c    +                                      PNS(LM1,LM2,IR,J)
c               IF (J==1) THEN
c                 WRITE(31,"((3I5),(4e17.9))") IR,LM1,LM2, 
c    +                             AMAT(LM1,LM2,IR),BMAT(LM1,LM2,IR)
c               END IF
   30         CONTINUE
   40       CONTINUE
   50     CONTINUE
   60   CONTINUE
c        CLOSE(30)
c        CLOSE(31)
c        CLOSE(32)
c         STOP " "
c-----------------------------------------------------------------------
c check convergence
      DO J = 1,NSRA
        DO IR = IRMIND,IRC1
          DO LM1 = 1,LMMAXD
            DO LM2 = 1,LMMAXD
              PNS0(LM1,LM2,IR,J) = PNS0(LM1,LM2,IR,J)-PNS(LM1,LM2,IR,J)
            END DO
          END DO
        END DO
      END DO

      ERR=0D0

      PNS1=0d0
      DO J=1,NSRA
        CALL CSOUT(PNS0(:,:,:,J),PNS1,LMMAXD**2,IRMIND,IRMD,IPAN,
     +           IRCUT)
        DO LM2=1,LMMAXD
          DO LM1=1,LMMAXD
            ERR=MAX(ERR,ABS(PNS1(LM1,LM2,IRC1)))
          END DO
        END DO
      END DO

      WRITE(*,*) 'Born_Fred',I,ERR
      IF(ERR .LT. 1D-7) EXIT
c      IF(I.EQ.ICST.AND.ERR.GT.1D-3)WRITE(*,*)'NOT CONVERGENT',ERR
      DO 280 J = 1,NSRA
      DO 280 IR = IRMIND,IRC1
      DO 280 LM2 = 1,LMMAXD
      DO 280 LM1 = 1,LMMAXD
 280  PNS0(LM1,LM2,IR,J) = PNS(LM1,LM2,IR,J)
c-----------------------------------------------------------------------
   70 CONTINUE
      ELSE
c-----------------------------------------------------------------------
c Volterra equation
      DO 200 I = 0,ICST
c      DO 200 I = 0,1
c       WRITE(160, *) "I",I
c       WRITE(34, *) "I",I
c      write(*,*)ek
c---> set up integrands for i-th born approximation
        IF (I.EQ.0) THEN

          CALL WFINT0(ADER,BDER,PZLM,QZEKDR,PZEKDR,VNSPLL,IRMIND,
     +                  IRMD,NSRA)
        ELSE
          CALL WFINT(PNS,ADER,BDER,QZEKDR,PZEKDR,VNSPLL,IRMIND,
     +                 IRMD,NSRA)
c          CALL WFINT(ADER,BDER,QZEKDR,PZEKDR,VNSPLL,IRMIND,
c     +                 IRMD,NSRA,PNS)
        END IF
c       DO IR = IRMIND,IRC1
c         DO LM1 = 1,LMMAXD
c           DO LM2 = 1,LMMAXD
c             WRITE(332,"((3I5),(4e17.9))") IR,LM2,LM1,
c    +                 ADER(LM1,LM2,IR),BDER(LM2,LM1,IR)
c           END DO
c         END DO
c       END DO
c       DO IR=IRMIND,IRC1
c         DO LM1=1,LMMAXD
c           DO LM2=1,LMMAXD
c             WRITE(34,"((3I5),(4e17.9))") IR,LM1,LM2,ADER(LM2,LM1,IR),
c    +                                     BDER(LM2,LM1,IR)
c           END DO
c         END DO
c       END DO
c---> call integration subroutines
        CALL CSOUT(ADER,AMAT,LMMAXD**2,IRMIND,IRMD,IPAN,IRCUT)
        CALL CSOUT(BDER,BMAT,LMMAXD**2,IRMIND,IRMD,IPAN,IRCUT)
        DO 150 IR = IRMIND,IRC1
          DO 140 LM1 = 1,LMMAXD
          DO 140 LM2 = 1,LMMAXD
            IF(LM1.EQ.LM2)THEN
            AMAT(LM1,LM2,IR) = CONE - AMAT(LM1,LM2,IR)
            ELSE
            AMAT(LM1,LM2,IR) = - AMAT(LM1,LM2,IR)
            ENDIF
c            WRITE(33,"((3I5),(4e17.9))") IR,LM1,LM2,AMAT(LM2,LM1,IR),
c     +                                   BMAT(LM2,LM1,IR)
  140     CONTINUE
  150   CONTINUE
c---> calculate non sph. wft. in i-th born approximation
        DO 190 J = 1,NSRA
          DO 180 IR = IRMIND,IRC1
            DO 170 LM2 = 1,LMMAXD
              DO 160 LM1 = 1,LMMAXD
                PNS(LM1,LM2,IR,J) =  AMAT(LM1,LM2,IR)*PZLM(LM1,IR,J)
     +                              +BMAT(LM1,LM2,IR)*QZLM(LM1,IR,J)
c                WRITE(160, "(2e17.9)") PNS(LM1,LM2,IR,J)
  160         CONTINUE
  170       CONTINUE
  180     CONTINUE
  190   CONTINUE
c-----------------------------------------------------------------------
c check convergence
       DO 290 J = 1,NSRA
       DO 290 IR = IRMIND,IRC1
       DO 290 LM2 = 1,LMMAXD
       DO 290 LM1 = 1,LMMAXD
  290  PNS0(LM1,LM2,IR,J) = PNS0(LM1,LM2,IR,J)-PNS(LM1,LM2,IR,J)
       ERR=0D0
       DO 300 J=1,NSRA
       CALL CSOUT(PNS0(1,1,IRMIND,J),PNS1,LMMAXD**2,IRMIND,IRMD,IPAN,
     +           IRCUT)
       DO 300 LM2=1,LMMAXD
       DO 300 LM1=1,LMMAXD
  300  ERR=MAX(ERR,ABS(PNS1(LM1,LM2,IRC1)))
       WRITE(*,*) 'Born Volterra',I,ERR 
       IF(ERR .LT. 1D-7) EXIT
c      IF(I.EQ.ICST.AND.ERR.GT.1D-3)WRITE(*,*)'NOT CONVERGENT',ERR
       DO 310 J = 1,NSRA
       DO 310 IR = IRMIND,IRC1
       DO 310 LM2 = 1,LMMAXD
       DO 310 LM1 = 1,LMMAXD
  310  PNS0(LM1,LM2,IR,J) = PNS(LM1,LM2,IR,J)
c-----------------------------------------------------------------------
  200 CONTINUE
c      STOP " " 
      CALL ZGEINV1(AMAT(1,1,IRC1),AR,BR,IPIV,LMMAXD)
      DO 220 LM2=1,LMMAXD
      DO 220 LM1=1,LMMAXD
      DO 210 IR=IRMIND,IRC1
        ADER(LM1,LM2,IR)=(0d0,0d0)
  210   BDER(LM1,LM2,IR)=(0d0,0d0)
      DO 220 LM3=1,LMMAXD
      DO 220 IR=IRMIND,IRC1
        ADER(LM1,LM2,IR)=ADER(LM1,LM2,IR)+AMAT(LM1,LM3,IR)*AR(LM3,LM2)
  220   BDER(LM1,LM2,IR)=BDER(LM1,LM2,IR)+BMAT(LM1,LM3,IR)*AR(LM3,LM2)
      DO 230 LM2=1,LMMAXD
      DO 230 LM1=1,LMMAXD
      DO 230 IR=IRMIND,IRC1
        AMAT(LM1,LM2,IR)=ADER(LM1,LM2,IR)
  230   BMAT(LM1,LM2,IR)=BDER(LM1,LM2,IR)
      DO 240 J = 1,NSRA
      DO 240 IR = IRMIND,IRC1
      DO 240 LM2 = 1,LMMAXD
      DO 240 LM1 = 1,LMMAXD
        PNS(LM1,LM2,IR,J) =  AMAT(LM1,LM2,IR)*PZLM(LM1,IR,J)
     +                      +BMAT(LM1,LM2,IR)*QZLM(LM1,IR,J)
  240 CONTINUE
c Volterra equation
c-----------------------------------------------------------------------
      ENDIF

c     DO IR = IRMIND,IRC1
c       DO LM2 = 1,LMMAXD
c         DO LM1 = 1,LMMAXD
c           WRITE(334,"((3I5),(4e17.9))") LM1,LM2,IR,
c    +          PNS(LM1,LM2,IR,1)
c         END DO
c       END DO
c     END DO

      DO 90 LM2 = 1,LMMAXD
        EFAC2 = EFAC(LM2)
c---> store alpha and t - matrix
        DO 80 LM1 = 1,LMMAXD
          EFAC1 = EFAC(LM1)
          AR(LM1,LM2) = AMAT(LM1,LM2,IRMIND)/EFAC1*EFAC2!*EK
c          AR(LM1,LM2) = AMAT(LM1,LM2,IRMIND)
c---> t-matrix
          BR(LM1,LM2) = BMAT(LM1,LM2,IRC1)*EFAC1*EFAC2/EK
   80   CONTINUE
   90 CONTINUE
c---> rescale with efac
      DO 130 J = 1,NSRA
        DO 120 LM2 = 1,LMMAXD
          EFAC2 = EFAC(LM2)
          DO 110 IR = IRMIND,IRC1
            DO 100 LM1 = 1,LMMAXD
              PNS(LM1,LM2,IR,J) = PNS(LM1,LM2,IR,J)*EFAC2
  100       CONTINUE
  110     CONTINUE
  120   CONTINUE
  130 CONTINUE

      END
c ************************************************************************
c     SUBROUTINE ZGEINV1(A,U,AUX,IPIV,DIM)
c ************************************************************************
c   - inverts a general double complex matrix A,
c   - the result is return in U,
c   - input matrix A is returned unchanged,
c   - AUX is a auxiliary matrix,
c   - A,U and AUX are of dimension (DIM,DIM),
c ------------------------------------------------------------------------
c     INTEGER DIM,IPIV(*)
c     DOUBLE COMPLEX A(DIM,*),AUX(DIM,*),U(DIM,*)
c
c     .. PARAMETER
c
c     DOUBLE COMPLEX CONE
c     PARAMETER(CONE=(1.D0,0.D0))
c
c     INTEGER LM1,INFO
c     EXTERNAL ZCOPY,ZGETRS,ZGETRF
c ------------------------------------------------------------------------
c     CALL CINIT(DIM*DIM,U)
c     DO 10 LM1=1,DIM
c       U(LM1,LM1) = CONE
c10   END DO

c     CALL ZCOPY(DIM*DIM,A,1,AUX,1)
c     CALL ZGETRF(DIM,DIM,AUX,DIM,IPIV,INFO)
c     CALL ZGETRS('N',DIM,DIM,AUX,DIM,IPIV,U,DIM,INFO)
c     
c     RETURN
c     END
