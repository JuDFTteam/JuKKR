      SUBROUTINE REGNS_SO(AR,BR,EFAC,PNS,VNSPLL,LSM,HSOFAC,ICST,IPAN,
     +             IRCUT,PZLM,QZLM,PZEKDR,QZEKDR,RM,
     +             EK,NSRA,IRMIN,IRMINSO,IRWS,IPAND,LMMAXD,LMMAXSO,
     +              Z,INS,NSPIN)
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
      DOUBLE COMPLEX      EK
      DOUBLE PRECISION    Z
      INTEGER             ICST,IPAN,IPAND,IRWS,IRMIN,LMMAXD,NSRA,INS,
     +                    IRMINSO,NSPIN,LMMAXSO
C     .. Array Arguments ..
      DOUBLE COMPLEX      AR(LMMAXSO,LMMAXSO),
     +                    BR(LMMAXSO,LMMAXSO),
     +                    EFAC(*),PNS(LMMAXSO,LMMAXSO,IRMINSO:IRWS,2),
     +                    PZEKDR(LMMAXD,IRMINSO:IRWS,2,NSPIN),
     +                    PZLM(LMMAXD,IRMINSO:IRWS,2,NSPIN),
     +                    QZEKDR(LMMAXD,IRMINSO:IRWS,2,NSPIN),
     +                    QZLM(LMMAXD,IRMINSO:IRWS,2,NSPIN),
     +                    LSM(LMMAXSO,LMMAXSO)
      DOUBLE PRECISION    VNSPLL(LMMAXD,LMMAXD,IRMIN:IRWS,NSPIN),
     +                    RM(IRWS),HSOFAC(IRMINSO:IRWS)

      INTEGER             IRCUT(0:IPAND)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX      EFAC1,EFAC2
      DOUBLE PRECISION    ERR,ERRSRA
      INTEGER             I,IR,IRC1,J,LM1,LM2,LM3,LM1ERR,LM2ERR,
     +                    LM1MOD,ISP
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX,ALLOCATABLE  ::  ADER(:,:,:),
     +                                AMAT(:,:,:),
     +                                BDER(:,:,:),
     +                                BMAT(:,:,:),
     +                                PNS0(:,:,:,:),
     +                                PNS1(:,:,:)
      DOUBLE PRECISION,ALLOCATABLE :: VNSPLL0(:,:,:)
      INTEGER IPIV(LMMAXSO)
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
c      DATA FRED/.true./

      IRC1 = IRCUT(IPAN)

      ALLOCATE(ADER(LMMAXSO,LMMAXSO,IRMINSO:IRWS))
      ALLOCATE(AMAT(LMMAXSO,LMMAXSO,IRMINSO:IRWS))   
      ALLOCATE(BDER(LMMAXSO,LMMAXSO,IRMINSO:IRWS))     
      ALLOCATE(BMAT(LMMAXSO,LMMAXSO,IRMINSO:IRWS))  
      ALLOCATE(PNS0(LMMAXSO,LMMAXSO,IRMINSO:IRWS,2))
      ALLOCATE(PNS1(LMMAXSO,LMMAXSO,IRMINSO:IRWS))
      ALLOCATE(VNSPLL0(LMMAXSO,LMMAXSO,IRMINSO:IRWS))
      
      VNSPLL0(:,:,:)=0d0

      IF (INS.GT.0 ) THEN
        VNSPLL0(1:LMMAXD,1:LMMAXD,IRMIN:IRWS)=
     +                                      VNSPLL(:,:,IRMIN:IRWS,1)
        IF (LMMAXD .NE. LMMAXSO) THEN
          VNSPLL0(LMMAXD+1:2*LMMAXD,LMMAXD+1:2*LMMAXD,IRMIN:IRWS)=
     +                                      VNSPLL(:,:,IRMIN:IRWS,NSPIN)
        END IF

      END IF

      PNS0=0d0
      DO 1 J = 1,NSRA
        DO 2 IR = IRMINSO,IRC1
          DO 3 LM2 = 1,LMMAXD
            DO 4 LM1 = 1,LMMAXD
              PNS0(LM1,LM2,IR,J)=0d0
              IF (LMMAXD .NE. LMMAXSO) THEN
                PNS0(LMMAXD+LM1,LMMAXD+LM2,IR,J)=0d0
              END IF
              IF(LM1.EQ.LM2)THEN
                PNS0(LM1,LM2,IR,J) = PZLM(LM1,IR,J,1)
                IF (LMMAXD .NE. LMMAXSO) THEN
                  PNS0(LMMAXD+LM1,LMMAXD+LM2,IR,J) =PZLM(LM1,IR,J,NSPIN)
                END IF
              ELSE
                PNS0(LM1,LM2,IR,J) = (0D0,0D0)
                IF (LMMAXD .NE. LMMAXSO) THEN
                  PNS0(LMMAXD+LM1,LMMAXD+LM2,IR,J) = PNS0(LM1,LM2,IR,J)
                ENDIF
              ENDIF
 4          CONTINUE
 3        CONTINUE
 2      CONTINUE
 1    CONTINUE

      IF(FRED)THEN

c     DO 70 I = 0,ICST

c---> set up integrands for i-th born approximation
c       IF (I.EQ.0) THEN


c         CALL WFINT0_SO(ADER,BDER,PZLM,QZEKDR,PZEKDR,VNSPLL0,
c    +                                   LSM,HSOFAC,IRMINSO,NSRA,LMMAXD)


c       ELSE
c         CALL WFINT_SO(PNS,ADER,BDER,QZEKDR,PZEKDR,VNSPLL0,
c    +                              LSM,HSOFAC,NSRA,IRMINSO,IRWS,LMMAXD)
c       END IF

c---> call integration subroutines
c         CALL CSINWD_SO(ADER,AMAT,2*LMMAXD,IRMINSO,IRWS,IPAN,IRCUT)
c         OPEN (unit=31, file="AMAT_SO", form="formatted")
c         DO LM2=1,2*LMMAXD
c           DO LM1=1,2*LMMAXD
c             DO IR=IRMINSO,IRWS
c               WRITE(31,"((3I5),(8e17.9))") LM2,LM1,IR,AMAT(LM1,LM2,IR)
c             END DO
c           END DO
c         END DO
c         CLOSE(31)
c       CALL CSOUT(BDER,BMAT,(2*LMMAXD)**2,IRMINSO,IRWS,IPAN,IRCUT)
c         OPEN (unit=31, file="BMAT_SO", form="formatted")
c         DO LM2=1,2*LMMAXD
c           DO LM1=1,2*LMMAXD
c             DO IR=IRMINSO,IRWS
c               WRITE(31,"((3I5),(8e17.9))") LM2,LM1,IR,BMAT(LM1,LM2,IR)
c             END DO
c           END DO
c         END DO
c         CLOSE(31)

c       DO 20 IR = IRMINSO,IRC1
c---> here you have to decide, in which direction the incoming wave is (spin-)polarized
c     up or down... 
c         DO 10 LM2 = 1,2*LMMAXD
c           AMAT(LM2,LM2,IR) = CONE + AMAT(LM2,LM2,IR)
c  10     CONTINUE
c  20   CONTINUE

c---> calculate non sph. wft. in i-th born approximation
c       DO 60 J = 1,NSRA
c         DO 50 IR = IRMINSO,IRC1
c           DO 40 LM1 = 1,2*LMMAXD
c             DO 30 LM2 = 1,2*LMMAXD
c               PNS(LM1,LM2,IR,J) = 0.99d0*PNS(LM1,LM2,IR,J)+
c    +                    0.01d0*(AMAT(LM1,LM2,IR)*PZLM(LM1,IR,J)+
c    +                              BMAT(LM1,LM2,IR)*QZLM(LM1,IR,J))
c  30         CONTINUE
c  40       CONTINUE
c  50     CONTINUE
c  60   CONTINUE
c-----------------------------------------------------------------------
c check convergence
c     DO 260 J = 1,NSRA
c     DO 260 IR = IRMINSO,IRC1
c     DO 260 LM1 = 1,2*LMMAXD
c     DO 260 LM2 = 1,2*LMMAXD
c260  PNS0(LM1,LM2,IR,J) = PNS0(LM1,LM2,IR,J)-PNS(LM1,LM2,IR,J)
c     ERR=0d0
c     DO J=1,NSRA
c       ERR=0d0
c       CALL CSOUT(PNS0(1,1,IRMINSO,J),PNS1,(2*LMMAXD)**2,IRMINSO,IRWS,
c    +           IPAN,IRCUT)
c       DO LM1=1,2*LMMAXD
c         DO LM2=1,2*LMMAXD
c           ERR=MAX(ERR,ABS(PNS1(LM1,LM2,IRC1)))
c           IF (ERR .EQ. ABS(PNS1(LM1,LM2,IRC1))) THEN
c             LM1ERR=LM1
c             LM2ERR=LM2
c           END IF
c         END DO
c       END DO
c       WRITE(*,*) 'Born_Fred,NSRA',I,J,ERR
c     END DO
c     
c     IF(ERR .LT. 1D-8) EXIT
c     IF(I.EQ.ICST.AND.ERR.GT.1D-5) THEN
c       WRITE(*,*)'NOT CONVERGENT',ERR
c       STOP " "
c     END IF

c     DO 280 J = 1,NSRA
c     DO 280 IR = IRMINSO,IRC1
c     DO 280 LM1 = 1,2*LMMAXD
c     DO 280 LM2 = 1,2*LMMAXD
c280  PNS0(LM1,LM2,IR,J) = PNS(LM1,LM2,IR,J)
c-----------------------------------------------------------------------
c   70 CONTINUE
      ELSE
c-----------------------------------------------------------------------
c Volterra equation
      DO 200 I = 0,ICST
c        WRITE(6,*) "I",I
c---> set up integrands for i-th born approximation
        IF (I.EQ.0) THEN
c          CALL WFINT0_SO(ADER,BDER,PZLM,QZEKDR,PZEKDR,VNSPLL0,
c     +                  LSM,HSOFAC,IRMINSO,NSRA,LMMAXD,LMMAXSO,NSPIN)
          CALL WFINT0_SO(ADER,BDER,PZLM,QZEKDR,PZEKDR,VNSPLL0,
     +                         LSM,HSOFAC,IRMINSO,NSRA,LMMAXD,NSPIN)

c         DO IR = IRMINSO,IRWS
c           DO LM2 = 1,2*LMMAXD
c             DO LM1 = 1,2*LMMAXD
c               WRITE(332,"((3I5),(4e17.9))") LM1,LM2,IR,
c    +             ADER(LM1,LM2,IR),BDER(LM1,LM2,IR)
c             END DO
c           END DO
c         END DO

c          STOP "ADER"
        ELSE
c          CALL WFINT_SO(PNS,ADER,BDER,QZEKDR,PZEKDR,VNSPLL0,
c     +               LSM,HSOFAC,NSRA,IRMINSO,IRWS,LMMAXD,LMMAXSO,NSPIN)
          CALL WFINT_SO(PNS,ADER,BDER,QZEKDR,PZEKDR,VNSPLL0,
     +                    LSM,HSOFAC,NSRA,IRMINSO,IRWS,LMMAXD,NSPIN)
        END IF

c---> call integration subroutines

c        AMAT=0d0
c        BMAT=0d0
c       WRITE(6,*) "IRMINSO",IRMINSO
c        WRITE(6,*) "IRWS",IRWS
c       WRITE(6,*) "(2*LMMAXD)**2",(2*LMMAXD)**2
c       WRITE(6,*) "IPAN",IPAN
c       WRITE(6,*) "IPAND",IPAND

        CALL CSOUT_SPLINE(ADER(1:LMMAXSO,1:LMMAXSO,IRMINSO:IRWS),
     +                    AMAT(1:LMMAXSO,1:LMMAXSO,IRMINSO:IRWS),
     +                    (LMMAXSO)**2,IRMINSO,IRWS,IPAN,IRCUT(0:IPAN))

c        CALL CSOUT_SPLINE(ADER(1:2*LMMAXD,1:2*LMMAXD,IRMINSO:IRWS),
c     +                    AMAT(1:2*LMMAXD,1:2*LMMAXD,IRMINSO:IRWS),
c     +                    (2*LMMAXD),IRMINSO,IRWS,IPAN,IRCUT(0:IPAN))

c       WRITE(6,*) "IRMINSO",IRMINSO
c       WRITE(6,*) "IRWS",IRWS
c       WRITE(6,*) "(2*LMMAXD)**2",(2*LMMAXD)**2
c       WRITE(6,*) "IPAN",IPAN
c       WRITE(6,*) "IPAND",IPAND

        CALL CSOUT_SPLINE(BDER(1:LMMAXSO,1:LMMAXSO,IRMINSO:IRWS),
     +                    BMAT(1:LMMAXSO,1:LMMAXSO,IRMINSO:IRWS),
     +                    (LMMAXSO)**2,IRMINSO,IRWS,IPAN,IRCUT(0:IPAN))

c        IF (I==0) THEN
        DO 150 IR = IRMINSO,IRWS
          DO LM2 = 1,LMMAXSO
            DO LM1 = 1,LMMAXSO
              IF(LM1.EQ.LM2)THEN
                AMAT(LM1,LM2,IR) = CONE - AMAT(LM1,LM2,IR)
              ELSE
                AMAT(LM1,LM2,IR) = - AMAT(LM1,LM2,IR)
              ENDIF
c              WRITE(333,"((3I5),(4e17.9))") LM1,LM2,IR,
c     +            AMAT(LM1,LM2,IR),BMAT(LM1,LM2,IR)
            END DO
          END DO
  150   CONTINUE
        
c         WRITE(47,*) "I",I
c       IF (I==0) THEN
c         DO IR=IRMINSO,IRWS
c           DO LM2=1,2*LMMAXD
c             DO LM1=1,2*LMMAXD
c               WRITE(47,"((3I5),(4e17.9))") IR,LM1,LM2,      
c    +            AMAT(LM1,LM2,IR),BMAT(LM1,LM2,IR)
c             END DO
c           END DO
c         END DO
cc          STOP "PNS " 
c       END IF
c         STOP " AMAT"
c---> calculate non sph. wft. in i-th born approximation
c       WRITE(6,*) "IRC1", IRC1
c       WRITE(6,*) "IRCUT"
c       DO J=0,IPAN
c         WRITE(6,"(2I5)") J, IRCUT(J)
c       END DO
c       DO 190 J = 1,NSRA
c         DO 180 IR = IRMINSO,IRC1
c           DO 170 LM2 = 1,2*LMMAXD
c             DO 160 LM1 = 1,LMMAXD
c                LM1MOD=MOD(LM1-1,LMMAXD)+1
c               PNS(LM1,LM2,IR,J) = AMAT(LM1,LM2,IR)*PZLM(LM1,IR,J,1)
c    +                            + BMAT(LM1,LM2,IR)*QZLM(LM1,IR,J,1)
c 160         CONTINUE
c             DO LM1 = LMMAXD+1,2*LMMAXD
c                LM1MOD=MOD(LM1-1,LMMAXD)+1
c               LM1MOD=LM1-LMMAXD
c               PNS(LM1,LM2,IR,J) = 
c    +                     AMAT(LM1,LM2,IR)*PZLM(LM1MOD,IR,J,NSPIN)
c    +                  +  BMAT(LM1,LM2,IR)*QZLM(LM1MOD,IR,J,NSPIN)
c             END DO   
c             DO LM1 = 1,2*LMMAXD
c               WRITE(334,"((4I5),(4e17.9))") J,LM1,LM2,IR,
c    +            PNS(LM1,LM2,IR,J),PNS0(LM1,LM2,IR,J)
c             END DO   
c 170       CONTINUE
c           DO LM1 = 1,LMMAXD
c             WRITE(335,"((3I5),(4e17.9))") J,LM2,IR,
c    +            PZLM(LM1,IR,J,1),QZLM(LM1,IR,J,1)
c           END DO
c 180     CONTINUE
c 190   CONTINUE
        DO J = 1,NSRA
          DO IR = IRMINSO,IRC1
            DO LM1 = 1,LMMAXSO
              IF (LM1 .GT. LMMAXD) THEN
                LM1MOD=LM1-LMMAXD
              ELSE
                LM1MOD=LM1
              END IF
              DO LM2 = 1,LMMAXSO
                PNS(LM1,LM2,IR,J) = AMAT(LM1,LM2,IR)*PZLM(LM1MOD,IR,J,1)
     +                            + BMAT(LM1,LM2,IR)*QZLM(LM1MOD,IR,J,1)
              END DO
            END DO   
          END DO   
        END DO   
c        STOP "PNS " 
c-----------------------------------------------------------------------
c  check convergence
c      DO 290 J = 1,NSRA
c      DO 290 IR = IRMINSO,IRC1
c      DO 290 LM2 = 1,2*LMMAXD
c      DO 290 LM1 = 1,2*LMMAXD
c  290  PNS0(LM1,LM2,IR,J) = PNS0(LM1,LM2,IR,J)-PNS(LM1,LM2,IR,J)
       DO J = 1,NSRA
         DO IR = IRMINSO,IRC1
           DO LM2 = 1,LMMAXSO
             DO LM1 = 1,LMMAXSO
               PNS0(LM1,LM2,IR,J) = PNS0(LM1,LM2,IR,J)-PNS(LM1,LM2,IR,J)
c               IF (I==20) THEN
c                 WRITE(335,"((4I5),(4e17.9))") J,LM1,LM2,IR,
c     +                  PNS(LM1,LM2,IR,J)
c                 IF (J==1) THEN
c                   WRITE(336,"((3I5),(4e17.9))") IR,LM1,LM2,
c     +                  AMAT(LM1,LM2,IR),BMAT(LM1,LM2,IR)
c                 END IF
c               END IF
             END DO
           END DO
         END DO
       END DO
c       IF (I==20)  STOP "PNS " 

       ERR=0D0

       DO J=1,NSRA
         CALL CSOUT_SPLINE(PNS0(1,1,IRMINSO,J),PNS1,(LMMAXSO)**2,
     +           IRMINSO,IRWS,IPAN,IRCUT)
         DO LM2=1,LMMAXSO
           DO LM1=1,LMMAXSO
             ERR=MAX(ERR,ABS(PNS1(LM1,LM2,IRC1)))
           END DO
         END DO
c         WRITE(*,*) 'Born_Volterra,NSRA',I,J,ERRSRA
       END DO
       WRITE(*,*) 'ERR Volterra',I,ERR
       IF(ERR.GT.1D+15) STOP 'NOT CONVERGENT'

       IF(ERR .LT. 1D-8) EXIT

c       IF(I.EQ.ICST.AND.ERR.GT.1D-08)WRITE(*,*)'NOT CONVERGENT',ERR
c       IF(I.EQ.ICST.AND.ERR.GT.1D-05) STOP "NOT CONVERGED"
       DO 310 J = 1,NSRA
       DO 310 IR = IRMINSO,IRC1
       DO 310 LM2 = 1,LMMAXSO
       DO 310 LM1 = 1,LMMAXSO
  310  PNS0(LM1,LM2,IR,J) = PNS(LM1,LM2,IR,J)
c-----------------------------------------------------------------------
  200 CONTINUE
c      DO J = 1,NSRA
c        DO IR = IRMINSO,IRC1
c          DO LM2 = 1,2*LMMAXD
c            DO LM1 = 1,2*LMMAXD
c                WRITE(335,"((4I5),(4e17.9))") J,LM1,LM2,IR,
c    +                  PNS(LM1,LM2,IR,J)
c                IF (J==1) THEN
c                  WRITE(336,"((3I5),(4e17.9))") IR,LM1,LM2,
c    +                  AMAT(LM1,LM2,IR),BMAT(LM1,LM2,IR)
c                END IF
c            END DO
c          END DO
c        END DO
c      END DO
c       STOP "PNS " 
c      WRITE(6,*) "IRC1",IRC1
c     WRITE(6,*) "LMMAXSO",LMMAXSO
c     WRITE(6,*) "LMMAXD",LMMAXD
c     DO LM2=1,2*LMMAXD
c       DO LM1=1,2*LMMAXD
c         WRITE(47,"((2I5),(4e17.9))") LM2,LM1,AMAT(LM1,LM2,IRC1),
c    +                                         BMAT(LM1,LM2,IRC1)  
c       END DO
c     END DO
      AR=0d0
      BR=0d0
      CALL ZGEINV1(AMAT(1,1,IRC1),AR,BR,IPIV,LMMAXSO)
c      WRITE(6,*) "after ZGEINV1"
      DO 220 LM2=1,LMMAXSO
      DO 220 LM1=1,LMMAXSO
      DO 210 IR=IRMINSO,IRC1
        ADER(LM1,LM2,IR)=(0d0,0d0)
  210   BDER(LM1,LM2,IR)=(0d0,0d0)
      DO 220 LM3=1,LMMAXSO
      DO 220 IR=IRMINSO,IRC1
        ADER(LM1,LM2,IR)=ADER(LM1,LM2,IR)+AMAT(LM1,LM3,IR)*AR(LM3,LM2)
  220   BDER(LM1,LM2,IR)=BDER(LM1,LM2,IR)+BMAT(LM1,LM3,IR)*AR(LM3,LM2)
      DO 230 LM2=1,2*LMMAXD
      DO 230 LM1=1,2*LMMAXD
      DO 230 IR=IRMINSO,IRC1
        AMAT(LM1,LM2,IR)=ADER(LM1,LM2,IR)
  230   BMAT(LM1,LM2,IR)=BDER(LM1,LM2,IR)
      DO J = 1,NSRA
        DO IR = IRMINSO,IRC1
          DO LM2 = 1,LMMAXSO
            DO LM1 = 1,LMMAXD
              PNS(LM1,LM2,IR,J) =  AMAT(LM1,LM2,IR)*PZLM(LM1,IR,J,1)
     +                            +BMAT(LM1,LM2,IR)*QZLM(LM1,IR,J,1)
            END DO  
            IF (LMMAXSO .NE. LMMAXD) THEN
              DO LM1 = LMMAXD+1,2*LMMAXD
                PNS(LM1,LM2,IR,J) = 
     +                      AMAT(LM1,LM2,IR)*PZLM(LM1-LMMAXD,IR,J,NSPIN)
     +                     +BMAT(LM1,LM2,IR)*QZLM(LM1-LMMAXD,IR,J,NSPIN)
              END DO
            END IF
          END DO
        END DO
      END DO  

      WRITE(6,*) "after Volterra"
c Volterra equation
c-----------------------------------------------------------------------
      ENDIF
c      WRITE(6,*) "IRC1", IRC1
      DO 90 LM2 = 1,LMMAXSO
        EFAC2 = EFAC(MOD(LM2-1,LMMAXD)+1)
c---> store alpha and t - matrix
        DO 80 LM1 = 1,LMMAXSO
          EFAC1 = EFAC(MOD(LM1-1,LMMAXD)+1)
          AR(LM1,LM2) = AMAT(LM1,LM2,IRMINSO)
c---> t-matrix
          BR(LM1,LM2) = BMAT(LM1,LM2,IRC1)*EFAC1*EFAC2/EK
   80   CONTINUE
   90 CONTINUE
c      write(86,*) ""
c---> rescale with efac
      DO 130 J = 1,NSRA
        DO 120 LM2 = 1,LMMAXSO
          EFAC2 = EFAC(MOD(LM2-1,LMMAXD)+1)
c         WRITE(6,*) "LM2,EFAC(LM)", LM2,EFAC(LM2)
          DO 110 IR = IRMINSO,IRC1
            DO 100 LM1 = 1,LMMAXSO
              PNS(LM1,LM2,IR,J) = PNS(LM1,LM2,IR,J)*EFAC2
  100       CONTINUE
  110     CONTINUE
  120   CONTINUE
  130 CONTINUE
c      WRITE(6,*) "after PNS"
      DEALLOCATE(ADER)
      DEALLOCATE(AMAT)   
      DEALLOCATE(BDER)     
      DEALLOCATE(BMAT)  
      DEALLOCATE(PNS0)
      DEALLOCATE(PNS1)
      DEALLOCATE(VNSPLL0)
      
      END
c ************************************************************************
      SUBROUTINE ZGEINV1(A,U,AUX,IPIV,DIM)
c ************************************************************************
C   - inverts a general double complex matrix A,
c   - the result is return in U,
c   - input matrix A is returned unchanged,
c   - AUX is a auxiliary matrix,
c   - A,U and AUX are of dimension (DIM,DIM),
c ------------------------------------------------------------------------
      INTEGER DIM,IPIV(*)
      DOUBLE COMPLEX A(DIM,*),AUX(DIM,*),U(DIM,*)
c
C     .. PARAMETER
C
      DOUBLE COMPLEX CONE
      PARAMETER(CONE=(1.D0,0.D0))
C
      INTEGER LM1,INFO,I1,I2
      EXTERNAL ZCOPY,ZGETRS,ZGETRF
c ------------------------------------------------------------------------
      CALL CINIT(DIM*DIM,U)
      CALL CINIT(DIM*DIM,AUX)

      DO 10 LM1=1,DIM
        U(LM1,LM1) = CONE
 10   END DO

c     write(6,*) "DIM",DIM
c     write(6,*) "DIM/2",DIM/2
c      do i1=1,DIM
c        do i2=1,DIM
c     do i1=1,DIM/2
c       do i2=1,DIM/2
c          write(24,"((2I5),(4e17.9))") i1,i2,A(i1,i2),
c    +                      A(i1+DIM/2,i2+DIM/2)
c       end do
c     end do

c     do i1=1,DIM/2
c       do i2=1,DIM/2
c          write(124,"((2I5),(4e17.9))") i1,i2,A(i1,i2+DIM/2),
c    +                      A(i1+DIM/2,i2)
c       end do
c     end do

      CALL ZCOPY(DIM*DIM,A,1,AUX,1)

c     do i1=1,DIM
c       do i2=1,DIM
c          write(224,"((2I5),(4e17.9))") i1,i2,AUX(i1,i2),
c    +                                    A(i1,i2)
c       end do
c     end do

      CALL ZGETRF(DIM,DIM,AUX,DIM,IPIV,INFO)
c      CALL ZGETRF(DIM/2,DIM/2,AUX(1:DIM/2,1:DIM/2),DIM/2,
c     +                           IPIV(1:DIM/2),INFO)

c      WRITE(6,*) "INFO",INFO
c      STOP "after ZGETRF"
c      CALL ZGETRF(DIM,DIM,AUX,DIM,IPIV,INFO)
c      write(324,*) "INFO",INFO
c     do i1=1,DIM
c       do i2=1,DIM
c          write(324,"((2I5),(4e17.9))") i1,i2,AUX(i1,i2)
c       end do
c     end do
      CALL ZGETRS('N',DIM,DIM,AUX(1:DIM,1:DIM),DIM,IPIV(1:DIM),
     +                      U(1:DIM,1:DIM),DIM,INFO)
c      WRITE(6,*) "INFO",INFO
      
      RETURN
      END
