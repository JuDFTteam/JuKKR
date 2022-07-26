      MODULE mod_REGNS
      CONTAINS

      SUBROUTINE REGNS(AR,BR,EFAC,PNS,VNSPLL,ICST,IPAN,IRCUT,PZLM,
     +                   QZLM,PZEKDR,QZEKDR,EK,ADER,AMAT,BDER,BMAT,NSRA,
     +                   IRMIND,IRMD,IPAND,LMMAXATOM)
      USE mod_zgeinv1
      USE mod_wfint0
      USE mod_wfint
      USE mod_csout
      USE mod_csinwd

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
      INTEGER ICST,IPAN,IPAND,IRMD,IRMIND,LMMAXATOM,NSRA
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX ADER(LMMAXATOM,LMMAXATOM,IRMIND:IRMD),
     +               AMAT(LMMAXATOM,LMMAXATOM,IRMIND:IRMD),
     +               AR(LMMAXATOM,LMMAXATOM),
     +               BDER(LMMAXATOM,LMMAXATOM,IRMIND:IRMD),
     +               BMAT(LMMAXATOM,LMMAXATOM,IRMIND:IRMD),
     +               BR(LMMAXATOM,LMMAXATOM),
     +               EFAC(:),PNS(LMMAXATOM,LMMAXATOM,IRMIND:IRMD,2),
     +               PZEKDR(LMMAXATOM,IRMIND:IRMD,2),
     +               PZLM(LMMAXATOM,IRMIND:IRMD,2),
     +               QZEKDR(LMMAXATOM,IRMIND:IRMD,2),
     +               QZLM(LMMAXATOM,IRMIND:IRMD,2)
      DOUBLE PRECISION VNSPLL(LMMAXATOM,LMMAXATOM,IRMIND:IRMD)
      INTEGER IRCUT(0:IPAND)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX EFAC1,EFAC2
      DOUBLE PRECISION ERR
      INTEGER I,IR,IRC1,J,LM1,LM2,LM3
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX,ALLOCATABLE :: PNS0(:,:,:,:), PNS1(:,:,:)
!       DOUBLE COMPLEX PNS0(LMMAXATOM,LMMAXATOM,IRMIND:IRMD,2),
!      +               PNS1(LMMAXATOM,LMMAXATOM,IRMIND:IRMD)
      INTEGER IPIV(LMMAXATOM)
C     ..
C     .. External Subroutines ..
!       EXTERNAL CSINWD !,!CSOUT,WFINT,WFINT0 !,ZGEINV1
C     ..
C     .. Parameters ..
      DOUBLE COMPLEX CONE,CZERO
      PARAMETER (CONE= (1.0D0,0.0D0), CZERO= (0.D0,0.D0))
C     ..
      LOGICAL FRED
      DATA FRED/.false./

      ALLOCATE( PNS0(LMMAXATOM,LMMAXATOM,IRMIND:IRMD,2),
     +         PNS1(LMMAXATOM,LMMAXATOM,IRMIND:IRMD) )


C     ..
c      write(*,*)ek
      IRC1 = IRCUT(IPAN)
c      DO 1 J = 1,NSRA
c        DO 2 IR = IRMIND,IRC1
c          DO 3 LM1 = 1,LMMAXATOM
c            DO 4 LM2 = 1,LMMAXATOM
c              IF(LM1.EQ.LM2)THEN
c              PNS0(LM1,LM2,IR,J) =  AMAT(LM1,LM2,IR)*PZLM(LM1,IR,J)
c              ELSE
c              PNS0(LM1,LM2,IR,J) = (0D0,0D0)
c              ENDIF
c 4          CONTINUE
c 3        CONTINUE
c 2      CONTINUE
c 1    CONTINUE
      IF(FRED)THEN
      DO 70 I = 0,ICST
c---> set up integrands for i-th born approximation
        IF (I.EQ.0) THEN
          CALL WFINT0(ADER,BDER,PZLM,QZEKDR,PZEKDR,VNSPLL,NSRA,IRMIND,
     +                  IRMD,LMMAXATOM)
        ELSE
          CALL WFINT(PNS,ADER,BDER,QZEKDR,PZEKDR,VNSPLL,NSRA,IRMIND,
     +                 IRMD,LMMAXATOM)
        END IF
c---> call integration subroutines
        CALL CSINWD(ADER,AMAT,LMMAXATOM**2,IRMIND,IRMD,IPAN,IRCUT)
        CALL CSOUT(BDER,BMAT,LMMAXATOM**2,IRMIND,IRMD,IPAN,IRCUT)
        DO 20 IR = IRMIND,IRC1
          DO 10 LM2 = 1,LMMAXATOM
            AMAT(LM2,LM2,IR) = CONE + AMAT(LM2,LM2,IR)
   10     CONTINUE
   20   CONTINUE
c---> calculate non sph. wft. in i-th born approximation
        DO 60 J = 1,NSRA
          DO 50 IR = IRMIND,IRC1
            DO 40 LM1 = 1,LMMAXATOM
              DO 30 LM2 = 1,LMMAXATOM
                PNS(LM1,LM2,IR,J) = (AMAT(LM1,LM2,IR)*PZLM(LM1,IR,J)+
     +                              BMAT(LM1,LM2,IR)*QZLM(LM1,IR,J))
   30         CONTINUE
   40       CONTINUE
   50     CONTINUE
   60   CONTINUE
c-----------------------------------------------------------------------
c check convergence
      DO 260 J = 1,NSRA
      DO 260 IR = IRMIND,IRC1
      DO 260 LM1 = 1,LMMAXATOM
      DO 260 LM2 = 1,LMMAXATOM
 260  PNS0(LM1,LM2,IR,J) = PNS0(LM1,LM2,IR,J)-PNS(LM1,LM2,IR,J)
      ERR=0D0
      DO 270 J=1,NSRA
      CALL CSOUT(PNS0(1,1,IRMIND,J),PNS1,LMMAXATOM**2,IRMIND,IRMD,IPAN,
     +           IRCUT)
      DO 270 LM1=1,LMMAXATOM
      DO 270 LM2=1,LMMAXATOM
 270  ERR=MAX(ERR,ABS(PNS1(LM1,LM2,IRC1)))
      WRITE(1337,*) 'Born_Fred',I,ERR
c      IF(I.EQ.ICST.AND.ERR.GT.1D-3)WRITE(*,*)'NOT CONVERGENT',ERR
      DO 280 J = 1,NSRA
      DO 280 IR = IRMIND,IRC1
      DO 280 LM1 = 1,LMMAXATOM
      DO 280 LM2 = 1,LMMAXATOM
 280  PNS0(LM1,LM2,IR,J) = PNS(LM1,LM2,IR,J)
c-----------------------------------------------------------------------
   70 CONTINUE
      ELSE
c-----------------------------------------------------------------------
c Volterra equation
      DO 200 I = 0,ICST
c---> set up integrands for i-th born approximation
        IF (I.EQ.0) THEN
          CALL WFINT0(ADER,BDER,PZLM,QZEKDR,PZEKDR,VNSPLL,NSRA,IRMIND,
     +                  IRMD,LMMAXATOM)
        ELSE
          CALL WFINT(PNS,ADER,BDER,QZEKDR,PZEKDR,VNSPLL,NSRA,IRMIND,
     +                 IRMD,LMMAXATOM)
        END IF
c---> call integration subroutines
        CALL CSOUT(ADER,AMAT,LMMAXATOM**2,IRMIND,IRMD,IPAN,IRCUT)
        CALL CSOUT(BDER,BMAT,LMMAXATOM**2,IRMIND,IRMD,IPAN,IRCUT)
        DO 150 IR = IRMIND,IRC1
          DO 140 LM2 = 1,LMMAXATOM
          DO 140 LM1 = 1,LMMAXATOM
            IF(LM1.EQ.LM2)THEN
            AMAT(LM1,LM2,IR) = CONE - AMAT(LM1,LM2,IR)
            ELSE
            AMAT(LM1,LM2,IR) = - AMAT(LM1,LM2,IR)
            ENDIF
  140     CONTINUE
  150   CONTINUE
c---> calculate non sph. wft. in i-th born approximation
        DO 190 J = 1,NSRA
          DO 180 IR = IRMIND,IRC1
            DO 170 LM2 = 1,LMMAXATOM
              DO 160 LM1 = 1,LMMAXATOM
                PNS(LM1,LM2,IR,J) =  AMAT(LM1,LM2,IR)*PZLM(LM1,IR,J)
     +                              +BMAT(LM1,LM2,IR)*QZLM(LM1,IR,J)
  160         CONTINUE
  170       CONTINUE
  180     CONTINUE
  190   CONTINUE
c-----------------------------------------------------------------------
c check convergence
       DO 290 J = 1,NSRA
       DO 290 IR = IRMIND,IRC1
       DO 290 LM2 = 1,LMMAXATOM
       DO 290 LM1 = 1,LMMAXATOM
  290  PNS0(LM1,LM2,IR,J) = PNS0(LM1,LM2,IR,J)-PNS(LM1,LM2,IR,J)
       ERR=0D0
       DO 300 J=1,NSRA
       CALL CSOUT(PNS0(1,1,IRMIND,J),PNS1,LMMAXATOM**2,IRMIND,IRMD,IPAN,
     +           IRCUT)
       DO 300 LM2=1,LMMAXATOM
       DO 300 LM1=1,LMMAXATOM
  300  ERR=MAX(ERR,ABS(PNS1(LM1,LM2,IRC1)))
       WRITE(1337,*) 'Born',I,ERR 
c      IF(I.EQ.ICST.AND.ERR.GT.1D-3)WRITE(*,*)'NOT CONVERGENT',ERR
       DO 310 J = 1,NSRA
       DO 310 IR = IRMIND,IRC1
       DO 310 LM2 = 1,LMMAXATOM
       DO 310 LM1 = 1,LMMAXATOM
  310  PNS0(LM1,LM2,IR,J) = PNS(LM1,LM2,IR,J)
c-----------------------------------------------------------------------
  200 CONTINUE
      CALL ZGEINV1(AMAT(:,:,IRC1),AR,BR,IPIV,LMMAXATOM)
      DO 210 IR=IRMIND,IRC1
      DO 210 LM2=1,LMMAXATOM
      DO 210 LM1=1,LMMAXATOM
        ADER(LM1,LM2,IR)= CZERO
  210   BDER(LM1,LM2,IR)= CZERO

      DO 220 IR=IRMIND,IRC1
      DO 220 LM2=1,LMMAXATOM
      DO 220 LM3=1,LMMAXATOM
      DO 220 LM1=1,LMMAXATOM
        ADER(LM1,LM2,IR)=ADER(LM1,LM2,IR)+AMAT(LM1,LM3,IR)*AR(LM3,LM2)
  220   BDER(LM1,LM2,IR)=BDER(LM1,LM2,IR)+BMAT(LM1,LM3,IR)*AR(LM3,LM2)

      DO 230 IR=IRMIND,IRC1
      DO 230 LM2=1,LMMAXATOM
      DO 230 LM1=1,LMMAXATOM
        AMAT(LM1,LM2,IR)=ADER(LM1,LM2,IR)
  230   BMAT(LM1,LM2,IR)=BDER(LM1,LM2,IR)

      DO 240 J = 1,NSRA
      DO 240 IR = IRMIND,IRC1
      DO 240 LM2 = 1,LMMAXATOM
      DO 240 LM1 = 1,LMMAXATOM
        PNS(LM1,LM2,IR,J) =  AMAT(LM1,LM2,IR)*PZLM(LM1,IR,J)
     +                      +BMAT(LM1,LM2,IR)*QZLM(LM1,IR,J)
  240 CONTINUE
c Volterra equation
c-----------------------------------------------------------------------
      ENDIF
      DO 90 LM2 = 1,LMMAXATOM
        EFAC2 = EFAC(LM2)
c---> store alpha and t - matrix
        DO 80 LM1 = 1,LMMAXATOM
          EFAC1 = EFAC(LM1)
          AR(LM1,LM2) = AMAT(LM1,LM2,IRMIND)
c---> t-matrix
          BR(LM1,LM2) = BMAT(LM1,LM2,IRC1)*EFAC1*EFAC2/EK
   80   CONTINUE
   90 CONTINUE
c---> rescale with efac
      DO 130 J = 1,NSRA
         DO 120 IR = IRMIND,IRC1
            DO 110 LM2 = 1,LMMAXATOM
               EFAC2 = EFAC(LM2)
               DO 100 LM1 = 1,LMMAXATOM
                  PNS(LM1,LM2,IR,J) = PNS(LM1,LM2,IR,J)*EFAC2
 100           CONTINUE
 110        CONTINUE
 120     CONTINUE
 130  CONTINUE
c      if(dreal(ek).gt.3.96d-2.and.dreal(ek).lt.3.97d-2)then
c      if(dimag(ek).gt.0.5018d0.and.dreal(ek).lt.0.5019d0)then
c      write(*,*)ek
c      do 250 lm1=5,7,2
c      write(*,*)'l=',lm1
c      do 250 ir=irmind,irc1
c 250  write(*,*)ir,dreal(pns(lm1,lm1,ir,1)),dimag(pns(lm1,lm1,ir,1))
c      endif
c      endif
      END SUBROUTINE
c ************************************************************************

      END MODULE mod_REGNS
