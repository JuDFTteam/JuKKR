C **********************************************************************
      SUBROUTINE IRWNS(CR,DR,EFAC,QNS,VNSPLL,ICST,IRMIN,IRWS,IPAN,IRCUT,
     +                 NSRA,PZLM,QZLM,PZEKDR,QZEKDR,CDER,CMAT,DDER,DMAT)
C ************************************************************************
c     determines the irregular non spherical wavefunctions in the n-th.
c       born approximation ( n given by input parameter icst ) .
c
c
c     using the wave functions pz and qz ( regular and irregular
c       solution ) of the spherically averaged potential , the ir-
c       regular wavefunction qns is determined by
c
c           qns(ir,lm1,lm2) = cr(ir,lm1,lm2)*pz(ir,l1)
c
c                                   + dr(ir,lm1,lm2)*qz(ir,l1)
c
c      the matrices cr and dr are determined by integral equations
c        containing qns and only the non spherical contributions of
c        the potential , stored in vinspll . these integral equations
c        are solved iteratively with born approximation up to given n.
c
c     the original way of writing the cr and dr matrices in the equa-
c        tion above caused numerical troubles . therefore here are used
c        rescaled cr and dr matrices (compare subroutine wftsca):
c
c              ~
c              cr(ir,lm1,lm2) = sqrt(e)**(l1+l2)
c
c                             * cr(ir,lm1,lm2)/((2*l1-1)!!*(2*l2-1)!!)
c
c              ~
c              dr(ir,lm1,lm2) = sqrt(e)**(l2-l1)
c
c                             * dr(ir,lm1,lm2)*((2*l1-1)!!/(2*l2-1)!!)
c
c
c
c     numerical tests showed that it is sufficient to integrate only
c        ca. 70 points inwards to get reliable results . then the
c        rescaled cr and dr matrices are nearly r independent . the
c        first point of the inwards integration is irmin
c
c     total energies are converged for first born approximation
c
c
c     attention :  the sign of the dr matrix is changed to reduce the
c     ===========  number of floating point operations
c
c     modified for the use of shape functions
c
c                              (see notes by b.drittler)
c
c                                b.drittler   mar.  1989
c-----------------------------------------------------------------------
c     modified by R. Zeller      Aug. 1994
C ************************************************************************
C     .. Parameters ..
      include 'inc.p'
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (LMAXD+1)**2)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
      DOUBLE COMPLEX CONE
      PARAMETER (CONE= (1.D0,0.D0))
C     ..
C     .. Scalar Arguments ..
      INTEGER ICST,IPAN,IRMIN,IRWS,NSRA
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX CDER(LMMAXD,LMMAXD,IRMIND:IRMD),
     +               CMAT(LMMAXD,LMMAXD,IRMIND:IRMD),CR(LMMAXD,LMMAXD),
     +               DDER(LMMAXD,LMMAXD,IRMIND:IRMD),
     +               DMAT(LMMAXD,LMMAXD,IRMIND:IRMD),DR(LMMAXD,LMMAXD),
     +               EFAC(LMMAXD),PZEKDR(LMMAXD,IRMIND:IRMD,2),
     +               PZLM(LMMAXD,IRMIND:IRMD,2),
     +               QNS(LMMAXD,LMMAXD,IRMIND:IRMD,*),
     +               QZEKDR(LMMAXD,IRMIND:IRMD,2),
     +               QZLM(LMMAXD,IRMIND:IRMD,2)
      DOUBLE PRECISION VNSPLL(LMMAXD,LMMAXD,IRMIND:IRMD)
      INTEGER IRCUT(0:IPAND)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX EFAC2
      INTEGER I,IQNS,IR,IRC1,J,LM1,LM2
C     ..
C     .. External Subroutines ..
      EXTERNAL CSINWD,WFINT,WFINT0
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX0
C     ..
      IRC1 = IRCUT(IPAN)
      DO 70 I = 0,ICST
c---> set up integrands for i-th born approximation
        IF (I.EQ.0) THEN
          CALL WFINT0(CDER,DDER,QZLM,QZEKDR,PZEKDR,VNSPLL,IRMIN,IRC1,
     +                NSRA)

        ELSE
          CALL WFINT(QNS,CDER,DDER,QZEKDR,PZEKDR,VNSPLL,IRMIN,IRC1,NSRA)
        END IF
c---> call integration subroutines
        CALL CSINWD(CDER,CMAT,IRMIN,IPAN,IRCUT)
        CALL CSINWD(DDER,DMAT,IRMIN,IPAN,IRCUT)
        DO 20 IR = IRMIN,IRC1
          DO 10 LM2 = 1,LMMAXD
            DMAT(LM2,LM2,IR) = DMAT(LM2,LM2,IR) - CONE
   10     CONTINUE
   20   CONTINUE
c---> calculate non sph. wft. in i-th born approximation
        DO 60 J = 1,NSRA
          DO 50 IR = IRMIN,IRC1
            DO 40 LM1 = 1,LMMAXD
              DO 30 LM2 = 1,LMMAXD
                QNS(LM1,LM2,IR,J) = CMAT(LM1,LM2,IR)*PZLM(LM1,IR,J) -
     +                              DMAT(LM1,LM2,IR)*QZLM(LM1,IR,J)
   30         CONTINUE
   40       CONTINUE
   50     CONTINUE
   60   CONTINUE
   70 CONTINUE
      DO 90 LM2 = 1,LMMAXD
c---> store c - and d - matrix
        DO 80 LM1 = 1,LMMAXD
          CR(LM1,LM2) = CMAT(LM1,LM2,IRMIN)
          DR(LM1,LM2) = -DMAT(LM1,LM2,IRMIN)
   80   CONTINUE
   90 CONTINUE
c---> in case of muffin tin calculation : fill qns
      DO 130 J = 1,NSRA
        DO 120 IR = IRC1 + 1,IRWS
          DO 110 LM1 = 1,LMMAXD
            DO 100 LM2 = 1,LMMAXD
              QNS(LM1,LM2,IR,J) = -DMAT(LM1,LM2,IRC1)*QZLM(LM1,IR,J)
  100       CONTINUE
  110     CONTINUE
  120   CONTINUE
  130 CONTINUE
c---> rescale with efac
      IQNS = MAX0(IRC1,IRWS)
      DO 170 J = 1,NSRA
        DO 160 LM2 = 1,LMMAXD
          EFAC2 = 1.D0/EFAC(LM2)
          DO 150 IR = IRMIN,IQNS
            DO 140 LM1 = 1,LMMAXD
              QNS(LM1,LM2,IR,J) = QNS(LM1,LM2,IR,J)*EFAC2
  140       CONTINUE
  150     CONTINUE
  160   CONTINUE
  170 CONTINUE
      END
