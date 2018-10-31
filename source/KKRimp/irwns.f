      MODULE MOD_IRWNS
      CONTAINS
      SUBROUTINE IRWNS(CR,DR,EFAC,QNS,VNSPLL,ICST,IPAN,IRCUT,NSRA,
     +                   PZLM,QZLM,PZEKDR,QZEKDR,CDER,CMAT,DDER,DMAT,
     +                   IRMIND,IRMD,IPAND,LMMAXD)
      USE MOD_WFINT0
      USE MOD_WFINT
      USE MOD_CSINWD
      IMPLICIT NONE
c-----------------------------------------------------------------------
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
c                             * cr(ir,lm1,lm2)/((2*l1-1)!!*(2*l2-1)!!)
c
c              ~
c              dr(ir,lm1,lm2) = sqrt(e)**(l2-l1)
c                             * dr(ir,lm1,lm2)*((2*l1-1)!!/(2*l2-1)!!)
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
c-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE COMPLEX CONE
      PARAMETER (CONE= (1.D0,0.D0))
C     ..
C     .. Scalar Arguments ..
      INTEGER ICST,IPAN,IPAND,IRMD,IRMIND,LMMAXD,NSRA
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX CDER(LMMAXD,LMMAXD,IRMIND:IRMD),
     +               CMAT(LMMAXD,LMMAXD,IRMIND:IRMD),CR(LMMAXD,LMMAXD),
     +               DDER(LMMAXD,LMMAXD,IRMIND:IRMD),
     +               DMAT(LMMAXD,LMMAXD,IRMIND:IRMD),DR(LMMAXD,LMMAXD),
     +               EFAC(LMMAXD),PZEKDR(LMMAXD,IRMIND:IRMD,2),
     +               PZLM(LMMAXD,IRMIND:IRMD,2),
     +               QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),
     +               QZEKDR(LMMAXD,IRMIND:IRMD,2),
     +               QZLM(LMMAXD,IRMIND:IRMD,2)
      DOUBLE PRECISION VNSPLL(LMMAXD,LMMAXD,IRMIND:IRMD)
      INTEGER IRCUT(0:IPAND)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX EFAC2
      INTEGER I,IR,IRC1,J,LM1,LM2
C     ..
C     .. External Subroutines ..
!       EXTERNAL CSINWD
C     ..
      IRC1 = IRCUT(IPAN)
      DO 70 I = 0,ICST
c---> set up integrands for i-th born approximation
        IF (I.EQ.0) THEN
          CALL WFINT0(CDER,DDER,QZLM,QZEKDR,PZEKDR,VNSPLL,NSRA,IRMIND,
     +                  IRMD,LMMAXD)
        ELSE
          CALL WFINT(QNS,CDER,DDER,QZEKDR,PZEKDR,VNSPLL,NSRA,IRMIND,
     +                 IRMD,LMMAXD)
        END IF
c---> call integration subroutines
        CALL CSINWD(CDER,CMAT,LMMAXD**2,IRMIND,IRMD,IPAN,IRCUT)
        CALL CSINWD(DDER,DMAT,LMMAXD**2,IRMIND,IRMD,IPAN,IRCUT)
        DO 20 IR = IRMIND,IRC1
          DO 10 LM2 = 1,LMMAXD
            DMAT(LM2,LM2,IR) = DMAT(LM2,LM2,IR) - CONE
   10     CONTINUE
   20   CONTINUE
c---> calculate non sph. wft. in i-th born approximation
        DO 60 J = 1,NSRA
          DO 50 IR = IRMIND,IRC1
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
          CR(LM1,LM2) = CMAT(LM1,LM2,IRMIND)
          DR(LM1,LM2) = -DMAT(LM1,LM2,IRMIND)
   80   CONTINUE
   90 CONTINUE
c---> rescale with efac
      DO 130 J = 1,NSRA
        DO 120 LM2 = 1,LMMAXD
          EFAC2 = 1.D0/EFAC(LM2)
          DO 110 IR = IRMIND,IRC1
            DO 100 LM1 = 1,LMMAXD
              QNS(LM1,LM2,IR,J) = QNS(LM1,LM2,IR,J)*EFAC2
  100       CONTINUE
  110     CONTINUE
  120   CONTINUE
  130 CONTINUE
      END SUBROUTINE IRWNS
      END MODULE MOD_IRWNS
