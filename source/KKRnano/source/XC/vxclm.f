      SUBROUTINE VXCLM(EXC,KTE,KXC,LPOT,NSPIN,RHO2NS,V,R,DRDI,
     +                 IRCUT,IPAN,GSH,ILM,IMAXSH,
     +                 IFUNM,THETAS,YR,WTYR,IJEND,LMSP,
C                      new input parameters after inc.p removal
     &                 irmd, irid, nfund, ngshd, ipand)
c-----------------------------------------------------------------------
c     add the exchange-correlation-potential to the given potential
c     and if total energies should be calculated (kte=1) the exchange-
c     correlation-energies are calculated .
c     use as input the charge density times r**2 (rho2ns(...,1)) and
c     in the spin-polarized case (nspin=2) the spin density times r**2
c     (rho2ns(...,2)) .
c     the density times 4 pi is generated at an angular mesh .
c     the exchange-correlation potential and the exchange-correlation
c     energy are calculated at those mesh points with a subroutine .
c     in the paramagnetic case the "spin-density" is set equal zero .
c     after that the exchange-correlation potential and in the case of
c     total energies (kte=1) the exchange-correlation energy are
c     expanded into spherical harmonics .
c     the ex.-cor. potential is added to the given potential .
c     the expansion into spherical harmonics uses the orthogonality
c     of these harmonics . - therefore a gauss-legendre integration
c     for "theta" and a gauss-tschebyscheff integration for "phi"
c     is used .
c     all needed values for the angular mesh and angular integration
c     are generate in the subroutine sphere .
c
c     the ex.-cor. potential is extrapolated to the origin only
c     for the lm=1 value .
c
c                               b.drittler   june 1987
c
c     modified for shape functions
c                                       b. drittler oct. 1989
c     simplified and modified for Paragon X/PS
c                                       R. Zeller Nov. 1993
c                            cor error 23/6/1996
c-----------------------------------------------------------------------
      use Quadrature_mod, only: simpson
      IMPLICIT NONE

      INTEGER irmd
      INTEGER irid
      INTEGER nfund
      INTEGER ngshd
      INTEGER ipand

C     .. Parameters ..
C     INTEGER LMPOTD
C     PARAMETER (LMPOTD= (LPOTD+1)**2)
C     INTEGER LMXSPD
C     PARAMETER (LMXSPD= (2*LPOTD+1)**2)

C     .. Scalar Arguments ..
      INTEGER IJEND,IPAN,KTE,KXC,LPOT,NSPIN
C     ..
C     .. Array Arguments ..
C     DOUBLE PRECISION DRDI(IRMD),EXC(0:LPOTD),GSH(*),R(IRMD),
C    +                 RHO2NS(IRMD,LMPOTD,2),THETAS(IRID,NFUND),
C    +                 V(IRMD,LMPOTD,2),WTYR(IJEND,*),YR(IJEND,*)
C     INTEGER IFUNM(LMXSPD)
C     INTEGER ILM(NGSHD,3),IMAXSH(0:LMPOTD),IRCUT(0:IPAND),
C    +        LMSP(LMXSPD)

      DOUBLE PRECISION DRDI(IRMD)
      DOUBLE PRECISION EXC(0:LPOT)
      DOUBLE PRECISION GSH(*)
      DOUBLE PRECISION R(IRMD)
      DOUBLE PRECISION RHO2NS(IRMD,(LPOT+1)**2,2)
      DOUBLE PRECISION THETAS(IRID,NFUND)
      DOUBLE PRECISION V(IRMD,(LPOT+1)**2,2)
      DOUBLE PRECISION WTYR(IJEND,*)
      DOUBLE PRECISION YR(IJEND,*)
      INTEGER IFUNM((2*LPOT+1)**2)
      INTEGER ILM(NGSHD,3)
      INTEGER IMAXSH(0:(LPOT+1)**2)
      INTEGER IRCUT(0:IPAND)
      INTEGER LMSP((2*LPOT+1)**2)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ELMXC,FPI,FPIPR2,VLMXC,VXC1,VXC2,VXC3
      INTEGER IFUN,IJ,IPOT,IR,IRC1,IRH,IRS1,IS,ISPIN,J,L,LM,LM2,LMMAX,M
C     ..
C     .. Local Arrays ..
C     DOUBLE PRECISION ER(IRMD,0:LPOTD),ESTOR(IRMD,LMPOTD),EXCIJ(IJEND),
C    +                 FPRHO(IJEND,2),VXC(IJEND,2),VXCR(2:3,2)
      DOUBLE PRECISION ER(IRMD,0:LPOT)
      DOUBLE PRECISION ESTOR(IRMD,(LPOT+1)**2)
      DOUBLE PRECISION EXCIJ(IJEND)
      DOUBLE PRECISION FPRHO(IJEND,2)
      DOUBLE PRECISION VXC(IJEND,2)
      DOUBLE PRECISION VXCR(2:3,2)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,SIMPK,VOSKO,VXCSPO
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN
C     ..
      FPI = 16.0D0*ATAN(1.0D0)
      LMMAX = (LPOT+1)* (LPOT+1)
c
c---> loop over given representive atoms
c
      IRC1 = IRCUT(IPAN)
      IRS1 = IRCUT(1)

      DO 10 ISPIN = 1,NSPIN
        VXCR(2,ISPIN) = 0.0D0
        VXCR(3,ISPIN) = 0.0D0
   10 CONTINUE
c
c---> initialize for ex.-cor. energy
c
      IF (KTE.EQ.1) THEN
        DO 30 L = 0,LPOT
          EXC(L) = 0.0D0
          DO 20 IR = 1,IRC1
            ER(IR,L) = 0.0D0
   20     CONTINUE
   30   CONTINUE
c
        DO 50 LM = 1,LMMAX
          DO 40 IR = 1,IRC1
            ESTOR(IR,LM) = 0.0D0
   40     CONTINUE
   50   CONTINUE
      END IF
c
c---> loop over radial mesh
c

      DO 140 IR = 2,IRC1

c
c---> generate the densities on an angular mesh
c
        DO 70 IS = 1,2
          DO 60 IJ = 1,IJEND
            FPRHO(IJ,IS) = 0.D0
   60     CONTINUE
   70   CONTINUE

        FPIPR2 = FPI/R(IR)**2
        DO 90 ISPIN = 1,NSPIN
          DO 80 LM = 1,LMMAX
            CALL DAXPY(IJEND,RHO2NS(IR,LM,ISPIN)*FPIPR2,YR(1,LM),1,
     +                 FPRHO(1,ISPIN),1)
   80     CONTINUE
   90   CONTINUE
c
c---> calculate the ex.-cor. potential
c
        IF (KXC.LE.1) THEN
          CALL VXCSPO(EXCIJ,FPRHO,VXC,KXC,IJEND,IJEND)
        ELSE
          CALL VOSKO(EXCIJ,FPRHO,VXC,IJEND,IJEND)
        END IF
c
c---> expand the ex.-cor. potential into spherical harmonics ,
c       using the orthogonality
c
        DO 110 ISPIN = 1,NSPIN
c
c---> determine the corresponding potential number
c
          IPOT = ISPIN
          DO 100 LM = 1,LMMAX
            VLMXC = DDOT(IJEND,VXC(1,ISPIN),1,WTYR(1,LM),1)
            V(IR,LM,IPOT) = V(IR,LM,IPOT) + VLMXC
c
c---> store the ex.-c. potential of ir=2 and =3 for the extrapolation
c
            IF (LM.EQ.1 .AND. (IR.EQ.2.OR.IR.EQ.3)) VXCR(IR,
     +          ISPIN) = VLMXC
  100     CONTINUE
  110   CONTINUE
c
c---> file er in case of total energies
c
        IF (KTE.EQ.1) THEN
c
c---> expand ex.-cor. energy into spherical harmonics
c       using the orthogonality
c
          DO 130 L = 0,LPOT
            DO 120 M = -L,L
              LM = L*L + L + M + 1
              ELMXC = DDOT(IJEND,EXCIJ,1,WTYR(1,LM),1)
c
c---> multiply the lm-component of the ex.-cor. energy with the same
c     lm-component of the charge density times r**2 and sum over lm
c     this corresponds to a integration over the angular .
c
              IF (IR.GT.IRS1) THEN
                ESTOR(IR,LM) = ELMXC

              ELSE

                ER(IR,L) = ER(IR,L) + RHO2NS(IR,LM,1)*ELMXC
              END IF

  120       CONTINUE

  130     CONTINUE

        END IF

  140 CONTINUE

c
c---> integrate er in case of total energies to get exc
c
      IF (KTE.EQ.1) THEN

        DO 190 L = 0,LPOT
          DO 180 M = -L,L
            LM = L*L + L + M + 1
c
c---> convolute with shape function
c
            DO 170 J = IMAXSH(LM-1) + 1,IMAXSH(LM)
              LM2 = ILM(J,2)
              IF (LMSP(ILM(J,3)).GT.0) THEN
                IFUN = IFUNM(ILM(J,3))
                DO 160 IR = IRS1 + 1,IRC1
                  IRH = IR - IRS1
                  ER(IR,L) = ER(IR,L) + RHO2NS(IR,LM,1)*GSH(J)*
     +                       THETAS(IRH,IFUN)*ESTOR(IR,LM2)
  160           CONTINUE
              END IF
  170       CONTINUE
  180     CONTINUE
c         CALL SIMPK(ER(1,L),EXC(L),IPAN,IRCUT,DRDI)
          EXC(L) = simpson(ER(1:,L), IPAN, IRCUT, DRDI)
  190     CONTINUE

      END IF
c
c---> extrapolate ex.-cor potential to the origin only for lm=1
c
      DO 200 ISPIN = 1,NSPIN
        IPOT = ISPIN
c
        VXC2 = VXCR(2,ISPIN)
        VXC3 = VXCR(3,ISPIN)
        VXC1 = VXC2 - R(2)* (VXC3-VXC2)/ (R(3)-R(2))
c
        V(1,1,IPOT) = V(1,1,IPOT) + VXC1
  200 CONTINUE
c
      END
