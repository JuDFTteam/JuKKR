      SUBROUTINE VXCGGA(EXC,KTE,KXC,LPOT,NSPIN,IATYP,RHO2NS,V,R,DRDI,A,
     +                  IRWS,IRCUT,IPAN,GSH,ILM,IMAXSH,
     +                  IFUNM,THETAS,YR,WTYR,IJEND,LMSP,THET,YLM,DYLMT1,
     +                  DYLMT2,DYLMF1,DYLMF2,DYLMTF,
C                       new input parameters after inc.p removal
     &                  irmd, irid, nfund, ngshd, ipand)
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
c     in the non-spin-polarized case the "spin-density" is
c     set equal zero .
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
c-----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER irmd
      INTEGER irid
      INTEGER nfund
      INTEGER ngshd
      INTEGER ipand

C     INTEGER LMPOTD
C     PARAMETER (LMPOTD= (LPOTD+1)**2)
C     INTEGER LMXSPD
C     PARAMETER (LMXSPD= (2*LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION A
      INTEGER IATYP,IJEND,IPAN,IRWS,KTE,KXC,LPOT,NSPIN
C     ..
C     .. Array Arguments ..
C     DOUBLE PRECISION DRDI(IRMD),R(IRMD),
C    +                 DYLMF1(IJEND,LMPOTD),DYLMF2(IJEND,LMPOTD),
C    +                 DYLMT1(IJEND,LMPOTD),DYLMT2(IJEND,LMPOTD),
C    +                 DYLMTF(IJEND,LMPOTD),EXC(0:LPOTD),GSH(*),
C    +                 RHO2NS(IRMD,LMPOTD,2),THET(IJEND),WTYR(IJEND,*),
C    +                 THETAS(IRID,NFUND),V(IRMD,LMPOTD,2),
C    +                 YLM(IJEND,LMPOTD),YR(IJEND,*)
C     INTEGER IFUNM(LMXSPD)
C     INTEGER ILM(NGSHD,3),IMAXSH(0:LMPOTD),IRCUT(0:IPAND),
C    +        LMSP(LMXSPD)

      DOUBLE PRECISION DRDI(IRMD)
      DOUBLE PRECISION R(IRMD)
      DOUBLE PRECISION DYLMF1(IJEND,(LPOT+1)**2)
      DOUBLE PRECISION DYLMF2(IJEND,(LPOT+1)**2)
      DOUBLE PRECISION DYLMT1(IJEND,(LPOT+1)**2)
      DOUBLE PRECISION DYLMT2(IJEND,(LPOT+1)**2)
      DOUBLE PRECISION DYLMTF(IJEND,(LPOT+1)**2)
      DOUBLE PRECISION EXC(0:(LPOT+1)**2)
      DOUBLE PRECISION GSH(*)
      DOUBLE PRECISION RHO2NS(IRMD,(LPOT+1)**2,2)
      DOUBLE PRECISION THET(IJEND)
      DOUBLE PRECISION WTYR(IJEND,*)
      DOUBLE PRECISION THETAS(IRID,NFUND)
      DOUBLE PRECISION V(IRMD,(LPOT+1)**2,2)
      DOUBLE PRECISION YLM(IJEND,(LPOT+1)**2)
      DOUBLE PRECISION YR(IJEND,*)

      INTEGER IFUNM((2*LPOT+1)**2)
      INTEGER ILM(NGSHD,3)
      INTEGER IMAXSH(0:(LPOT+1)**2)
      INTEGER IRCUT(0:IPAND)
      INTEGER LMSP((2*LPOT+1)**2)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CHGDEN,DX,ELMXC,FPI,R1,R2,RPOINT,SPIDEN,VLMXC,
     +                 VXC1,VXC2,VXC3,ZERO,ZERO1
      INTEGER IFUN,IPAN1,IPOT,IR,IRC0,IRC1,IRH,IRS1,ISPIN,J,L,L1MAX,LM,
     +        LM2,LMMAX,M,MESH,NSPIN2
C     ..
C     .. Local Arrays ..
C     DOUBLE PRECISION DDRRL(IRMD,LMPOTD),DDRRUL(IRMD,LMPOTD),
C    +                 DRRL(IRMD,LMPOTD),DRRUL(IRMD,LMPOTD),
C    +                 ER(IRMD,0:LPOTD),ESTOR(IRMD,LMPOTD),EXCIJ(IJEND),
C    +                 RHOL(IRMD,2,LMPOTD),RHOLM(LMPOTD,2),VXC(IJEND,2),
C    +                 VXCR(2:3,2)

      DOUBLE PRECISION DDRRL(IRMD,(LPOT+1)**2)
      DOUBLE PRECISION DDRRUL(IRMD,(LPOT+1)**2)
      DOUBLE PRECISION DRRL(IRMD,(LPOT+1)**2)
      DOUBLE PRECISION DRRUL(IRMD,(LPOT+1)**2)
      DOUBLE PRECISION ER(IRMD,0:(LPOT+1)**2)
      DOUBLE PRECISION ESTOR(IRMD,(LPOT+1)**2)
      DOUBLE PRECISION EXCIJ(IJEND)
      DOUBLE PRECISION RHOL(IRMD,2,(LPOT+1)**2)
      DOUBLE PRECISION RHOLM((LPOT+1)**2,2)
      DOUBLE PRECISION VXC(IJEND,2)
      DOUBLE PRECISION VXCR(2:3,2)
C     ..

C     .. External Functions ..
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT
C     ..
C     .. External Subroutines ..
      EXTERNAL GRADRL,MKXCPE,SIMPK
C     ..


C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,MOD
C     ..

C     .. Data statements ..
      DATA ZERO,ZERO1/0.d0,1.d-12/

      INTEGER LMPOTD
      LMPOTD= (LPOT+1)**2

C     ..
      WRITE (6,FMT=*) ' GGA CALCULATION '
      FPI = 16.0D0*ATAN(1.0D0)
      LMMAX = (LPOT+1)* (LPOT+1)
c
c---> loop over given representive atoms
c
      IPAN1 = IPAN
      IRC1 = IRCUT(IPAN)
      IRS1 = IRCUT(1)
      IRC0 = 2

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
C
      L1MAX = LPOT + 1
      MESH = IRWS
      DX = A
C
      IF (NSPIN.EQ.2) THEN
        DO 70 LM = 1,LMMAX
          DO 60 IR = 2,MESH
            R1 = R(IR)
            R2 = R1*R1
            CHGDEN = RHO2NS(IR,LM,1)/R2
            SPIDEN = RHO2NS(IR,LM,2)/R2
            IF (ABS(CHGDEN).LE.ZERO1) CHGDEN = ZERO
            IF (ABS(SPIDEN).LE.ZERO1) SPIDEN = ZERO
            RHOL(IR,2,LM) = (CHGDEN+SPIDEN)/2.d0
            RHOL(IR,1,LM) = (CHGDEN-SPIDEN)/2.d0
   60     CONTINUE
c
c       extrapolate
c
          RHOL(1,1,LM) = RHOL(2,1,LM)
          RHOL(1,2,LM) = RHOL(2,2,LM)
   70   CONTINUE
C
      ELSE
C
        DO 90 LM = 1,LMMAX
          DO 80 IR = 2,MESH
            R1 = R(IR)
            R2 = R1*R1
c
            CHGDEN = RHO2NS(IR,LM,1)/R2
            IF (ABS(CHGDEN).LE.ZERO1) CHGDEN = ZERO
            RHOL(IR,1,LM) = CHGDEN/2.d0
            RHOL(IR,2,LM) = CHGDEN/2.d0
   80     CONTINUE
c
c       extrapolate
          RHOL(1,1,LM) = RHOL(2,1,LM)
          RHOL(1,2,LM) = RHOL(2,2,LM)
   90   CONTINUE
      END IF
c

      CALL GRADRL(NSPIN,MESH,L1MAX,DX,RHOL,R,DRDI,IPAN1,IPAND,IRCUT,
     +            DRRL,DDRRL,DRRUL,DDRRUL,IRMD,LMPOTD)

c
c---> loop over radial mesh
c

      DO 160 IR = IRC0,IRC1
        RPOINT = R(IR)
c
c---> calculate the ex.-cor. potential
c
        NSPIN2 = 2

        DO 110 ISPIN = 1,NSPIN2
          DO 100 LM = 1,LMMAX
            RHOLM(LM,ISPIN) = RHOL(IR,ISPIN,LM)
  100     CONTINUE
  110   CONTINUE
c
c    only for spin-polarized
c
c        write(6,*) ' before mkxcpe '
        CALL MKXCPE(NSPIN2,IR,IJEND,L1MAX,RPOINT,RHOLM,VXC,EXCIJ,THET,
     +              YLM,DYLMT1,DYLMT2,DYLMF1,DYLMF2,DYLMTF,DRRL,DDRRL,
     +              DRRUL,DDRRUL,IRMD,LMPOTD)
c
c
c
c
c---> expand the ex.-cor. potential into spherical harmonics ,
c       using the orthogonality
c
        DO 130 ISPIN = 1,NSPIN
c
c---> determine the corresponding potential number
c
          IPOT = ISPIN
          DO 120 LM = 1,LMMAX
            VLMXC = DDOT(IJEND,VXC(1,ISPIN),1,WTYR(1,LM),1)
            V(IR,LM,IPOT) = V(IR,LM,IPOT) + VLMXC
c
c---> store the ex.-c. potential of ir=2 and =3 for the extrapolation
c
            IF (LM.EQ.1 .AND. (IR.EQ.2.OR.IR.EQ.3)) VXCR(IR,
     +          ISPIN) = VLMXC
  120     CONTINUE
  130   CONTINUE
c
c---> file er in case of total energies
c
        IF (KTE.EQ.1) THEN
c
c---> expand ex.-cor. energy into spherical harmonics
c       using the orthogonality
c
          DO 150 L = 0,LPOT
            DO 140 M = -L,L
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

  140       CONTINUE

  150     CONTINUE

        END IF

  160 CONTINUE
c
c---> integrate er in case of total energies to get exc
c
      IF (KTE.EQ.1) THEN

        DO 210 L = 0,LPOT
          DO 200 M = -L,L
            LM = L*L + L + M + 1
c
c---> convolute with shape function
c
            DO 190 J = IMAXSH(LM-1) + 1,IMAXSH(LM)
              LM2 = ILM(J,2)
              IF (LMSP(ILM(J,3)).GT.0) THEN
                IFUN = IFUNM(ILM(J,3))
                DO 180 IR = IRS1 + 1,IRC1
                  IRH = IR - IRS1
                  ER(IR,L) = ER(IR,L) + RHO2NS(IR,LM,1)*GSH(J)*
     +                         THETAS(IRH,IFUN)*ESTOR(IR,LM2)
  180           CONTINUE
              END IF
  190       CONTINUE
  200     CONTINUE
          CALL SIMPK(ER(1,L),EXC(L),IPAN1,IRCUT,DRDI)
  210     CONTINUE

      END IF
c
c---> extrapolate ex.-cor potential to the origin only for lm=1
c
      DO 220 ISPIN = 1,NSPIN
        IPOT = ISPIN
c
        VXC2 = VXCR(2,ISPIN)
        VXC3 = VXCR(3,ISPIN)
        VXC1 = VXC2 - R(2)* (VXC3-VXC2)/ (R(3)-R(2))
c
        V(1,1,IPOT) = V(1,1,IPOT) + VXC1
  220 CONTINUE
c
      END
