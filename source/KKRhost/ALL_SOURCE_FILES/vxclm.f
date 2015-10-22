      SUBROUTINE VXCLM(EXC,KTE,KXC,LMAX,NSPIN,IATYP,RHO2NS,V,R,DRDI,
     +                 IRWS,IRCUT,IPAN,KSHAPE,GSH,ILM,IMAXSH,
     +                 IFUNM,THETAS,YR,WTYR,IJEND,LMSP)
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
      INCLUDE 'inc.p'
C     .. Parameters ..
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER LMXSPD
      PARAMETER (LMXSPD= (2*LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER IATYP,IJEND,IPAN,IRWS,KSHAPE,KTE,KXC,LMAX,NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DRDI(IRMD),EXC(0:LPOTD,*),GSH(*),R(IRMD),
     +                 RHO2NS(IRMD,LMPOTD,2),THETAS(IRID,NFUND),
     +                 V(IRMD,LMPOTD,2),WTYR(IJEND,*),YR(IJEND,*)
      INTEGER IFUNM(LMXSPD)
      INTEGER ILM(NGSHD,3),IMAXSH(0:LMPOTD),IRCUT(0:IPAND),
     +        LMSP(LMXSPD)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ELMXC,FPI,FPIPR2,VLMXC,VXC1,VXC2,VXC3,factor
      INTEGER IFUN,IJ,IPOT,IR,IRC1,IRH,IRS1,IS,ISPIN,J,L,LM,LM2,LMMAX,M
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION ER(IRMD,0:LPOTD),ESTOR(IRMD,LMPOTD),EXCIJ(IJEND),
     +                 FPRHO(IJEND,2),VXC(IJEND,2),VXCR(2:3,2)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,SIMP3,SIMPK,VOSKO,VXCSPO
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN
C     ..
      WRITE(1337,*) 'Including cutoff of vxc for small density'
      FPI = 16.0D0*ATAN(1.0D0)
      LMMAX = (LMAX+1)* (LMAX+1)

c
c---> loop over given representive atoms
c
      IF (KSHAPE.NE.0) THEN
        IRC1 = IRCUT(IPAN)
        IRS1 = IRCUT(1)

      ELSE

        IRC1 = IRWS
        IRS1 = IRC1
      END IF

      DO 10 ISPIN = 1,NSPIN
        VXCR(2,ISPIN) = 0.0D0
        VXCR(3,ISPIN) = 0.0D0
   10 CONTINUE
c
c---> initialize for ex.-cor. energy
c
      IF (KTE.EQ.1) THEN
        DO 30 L = 0,LMAX
          EXC(L,IATYP) = 0.0D0
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

          do ij=1,ijend
          factor = (1.d0-dexp(-dabs(fprho(ij,1))*1000.d0))
          do ispin=1,nspin
          vxc(ij,ispin) =
     &      vxc(ij,ispin) * factor  !cutoff
          enddo
          enddo

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
          DO 130 L = 0,LMAX
            DO 120 M = -L,L
              LM = L*L + L + M + 1
              ELMXC = DDOT(IJEND,EXCIJ,1,WTYR(1,LM),1)
c
c---> multiply the lm-component of the ex.-cor. energy with the same
c     lm-component of the charge density times r**2 and sum over lm
c     this corresponds to a integration over the angular .
c
              IF ((KSHAPE.NE.0) .AND. (IR.GT.IRS1)) THEN
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
        IF (KSHAPE.EQ.0) THEN
          DO 150 L = 0,LMAX
            CALL SIMP3(ER(1,L),EXC(L,IATYP),1,IRS1,DRDI)
  150     CONTINUE

        ELSE

          DO 190 L = 0,LMAX
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
     +                         THETAS(IRH,IFUN)*ESTOR(IR,LM2)
  160             CONTINUE
                END IF
  170         CONTINUE
  180       CONTINUE
            CALL SIMPK(ER(1,L),EXC(L,IATYP),IPAN,IRCUT,DRDI)
  190     CONTINUE
        END IF

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
