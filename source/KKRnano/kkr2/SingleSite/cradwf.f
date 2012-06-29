      SUBROUTINE CRADWF(E,EK,NSRA,ALPHA,IPAN,IRCUT,CVLIGHT,RS,S,PZ,FZ,
     +                  QZ,SZ,TMAT,VM2Z,DRDI,R,Z,
     >                  LDAU,NLDAU,LLDAU,WMLDAUAV,LDAUCUT,
C                       new input parameters after inc.p removal
     &                  lmaxd, irmd, ipand)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c  subroutine for radial wave functions of spherical potentials
c
c             the generalized phase shifts are calculated by
c             a wronski relation :
c
c                 alpha(z,l) =-sqrt(z)*wronski{hl(r;z),rl(r;z)}; r->0
c
c             where hl is the free hankel function and rl the regular
c             solution . Using the analytical behaviour of rl at the
c             origin (rl = alphal * r**(l+1)  ; r->0),
c             the generalized phase shifts can be calculated
c             directly with the renormalization alphal .
c                                           b.drittler nov.1987
c-----------------------------------------------------------------------

      INTEGER lmaxd
      INTEGER irmd
      INTEGER ipand

      DOUBLE COMPLEX      CI,CZERO
      PARAMETER          (CI= (0.D0,1.D0),CZERO= (0.0D0,0.0D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX     E,EK
      DOUBLE PRECISION   CVLIGHT,Z
      INTEGER            IPAN,NSRA,NLDAU
      LOGICAL            LDAU
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX ALPHA(0:LMAXD),FZ(IRMD,0:LMAXD),PZ(IRMD,0:LMAXD),
     +               QZ(IRMD,0:LMAXD),SZ(IRMD,0:LMAXD),TMAT(0:LMAXD)
      DOUBLE PRECISION   DRDI(IRMD),R(IRMD),
     +                   RS(IRMD,0:LMAXD),S(0:LMAXD),
     +                   VM2Z(IRMD),
     +                   LDAUCUT(IRMD),
     +                   WMLDAUAV(LMAXD + 1)
      INTEGER            IRCUT(0:IPAND),
     +                   LLDAU(LMAXD + 1)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX ALPHAL,ARG,BL,EKLFAC,HL,PN,QF,SLOPE,TL,TLSQEZ,
     +               VALUE,W,X,Y
      DOUBLE PRECISION RIRC,RSIRC,S1
      INTEGER I,IR,IRC1,L,N
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX BESSJW(0:LMAXD+1),BESSYW(0:LMAXD+1),DLOGDP(0:LMAXD)
      DOUBLE COMPLEX HAMF(IRMD,0:LMAXD),HANKWS(0:LMAXD+1),MASS(IRMD)
      DOUBLE PRECISION DROR(IRMD)
C     ..
C     .. External Subroutines ..
C      EXTERNAL BESHAN,IRWSOL,REGSOL
C     ..
C     .. Save statement ..
C      SAVE
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE

      INTEGER LMAXP1

      LMAXP1=LMAXD+1

C INITIALISATIONS
      DROR = 0.0d0
      HAMF = CZERO
      MASS = CZERO
C INITIALISATIONS

      IRC1 = IRCUT(IPAN)
      DO IR = 2,IRC1
        DROR(IR) = DRDI(IR)/R(IR)
      END DO
      RIRC = R(IRC1)
      ARG = RIRC*EK
      CALL BESHAN(HANKWS,BESSJW,BESSYW,ARG,LMAXP1)
c
c    attention : contrary to abramowitz and stegun and
c                the bessel functions of third kind ( hankel functions)
c                are defined as:      hl(l) = nl(l) - i * jl(l)
c
      DO L = 0,LMAXP1
        HANKWS(L) = BESSYW(L) - CI*BESSJW(L)
      END DO
c
c---> calculate regular wavefunctions
c
      CALL REGSOL(CVLIGHT,E,NSRA,DLOGDP,FZ,HAMF,MASS,PZ,DROR,R,S,VM2Z,
     +              Z,IPAN,IRCUT,IRMD,IPAND,LMAXD,
     +              LDAU,NLDAU,LLDAU,WMLDAUAV,LDAUCUT)
c
      EKLFAC = EK
c
      DO 20 L = 0,LMAXD
        S1 = S(L)
        RSIRC = RS(IRC1,L)
        EKLFAC = EKLFAC/EK*DBLE(2*L+1)
c
c---> determine t - matrix
c
        PN = PZ(IRC1,L)*RSIRC
        N = L + 1
        QF = DBLE(L)/RIRC
        HL = HANKWS(L)
        BL = BESSJW(L)
        X = QF*HL - EK*HANKWS(N)
        Y = QF*BL - EK*BESSJW(N)
        W = DLOGDP(L)
        TLSQEZ = (BL*W-Y)/ (X-HL*W)
        TL = TLSQEZ/EK
        TMAT(L) = TL
c
c---> determine the renormalization
c
        ALPHAL = (BL+HL*TLSQEZ)*RIRC/PN
c
c---> determine the alpha matrix
c
        ALPHA(L) = ALPHAL*EKLFAC

        DO 10 I = 2,IRC1
          PZ(I,L) = PZ(I,L)*ALPHAL
          FZ(I,L) = FZ(I,L)*ALPHAL
   10   CONTINUE
c
        VALUE = HL*RIRC*RSIRC
        SLOPE = DBLE(L+1)*HL - RIRC*EK*HANKWS(L+1)
        SLOPE = (SLOPE*RSIRC+S1/RIRC*VALUE)
        QZ(IRC1,L) = VALUE
        SZ(IRC1,L) = (SLOPE*RIRC- (S1+1.0D0)*VALUE)/MASS(IRC1)*
     +               DROR(IRC1)
   20 CONTINUE
c
c---> calculate irregular wavefunctions
c
      CALL IRWSOL(EK,FZ,HAMF,MASS,PZ,QZ,SZ,DROR,S,IPAN,IRCUT,IRMD,
     +              IPAND,LMAXD)

      DO 50 L = 0,LMAXD
        IF (NSRA.EQ.2) THEN
          DO 30 I = 2,IRC1
            PZ(I,L) = PZ(I,L)*RS(I,L)
            QZ(I,L) = QZ(I,L)/RS(I,L)
            FZ(I,L) = FZ(I,L)*RS(I,L)/CVLIGHT
            SZ(I,L) = SZ(I,L)/RS(I,L)/CVLIGHT
   30     CONTINUE

        ELSE
          DO 40 I = 2,IRC1
            PZ(I,L) = PZ(I,L)*RS(I,L)
            QZ(I,L) = QZ(I,L)/RS(I,L)
            FZ(I,L) = CZERO
            SZ(I,L) = CZERO
   40     CONTINUE
        END IF
   50 CONTINUE

      END
