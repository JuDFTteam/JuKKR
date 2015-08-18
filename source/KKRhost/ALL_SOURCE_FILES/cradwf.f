      SUBROUTINE CRADWF(ERYD,EK,NSRA,ALPHA,IPAN,IRCUT,CVLIGHT,RS,SL,
     +                  PZ,FZ,QZ,SZ,TMAT,VM2Z,DRDI,RMESH,ZAT,LIRRSOL,
     +                  IDOLDAU,LOPT,WLDAUAV,CUTOFF)
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
c
c   LDA+U added, March 2003 - Dec 2004, Munich/Juelich
c-----------------------------------------------------------------------
C     .. Parameters ..
      INCLUDE 'inc.p'
      INTEGER LMAXP1
      PARAMETER (LMAXP1=LMAXD+1)
      DOUBLE COMPLEX CI,CZERO
      PARAMETER (CI= (0.D0,1.D0),CZERO= (0.0D0,0.0D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX ERYD,EK
      DOUBLE PRECISION CVLIGHT,ZAT
      DOUBLE PRECISION WLDAUAV
      INTEGER IPAN,NSRA,IDOLDAU,LOPT
      LOGICAL LIRRSOL
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX ALPHA(0:LMAXD),FZ(IRMD,0:LMAXD),PZ(IRMD,0:LMAXD),
     +               QZ(IRMD,0:LMAXD),SZ(IRMD,0:LMAXD),TMAT(0:LMAXD)
      DOUBLE PRECISION DRDI(IRMD),RMESH(IRMD),RS(IRMD,0:LMAXD),
     +                 SL(0:LMAXD),VM2Z(IRMD),CUTOFF(IRMD)
      INTEGER IRCUT(0:IPAND)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX ALPHAL,ARG,BL,EKLFAC,HL,PN,QF,SLOPE,TLSQEZ,VALUE
      DOUBLE PRECISION RIRC,RIRC1,RSIRC,S1
      INTEGER I,IR,IRC1,L1
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX BESSJW(0:LMAXP1),BESSYW(0:LMAXP1),DLOGDP(0:LMAXD),
     +               HAMF(IRMD,0:LMAXD),HANKWS(0:LMAXP1),MASS(IRMD)
      DOUBLE PRECISION DROR(IRMD)
C     ..
C     .. External Subroutines ..
      EXTERNAL BESHAN,IRWSOL,REGSOL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE
C     ..
      IRC1 = IRCUT(IPAN)
      DO IR = 2,IRC1
         DROR(IR) = DRDI(IR)/RMESH(IR)
      END DO
      RIRC = RMESH(IRC1)
      RIRC1 = 1D0/RIRC
      ARG = RIRC*EK
      CALL BESHAN(HANKWS,BESSJW,BESSYW,ARG,LMAXP1)
c
c---> calculate regular wavefunctions
c
      CALL REGSOL(CVLIGHT,ERYD,NSRA,DLOGDP,FZ,HAMF,MASS,PZ,DROR,RMESH,
     +            SL,VM2Z,ZAT,IPAN,IRCUT,IDOLDAU,LOPT,WLDAUAV,CUTOFF,
     +            IRMD,IPAND,LMAXD)
c
      EKLFAC = EK
C ======================================================================
      DO L1 = 0,LMAXD
c
c---> determine t - matrix
c
        QF = DBLE(L1)*RIRC1
        HL = HANKWS(L1) * DLOGDP(L1)
        BL = BESSJW(L1) * DLOGDP(L1)
        HL = QF*HANKWS(L1) - EK*HANKWS(L1+1) - HL
        BL = BL - QF*BESSJW(L1) + EK*BESSJW(L1+1)
        HL = HL * EK
        TMAT(L1) = CI * BL/HL
c     
c---> determine the renormalization
c
        TLSQEZ = TMAT(L1) * EK
        S1 = SL(L1)
        RSIRC = RS(IRC1,L1)
        EKLFAC = EKLFAC/EK*DBLE(2*L1+1)
        PN = PZ(IRC1,L1)*RSIRC
        ALPHAL = (BESSJW(L1) - CI*HANKWS(L1)*TLSQEZ)*RIRC/PN
c
c---> determine the alpha matrix
c
        ALPHA(L1) = ALPHAL*EKLFAC

        DO I = 2,IRC1
           PZ(I,L1) = PZ(I,L1)*ALPHAL
           FZ(I,L1) = FZ(I,L1)*ALPHAL
        END DO
c
        VALUE = -CI*HANKWS(L1)*RIRC*RSIRC
        SLOPE = DBLE(L1+1)*HANKWS(L1) - RIRC*EK*HANKWS(L1+1)
        SLOPE = (-CI*SLOPE*RSIRC+S1/RIRC*VALUE)
        QZ(IRC1,L1) = VALUE
        SZ(IRC1,L1) = (SLOPE*RIRC - (S1+1.0D0)*VALUE)/MASS(IRC1)
     &               * DROR(IRC1)
      END DO
C ======================================================================
C
C -> calculate irregular wavefunctions
C
      IF ( LIRRSOL ) CALL IRWSOL(EK,FZ,HAMF,MASS,PZ,QZ,SZ,DROR,SL,
     &                           IPAN,IRCUT,IRMD,IPAND,LMAXD)
C ======================================================================
      DO L1 = 0,LMAXD
         IF (NSRA.EQ.2) THEN
            DO I = 2,IRC1
               PZ(I,L1) = PZ(I,L1)*RS(I,L1)
               QZ(I,L1) = QZ(I,L1)/RS(I,L1)
               FZ(I,L1) = FZ(I,L1)*RS(I,L1)/CVLIGHT
               SZ(I,L1) = SZ(I,L1)/RS(I,L1)/CVLIGHT
            END DO
         ELSE
            DO I = 2,IRC1
               PZ(I,L1) = PZ(I,L1)*RS(I,L1)
               QZ(I,L1) = QZ(I,L1)/RS(I,L1)
               FZ(I,L1) = CZERO
               SZ(I,L1) = CZERO
            END DO
         END IF
      END DO
C ======================================================================
      END
