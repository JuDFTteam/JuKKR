      MODULE mod_cradwf
      CONTAINS
      SUBROUTINE CRADWF(E,EK,NSRA,ALPHA,NPAN,NRCUT,CVLIGHT,RS,S,PZ,FZ,
     +                  QZ,SZ,TMAT,VM2Z,DRDI,R,Z,LIRRSOL,IDOLDAU,LOPT,
     +                  WLDAUAV,CUTOFF,LMAXATOM,LMAXP1,NRMAX)
      USE mod_regsol
      USE mod_beshan
      USE mod_irwsol
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
!       INCLUDE 'inc.p'
      INTEGER LMAXATOM
      INTEGER NRMAX
      INTEGER LMAXP1
!       PARAMETER (LMAXP1=LMAXATOM+1)
      DOUBLE COMPLEX CI,CZERO
      PARAMETER (CI= (0.D0,1.D0),CZERO= (0.0D0,0.0D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX E,EK
      DOUBLE PRECISION CVLIGHT,Z
      DOUBLE PRECISION WLDAUAV
      INTEGER NPAN,NSRA,IDOLDAU,LOPT
      LOGICAL LIRRSOL
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX ALPHA(0:LMAXATOM),FZ(NRMAX,0:LMAXATOM),
     +               PZ(NRMAX,0:LMAXATOM),
     +               QZ(NRMAX,0:LMAXATOM),SZ(NRMAX,0:LMAXATOM),
     +               TMAT(0:LMAXATOM)
      DOUBLE PRECISION DRDI(NRMAX),R(NRMAX),RS(NRMAX,0:LMAXATOM),
     +                 S(0:LMAXATOM),VM2Z(NRMAX),CUTOFF(NRMAX)
      INTEGER NRCUT(0:NPAN)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX ALPHAL,ARG,BL,EKLFAC,HL,PN,QF,SLOPE,TLSQEZ,VALUE
      DOUBLE PRECISION RIRC,RIRC1,RSIRC,S1
      INTEGER I,IR,IRC1,L
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX BESSJW(0:LMAXP1),BESSYW(0:LMAXP1),
     +               DLOGDP(0:LMAXATOM),
     +               HAMF(NRMAX,0:LMAXATOM),HANKWS(0:LMAXP1),MASS(NRMAX)
      DOUBLE PRECISION DROR(NRMAX)
C     ..
C     .. External Subroutines ..
  !    EXTERNAL BESHAN,IRWSOL !,REGSOL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE

C     ..
!        write(*,*) NRCUT(NPAN)
      IRC1 = NRCUT(NPAN)
      DO IR = 2,IRC1
         DROR(IR) = DRDI(IR)/R(IR)
      END DO
      RIRC = R(IRC1)
      RIRC1 = 1D0/RIRC
      ARG = RIRC*EK

      CALL BESHAN(HANKWS,BESSJW,BESSYW,ARG,LMAXP1)
c
c---> calculate regular wavefunctions
c



        
      CALL REGSOL(CVLIGHT,E,NSRA,DLOGDP,FZ,HAMF,MASS,PZ,DROR,R,S,VM2Z,
     +            Z,NPAN,NRCUT,IDOLDAU,LOPT,WLDAUAV,CUTOFF,
     +            NRMAX,NPAN,LMAXATOM)
c

      EKLFAC = EK
C ======================================================================
      DO L = 0,LMAXATOM
c
c---> determine t - matrix
c
        QF = DBLE(L)*RIRC1
        HL = HANKWS(L) * DLOGDP(L)
        BL = BESSJW(L) * DLOGDP(L)
        HL = QF*HANKWS(L) - EK*HANKWS(L+1) - HL
        BL = BL - QF*BESSJW(L) + EK*BESSJW(L+1)
        HL = HL * EK
        TMAT(L) = CI * BL/HL
c     
c---> determine the renormalization
c
        TLSQEZ = TMAT(L) * EK
        S1 = S(L)
        RSIRC = RS(IRC1,L)
        EKLFAC = EKLFAC/EK*DBLE(2*L+1)
        PN = PZ(IRC1,L)*RSIRC
        ALPHAL = (BESSJW(L) - CI*HANKWS(L)*TLSQEZ)*RIRC/PN
c
c---> determine the alpha matrix
c
        ALPHA(L) = ALPHAL*EKLFAC

        DO I = 2,IRC1
           PZ(I,L) = PZ(I,L)*ALPHAL
           FZ(I,L) = FZ(I,L)*ALPHAL
        END DO
c
        VALUE = -CI*HANKWS(L)*RIRC*RSIRC
        SLOPE = DBLE(L+1)*HANKWS(L) - RIRC*EK*HANKWS(L+1)
        SLOPE = (-CI*SLOPE*RSIRC+S1/RIRC*VALUE)
        QZ(IRC1,L) = VALUE
        SZ(IRC1,L) = (SLOPE*RIRC - (S1+1.0D0)*VALUE)/MASS(IRC1)
     &               * DROR(IRC1)
      END DO
C ======================================================================
C
C -> calculate irregular wavefunctions
C
      IF ( LIRRSOL ) CALL IRWSOL(EK,FZ,HAMF,MASS,PZ,QZ,SZ,DROR,S,
     &                           NPAN,NRCUT,NRMAX,NPAN,LMAXATOM)
C ======================================================================
      DO L = 0,LMAXATOM
         IF (NSRA.EQ.2) THEN
            DO I = 2,IRC1
               PZ(I,L) = PZ(I,L)*RS(I,L)
               QZ(I,L) = QZ(I,L)/RS(I,L)
               FZ(I,L) = FZ(I,L)*RS(I,L)/CVLIGHT
               SZ(I,L) = SZ(I,L)/RS(I,L)/CVLIGHT
            END DO
         ELSE
            DO I = 2,IRC1
               PZ(I,L) = PZ(I,L)*RS(I,L)
               QZ(I,L) = QZ(I,L)/RS(I,L)
               FZ(I,L) = CZERO
               SZ(I,L) = CZERO
            END DO
         END IF
      END DO
C ======================================================================
      END SUBROUTINE
      END MODULE mod_cradwf
