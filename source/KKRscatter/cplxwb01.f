      SUBROUTINE CPLXWB01(PZ,FZ,QZ,SZ,TMAT,ALPHA,E,EK,KSRA,IPAN,IRCUT,
     +     IRWS,VM2Z,DRDI,R,Z,LMAX,C,DROR,RS,S,IE,IELAST,EF,
     +     REALFLAG,PHASE_SHIFT)
C ************************************************************************
c  on output: the regular radial wavefunctions (in array pz)
c             the non-regular ones corresponding to the hankel-functions
c             in free space (in array qz)
c             the alpha-matrices in array alpha
c             the t-matrices in array tmat.
c
c             the generalized phase shifts can be calculated by using
c             a wronski relation :
c
c                 alpha(z,l) =-sqrt(z)*wronski{hl(r;z),rl(r;z)}; r->0
c
c             where hl is the free hankel function and rl the regular
c             solution . using the analytical behaviour of rl at the
c             origin : rl = alphal * r**(l+1)  ; r->0
c             therefore the generalized phase shifts can be calcu-
c             lated directly with the renormalization alphal .
c
c             modified for bandstructure code
c
c                                           b.drittler nov.1989
c             valerio change 14.6.99
C ************************************************************************
C     .. Parameters ..
      include 'inc.p'
      INTEGER LMX
      PARAMETER (LMX=LMAXD+1)
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX E,EK,EF
      DOUBLE PRECISION C,Z
      INTEGER IPAN,IRWS,KSRA,KVSRA,LMAX,IE,IELAST
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX ALPHA(0:LMAXD),FZ(IRMD,0:LMAXD),PZ(IRMD,0:LMAXD),
     +               QZ(IRMD,0:LMAXD),SZ(IRMD,0:LMAXD),TMAT(0:LMAXD),
     +               PHASE_SHIFT(0:LMAXD)
      DOUBLE PRECISION DRDI(IRMD),DROR(IRMD),R(IRMD),RS(IRMD,0:LMAXD),
     +                 S(0:LMAXD),VM2Z(IRMD)
      INTEGER IRCUT(0:IPAND)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX ALPHAL,ARG,BL,CZERO,EKLFAC,FAC,HL,PN,QF,SLOPE,
     +               TLSQEZ,VALUE,W,X,Y
      DOUBLE PRECISION R1,RIRC,RSIRC,S1,PI
      INTEGER I,IRC1,L,LMAXP1,N
      LOGICAL LCALL,REALFLAG
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX BESSJW(0:LMX),BESSYW(0:LMX),DLOGDP(0:LMAXD),
     +               HAMF(IRMD,0:LMAXD),HANKWS(0:LMX),MASS(IRMD)
C     ..
C     .. External Subroutines ..
      EXTERNAL BESSEL,IRWSOL,REGSOL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC REAL,SQRT
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Data statements ..
      DATA CZERO/ (0.0D0,0.0D0)/
      DATA LCALL/.false./
C     ..
      LMAXP1 = LMAX + 1
      IF (KSRA.GE.1) THEN

         KVSRA=1

        EK = SQRT(E+E*E/ (C*C))

      ELSE

         KVSRA=0

        EK = SQRT(E)
      END IF
      ! write(6,*) 'cplx01'

      IRC1 = IRCUT(IPAN)
      RIRC = R(IRC1)
      CALL BESSEL(BESSJW,BESSYW,HANKWS,RIRC*EK,LMX,LMAXP1,.true.,.true.,
     +            .true.,LCALL)
c
c---> calculate regular wavefunctions
c
      CALL REGSOL(C,E,KVSRA,LMAX,DLOGDP,FZ,HAMF,MASS,PZ,DRDI,DROR,R,S,
     +            VM2Z,Z,IPAN,IRCUT)
      ! write(6,*) 'after regsol'
c
      EKLFAC = EK
c

      PI = 4.0D0*DATAN(1.0D0)

      DO 20 L = 0,LMAX

          S1 = S(L)
          RSIRC = RS(IRC1,L)
          EKLFAC = EKLFAC/EK*REAL(2*L+1)
c
c---> determine t - matrix
c
          PN = PZ(IRC1,L)*RSIRC
          N = L + 1
          QF = REAL(L)/RIRC
          HL = HANKWS(L)
          BL = BESSJW(L)
          X = QF*HL - EK*HANKWS(N)
          Y = QF*BL - EK*BESSJW(N)
          W = DLOGDP(L)
          TLSQEZ = (BL*W-Y)/ (X-HL*W)
          TMAT(L) = TLSQEZ/EK

c write out the phase shifts

          IF(REALFLAG == .TRUE.) then
            PHASE_SHIFT(L)=
     +          1d0/2d0*(0d0,-1d0)*log(-2*(0,1)*TLSQEZ+1d0)
c           IF (REAL(PHASE_SHIFT(L)) .LT. 0.d0) THEN
c             PHASE_SHIFT(L)=PHASE_SHIFT(L)+2d0*PI
c           END IF
c            write(601,"((2e17.9),(I5),(2e17.9))") E,L,PHASE_SHIFT(L)
          END IF
 
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
        SLOPE = REAL(L+1)*HL - RIRC*EK*HANKWS(L+1)
        SLOPE = (SLOPE*RSIRC+S1/RIRC*VALUE)
        QZ(IRC1,L) = VALUE
        SZ(IRC1,L) = (SLOPE*RIRC- (S1+1.0D0)*VALUE)/MASS(IRC1)*
     +               DROR(IRC1)
   20 CONTINUE
c
c---> calculate irregular wavefunctions
c
      CALL IRWSOL(EK,LMAX,FZ,HAMF,MASS,PZ,QZ,SZ,DROR,S,IPAN,IRCUT)

      DO 50 L = 0,LMAX
        IF (KVSRA.EQ.1) THEN
          DO 30 I = 2,IRC1
            PZ(I,L) = PZ(I,L)*RS(I,L)
            QZ(I,L) = QZ(I,L)/RS(I,L)
            FZ(I,L) = FZ(I,L)*RS(I,L)/C
            SZ(I,L) = SZ(I,L)/RS(I,L)/C
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

c
c---> in case of mt-calculation : calculate the wavefunctions between
c             mt-radius and ws-radius analytically
c
      IF (IRC1.LT.IRWS) THEN
        DO 80 I = IRC1 + 1,IRWS
          R1 = R(I)
          ARG = R1*EK
c
          CALL BESSEL(BESSJW,BESSYW,HANKWS,ARG,LMX,LMAXP1,.true.,.true.,
     +                .true.,LCALL)
c
          DO 60 L = 0,LMAX
            TLSQEZ = TMAT(L)*EK
            PZ(I,L) = (BESSJW(L)+TLSQEZ*HANKWS(L))*R1
            QZ(I,L) = HANKWS(L)*R1
            FZ(I,L) = CZERO
            SZ(I,L) = CZERO
   60     CONTINUE
c
c---> calculate small component in case of sra
c
          IF (KVSRA.EQ.1) THEN
            FAC = C/ (E+C*C)
            DO 70 L = 0,LMAX
              N = L + 1
              TLSQEZ = TMAT(L)*EK
              FZ(I,L) = (REAL(L)* (BESSJW(L)+TLSQEZ*HANKWS(L))-
     +                  ARG* (BESSJW(N)+TLSQEZ*HANKWS(N)))*FAC
              SZ(I,L) = (REAL(L)*HANKWS(L)-ARG*HANKWS(N))*FAC
   70       CONTINUE

          END IF

   80   CONTINUE

      END IF



      END




