C **********************************************************************
      SUBROUTINE BESHAN(HL,JL,NL,Z,LMIN,LMAX)
C **********************************************************************
C  CALCULATES SPHERICAL BESSEL, HANKEL AND NEUMANN FUNCTIONS
C  FOR THE ORDERS 0 .LE. L .LE. LMAX.
C  FOR ARGUMENTS Z .LT. 1 THE TAYLOR EXPANSIONS OF JL AND NL ARE USED.
C  FOR ARGUMENTS Z .GE. 1 THE EXPLICIT EXPRESSIONS FOR HL(+), HL(-) ARE
C  USED.
C     .. Parameters ..
      DOUBLE COMPLEX CI
      PARAMETER (CI= (0.0D0,1.0D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX Z
      INTEGER LMAX,LMIN
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX HL(0:LMAX),JL(0:LMAX),NL(0:LMAX)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX TERMJ,TERMN,Z2,ZJ,ZN
      DOUBLE PRECISION RL,RN,RNM
      INTEGER L,M,N
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC CDABS,CDEXP
C     ..
      IF (CDABS(Z).LT.1.D0) THEN
        ZJ = 1.D0
        ZN = 1.D0
        Z2 = Z*Z
        DO 20 L = 0,LMAX
          RL = L + L
          TERMJ = -0.5D0/ (RL+3.D0)*Z2
          TERMN = 0.5D0/ (RL-1.D0)*Z2
          JL(L) = 1.D0
          NL(L) = 1.D0
          N = 1
          DO 10 WHILE ((CDABS(TERMJ) +CDABS(TERMN)).GE.1.0D-20)
             N = N + 1 
             JL(L) = JL(L) + TERMJ
             NL(L) = NL(L) + TERMN
             RN = N + N
             TERMJ = -TERMJ/ (RL+RN+1.D0)/RN*Z2
             TERMN = TERMN/ (RL-RN+1.D0)/RN*Z2
 10       CONTINUE
          JL(L) = JL(L)*ZJ
          NL(L) = -NL(L)*ZN/Z
          HL(L) = JL(L) + NL(L)*CI

          ZJ = ZJ*Z/ (RL+3.D0)
          ZN = ZN/Z* (RL+1.D0)
   20   CONTINUE

      ELSE
        DO 40 L = 0,LMAX
          HL(L) = 0.D0
          NL(L) = 0.D0
          RNM = 1.D0
          DO 30 M = 0,L
            HL(L) = HL(L) + RNM/ (-CI* (Z+Z))**M
            NL(L) = NL(L) + RNM/ (CI* (Z+Z))**M
            RNM = RNM* (L*L+L-M*M-M)/ (M+1.D0)
   30     CONTINUE
          HL(L) = HL(L)* (-CI)**L*CDEXP(CI*Z)/ (CI*Z)   ! HL +
          NL(L) = NL(L)*CI**L*CDEXP(-CI*Z)/ (-CI*Z)     ! HL -
          JL(L) = (HL(L)+NL(L))*0.5D0
          NL(L) = (HL(L)-JL(L))/CI
   40   CONTINUE
      END IF

      RETURN

      END
C **********************************************************************
      SUBROUTINE BESHAN1(HL,JL,NL,Z,LMAX)
C **********************************************************************
C  CALCULATES BESSEL, HANKEL AND NEUMANN FUNCTIONS OF INTEGER ORDER
C  FOR THE ORDERS 0 .LE. L .LE. LMAX.
C  THE TAYLOR EXPANSIONS OF JL AND NL ARE USED.
C  FROM M.ABRAMOWITSCH, I.A. STEGUN, 'HANDBOOK OF MATH. FUNCTIONS'   
c ------------------------------------------------------------------------
C
C     .. PARAMETERS ..
      DOUBLE COMPLEX CI
      PARAMETER (CI= (0.0D0,1.0D0))
C     ..
C     .. SCALAR ARGUMENTS ..
      DOUBLE COMPLEX Z
      INTEGER LMAX
C     ..
C     .. ARRAY ARGUMENTS ..
      DOUBLE COMPLEX HL(0:LMAX),JL(0:LMAX),NL(0:LMAX)
C     ..
C     .. LOCAL SCALARS ..
      DOUBLE COMPLEX TERMJ,TERMN,TERMN1,Z2,ZJ,ZN
      DOUBLE PRECISION L0,LN,LLN,RL,RN,RNM,PI
      INTEGER L,M,N
C     ..
C     .. INTRINSIC FUNCTIONS ..
      INTRINSIC CDABS,CDEXP
C     ..
      PI = 4.0D0 *DATAN(1.0D0)
      ZJ = 1.D0
      ZN = -2.0D0*(CDLOG(Z/2.0D0) +0.5772156649D0)
      Z2 = Z*Z
      DO 20 L = 0,LMAX
         RL = L + L
         L0 = L
         TERMJ = ZJ
         TERMN1= ZN
         JL(L) = ZJ
         IF (L.EQ.0) THEN 
            NL(L) = ZJ*ZN
            TERMN = 1.0D0 /ZJ
         ELSE 
            NL(L) = 1.0D0 /ZJ /L0 + ZJ*ZN
            TERMN = 1.0D0 /ZJ /L0
         END IF   
         N=0
         DO 10 WHILE ((CDABS(TERMJ)+CDABS(TERMJ*TERMN1)).GE.1.0D-20)
            N  = N + 1
            RN = N + N
            LN = N + L
            IF (N.LT.L) THEN
               LLN = L - N
               TERMN = 0.5D0 *TERMN /RN /LLN *Z2 
               NL(L) = NL(L) + TERMN
            END IF   
C     *********************************************************
C     
C     TERMJ= (-Z*Z/4)**N L! /N! /(N+L)!
C     
C                    Z              1         1
C     TERMN1 = -2(LN(-)-EXP(-1)) + (- + ... + -) +
C                    2              1         L
C     
C                1        1     1        1
C             +( - + .. + - )+( - + .. + - )
C               L+1      L+N    1        N
C     
C              WRITE(6,9000) TERMJ
C 9000         FORMAT(2F35.20)
C     *********************************************************
            TERMJ  = - 0.5D0 *TERMJ /RN /LN *Z2 
            TERMN1 = TERMN1 + 2.0D0 /RN + 1.0D0 / LN
            JL(L) = JL(L) + TERMJ
            NL(L) = NL(L) + TERMN1 * TERMJ
 10      CONTINUE                       ! WHILE (ABS() > ..)
         NL(L) = -NL(L)/PI
         HL(L) = JL(L) + NL(L)*CI
         
         ZJ = ZJ *Z /(RL+2.D0)
         ZN = ZN +2.0D0 /(RL+2.0D0)
C        WRITE(6,9001) L,ZJ
C 9001   FORMAT(I10,2F35.20)
 20   CONTINUE                          !  L = 0,LMAX
      
      RETURN

      END
c *******************************************************************
      SUBROUTINE BESSEL(BJ,Y,H,ARG,LMX,LMAX,LJ,LY,LH,LCALL)
c *******************************************************************
c   * this subroutine computes the spherical bessel functions of
c   * first, second and third kind using a chebychev expansion
c   * given by y.l.luke ,algorithms for the computation of
c   * mathematical functions, academic press,london 1977
c   *
c   *
c   * using subroutine cnwf01
c   *
c   *
c   * description of variables
c   *
c   * arg   -input  - argument of the bessel functions
c   *
c   * lmax  -input  - max. order of the bessel functions
c   *                (limited up to 25 in that version)
c   *
c   * lj    -input  - logical : if lj is true the spherical bessel
c   *                 functions of the first kind are calculated
c   *                 up to lmax
c   *
c   * ly    -input  - logical : if ly is true the spherical bessel
c   *                 functions of the second kind are calculated
c   *                 up to lmax
c   *
c   * lh    -input  - logical : if lh is true the spherical bessel
c   *                 functions of the third kind are calculated
c   *                 up to lmax
c   *
c   * lcall -input  - logical : if lh is false the chebychev coefficients
c   *                 are calculated - this part has to be called once
c   *
c   *
c   * bj    -output - an array containing the bessel functions of the
c   *                 first kind up to lmax if lj is true . remember ,
c   *                 that bj(1) contains the function of l=0 and so on.
c   *
c   * y     -output - an array containing the bessel functions of the
c   *                 second kind up to lmax if ly is true . remember ,
c   *                 that y(1) contains the function of l=0 and so on.
c   *
c   * h     -output - an array containing the bessel functions of the
c   *                 third kind up to lmax if lh is true . remember ,
c   *                 that h(1) contains the function of l=0 and so on.
c   *
c   *        attention : contrary to abramowitz and stegun the bessel
c   *                    functions of third kind ( hankel functions)
c   *                    are definied as:
c   *                            h(l) = y(l) - i * bj(l)
c   *
c   *
c   * all other variables are for internal use
c **********************************************************************
C     .. Parameters ..
      INTEGER NDIM
      PARAMETER (NDIM=24)
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX ARG
      INTEGER LMAX,LMX
      LOGICAL LCALL,LH,LJ,LY
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX BJ(LMX+1),H(LMX+1),Y(LMX+1)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX CI,CONE,CTWO,CZERO,RES,T0,T1,T2,TARG,TT1
      DOUBLE PRECISION CPJ,CPY,FJ,FY,ONE,SUM,W,W2
      INTEGER L,LMSAVE,N,NCNW,NMAX
C     ..
      PARAMETER (NMAX = 20)
C     .. Local Arrays ..
      DOUBLE COMPLEX TN(NMAX),ZN(NDIM+2)
      DOUBLE PRECISION CNWJ(NMAX,NDIM+1),CNWY(NMAX,NDIM+1)
C     ..
C     .. External Subroutines ..
      EXTERNAL CNWF01,RCSTOP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     ..
C     .. Data statements ..
      DATA CZERO,CONE,CTWO / (0.D0,0.D0), (1.D0,0.D0), (2.D0,0.D0)/
      DATA CI              / (0.D0,1.D0) /
      DATA W2,ONE          / 10.D0,1.D0  /
      DATA LMSAVE          / 0 /
C     ..
C     .. Save statement ..
      SAVE
c
c ------------------------------------------------------------------------
      IF (.NOT.LCALL) THEN
        IF (LMAX.GT.NDIM) THEN
          WRITE (6,FMT=9000) LMAX
          CALL RCSTOP('27      ')
c
        ELSE
          LMSAVE = LMAX + 1
          NCNW = NMAX - 2
          W = -W2*W2*.25D0
          FJ = ONE
          FY = -ONE
          DO 20 L = 1,LMSAVE
            CPJ = 0.5D0 + L
            CPY = 1.5D0 - L
            CALL CNWF01(CPJ,W,NCNW,CNWJ(1,L),SUM)
            CALL CNWF01(CPY,W,NCNW,CNWY(1,L),SUM)
            DO 10 N = 1,NMAX
              CNWJ(N,L) = CNWJ(N,L)/FJ
              CNWY(N,L) = CNWY(N,L)*FY
   10       CONTINUE
            FJ = FJ* (L+L+1)
            FY = FY* (L+L-1)
   20     CONTINUE
          LCALL = .true.
        END IF
c
      END IF
c
      IF (LMAX.GT. (LMSAVE-1) .OR. ABS(ARG).GT.W2) THEN
        WRITE (6,FMT=9010) LMAX,LMSAVE-1,ABS(ARG),W2
        CALL RCSTOP('28      ')
c
      ELSE
c
c--->   calculate arg**n and tn*(-arg**2/4)
c
        ZN(1) = CONE
        DO 30 L = 2,LMSAVE+1
          ZN(L) = ZN(L-1)*ARG
   30   CONTINUE
c
        T0 = CONE
        TARG = -ZN(3)*0.25D0/W
        T1 = CTWO*TARG - CONE
        TN(1) = T0
        TN(2) = T1
        TT1 = T1 + T1
        DO 40 N = 3,NMAX
          T2 = TT1*T1 - T0
          TN(N) = T2
          T0 = T1
          T1 = T2
   40   CONTINUE
c
        IF (LJ .OR. LH) THEN
          DO 60 L = 1,LMSAVE
            RES = CZERO
            DO 50 N = 1,NMAX
              RES = RES + TN(N)*CNWJ(N,L)
   50       CONTINUE
            BJ(L) = RES*ZN(L)
   60     CONTINUE
          IF (.NOT. (LY.OR.LH)) RETURN
        END IF
c
        IF (LY .OR. LH) THEN
          DO 80 L = 1,LMSAVE
            RES = CZERO
            DO 70 N = 1,NMAX
              RES = RES + TN(N)*CNWY(N,L)
   70       CONTINUE
            Y(L) = RES/ZN(L+1)
   80     CONTINUE
          IF (LH) THEN
            DO 90 L = 1,LMSAVE
c
c--->         changed !   h(l) = (-i)*( bj(l) + ci*y(l) )
c
              H(L) = Y(L) - CI*BJ(L)
   90       CONTINUE
          END IF
c
        ELSE
          WRITE (6,FMT=9020)
        END IF
c
      END IF
c
c
 9000 FORMAT (' lmax is too high, error stop in bessel *******',/,
     +       ' change NDIM to ',i4,/)
 9010 FORMAT (' lmax is higher than previously given',/,
     +     ' or argument too high, error stop in bessel ************'
     +     ,/,
     +     ' lmax     : ',i5,/,
     +     ' lmsave-1 : ',i5,/,
     +     ' abs(arg) : ',f12.6,/,
     +     ' w2       : ',f12.6)
 9020 FORMAT (' *******  error warning from bessel: no output',
     +       ' required   ********',/,/)
      END
C ************************************************************************
      SUBROUTINE GAUNT(LMAX,LPOT,W,YR,CLEB,LOFLM,ICLEB,IEND,JEND)
C ************************************************************************
c
c   - fills the array cleb with the gaunt coeffients ,i.e.
c      the integral of y(l1,m1)*y(l2,m2)*y(l3,m3)
c      but only for lm2.le.lm1 and lm3>1
c   - calculate the pointer array jend  to project the indices
c      array cleb with the same lm3,l1,l2 values - because of
c      the special ordering of array cleb only the last index
c      has to be determined .
c     (the parameter n has to be chosen that l1+l2+l3 .lt. 2*n)
c     using gaussian quadrature as given by
c     m. abramowitz and i.a. stegun, handbook of mathematical functions,
c     nbs applied mathematics series 55 (1968), pages 887 and 916
c     m. weinert and e. wimmer
c     northwestern university march 1980
c
c     an index array -icleb- is used to save storage place .
c     fills the array loflm which is used to determine the
c     l-value of a given lm-value .
c     this subroutine has to be called only once !
c
c                               b.drittler   november 1987
c
c     modified gaunt coefficients are als calculated defined by
c     the integral of y(l1,m1)*y(l2,m2)*y(l3,m3)*i**(l2-l1+l3)
c-----------------------------------------------------------------------
c
c---> attention : ncleb is an empirical factor - it has to be optimized
c
C     .. Parameters ..
      include 'inc.p'
C     ..
      DOUBLE COMPLEX CI
      PARAMETER (CI= (0.0D0,1.0D0))
C     ..
      INTEGER LMMAXD,LMPOTD
      PARAMETER (LMMAXD= (LMAXD+1)**2,LMPOTD= (LPOTD+1)**2)
      INTEGER N
      PARAMETER (N=4*LMAXD)
c      INTEGER NCLEB
c      PARAMETER (NCLEB=LMPOTD*LMMAXD)
C     ..
C     .. Scalar Arguments ..
      INTEGER IEND,LMAX,LPOT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CLEB(NCLEB,2),W(*),YR(N,0:N,0:N)
      INTEGER ICLEB(NCLEB,4),JEND(LMPOTD,0:LMAXD,0:LMAXD),LOFLM(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CLECG,FACTOR,FCI,S
      INTEGER I,J,L,L1,L1P,L2,L2P,L3,LM1,LM2,LM3,LM3P,LMPOT,M,M1,M1A,
     +        M1S,M2,M2A,M2S,M3,M3A,M3S
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DREAL,MOD,REAL,SIGN
C     ..
C     .. External Subroutines ..
      EXTERNAL RCSTOP
C     ..
      I = 1
      DO 20 L = 0,2*LMAX
        DO 10 M = -L,L
          LOFLM(I) = L
          I = I + 1
   10   CONTINUE
   20 CONTINUE
c
      IF (LPOT.EQ.0) THEN
        IEND = 1
        ICLEB(1,1) = (LMAX+1)**2
        ICLEB(1,3) = 1
      END IF
c
      IF (LPOT.NE.0) THEN
c
c---> set up of the gaunt coefficients with an index field
c
        I = 1
        DO 90 L3 = 1,LPOT
          DO 80 M3 = -L3,L3
c
            DO 70 L1 = 0,LMAX
              DO 60 L2 = 0,L1
c
                IF (MOD((L1+L2+L3),2).NE.1 .AND. (L1+L2-L3).GE.0 .AND.
     +              (L1-L2+L3).GE.0 .AND. (L2-L1+L3).GE.0) THEN
c
                  FCI = DREAL(CI** (L2-L1+L3))
                  DO 50 M1 = -L1,L1
                    DO 40 M2 = -L2,L2
c
c---> store only gaunt coeffients for lm2.le.lm1
c
                      LM1 = L1* (L1+1) + M1 + 1
                      LM2 = L2* (L2+1) + M2 + 1
                      IF (LM2.LE.LM1) THEN
c
                        M1S = SIGN(1,M1)
                        M2S = SIGN(1,M2)
                        M3S = SIGN(1,M3)
c
                        IF (M1S*M2S*M3S.GE.0) THEN
c
                          M1A = ABS(M1)
                          M2A = ABS(M2)
                          M3A = ABS(M3)
c
                          FACTOR = 0.0
c
                          IF (M1A+M2A.EQ.M3A) FACTOR = FACTOR +
     +                        REAL(3*M3S+SIGN(1,-M3))/8.0d0
                          IF (M1A-M2A.EQ.M3A) FACTOR = FACTOR +
     +                        REAL(M1S)/4.0d0
                          IF (M2A-M1A.EQ.M3A) FACTOR = FACTOR +
     +                        REAL(M2S)/4.0d0
c
                          IF (FACTOR.NE.0.0) THEN
c
                            IF (M1S*M2S.NE.1 .OR. M2S*M3S.NE.1 .OR.
     +                          M1S*M3S.NE.1) FACTOR = -FACTOR
c
                            S = 0.0
                            DO 30 J = 1,N
                              S = S + W(J)*YR(J,L1,M1A)*YR(J,L2,M2A)*
     +                            YR(J,L3,M3A)
   30                       CONTINUE
                            CLECG = S*FACTOR
                            IF (ABS(CLECG).GT.1.D-10) THEN
                              CLEB(I,1) = CLECG
                              CLEB(I,2) = FCI*CLECG
                              ICLEB(I,1) = LM1
                              ICLEB(I,2) = LM2
                              ICLEB(I,3) = L3* (L3+1) + M3 + 1
                              ICLEB(I,4) = LM2*LMMAXD -
     +                                     (LM2*LM2-LM2)/2 + LM1 -
     +                                     LMMAXD
                              I = I + 1
                            END IF

                          END IF

                        END IF

                      END IF

   40               CONTINUE
   50             CONTINUE
                END IF

   60         CONTINUE
   70       CONTINUE
   80     CONTINUE
   90   CONTINUE
        IEND = I - 1
        IF (NCLEB.LT.IEND) THEN
          WRITE (6,FMT=9000) NCLEB,IEND
          CALL RCSTOP('33      ')

        ELSE
c
c---> set up of the pointer array jend,use explicitly
c     the ordering of the gaunt coeffients
c
          LMPOT = (LPOT+1)* (LPOT+1)
          DO 120 L1 = 0,LMAX
            DO 110 L2 = 0,L1
              DO 100 LM3 = 2,LMPOT
                JEND(LM3,L1,L2) = 0
  100         CONTINUE
  110       CONTINUE
  120     CONTINUE
c
          LM3 = ICLEB(1,3)
          L1 = LOFLM(ICLEB(1,1))
          L2 = LOFLM(ICLEB(1,2))
c
          DO 130 J = 2,IEND
            LM3P = ICLEB(J,3)
            L1P = LOFLM(ICLEB(J,1))
            L2P = LOFLM(ICLEB(J,2))
c
            IF (LM3.NE.LM3P .OR. L1.NE.L1P .OR. L2.NE.L2P) THEN
              JEND(LM3,L1,L2) = J - 1
              LM3 = LM3P
              L1 = L1P
              L2 = L2P
            END IF

  130     CONTINUE
          JEND(LM3,L1,L2) = IEND
c
c
        END IF

      END IF
c
c

 9000 FORMAT (13x,'error stop in gaunt : dimension of NCLEB = ',i10,
     +       ' too small ',/,
     +       13x,'change NCLEB to ',i6)
      END
C***********************************************************************
      SUBROUTINE GAUNT2(W,YR)
C ************************************************************************
c     sets up values needed for gaunt
c        m. weinert  january 1982
c
c     changed for calculating with real spherical harmonics
c                                           b.drittler  july 1987
c
c     W(N)        integration weights on 4*LMAXD points in the intervall
c                 (-1,0) (from routine GRULE)
c
c     YR(N,L,M)   spherical harmonics on 4*LMAXD points to angular 
c                 momentum indices (l,m) scaled with a factor 
c                 of RF=(4*pi)**(1/3)
c
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.p'
      INTEGER N
      PARAMETER (N=4*LMAXD)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,CD,CTH,FAC,FPI,RF,STH,T
      INTEGER K,L,LOMAX,M
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION P(0:N+1,0:N),X(N)
C     ..
C     .. External Subroutines ..
      EXTERNAL GRULE
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,SQRT
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION W(*),YR(N,0:N,0:N)
C     ..
      FPI = 16.*ATAN(1.0)
      RF = FPI** (1.0/3.0)
      LOMAX = N
c
c--->    obtain gauss-legendre points and weights
c
      CALL GRULE(2*N,X,W)
c
c--->    generate associated legendre functions for m.ge.0
c
      DO 50 K = 1,N
        CTH = X(K)
        STH = SQRT(1.0-CTH*CTH)
        FAC = 1.0
c
c--->    loop over m values
c
        DO 20 M = 0,LOMAX
          FAC = - (2*M-1)*FAC
          P(M,M) = FAC
          P(M+1,M) = (2*M+1)*CTH*FAC
c
c--->    recurse upward in l
c
          DO 10 L = M + 2,LOMAX
            P(L,M) = ((2*L-1)*CTH*P(L-1,M)- (L+M-1)*P(L-2,M))/ (L-M)
   10     CONTINUE
          FAC = FAC*STH
   20   CONTINUE
c
c--->    multiply in the normalization factors
c
        DO 40 L = 0,LOMAX
          A = RF*SQRT((2*L+1)/FPI)
          CD = 1
          YR(K,L,0) = A*P(L,0)
          DO 30 M = 1,L
            T = (L+1-M)* (L+M)
            CD = CD/T
            YR(K,L,M) = A*SQRT(2.0*CD)*P(L,M)
   30     CONTINUE
   40   CONTINUE
   50 CONTINUE
      END
C***********************************************************************
      SUBROUTINE GAUNT3(W,YR,CTH,P) ! test of Legendre polynoms
C ************************************************************************
c     sets up values needed for gaunt
c        m. weinert  january 1982
c
c     changed for calculating with real spherical harmonics
c                                           b.drittler  july 1987
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.p'
      INTEGER N
      PARAMETER (N=4*LMAXD)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,CD,CTH,FAC,FPI,RF,STH,T
      INTEGER I,J,K,L,LOMAX,M
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION P(0:N+1,0:N),X(N)
C     ..
C     .. External Subroutines ..
      EXTERNAL GRULE
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,SQRT
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION W(*),YR(N,0:N,0:N)
C     ..
      LOMAX = N
      FPI = 16.d0*ATAN(1.0d0)
      RF = FPI** (1.0d0/3.0d0)
c
c--->    obtain gauss-legendre points and weights
c
c      CALL GRULE(2*N,X,W)
c
c--->    generate associated legendre functions for m.ge.0
c
c      DO 50 K = 1,N
      K = 1
c        CTH = X(K)
        STH = SQRT(1.0-CTH*CTH)
        FAC = 1.0
c
c--->    loop over m values
c
        DO 20 M = 0,LOMAX
          FAC = - (2*M-1)*FAC
          P(M,M) = FAC
          P(M+1,M) = (2*M+1)*CTH*FAC
c
c--->    recurse upward in l
c
          DO 10 L = M + 2,LOMAX
            P(L,M) = ((2*L-1)*CTH*P(L-1,M)- (L+M-1)*P(L-2,M))/ (L-M)
   10     CONTINUE
          FAC = FAC*STH
   20   CONTINUE
c
c--->    multiply in the normalization factors
c
        DO 40 L = 0,LOMAX
          A = RF*SQRT((2*L+1)/FPI)
          CD = 1
          YR(K,L,0) = A*P(L,0)
          DO 30 M = 1,L
            T = (L+1-M)* (L+M)
            CD = CD/T
            YR(K,L,M) = A*SQRT(2.0*CD)*P(L,M)
   30     CONTINUE
   40   CONTINUE
c   50 CONTINUE
      END
c ************************************************************************
      SUBROUTINE HANKEL(H,L,ARG)
c ************************************************************************
c  this subroutine uses the explicit formulas for the hankel
c  functions. for higher l-values these formulas may lead to
c  loss of significant figures.
c ------------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE COMPLEX ARG
      INTEGER L
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX H(*)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX A1,A2,A3,A4,I
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC EXP
C     ..
C     .. Save statement ..
      SAVE I
C     ..
C     .. External Subroutines ..
      EXTERNAL RCSTOP
C     ..
C     .. Data statements ..
      DATA I/ (0.D0,1.D0)/
C     ..

      H(1) = -EXP(ARG*I)/ARG
      IF (L.NE.1) THEN
        A1 = (1.D0,0.D0) - ARG*I
        H(2) = H(1)*A1/ARG
        IF (L.NE.2) THEN
          A1 = 3.D0*A1
          A2 = ARG*ARG
          H(3) = H(1)* (A1-A2)/A2
          IF (L.NE.3) THEN
            A1 = 5.D0*A1
            A3 = A2*ARG*I
            A4 = A2*ARG
            A2 = 6.D0*A2
            H(4) = H(1)* (A1-A2+A3)/A4
            IF (L.NE.4) THEN
              A1 = 7.D0*A1
              A2 = 7.5D0*A2
              A3 = 10.D0*A3
              A4 = A4*ARG
              H(5) = H(1)* (A1-A2+A3+A4)/A4
              IF (L.NE.5) THEN
                H(6) = (9.0D0,0.0D0)*H(5)/ARG - H(4)
                IF (L.NE.6) THEN
                  WRITE (6,FMT=9000) L
                  CALL RCSTOP('34      ')

                END IF

              END IF

            END IF

          END IF

        END IF

      END IF

      RETURN

 9000 FORMAT (2x,' hankel :  l=',i2,' is too large')
      END
c ************************************************************************
      SUBROUTINE SHAPE(LPOT,NATYP,GSH,ILM,IMAXSH,LMSP,NTCELL,W,YR)
c ************************************************************************
c   - prepares shape corrections
c     (the parameter n has to be chosen that l1+l2+l3 .le. 2*n)
c     using gaussian quadrature as given by
c     m. abramowitz and i.a. stegun, handbook of mathematical functions,
c     nbs applied mathematics series 55 (1968), pages 887 and 916
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.p'
c      INTEGER NATYPD
c      PARAMETER (NATYPD=1)
c      INTEGER LMAXD,LPOTD
c      PARAMETER (LMAXD=4,LPOTD=8)
c      INTEGER NGSHD
c      PARAMETER (NGSHD=3079)
      INTEGER N,LASSLD,LMPOTD
      PARAMETER (N=4*LMAXD,LASSLD=N,LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER LPOT,NATYP
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION GSH(*),W(N),YR(N,0:LASSLD,0:LASSLD)
      INTEGER ILM(NGSHD,3),IMAXSH(0:LMPOTD),LMSP(NATYPD,*),NTCELL(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION FACTOR,GAUNT,S
      INTEGER I,IAT,ICELL,ISUM,J,L1,L2,L3,LM1,LM2,LM3,M1,M1A,M1S,M2,M2A,
     +        M2S,M3,M3A,M3S
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MOD,REAL,SIGN
C     ..
C     .. External Subroutines ..
      EXTERNAL RCSTOP
C     ..
c
c---> set up of the gaunt coefficients with an index field
c     so c(lm,lm',lm'') is mapped to c(i)
      I = 1
      DO 80 L1 = 0,LPOT
        DO 70 M1 = -L1,L1
          LM1 = L1*L1 + L1 + M1 + 1
          IMAXSH(LM1-1) = I - 1
          DO 60 L3 = 0,LPOT*2
            DO 50 M3 = -L3,L3
              LM3 = L3*L3 + L3 + M3 + 1
              ISUM = 0
              DO 10 IAT = 1,NATYP
                ICELL = NTCELL(IAT)
                ISUM = ISUM + LMSP(ICELL,LM3)
   10         CONTINUE
              IF (ISUM.GT.0) THEN
                DO 40 L2 = 0,LPOT
                  IF (MOD((L1+L2+L3),2).NE.1 .AND. (L1+L2-L3).GE.0 .AND.
     +                (L1-L2+L3).GE.0 .AND. (L2-L1+L3).GE.0) THEN
                    DO 30 M2 = -L2,L2
                      LM2 = L2* (L2+1) + M2 + 1
c---> use the m-conditions for the gaunt coefficients not to be 0
                      M1S = SIGN(1,M1)
                      M2S = SIGN(1,M2)
                      M3S = SIGN(1,M3)
c
                      IF (M1S*M2S*M3S.GE.0) THEN
                        M1A = ABS(M1)
                        M2A = ABS(M2)
                        M3A = ABS(M3)
c
                        FACTOR = 0.0D0
c
                        IF (M1A+M2A.EQ.M3A) FACTOR = FACTOR +
     +                      REAL(3*M3S+SIGN(1,-M3))/8.0D0
                        IF (M1A-M2A.EQ.M3A) FACTOR = FACTOR +
     +                      REAL(M1S)/4.0D0
                        IF (M2A-M1A.EQ.M3A) FACTOR = FACTOR +
     +                      REAL(M2S)/4.0D0
c
                        IF (FACTOR.NE.0.0D0) THEN
c
                          IF (M1S*M2S.NE.1 .OR. M2S*M3S.NE.1 .OR.
     +                        M1S*M3S.NE.1) FACTOR = -FACTOR
                          S = 0.0D0
                          DO 20 J = 1,N
                            S = S + W(J)*YR(J,L1,M1A)*YR(J,L2,M2A)*
     +                          YR(J,L3,M3A)
   20                     CONTINUE
                          GAUNT = S*FACTOR
                          IF (ABS(GAUNT).GT.1.D-10) THEN
                            GSH(I) = GAUNT
                            ILM(I,1) = LM1
                            ILM(I,2) = LM2
                            ILM(I,3) = LM3
                            I = I + 1
                          END IF

                        END IF


                      END IF

   30               CONTINUE
                  END IF

   40           CONTINUE



              END IF

   50       CONTINUE

   60     CONTINUE

   70   CONTINUE

   80 CONTINUE
      IMAXSH(LM1) = I - 1

      WRITE (*,FMT=9000) IMAXSH(LM1),NGSHD
      IF (IMAXSH(LM1).GT.NGSHD) CALL RCSTOP('SHAPE   ')
c
 9000 FORMAT(' >>> SHAPE : IMAXSH(',I4,'),NGSHD :',2I6)
c
      END
c **********************************************************************
      SUBROUTINE SPHERE_NOGGA(LMAX,YR,WTYR,RIJ,IJEND,IJD,W) 
c **********************************************************************
C       IJD=434
c-----------------------------------------------------------------------
c     generate an angular mesh and spherical harmonics at those
c     mesh points. For an angular integration the weights are ge-
c     rated .
c
c     R. Zeller      Feb. 1996
c-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER IJD,IJEND,LMAX
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION PI,R,R1,R2,R3
      INTEGER IJ,LM1
C     ..
C     .. External Subroutines ..
      EXTERNAL YMY
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION RIJ(IJD,3),WTYR(IJD,*),YR(IJD,*)
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION W(1000),Y(1000)
C     ..
C     .. Data statements ..
      DATA PI/3.14159265359D0/
C     ..
c      write (6,*) 'SPHERE : read lebedev.mesh' 
       write (6,*) 'SPHERE : calculate lebedev.mesh' 
c      READ (89,*) IJEND
       IJEND=IJD
c      IF (IJEND.GT.IJD .OR. IJEND.GT.1000) STOP 'SPHERE'
      IF (IJD.GT.1000) STOP 'SPHERE'
c
c
      DO 30 IJ = 1,IJD
        CALL LEBEDEV(IJ,R1,R2,R3,W(IJ))
c        READ (89,*) R1,R2,R3,W(IJ)
c        write (189,"(3(e15.7,2x))") R1,R2,R3
        RIJ(IJ,1) = R1
        RIJ(IJ,2) = R2
        RIJ(IJ,3) = R3
        CALL YMY(R1,R2,R3,R,Y,LMAX)
        DO 10 LM1 = 1, (LMAX+1)**2
          YR(IJ,LM1) = Y(LM1)
   10   CONTINUE
c
c--->   multiply the spherical harmonics with the weights
c
        DO 20 LM1 = 1, (LMAX+1)**2
          WTYR(IJ,LM1) = YR(IJ,LM1)*W(IJ)*PI*4.D0
   20   CONTINUE
   30 CONTINUE

      END
c ************************************************************************
      SUBROUTINE SPHERE1(LMAX,YR,WTYR,RIJ,IJEND,IJDD)
c ************************************************************************
c     generate a angular mesh and the spherical harmonics at those
c     mesh points . for an angular integration the weights are ge-
c     rated .
c     the spherical harmonics and the spherical harmonics times the
c     weights generated on the unit sphere , the points and the
c     weights for the gauss-legendre integration are stored in the
c     common block rsphere .
c     for a better vectorization combined indices are introduced .
c
c                                  b.drittler    may 1987
c
c
c     using gaussian quadrature as given by
c     m. abramowitz and i.a. stegun, handbook of mathematical functions,
c     nbs applied mathematics series 55 (1968), pages 887 and 916
c     m. weinert and e. wimmer
c     northwestern university march 1980
c
c
c     modified to use calculated points and weights
c     to make it dynamic.   (m.w.  jan. 1982)
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.p'
c      INTEGER LPOTD
c      PARAMETER (LPOTD=8)
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (LPOTD+1)**2)
      INTEGER N
      PARAMETER (N=2* (LPOTD+1))
C     ..
C     .. Scalar Arguments ..
      INTEGER IJDD,IJEND,LMAX
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION FAC,HW,PI,R,R1,R2,R3,WJ
      INTEGER I,IJ,IS,J,LM1,NN
C     ..
C     .. External Subroutines ..
      EXTERNAL GRULE,YMY
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC COS,REAL,SIN,SQRT
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION RIJ(IJD,3),WTYR(IJD,*),YR(IJD,*)
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION W(N),X(N),Y(LMMAXD)
C     ..
C     .. Data statements ..
      DATA PI/ 3.14159265358979312d0/
C     ..
      write(6,*) '>>> SPHERE: angular mesh and spherical harmonics'
c
c---> obtain gauss-legendre points and weights
c
      NN = 2*N
      CALL GRULE(N,X,W)
c
c---> obtain gauss-tschebyscheff weights
c
      WJ = PI/REAL(N)
      HW = WJ/2.0D0
c
c---> set up of the real spherical harmonics (i=theta ,j=phi)
c     a vector on the unit sphere is given by     x=sin(theta)cos(phi)
c                                                 y=sin(theta)sin(phi)
c                                                 z=    cos(theta)
c
      IJ = 0
c
c---> is=0 :give the points with positive theta
c     is=1 :give the points with negative theta
c
      FAC = 1.0D0
      DO 50 IS = 0,1
        DO 40 I = 1,N/2
          DO 30 J = 1,NN
            IJ = IJ + 1
            R1 = SQRT(1.D0-X(I)*X(I))*COS(REAL(2*J-1)*HW)
            R2 = SQRT(1.D0-X(I)*X(I))*SIN(REAL(2*J-1)*HW)
            R3 = FAC*X(I)
c
            RIJ(IJ,1) = R1
            RIJ(IJ,2) = R2
            RIJ(IJ,3) = R3
c
            CALL YMY(R1,R2,R3,R,Y,LMAX)
            DO 10 LM1 = 1, (LMAX+1)**2
              YR(IJ,LM1) = Y(LM1)
   10       CONTINUE
c
c---> multiply the spherical harmonics with the weights
c
            DO 20 LM1 = 1, (LMAX+1)**2
              WTYR(IJ,LM1) = YR(IJ,LM1)*W(I)*WJ
   20       CONTINUE
   30     CONTINUE
   40   CONTINUE
        FAC = -FAC
   50 CONTINUE
      IJEND = IJ

      END
c ************************************************************************
      SUBROUTINE YMY(V1,V2,V3,R,YLM,LMAX)
c ************************************************************************
c    this subroutine calculates real spherical harmonics with the
c     normalization : <y|y> =1
c    returns also r = length of vector v
c
c     generate the complex spherical harmonics for the vector v
c     using a stable upward recursion in l.  (see notes
c     by m. weinert.)
c                                  m.weinert  1982
c
c     converted to real spherical harmonics .
c                                  b.drittler 1987
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.p'
c      INTEGER LMAXD
c      PARAMETER (LMAXD=4)
      INTEGER L4MAXD
      PARAMETER (L4MAXD=4*LMAXD)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION R,V1,V2,V3
      INTEGER LMAX
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION YLM(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,CD,CPH,CTH,FAC,FPI,PI,RTWO,SGM,SNULL,SPH,STH,T,
     +                 XY,XYZ
      INTEGER I,L,M
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION C(0:L4MAXD),P(0:L4MAXD,0:L4MAXD),S(0:L4MAXD)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DATAN,SQRT
C     ..
C     .. Save statement ..
      SAVE SNULL
C     ..
C     .. External Subroutines ..
      EXTERNAL RCSTOP
C     ..
C     .. Data statements ..
      DATA SNULL/1.0D-20/
C     ..
      PI = 4.D0*DATAN(1.D0)
      FPI = 4.D0*PI
      RTWO = SQRT(2.0D0)
c
      IF (LMAX.GT.L4MAXD) THEN
        CALL RCSTOP('ylm3    ')

      ELSE

c
c--->    calculate sin and cos of theta and phi
c
        XY = V1**2 + V2**2
        XYZ = XY + V3**2
c        write(6,*) "XYZ", XYZ
c
        R = SQRT(XYZ)
        IF (XYZ.LE.0.0D0) THEN
          CALL RCSTOP('ylm=0   ')

        ELSE

          IF (XY.GT.SNULL*XYZ) THEN
            XY = SQRT(XY)
            XYZ = SQRT(XYZ)
            CTH = V3/XYZ
            STH = XY/XYZ
            CPH = V1/XY
            SPH = V2/XY

          ELSE

            STH = 0.0D0
            CTH = 1.0D0
            IF (V3.LT.0) CTH = -1.0D0
            CPH = 1.0D0
            SPH = 0.0D0
          END IF
c
c--->    generate associated legendre functions for m.ge.0
c        loop over m values
c
          FAC = 1.0D0
          DO 20 M = 0,LMAX - 1
            FAC = - (2*M-1)*FAC
            P(M,M) = FAC
            P(M+1,M) = (2*M+1)*CTH*FAC
c
c--->    recurse upward in l
c
            DO 10 L = M + 2,LMAX
              P(L,M) = ((2*L-1)*CTH*P(L-1,M)- (L+M-1)*P(L-2,M))/ (L-M)
   10       CONTINUE
            FAC = FAC*STH
   20     CONTINUE
          P(LMAX,LMAX) = - (2*LMAX-1)*FAC
c
c--->    determine sin and cos of phi
c
          S(0) = 0.0D0
          S(1) = SPH
          C(0) = 1.0D0
          C(1) = CPH
          DO 30 M = 2,LMAX
            S(M) = 2*CPH*S(M-1) - S(M-2)
            C(M) = 2*CPH*C(M-1) - C(M-2)
   30     CONTINUE
c
c--->    multiply in the normalization factors
c
          I = 0
          DO 50 L = 0,LMAX
            I = I + L + 1
            A = SQRT((2*L+1)/FPI)
            CD = 1
            YLM(I) = A*P(L,0)
            SGM = -RTWO
            DO 40 M = 1,L
              T = (L+1-M)* (L+M)
              CD = CD/T
              T = A*SQRT(CD)
              YLM(I+M) = SGM*T*P(L,M)*C(M)
              YLM(I-M) = SGM*T*P(L,M)*S(M)
              SGM = -SGM
   40       CONTINUE
            I = I + L
   50     CONTINUE

        END IF

      END IF

      RETURN

      END
C *********************************************************************
      SUBROUTINE YSHYSH(X,Y,Z,R,YREALY)
C *********************************************************************
C     .. Parameters ..
      include 'inc.p'
      INTEGER LMAX
      PARAMETER (LMAX=LMAXD)
      INTEGER LMAX2
      PARAMETER (LMAX2=2*LMAX)
      INTEGER LMAX2P,LMXP
      PARAMETER (LMAX2P=LMAX2+1,LMXP=LMAX2P* (LMAX2P+1)/2)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION R,X,Y,Z
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION YREALY(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,ARG,B,C,COSPHI,COSTHE,EAT,P,PSQ,RSSQ,SAVE,
     +                 SINPHI,SINTHE,TAVE,TENT,W,XA,YA
      INTEGER I,IA,IB,IC,ISTOPZ,ISUZY,J,JC,K,KOV2,L,LAVE,LSUZY,LTWOQP,M,
     +        MAVE,MP1,MSUZY,N
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION COSMPH(LMAX2P),FACTOR(50),PLM(LMXP),
     +                 SINMPH(LMAX2P)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS,DSIGN,DSQRT
C     ..
      FACTOR(1) = 1.0D00
      DO 10 I = 2,50
        XA = I - 1
        FACTOR(I) = XA*FACTOR(I-1)
   10 CONTINUE
      PSQ = X*X + Y*Y
      RSSQ = PSQ + Z*Z
      IF (RSSQ-1.0D-10) 20,20,30
   20 SINTHE = 0.0D00
      COSTHE = 1.0D00
      COSPHI = 1.0D00
      SINPHI = 0.0D00
      R = 0.0D00
      GO TO 60

   30 IF (PSQ-1.0D-10) 40,40,50
   40 R = DSQRT(RSSQ)
      SINTHE = 0.0D00
      COSTHE = DSIGN(1.0D00,Z)
      SINPHI = 0.0D00
      COSPHI = 1.0D00
      GO TO 60

   50 R = DSQRT(RSSQ)
      P = DSQRT(PSQ)
      SINTHE = P/R
      COSTHE = Z/R
      SINPHI = Y/P
      COSPHI = X/P
C
   60 XA = DABS(COSTHE)
      YA = DABS(SINTHE)
c      write(6,*) costhe,sinthe
C
      IF (XA-1.0D-08) 70,70,150
   70 L = 0
      J = 0
      TAVE = 1.0D00
      LSUZY = 1
   80 M = 0
      MSUZY = 1
   90 J = J + 1
      ISUZY = LSUZY*MSUZY
      IF (ISUZY.GT.0) GO TO 100
      PLM(J) = 0.0D00
      GO TO 110

  100 K = L + M
      KOV2 = K/2
      IA = K + 1
      IB = KOV2 + 1
      JC = KOV2 - M
      IC = JC + 1
      PLM(J) = (((-1)**JC)*FACTOR(IA))/ (TAVE*FACTOR(IB)*FACTOR(IC))
  110 IF (M-L) 120,130,130
  120 M = M + 1
      MSUZY = -MSUZY
      GO TO 90

  130 IF (L-LMAX2) 140,370,370
  140 L = L + 1
      LSUZY = -LSUZY
      TAVE = 2.0D00*TAVE
      GO TO 80
C
  150 IF (XA-0.99999999D00) 250,160,160
  160 PLM(1) = 1.0D00
      PLM(2) = COSTHE
      L = 2
      J = 2
  170 J = J + L
      A = L
      LTWOQP = L + L
      B = LTWOQP - 1
      C = L - 1
      K = J - L
      M = J - LTWOQP + 1
      PLM(J) = (B*COSTHE*PLM(K)-C*PLM(M))/A
      IF (L-LMAX2) 180,190,190
  180 L = L + 1
      GO TO 170
C
  190 L = 1
      LAVE = 1
  200 M = 1
      LAVE = LAVE + L
  210 J = LAVE + M
      PLM(J) = 0.0D00
      IF (M-L) 220,230,230
  220 M = M + 1
      GO TO 210

  230 IF (L-LMAX2) 240,370,370
  240 L = L + 1
      GO TO 200
C
  250 TENT = (2.0D00*COSTHE)/YA
      PLM(1) = 1.0D00
      PLM(2) = COSTHE
      PLM(3) = YA
      PLM(5) = 3.0D00*YA*COSTHE
      L = 2
      J = 2
  260 J = J + L
      A = L
      LTWOQP = L + L
      B = LTWOQP - 1
      C = L - 1
      K = J - L
      M = J - LTWOQP + 1
      PLM(J) = (B*COSTHE*PLM(K)-C*PLM(M))/A
      IF (L-LMAX2) 270,280,280
  270 L = L + 1
      GO TO 260

  280 L = 3
      J = 5
  290 J = J + L
      A = L - 1
      LTWOQP = L + L
      B = LTWOQP - 1
      C = L
      K = J - L
      M = J - LTWOQP + 1
      PLM(J) = (B*COSTHE*PLM(K)-C*PLM(M))/A
      IF (L-LMAX2) 300,310,310
  300 L = L + 1
      GO TO 290

  310 L = 2
      LAVE = 3
  320 LAVE = LAVE + L
      M = 1
  330 J = LAVE + M
      K = J - 1
      N = K - 1
      EAT = M
      A = TENT*EAT
      B = (M+L)* (L-M+1)
      PLM(J) = A*PLM(K) - B*PLM(N)
      IF (M+1-L) 340,350,350
  340 M = M + 1
      GO TO 330

  350 IF (L-LMAX2) 360,370,370
  360 L = L + 1
      GO TO 320
C
  370 SINMPH(1) = 0.0D00
      COSMPH(1) = 1.0D00
      ISTOPZ = LMAX2 + 1
      DO 380 I = 2,ISTOPZ
        J = I - 1
        SINMPH(I) = SINPHI*COSMPH(J) + COSPHI*SINMPH(J)
        COSMPH(I) = COSPHI*COSMPH(J) - SINPHI*SINMPH(J)
  380 CONTINUE
C
      L = 0
  390 M = 0
      LAVE = L* (L+1) + 1
      MAVE = ((L* (L+1))/2) + 1
      SAVE = 2*L + 1
  400 IF (M.NE.0) GO TO 410
      ARG = SAVE/12.5663706144D00
      W = DSQRT(ARG)
      YREALY(LAVE) = W*PLM(MAVE)
      GO TO 420

  410 IA = L - M + 1
      IB = L + M + 1
      ARG = (SAVE*FACTOR(IA))/ (6.28318530718D00*FACTOR(IB))
      MP1 = M + 1
      W = DSQRT(ARG)
      I = LAVE + M
      J = MAVE + M
      YREALY(I) = W*PLM(J)*COSMPH(MP1)
      I = LAVE - M
      YREALY(I) = W*PLM(J)*SINMPH(MP1)
  420 IF (M.GE.L) GO TO 430
      M = M + 1
      GO TO 400

  430 IF (L.GE.LMAX2) GO TO 440
      L = L + 1
      GO TO 390
C
  440 RETURN
C
C
C
      END
c ************************************************************************
c EOF gaunt.f
c ************************************************************************
