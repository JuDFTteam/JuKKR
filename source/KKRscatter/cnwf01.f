C ************************************************************************
      SUBROUTINE CNWF01(CP,W,N,C,SUM)
c   *****************************************************************
c   *this subroutine computes the coefficients in the chebychev
c   *expansion of 0f1(;c;z). p.80 in algorithms for the computation
c   *of mathematical functions, y.l.luke, academic press,london 1977
c   *
c   *
c   *description of variables
c   *
c   *cp    -input  - parameter c in 0f1(;c;z)
c   *
c   *w     -input  - this is a preselected scale factor such that
c   *                0.le.(z/w).le.1
c   *
c   *n     -input  - two less than the number of coefficients to be
c   *                generated
c   *
c   *c     -output - a vector containing the n+2 chebychev coefficients
c   *
c   *sum   -output - the sum of the coefficients
c   *
c   *all other variables are for internal use
c**********************************************************************
C     .. Scalar Arguments ..
      DOUBLE PRECISION CP,SUM,W
      INTEGER N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION C(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A1,A2,A3,C1,DIVFAC,FOUR,ONE,P,RHO,START,TWO,X1,
     +                 Z1,ZERO
      INTEGER I,K,L,N1,N2,NCOUNT
C     ..
C     .. Save statement ..
      SAVE ZERO,ONE,TWO,FOUR,START
C     ..
C     .. Data statements ..
      DATA ZERO,ONE,TWO,FOUR,START/0.D0,1.D0,2.D0,4.D0,1.D-20/
C     ..

      N1 = N + 1
      N2 = N + 2
c
c     ------- start computing coefficients by means of    --------------
c     ------- backward recurrence scheme                  --------------
      A3 = ZERO
      A2 = ZERO
      A1 = START
      Z1 = FOUR/W
      NCOUNT = N2
      C(NCOUNT) = START
      X1 = N2
      C1 = ONE - CP
      DO 10 K = 1,N1
        DIVFAC = ONE/X1
        X1 = X1 - ONE
        NCOUNT = NCOUNT - 1
        C(NCOUNT) = X1* ((DIVFAC+Z1* (X1-C1))*A1+
     +              (ONE/X1+Z1* (X1+C1+ONE))*A2-DIVFAC*A3)
        A3 = A2
        A2 = A1
        A1 = C(NCOUNT)
   10 CONTINUE
      C(1) = C(1)/TWO
c
c     --------- compute scale factor                      --------------
      RHO = C(1)
      SUM = RHO
      P = ONE
      DO 20 I = 2,N2
        RHO = RHO - P*C(I)
        SUM = SUM + C(I)
        P = -P
   20 CONTINUE
c
c     ---------- scale coefficients                         ------------
      DO 30 L = 1,N2
        C(L) = C(L)/RHO
   30 CONTINUE
      SUM = SUM/RHO
      END
