!-------------------------------------------------------------------------------
!> Summary: Computes the coefficients in the Chebychev expansion
!> Author: 
!> Deprecated: False 
!-------------------------------------------------------------------------------
!> see: p.80 in algorithms for the computation
!> of mathematical functions, y.l.luke, academic press,london 1977
!>
!> Computes the coefficients in the chebychev expansion of 0f1(;c;z)
!>
!> Description of variables:
!> cp  input  - parameter c in 0f1(;c;z)
!> w   input  - this is a preselected scale factor such that \(0 \leq z/w \leq 1\)
!> n   input  - two less than the number of coefficients to be generated
!> c   output - a vector containing the n+2 Chebychev coefficients
!> sum output - the sum of the coefficients
!> all other variables are for internal use
!-------------------------------------------------------------------------------
      MODULE MOD_CNWF011
      CONTAINS
!-------------------------------------------------------------------------------
!> Summary: Computes the coefficients in the Chebychev expansion of 0f1(;c;z)
!> Author: 
!> Date: 
!> Category: KKRimp, special-functions, numerical-tools
!> Deprecated: False 
!-------------------------------------------------------------------------------
      SUBROUTINE CNWF011(CP,W,N,C,SUM)

c**********************************************************************
C     .. Scalar Arguments ..
      REAL*8 CP             !! Parameter c in 0f1(;c;z)
      REAL*8 SUM            !! Sum of the coefficients
      REAL*8 W              !! Preselected scale factor such that \(0 \leq z/w \leq 1\)
      INTEGER N             !! Two less than the number of coefficients to be generated
C     ..
C     .. Array Arguments ..
      REAL*8 C(*)           !! Vector containing the n+2 Chebychev coefficients
C     ..
C     .. Local Scalars ..
      REAL*8 A1,A2,A3,C1,DIVFAC,FOUR,ONE,P,RHO,START,TWO,X1,
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
      END SUBROUTINE
      END MODULE MOD_CNWF011
