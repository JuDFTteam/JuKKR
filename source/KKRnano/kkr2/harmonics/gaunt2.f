c>    sets up values needed for gaunt.

c>       m. weinert  january 1982
c>
c>    changed for calculating with real spherical harmonics
c>                                          b.drittler  july 1987
c>
c>    W(N)        integration weights on 4*LMAX points in the intervall
c>                (-1,0) (from routine GRULE)
c>
c>    YR(N,L,M)   spherical harmonics on 4*LMAX points to angular
c>                momentum indices (l,m) scaled with a factor
c>                of RF=(4*pi)**(1/3)
c>
c>    LMAX        l-cutoff, added after inc.p remove
C***********************************************************************
      SUBROUTINE GAUNT2(W,YR,LMAX)
C ************************************************************************
c     sets up values needed for gaunt.

c        m. weinert  january 1982
c
c     changed for calculating with real spherical harmonics
c                                           b.drittler  july 1987
c
c     W(N)        integration weights on 4*LMAX points in the intervall
c                 (-1,0) (from routine GRULE)
c
c     YR(N,L,M)   spherical harmonics on 4*LMAX points to angular
c                 momentum indices (l,m) scaled with a factor
c                 of RF=(4*pi)**(1/3)
c
c     LMAX        l-cutoff, added after inc.p remove
c
c-----------------------------------------------------------------------
C     .. Parameters ..
C      include 'inc.p'
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: LMAX
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,CD,CTH,FAC,FPI,RF,STH,T
      INTEGER K,L,LOMAX,M
C     ..
C     ..
C     .. External Subroutines ..
      EXTERNAL GRULE
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,SQRT
C     ..
C     .. Save statement ..
C     E.Rabel: is it really needed?
      SAVE
C     ..
C     .. Array Arguments ..
C      DOUBLE PRECISION W(*),YR(N,0:N,0:N)
C      N = 4 * LMAX
      DOUBLE PRECISION W(*),YR(4*LMAX,0:4*LMAX,0:4*LMAX)
C     ..

C     .. Local Arrays ..
      DOUBLE PRECISION P(0:4*LMAX+1,0:4*LMAX),X(4*LMAX)

      INTEGER N
      N=4*LMAX

      FPI = 16.D0*ATAN(1.0D0)
      RF = FPI** (1.0D0/3.0D0)
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
        STH = SQRT(1.0D0-CTH*CTH)
        FAC = 1.0D0
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
            YR(K,L,M) = A*SQRT(2.0D0*CD)*P(L,M)
   30     CONTINUE
   40   CONTINUE
   50 CONTINUE
      END
