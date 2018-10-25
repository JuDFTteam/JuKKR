C***********************************************************************
      SUBROUTINE GAUNT2(W,YR,N)
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
      IMPLICIT NONE
C     .. Arguments
      INTEGER N
      DOUBLE PRECISION W(*),YR(N,0:N,0:N)
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
      FPI = 16D0*ATAN(1D0)
      RF = FPI**(1D0/3D0)
      LOMAX = N
c
c--->    obtain gauss-legendre points and weights
c
      CALL GRULE(2*N,X,W)
c
c--->    generate associated legendre functions for m.ge.0
c
      DO K = 1,N
         CTH = X(K)
         STH = SQRT(1.D0-CTH*CTH)
         FAC = 1.D0
c
c--->    loop over m values
c
         DO M = 0,LOMAX
            FAC = - DBLE(2*M-1)*FAC
            P(M,M) = FAC
            P(M+1,M) = DBLE(2*M+1)*CTH*FAC
c
c--->    recurse upward in l
c
            DO L = M + 2,LOMAX
               P(L,M) = ( DBLE(2*L-1)*CTH*P(L-1,M) 
     &                  - DBLE(L+M-1)    *P(L-2,M) ) / DBLE(L-M)
            END DO
C
            FAC = FAC*STH
         END DO
c
c--->    multiply in the normalization factors
c
         DO L = 0,LOMAX
            A = RF*SQRT((2*L+1)/FPI)
            CD = 1.D0
            YR(K,L,0) = A*P(L,0)
C
            DO M = 1,L
               T = DBLE( (L+1-M)* (L+M))
               CD = CD/T
               YR(K,L,M) = A*SQRT(2.D0*CD)*P(L,M)
            END DO
         END DO
      END DO
      END
