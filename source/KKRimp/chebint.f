      MODULE MOD_CHEBINT
!-------------------------------------------------------------------------------
!> Summary: Calculates the matrices for the Chebyshev integration
!> Author:
!> Category: KKRimp, radial-grid
!>           
!-------------------------------------------------------------------------------
        CONTAINS
!-------------------------------------------------------------------------------
!> Summary: Calculates the matrices for the Chebyshev integration
!> Author:
!> Category: KKRimp, radial-grid
!>           
!-------------------------------------------------------------------------------
      SUBROUTINE CHEBINT(CSLC1,CSRC1,SLC1SUM,C1,N)
C---------------------------------------------------------------------
C this subroutine calculates the matrices for the Chebyshev integration
C as defined on page 141 and 142 of the article:
C Integral Equation Method for the Continuous Spectrum Radial
C Schroedinger Equation by R. A. Gonzales et al
C in Journal of computational physics 134, 134-149 (1997) 
C
C the matrix C is the discrete cosine transform matrix
C the matrix C1 is the inverse of C 
C the matrix SL is the left spectral integration matrix
C the matrix SR is the right spectral integration matrix
C the matrix CSLC1 is the product of C, SL and C1
C the matrix CSRC1 is the product of C, SR and C1
C---------------------------------------------------------------------
      IMPLICIT NONE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CSLC1(0:N,0:N),CSRC1(0:N,0:N),SLC1SUM(0:N)
C     ..
C     .. Scalar Arguments ..
      INTEGER N
C     .. Local Scalars ..
      DOUBLE PRECISION PI
      INTEGER J,K
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION C(0:N,0:N),C1(0:N,0:N),S1(0:N,0:N),S2(0:N,0:N),
     +                 SL(0:N,0:N),SLC1(0:N,0:N),SR(0:N,0:N),
     +                 SRC1(0:N,0:N)
C     ..
C     .. External Subroutines ..
      EXTERNAL DGEMM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,COS
C     ..
      PI = 4.D0*ATAN(1.D0)
C---------------------------------------------------------------------
C determine the discrete cosine transform matrix from the zeros of the
C Chebyshev polynomials
      DO J = 0,N
        DO K = 0,N
          C(K,J) = COS(((2*K+1)*J*PI)/ (2* (N+1)))
        END DO
      END DO
C---------------------------------------------------------------------
C determine the inverse of the discrete cosine transform matrix from
C the transpose of the discrete cosine transform matrix
      DO J = 0,N
        DO K = 0,N
          C1(K,J) = C(J,K)*2.D0/ (N+1)
        END DO
        C1(0,J) = C1(0,J)*0.5D0
      END DO
C---------------------------------------------------------------------
C next to statements can be used to check the products CT*C and C1*C
      CALL DGEMM('T','N',N+1,N+1,N+1,1.D0,C,N+1,C,N+1,0.D0,SR,N+1)
      CALL DGEMM('N','N',N+1,N+1,N+1,1.D0,C1,N+1,C,N+1,0.D0,SR,N+1)
C---------------------------------------------------------------------
C preparation of the left and right
C spectral integration matrices SL and SR
      DO J = 0,N
        DO K = 0,N
          S1(K,J) = 0.0d0
          S2(K,J) = 0.0d0
        END DO
      END DO
      DO J = 0,N
        S1(0,J) = (-1.D0)** (J+1)
        S1(J,J) = 1.D0
      END DO
      DO J = 2,N - 1
        S2(J,J-1) = 0.5D0/J
        S2(J,J+1) = -0.5D0/J
      END DO
      S2(N,N-1) = 0.5D0/N
      S2(1,0) = 1.D0
      S2(1,2) = -0.5D0
      CALL DGEMM('N','N',N+1,N+1,N+1,1.D0,S1,N+1,S2,N+1,0.D0,SL,N+1)
      DO J = 0,N
        DO K = 0,N
          S1(K,J) = 0.0d0
        END DO
      END DO
      DO J = 0,N
        S1(J,J) = -1.D0
        S1(0,J) = 1.D0
      END DO
      CALL DGEMM('N','N',N+1,N+1,N+1,1.D0,S1,N+1,S2,N+1,0.D0,SR,N+1)
C---------------------------------------------------------------------
C determination of the products C*SL*C1 and C*SR*C1
      CALL DGEMM('N','N',N+1,N+1,N+1,1.D0,SL,N+1,C1,N+1,0.D0,SLC1,N+1)
      CALL DGEMM('N','N',N+1,N+1,N+1,1.D0,C,N+1,SLC1,N+1,0.D0,CSLC1,N+1)
      CALL DGEMM('N','N',N+1,N+1,N+1,1.D0,SR,N+1,C1,N+1,0.D0,SRC1,N+1)
      CALL DGEMM('N','N',N+1,N+1,N+1,1.D0,C,N+1,SRC1,N+1,0.D0,CSRC1,N+1)
C---------------------------------------------------------------------
      DO K = 0,N
        SLC1SUM(K) = 0.0D0
        DO J = 0,N
          SLC1SUM(K) = SLC1SUM(K) + SLC1(J,K)
        END DO
      END DO
      RETURN
      END SUBROUTINE
      END MODULE MOD_CHEBINT
