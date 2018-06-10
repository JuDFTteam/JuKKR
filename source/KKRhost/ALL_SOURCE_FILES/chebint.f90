SUBROUTINE chebint(cslc1,csrc1,slc1sum,c1,n)
!---------------------------------------------------------------------
! this subroutine calculates the matrices for the Chebyshev integration
! as defined on page 141 and 142 of the article:
! Integral Equation Method for the Continuous Spectrum Radial
! Schroedinger Equation by R. A. Gonzales et al
! in Journal of computational physics 134, 134-149 (1997)

! the matrix C is the discrete cosine transform matrix
! the matrix C1 is the inverse of C
! the matrix SL is the left spectral integration matrix
! the matrix SR is the right spectral integration matrix
! the matrix CSLC1 is the product of C, SL and C1
! the matrix CSRC1 is the product of C, SR and C1
!---------------------------------------------------------------------
IMPLICIT NONE
!     .. Local Scalars ..
DOUBLE PRECISION :: pi
INTEGER :: j,k
!     ..
!     .. Local Arrays ..
DOUBLE PRECISION :: c(0:n,0:n),c1(0:n,0:n),s1(0:n,0:n),s2(0:n,0:n),  &
    sl(0:n,0:n),slc1(0:n,0:n),sr(0:n,0:n), src1(0:n,0:n)
!     ..
!     .. External Subroutines ..
EXTERNAL dgemm
!     ..
!     .. Intrinsic Functions ..
INTRINSIC ATAN,COS
!     ..
!     .. Array Arguments ..
DOUBLE PRECISION :: cslc1(0:n,0:n),csrc1(0:n,0:n),slc1sum(0:n)
!     ..
!     .. Scalar Arguments ..
INTEGER, intent(in) :: n
!     ..
pi = 4.d0*ATAN(1.d0)
!---------------------------------------------------------------------
! determine the discrete cosine transform matrix from the zeros of the
! Chebyshev polynomials
DO j = 0,n
  DO k = 0,n
    c(k,j) = COS(((2*k+1)*j*pi)/ (2* (n+1)))
  END DO
END DO
!---------------------------------------------------------------------
! determine the inverse of the discrete cosine transform matrix from
! the transpose of the discrete cosine transform matrix
DO j = 0,n
  DO k = 0,n
    c1(k,j) = c(j,k)*2.d0/ (n+1)
  END DO
  c1(0,j) = c1(0,j)*0.5D0
END DO
!---------------------------------------------------------------------
! next to statements can be used to check the products CT*C and C1*C
CALL dgemm('T','N',n+1,n+1,n+1,1.d0,c,n+1,c,n+1,0.d0,sr,n+1)
CALL dgemm('N','N',n+1,n+1,n+1,1.d0,c1,n+1,c,n+1,0.d0,sr,n+1)
!---------------------------------------------------------------------
! preparation of the left and right
! spectral integration matrices SL and SR
DO j = 0,n
  DO k = 0,n
    s1(k,j) = 0.0D0
    s2(k,j) = 0.0D0
  END DO
END DO
DO j = 0,n
  s1(0,j) = (-1.d0)** (j+1)
  s1(j,j) = 1.d0
END DO
DO j = 2,n - 1
  s2(j,j-1) = 0.5D0/j
  s2(j,j+1) = -0.5D0/j
END DO
s2(n,n-1) = 0.5D0/n
s2(1,0) = 1.d0
s2(1,2) = -0.5D0
CALL dgemm('N','N',n+1,n+1,n+1,1.d0,s1,n+1,s2,n+1,0.d0,sl,n+1)
DO j = 0,n
  DO k = 0,n
    s1(k,j) = 0.0D0
  END DO
END DO
DO j = 0,n
  s1(j,j) = -1.d0
  s1(0,j) = 1.d0
END DO
CALL dgemm('N','N',n+1,n+1,n+1,1.d0,s1,n+1,s2,n+1,0.d0,sr,n+1)
!---------------------------------------------------------------------
! determination of the products C*SL*C1 and C*SR*C1
CALL dgemm('N','N',n+1,n+1,n+1,1.d0,sl,n+1,c1,n+1,0.d0,slc1,n+1)
CALL dgemm('N','N',n+1,n+1,n+1,1.d0,c,n+1,slc1,n+1,0.d0,cslc1,n+1)
CALL dgemm('N','N',n+1,n+1,n+1,1.d0,sr,n+1,c1,n+1,0.d0,src1,n+1)
CALL dgemm('N','N',n+1,n+1,n+1,1.d0,c,n+1,src1,n+1,0.d0,csrc1,n+1)
!---------------------------------------------------------------------
DO k = 0,n
  slc1sum(k) = 0.0D0
  DO j = 0,n
    slc1sum(k) = slc1sum(k) + slc1(j,k)
  END DO
END DO
RETURN
END SUBROUTINE
