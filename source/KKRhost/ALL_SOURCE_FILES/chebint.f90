    Subroutine chebint(cslc1, csrc1, slc1sum, c1, n)
      Use mod_datatypes, Only: dp
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
      Implicit None
!     .. Scalar Arguments ..
      Integer, Intent (In) :: n
!     .. Array Arguments ..
      Real (Kind=dp) :: cslc1(0:n, 0:n), csrc1(0:n, 0:n), slc1sum(0:n)

!     .. Local Scalars ..
      Real (Kind=dp) :: pi
      Integer :: j, k
!     .. Local Arrays ..
      Real (Kind=dp) :: c(0:n, 0:n), c1(0:n, 0:n), s1(0:n, 0:n), s2(0:n, 0:n), &
        sl(0:n, 0:n), slc1(0:n, 0:n), sr(0:n, 0:n), src1(0:n, 0:n)
!     .. External Subroutines ..
      External :: dgemm
!     .. Intrinsic Functions ..
      Intrinsic :: atan, cos
!     ..
      pi = 4.E0_dp*atan(1.E0_dp)
!---------------------------------------------------------------------
! determine the discrete cosine transform matrix from the zeros of the
! Chebyshev polynomials
      Do j = 0, n
        Do k = 0, n
          c(k, j) = cos(((2*k+1)*j*pi)/(2*(n+1)))
        End Do
      End Do
!---------------------------------------------------------------------
! determine the inverse of the discrete cosine transform matrix from
! the transpose of the discrete cosine transform matrix
      Do j = 0, n
        Do k = 0, n
          c1(k, j) = c(j, k)*2.E0_dp/(n+1)
        End Do
        c1(0, j) = c1(0, j)*0.5E0_dp
      End Do
!---------------------------------------------------------------------
! next to statements can be used to check the products CT*C and C1*C
      Call dgemm('T', 'N', n+1, n+1, n+1, 1.E0_dp, c, n+1, c, n+1, 0.E0_dp, &
        sr, n+1)
      Call dgemm('N', 'N', n+1, n+1, n+1, 1.E0_dp, c1, n+1, c, n+1, 0.E0_dp, &
        sr, n+1)
!---------------------------------------------------------------------
! preparation of the left and right
! spectral integration matrices SL and SR
      Do j = 0, n
        Do k = 0, n
          s1(k, j) = 0.0E0_dp
          s2(k, j) = 0.0E0_dp
        End Do
      End Do
      Do j = 0, n
        s1(0, j) = (-1.E0_dp)**(j+1)
        s1(j, j) = 1.E0_dp
      End Do
      Do j = 2, n - 1
        s2(j, j-1) = 0.5E0_dp/j
        s2(j, j+1) = -0.5E0_dp/j
      End Do
      s2(n, n-1) = 0.5E0_dp/n
      s2(1, 0) = 1.E0_dp
      s2(1, 2) = -0.5E0_dp
      Call dgemm('N', 'N', n+1, n+1, n+1, 1.E0_dp, s1, n+1, s2, n+1, 0.E0_dp, &
        sl, n+1)
      Do j = 0, n
        Do k = 0, n
          s1(k, j) = 0.0E0_dp
        End Do
      End Do
      Do j = 0, n
        s1(j, j) = -1.E0_dp
        s1(0, j) = 1.E0_dp
      End Do
      Call dgemm('N', 'N', n+1, n+1, n+1, 1.E0_dp, s1, n+1, s2, n+1, 0.E0_dp, &
        sr, n+1)
!---------------------------------------------------------------------
! determination of the products C*SL*C1 and C*SR*C1
      Call dgemm('N', 'N', n+1, n+1, n+1, 1.E0_dp, sl, n+1, c1, n+1, 0.E0_dp, &
        slc1, n+1)
      Call dgemm('N', 'N', n+1, n+1, n+1, 1.E0_dp, c, n+1, slc1, n+1, 0.E0_dp, &
        cslc1, n+1)
      Call dgemm('N', 'N', n+1, n+1, n+1, 1.E0_dp, sr, n+1, c1, n+1, 0.E0_dp, &
        src1, n+1)
      Call dgemm('N', 'N', n+1, n+1, n+1, 1.E0_dp, c, n+1, src1, n+1, 0.E0_dp, &
        csrc1, n+1)
!---------------------------------------------------------------------
      Do k = 0, n
        slc1sum(k) = 0.0E0_dp
        Do j = 0, n
          slc1sum(k) = slc1sum(k) + slc1(j, k)
        End Do
      End Do
      Return
    End Subroutine
