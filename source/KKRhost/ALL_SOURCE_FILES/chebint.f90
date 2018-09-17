module mod_chebint
  use :: mod_datatypes, only: dp
  private :: dp

contains

  subroutine chebint(cslc1, csrc1, slc1sum, c1, n)
    ! ---------------------------------------------------------------------
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
    ! ---------------------------------------------------------------------
    implicit none
    ! .. Scalar Arguments ..
    integer, intent (in) :: n
    ! .. Array Arguments ..
    real (kind=dp) :: cslc1(0:n, 0:n), csrc1(0:n, 0:n), slc1sum(0:n)

    ! .. Local Scalars ..
    real (kind=dp) :: pi
    integer :: j, k
    ! .. Local Arrays ..
    real (kind=dp) :: c(0:n, 0:n), c1(0:n, 0:n), s1(0:n, 0:n), s2(0:n, 0:n), sl(0:n, 0:n), slc1(0:n, 0:n), sr(0:n, 0:n), src1(0:n, 0:n)
    ! ..
    pi = 4.e0_dp*atan(1.e0_dp)
    ! ---------------------------------------------------------------------
    ! determine the discrete cosine transform matrix from the zeros of the
    ! Chebyshev polynomials
    do j = 0, n
      do k = 0, n
        c(k, j) = cos(((2*k+1)*j*pi)/(2*(n+1)))
      end do
    end do
    ! ---------------------------------------------------------------------
    ! determine the inverse of the discrete cosine transform matrix from
    ! the transpose of the discrete cosine transform matrix
    do j = 0, n
      do k = 0, n
        c1(k, j) = c(j, k)*2.e0_dp/(n+1)
      end do
      c1(0, j) = c1(0, j)*0.5e0_dp
    end do
    ! ---------------------------------------------------------------------
    ! next to statements can be used to check the products CT*C and C1*C
    call dgemm('T', 'N', n+1, n+1, n+1, 1.e0_dp, c, n+1, c, n+1, 0.e0_dp, sr, n+1)
    call dgemm('N', 'N', n+1, n+1, n+1, 1.e0_dp, c1, n+1, c, n+1, 0.e0_dp, sr, n+1)
    ! ---------------------------------------------------------------------
    ! preparation of the left and right
    ! spectral integration matrices SL and SR
    do j = 0, n
      do k = 0, n
        s1(k, j) = 0.0e0_dp
        s2(k, j) = 0.0e0_dp
      end do
    end do
    do j = 0, n
      s1(0, j) = (-1.e0_dp)**(j+1)
      s1(j, j) = 1.e0_dp
    end do
    do j = 2, n - 1
      s2(j, j-1) = 0.5e0_dp/j
      s2(j, j+1) = -0.5e0_dp/j
    end do
    s2(n, n-1) = 0.5e0_dp/n
    s2(1, 0) = 1.e0_dp
    s2(1, 2) = -0.5e0_dp
    call dgemm('N', 'N', n+1, n+1, n+1, 1.e0_dp, s1, n+1, s2, n+1, 0.e0_dp, sl, n+1)
    do j = 0, n
      do k = 0, n
        s1(k, j) = 0.0e0_dp
      end do
    end do
    do j = 0, n
      s1(j, j) = -1.e0_dp
      s1(0, j) = 1.e0_dp
    end do
    call dgemm('N', 'N', n+1, n+1, n+1, 1.e0_dp, s1, n+1, s2, n+1, 0.e0_dp, sr, n+1)
    ! ---------------------------------------------------------------------
    ! determination of the products C*SL*C1 and C*SR*C1
    call dgemm('N', 'N', n+1, n+1, n+1, 1.e0_dp, sl, n+1, c1, n+1, 0.e0_dp, slc1, n+1)
    call dgemm('N', 'N', n+1, n+1, n+1, 1.e0_dp, c, n+1, slc1, n+1, 0.e0_dp, cslc1, n+1)
    call dgemm('N', 'N', n+1, n+1, n+1, 1.e0_dp, sr, n+1, c1, n+1, 0.e0_dp, src1, n+1)
    call dgemm('N', 'N', n+1, n+1, n+1, 1.e0_dp, c, n+1, src1, n+1, 0.e0_dp, csrc1, n+1)
    ! ---------------------------------------------------------------------
    do k = 0, n
      slc1sum(k) = 0.0e0_dp
      do j = 0, n
        slc1sum(k) = slc1sum(k) + slc1(j, k)
      end do
    end do
    return
  end subroutine chebint

end module mod_chebint
