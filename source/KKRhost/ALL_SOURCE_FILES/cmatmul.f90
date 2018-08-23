module mod_cmatmul

contains

subroutine cmatmul(n, m, a, b, c)
  use :: mod_datatypes, only: dp
  ! ********************************************************************
  ! *                                                                  *
  ! *   perform  the matrix-matrix operation           C = A * B       *
  ! *                                                                  *
  ! *   A,B,C   complex  SQUARE  N x N - matrices                      *
  ! *   N       dimension of A, B and C                                *
  ! *   M       array size of A, B, C with M >= N                      *
  ! *                                                                  *
  ! ********************************************************************

  implicit none

  ! PARAMETER definitions
  real (kind=dp), parameter :: eps=1.0D-12
  complex (kind=dp) :: c0
  parameter (c0=(0.0e0_dp,0.0e0_dp))

  ! Dummy arguments
  integer :: m, n
  complex (kind=dp) :: a(m, m), b(m, m), c(m, m)

  ! Local variables
  complex (kind=dp) :: blj
  integer :: i, j, l

  do j = 1, n
    do i = 1, n
      c(i, j) = c0
    end do
  end do

  do j = 1, n
    do l = 1, n
      blj = b(l, j)
      if (abs(blj)>eps) then
        do i = 1, n
          c(i, j) = c(i, j) + a(i, l)*blj
        end do
      end if
    end do
  end do

end subroutine cmatmul

end module mod_cmatmul
