subroutine cmatmul(n, m, a, b, c)
!   ********************************************************************
!   *                                                                  *
!   *   perform  the matrix-matrix operation           C = A * B       *
!   *                                                                  *
!   *   A,B,C   complex  SQUARE  N x N - matrices                      *
!   *   N       dimension of A, B and C                                *
!   *   M       array size of A, B, C with M >= N                      *
!   *                                                                  *
!   ********************************************************************

  implicit none

! PARAMETER definitions
  double complex :: c0
  parameter (c0=(0.0d0,0.0d0))

! Dummy arguments
  integer :: m, n
  double complex :: a(m, m), b(m, m), c(m, m)

! Local variables
  double complex :: blj
  integer :: i, j, l

  do j = 1, n
    do i = 1, n
      c(i, j) = c0
    end do
  end do

  do j = 1, n
    do l = 1, n
      blj = b(l, j)
      if (blj/=c0) then
        do i = 1, n
          c(i, j) = c(i, j) + a(i, l)*blj
        end do
      end if
    end do
  end do

end subroutine
