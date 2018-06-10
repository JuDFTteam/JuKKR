subroutine sumupint(sum, vg, g, wg, vf, f, wf, n)
!   ********************************************************************
!   *                                                                  *
!   ********************************************************************
  implicit none


! Dummy arguments
  integer :: n
  complex *16 :: sum
  real *8 :: vf, vg
  complex *16 :: f(2, 2), g(2, 2)
  real *8 :: wf(2, 2), wg(2, 2)

! Local variables
  integer :: i, j

  sum = 0.0d0
  do j = 1, n
    do i = 1, n
      sum = sum + vg*g(i, j)*wg(i, j) + vf*f(i, j)*wf(i, j)
    end do
  end do

end subroutine
