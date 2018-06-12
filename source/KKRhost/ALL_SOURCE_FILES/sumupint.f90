subroutine sumupint(sum, vg, g, wg, vf, f, wf, n)
  use :: mod_datatypes, only: dp
  ! ********************************************************************
  ! *                                                                  *
  ! ********************************************************************
  implicit none


  ! Dummy arguments
  integer :: n
  complex (kind=dp) :: sum
  real (kind=dp) :: vf, vg
  complex (kind=dp) :: f(2, 2), g(2, 2)
  real (kind=dp) :: wf(2, 2), wg(2, 2)

  ! Local variables
  integer :: i, j

  sum = 0.0e0_dp
  do j = 1, n
    do i = 1, n
      sum = sum + vg*g(i, j)*wg(i, j) + vf*f(i, j)*wf(i, j)
    end do
  end do

end subroutine sumupint
