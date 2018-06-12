subroutine nrmliz(n, r, rn)
  use :: mod_datatypes, only: dp
  ! -  normalizes a vector
  ! ----------------------------------------------------------------------
  ! i Inputs
  ! i   n     :number of vectors
  ! i   r     :vector
  ! o Outputs:
  ! o   rn    :normalized vector
  ! ----------------------------------------------------------------------
  implicit none
  real (kind=dp), parameter :: eps = 1e-14_dp
  ! Passed parameters:
  integer :: n
  real (kind=dp) :: r(3, *), rn(3, *)
  ! Local parameters
  integer :: i
  real (kind=dp) :: d, d2
  ! External calls
  external :: dcopy, dscal

  call dcopy(3*n, r, 1, rn, 1)
  do i = 1, n
    d2 = r(1, i)*r(1, i) + r(2, i)*r(2, i) + r(3, i)*r(3, i)
    d = sqrt(d2)
    if (abs(d)<eps) call dscal(3, 1.e0_dp/d, rn(1,i), 1)
  end do

end subroutine nrmliz
