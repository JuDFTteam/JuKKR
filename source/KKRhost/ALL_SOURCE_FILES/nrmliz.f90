!------------------------------------------------------------------------------------
!> Summary: Normalizes vectors
!> Author: 
!> Normalizes vectors
!------------------------------------------------------------------------------------
module mod_nrmliz
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Normalizes vectors
  !> Author: 
  !> Category: numerical-tools, KKRhost
  !> Deprecated: False 
  !> Normalizes vectors
  !-------------------------------------------------------------------------------
  subroutine nrmliz(n, r, rn)

    implicit none

    real (kind=dp), parameter :: eps = 1e-14_dp
    ! Passed parameters:
    integer, intent(in) :: n !! number of vectors
    real (kind=dp), dimension(3,*), intent(in) :: r   !! Input vector
    real (kind=dp), dimension(3,*), intent(out) :: rn !! Normalized vector
    ! Local parameters
    integer :: i
    real (kind=dp) :: d, d2

    call dcopy(3*n, r, 1, rn, 1)
    do i = 1, n
      d2 = r(1, i)*r(1, i) + r(2, i)*r(2, i) + r(3, i)*r(3, i)
      d = sqrt(d2)
      if (abs(d)<eps) call dscal(3, 1.e0_dp/d, rn(1,i), 1)
    end do

  end subroutine nrmliz

end module mod_nrmliz
