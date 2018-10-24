!------------------------------------------------------------------------------------
!> Summary: 
!> Author: 
!> 
!------------------------------------------------------------------------------------
module mod_sumupint
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: 
  !> Author: 
  !> Category: dirac, undefined, KKRhost
  !> Deprecated: False 
  !> 
  !-------------------------------------------------------------------------------
  subroutine sumupint(sum, vg, g, wg, vf, f, wf, n)

    implicit none

    ! .. Input variables
    integer, intent(in) :: n
    real (kind=dp), intent(in) :: vf
    real (kind=dp), intent(in) :: vg
    real (kind=dp), dimension(2,2), intent(in) :: wf
    real (kind=dp), dimension(2,2), intent(in) :: wg
    complex (kind=dp), dimension(2,2), intent(in) :: f
    complex (kind=dp), dimension(2,2), intent(in) :: g
    ! .. Output variables
    complex (kind=dp), intent(out) :: sum
    ! Local variables
    integer :: i, j

    sum = 0.0e0_dp
    do j = 1, n
      do i = 1, n
        sum = sum + vg*g(i, j)*wg(i, j) + vf*f(i, j)*wf(i, j)
      end do
    end do

  end subroutine sumupint

end module mod_sumupint
