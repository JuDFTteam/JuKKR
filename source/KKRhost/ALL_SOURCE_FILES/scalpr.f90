!------------------------------------------------------------------------------------
!> Summary: Calculates the scalar product between two vectors 
!> Author: 
!> Calculates the scalar product between two vectors 
!------------------------------------------------------------------------------------
!> @note Jonathan Chico: This can be generalized to vectors of any dimmension by using the 
!> call `sum(x(:)*y(:))`. This should be more efficient for `ifort` as it vectorizes
!> the `sum` call automatically
!> @endnote
!------------------------------------------------------------------------------------
module mod_scalpr
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculates the scalar product between two vectors 
  !> Author: 
  !> Category: numerical-tools, KKRhost
  !> Deprecated: False 
  !> Calculates the scalar product between two vectors 
  !-------------------------------------------------------------------------------
  !> @note Jonathan Chico: This can be generalized to vectors of any dimmension by using the 
  !> call `sum(x(:)*y(:))`. This should be more efficient for `ifort` as it vectorizes
  !> the `sum` call automatically
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine scalpr(x, y, z)

    implicit none

    real (kind=dp), dimension(*), intent (in) :: x
    real (kind=dp), dimension(*), intent (in) :: y
    real (kind=dp), intent (out) :: z

    z = x(1)*y(1) + x(2)*y(2) + x(3)*y(3)
    return
  end subroutine scalpr

end module mod_scalpr
