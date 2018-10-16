module mod_ddet33

  use :: mod_datatypes, only: dp
  private :: dp
  public :: ddet33

contains

  !-------------------------------------------------------------------------------
  !> Summary: Determinant of 3x3 double precision matrix
  !> Author: 
  !> Category: KKRhost, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Calculates the determinant of a 3X3 matrix
  !> 
  !> Inputs:
  !>   matrix:input matrix
  !> Outputs:
  !>   ddet33: determinant
  !-------------------------------------------------------------------------------
  real(kind=dp) function ddet33(matrix)
    use :: mod_cross, only: cross
    implicit none
    ! Passed parameters:
    real (kind=dp), intent(in) :: matrix(*) !! ipnut matrix
    ! Local parameters:
    real (kind=dp) :: m1cm2(3) !! temporary value of cross product
    ! external calls:
    real (kind=dp), external :: ddot !! ddot is a LAPACK function

    call cross(matrix(4), matrix(7), m1cm2)
    ddet33 = ddot(3, matrix(1), 1, m1cm2, 1)

  end function ddet33

end module mod_ddet33
