module mod_ddet33

contains

function ddet33(matrix)
  ! - calculates the determinant of a 3X3 matrix
  ! ----------------------------------------------------------------------
  ! i Inputs:
  ! i   matrix:input matrix
  ! o Outputs:
  ! o   ddet33:determinant
  ! ----------------------------------------------------------------------
  use :: mod_datatypes, only: dp
   use mod_cross
  implicit none
  real (kind=dp) :: ddet33
  ! Passed parameters:
  real (kind=dp) :: matrix(*)
  ! Local parameters:
  real (kind=dp) :: ddot, m1cm2(3)
  ! external calls:
  external :: cross, ddot

  call cross(matrix(4), matrix(7), m1cm2)
  ddet33 = ddot(3, matrix(1), 1, m1cm2, 1)

end function ddet33

end module mod_ddet33
