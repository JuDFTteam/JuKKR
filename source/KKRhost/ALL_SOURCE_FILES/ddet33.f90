DOUBLE PRECISION FUNCTION ddet33(matrix)
!- calculates the determinant of a 3X3 matrix
! ----------------------------------------------------------------------
!i Inputs:
!i   matrix:input matrix
!o Outputs:
!o   ddet33:determinant
! ----------------------------------------------------------------------
      implicit none
! Passed parameters:
      double precision matrix(*)
! Local parameters:
      double precision ddot,m1cm2(3)
! external calls:
      external cross,ddot

CALL cross(matrix(4),matrix(7),m1cm2)
ddet33 = ddot(3,matrix(1),1,m1cm2,1)

END FUNCTION ddet33
