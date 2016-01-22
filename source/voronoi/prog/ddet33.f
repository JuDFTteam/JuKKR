      REAL*8         function ddet33(matrix)
!- calculates the determinant of a 3X3 matrix
! ----------------------------------------------------------------------
!i Inputs:
!i   matrix:input matrix
!o Outputs:
!o   ddet33:determinant
! ----------------------------------------------------------------------
      implicit none
! Passed parameters:
      REAL*8         matrix(*)
! Local parameters:
      REAL*8         ddot,m1cm2(3),mapm(3)
      INTEGER   I
! external calls:
      external cross,ddot

      call cross(matrix(4),matrix(7),m1cm2)
      do i=1,3
         mapm(i) = matrix(i)
      end do
      ddet33 = ddot(3,mapm,1,m1cm2,1)

      end
