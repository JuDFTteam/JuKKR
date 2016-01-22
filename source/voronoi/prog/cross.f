      subroutine cross(x1,x2,cr12)
!- Cross product cr12 = X1 cross X2
! ----------------------------------------------------------------------
!i Inputs:
!i   x1    :first vector to multiply
!i   x2    :second vector to multiply
!o Outputs:
!o   cr12  :cross product
! ----------------------------------------------------------------------
      implicit none
! Passed parameters:
      REAL*8           x1(3),x2(3),cr12(3)

      cr12(1) = x1(2)*x2(3) - x1(3)*x2(2)
      cr12(2) = x1(3)*x2(1) - x1(1)*x2(3)
      cr12(3) = x1(1)*x2(2) - x1(2)*x2(1)
      return
      end
