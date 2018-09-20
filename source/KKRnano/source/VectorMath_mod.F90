module VectorMath_mod
!-------------------------------------------------------------------------------
!> Summary: Spectial operations on 3D vectors and 3x3 matrices
!> Author: Elias Rabel, Alexander R Thiess, Paul F Baumeister
!> Category: KKRnano
!-------------------------------------------------------------------------------
  implicit none
  private
  
  public :: ddet33, cross, spatpr
  
  contains

  function cross(x1, x2) result(x3)
!- Cross product x3 = x1 cross x2
! ----------------------------------------------------------------------
!i Inputs:
!i   x1    :first vector to multiply
!i   x2    :second vector to multiply
!o Outputs:
!o   x3  :cross product
! ----------------------------------------------------------------------
    double precision, intent(in) :: x1(3), x2(3)
    double precision :: x3(3) ! result

    x3(1) = x1(2)*x2(3) - x1(3)*x2(2)
    x3(2) = x1(3)*x2(1) - x1(1)*x2(3)
    x3(3) = x1(1)*x2(2) - x1(2)*x2(1)
  endfunction cross
      
      
  double precision function ddet33(matrix)
!- calculates the determinant of a 3X3 matrix
! ----------------------------------------------------------------------
!i Inputs:
!i   matrix:input matrix
!o Outputs:
!o   ddet33:determinant
! ----------------------------------------------------------------------
    double precision, intent(in) :: matrix(*)

    double precision :: m1cm2(3)
    double precision, external :: ddot

    m1cm2 = cross(matrix(4:6), matrix(7:9))
    ddet33 = ddot(3, matrix(1), 1, m1cm2, 1) ! = dot_product(matrix(1:3), m1cm2(1:3)) ! todo remove call to ddot

  endfunction ! ddet33


  double precision function spatpr(a, b, c) result(v)
! spatpr computes the spatial product of three vectors a,b and c returning it into v: v=axb.c
    double precision, intent(in) :: a(3), b(3), c(3)
    v = c(1)*(a(2)*b(3) - a(3)*b(2)) &
      + c(2)*(a(3)*b(1) - a(1)*b(3)) &
      + c(3)*(a(1)*b(2) - a(2)*b(1))
  endfunction spatpr

endmodule VectorMath_mod     
