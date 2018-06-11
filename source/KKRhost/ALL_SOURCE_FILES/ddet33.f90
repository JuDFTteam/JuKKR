    Function ddet33(matrix)
      Use mod_datatypes, Only: dp
!- calculates the determinant of a 3X3 matrix
! ----------------------------------------------------------------------
!i Inputs:
!i   matrix:input matrix
!o Outputs:
!o   ddet33:determinant
! ----------------------------------------------------------------------
      Implicit None
      Real (Kind=dp) :: ddet33
! Passed parameters:
      Real (Kind=dp) :: matrix(*)
! Local parameters:
      Real (Kind=dp) :: ddot, m1cm2(3)
! external calls:
      External :: cross, ddot

      Call cross(matrix(4), matrix(7), m1cm2)
      ddet33 = ddot(3, matrix(1), 1, m1cm2, 1)

    End Function
