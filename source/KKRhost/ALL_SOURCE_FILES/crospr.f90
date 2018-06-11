! ************************************************************************
    Subroutine crospr(x, y, z)
      Use mod_datatypes, Only: dp
! ************************************************************************
!     CROSP COMPUTES THE CROSS PRODUCT OF X AND Y RETURNING
!     IT INTO Z.
! ------------------------------------------------------------------------
      Implicit None
      Real (Kind=dp), Intent (In) :: x(*)
      Real (Kind=dp), Intent (In) :: y(*)
      Real (Kind=dp), Intent (Out) :: z(*)

      z(1) = x(2)*y(3) - x(3)*y(2)
      z(2) = x(3)*y(1) - x(1)*y(3)
      z(3) = x(1)*y(2) - x(2)*y(1)
      Return
    End Subroutine
