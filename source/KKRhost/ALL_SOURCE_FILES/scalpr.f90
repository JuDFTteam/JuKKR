! ************************************************************************
    Subroutine scalpr(x, y, z)
      Use mod_datatypes, Only: dp
! ************************************************************************
!     SCALSP COMPUTES THE scalar PRODUCT OF X AND Y RETURNING
!     IT INTO Z.
! ------------------------------------------------------------------------


      Real (Kind=dp), Intent (In) :: x(*)
      Real (Kind=dp), Intent (In) :: y(*)
      Real (Kind=dp), Intent (Out) :: z

      z = x(1)*y(1) + x(2)*y(2) + x(3)*y(3)
      Return
    End Subroutine
