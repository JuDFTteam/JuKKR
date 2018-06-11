! ************************************************************************
    Subroutine spatpr(a, b, c, v)
      Use mod_datatypes, Only: dp
! ************************************************************************
! SPATPR COMPUTES THE SPATIAL PRODUCT OF THREE VECTORS A,B AND C
! RETURNING IT INTO V: V=AXB.C.
! ------------------------------------------------------------------------


      Real (Kind=dp), Intent (In) :: a(*)
      Real (Kind=dp), Intent (In) :: b(*)
      Real (Kind=dp), Intent (In) :: c(*)
      Real (Kind=dp), Intent (Out) :: v

      v = 0.0E0_dp
      v = v + c(1)*(a(2)*b(3)-a(3)*b(2))
      v = v + c(2)*(a(3)*b(1)-a(1)*b(3))
      v = v + c(3)*(a(1)*b(2)-a(2)*b(1))
      Return
    End Subroutine
