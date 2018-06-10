! ************************************************************************
subroutine spatpr(a, b, c, v)
! ************************************************************************
! SPATPR COMPUTES THE SPATIAL PRODUCT OF THREE VECTORS A,B AND C
! RETURNING IT INTO V: V=AXB.C.
! ------------------------------------------------------------------------


  double precision, intent (in) :: a(*)
  double precision, intent (in) :: b(*)
  double precision, intent (in) :: c(*)
  double precision, intent (out) :: v

  v = 0.0d0
  v = v + c(1)*(a(2)*b(3)-a(3)*b(2))
  v = v + c(2)*(a(3)*b(1)-a(1)*b(3))
  v = v + c(3)*(a(1)*b(2)-a(2)*b(1))
  return
end subroutine
