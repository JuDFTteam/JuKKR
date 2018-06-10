! ************************************************************************
subroutine scalpr(x, y, z)
! ************************************************************************
!     SCALSP COMPUTES THE scalar PRODUCT OF X AND Y RETURNING
!     IT INTO Z.
! ------------------------------------------------------------------------


  double precision, intent (in) :: x(*)
  double precision, intent (in) :: y(*)
  double precision, intent (out) :: z

  z = x(1)*y(1) + x(2)*y(2) + x(3)*y(3)
  return
end subroutine
