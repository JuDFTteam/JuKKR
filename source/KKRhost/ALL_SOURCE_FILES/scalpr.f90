! ************************************************************************
subroutine scalpr(x, y, z)
  use :: mod_datatypes, only: dp
  ! ************************************************************************
  ! SCALSP COMPUTES THE scalar PRODUCT OF X AND Y RETURNING
  ! IT INTO Z.
  ! ------------------------------------------------------------------------


  real (kind=dp), intent (in) :: x(*)
  real (kind=dp), intent (in) :: y(*)
  real (kind=dp), intent (out) :: z

  z = x(1)*y(1) + x(2)*y(2) + x(3)*y(3)
  return
end subroutine scalpr
