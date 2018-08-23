module mod_crospr

contains

! ************************************************************************
subroutine crospr(x, y, z)
  use :: mod_datatypes, only: dp
  ! ************************************************************************
  ! CROSP COMPUTES THE CROSS PRODUCT OF X AND Y RETURNING
  ! IT INTO Z.
  ! ------------------------------------------------------------------------
  implicit none
  real (kind=dp), intent (in) :: x(*)
  real (kind=dp), intent (in) :: y(*)
  real (kind=dp), intent (out) :: z(*)

  z(1) = x(2)*y(3) - x(3)*y(2)
  z(2) = x(3)*y(1) - x(1)*y(3)
  z(3) = x(1)*y(2) - x(2)*y(1)
  return
end subroutine crospr

end module mod_crospr
