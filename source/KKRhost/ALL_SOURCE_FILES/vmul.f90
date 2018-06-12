! ************************************************************************
subroutine vmul(a, b, c)
  use :: mod_datatypes, only: dp
  ! ************************************************************************

  real (kind=dp), intent (in) :: a(*)
  real (kind=dp), intent (in) :: b
  real (kind=dp), intent (out) :: c(*)

  integer :: i

  do i = 1, 3
    c(i) = b*a(i)
  end do
  return
end subroutine vmul
