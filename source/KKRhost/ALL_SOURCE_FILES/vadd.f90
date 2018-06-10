! ************************************************************************
subroutine vadd(a, b, c)
! ************************************************************************

  double precision, intent (in) :: a(*)
  double precision, intent (in) :: b(*)
  double precision, intent (out) :: c(*)

  integer :: i

  do i = 1, 3
    c(i) = a(i) + b(i)
  end do
  return
end subroutine
