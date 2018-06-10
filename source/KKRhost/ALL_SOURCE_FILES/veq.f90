! ************************************************************************
subroutine veq(a, b)
! ************************************************************************
  implicit none

  double precision :: a(*)
  double precision :: b(*)

  integer :: i

  do i = 1, 3
    b(i) = a(i)
  end do
end subroutine
