! ************************************************************************
subroutine veq(a, b)
  use :: mod_datatypes, only: dp
  ! ************************************************************************
  implicit none

  real (kind=dp) :: a(*)
  real (kind=dp) :: b(*)

  integer :: i

  do i = 1, 3
    b(i) = a(i)
  end do
end subroutine veq
