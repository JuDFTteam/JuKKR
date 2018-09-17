module mod_veq
  use :: mod_datatypes, only: dp
  private :: dp

contains

  ! ************************************************************************
  subroutine veq(a, b)
    ! ************************************************************************
    implicit none

    real (kind=dp) :: a(*)
    real (kind=dp) :: b(*)

    integer :: i

    do i = 1, 3
      b(i) = a(i)
    end do
  end subroutine veq

end module mod_veq
