module mod_vadd
  use :: mod_datatypes, only: dp
  private :: dp

contains

  ! ************************************************************************
  subroutine vadd(a, b, c)
    ! ************************************************************************

    real (kind=dp), intent (in) :: a(*)
    real (kind=dp), intent (in) :: b(*)
    real (kind=dp), intent (out) :: c(*)

    integer :: i

    do i = 1, 3
      c(i) = a(i) + b(i)
    end do
    return
  end subroutine vadd

end module mod_vadd
