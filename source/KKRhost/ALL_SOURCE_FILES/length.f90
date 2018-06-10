! ************************************************************************
integer function length(s, max)
! ************************************************************************

  character (len=1), intent (inout) :: s(*)
  integer, intent (in) :: max

  integer :: i
! ------------------------------------------------------------------------
  i = max

  do while (s(i)==' ')
    i = i - 1
  end do

  length = i

  return
end function
