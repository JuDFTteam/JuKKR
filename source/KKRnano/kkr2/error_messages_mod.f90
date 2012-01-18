!> Module with routines that display error messages.
module ErrorMessages_mod
  implicit none

  contains

  !> Display memory error message, stop program.
  !> @param msg An optional error message for display
  subroutine fatalMemoryError(msg)
    implicit none
    character(len=*), intent(in), optional :: msg

    write(*,*) "FATAL error, failure to (de)allocate memory."
    if (present(msg)) then
      write(*,*) msg
    end if

    stop

  end subroutine

end module ErrorMessages_mod
