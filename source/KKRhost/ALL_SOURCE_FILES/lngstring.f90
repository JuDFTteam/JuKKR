module mod_lngstring

contains

  integer function lngstring(string, lstrmax)
    ! ********************************************************************
    ! *                                                                  *
    ! *  find position of last non-blank character in STRING(1:LSTRMAX)  *
    ! *                                                                  *
    ! ********************************************************************
    implicit none

    ! Dummy arguments
    integer :: lstrmax
    character (len=*) :: string

    ! Local variables
    character :: c
    integer :: i, ichar

    lngstring = 0
    do i = lstrmax, 1, -1
      c = string(i:i)
      if (c/=' ' .and. ichar(c)>0) then
        lngstring = i
        return
      end if
    end do
  end function lngstring

end module mod_lngstring
