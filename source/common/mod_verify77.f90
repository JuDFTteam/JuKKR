module mod_verify77

  contains
    !-------------------------------------------------------------------------------
  !> Summary: Return position of first whitespace and letter
  !> Author: 
  !> Category: KKRhost, input-output
  !> Deprecated: True ! This needs to be set to True for deprecated subroutines
  !>
  !> This sub returns the position of the first space character
  !> in ipos2, and the position of the first letter in the string
  !> STR1
  !-------------------------------------------------------------------------------
  subroutine verify77(nabc, abc, nchar, str1, ipos1, ipos2)

    implicit none
    integer, intent(in) :: nchar !! length of the string 
    integer, intent(in) :: nabc !! length of the second string
    character (len=nchar) :: str1
    character (len=nabc) :: abc
    character (len=1) :: char
    integer :: ipos, ipos1, ipos2, i, j

    ipos2 = 0

    ipos1 = index(str1, ' ')
    do j = 1,nchar-1
      char = str1(j:j+1)
      ipos = 0
      do i = 1, nchar
        ipos = index(char, abc(i:i))
        if (ipos>0) then
          ipos2 = j
          return
        end if
      end do
    end do
    return
  end subroutine verify77

end module mod_verify77