!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_ioinput

  private
  public :: ioinput

contains

  !-------------------------------------------------------------------------------
  !> Summary: Read value of keyword from inputcard
  !> Author: 
  !> Category: KKRhost, input-output
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> This subroutine is responsible for the I/O with the input file.
  !>
  !> GIVEN a KEYWORD: CHARKEY it positions the
  !> reading just after the CHARKEY if this
  !> includes a '=', or ILINE lines after the
  !> occurence of THE CHARKEY.
  !> USAGE :
  !> To read lmax include in the inputcard (ifile)
  !>
  !>     LMAX= 3      CORRECT!
  !>
  !>     or
  !>
  !>     LMAX         CORRECT!     (iline=1)
  !>      3
  !>   (without  the '=' )
  !>     LMAX
  !>  ---------                    (iline=2),etc
  !>      3
  !> be carefull in this case to put the value after the
  !> keyword example:
  !>
  !>    LMAX
  !>  3               WRONG!
  !>
  !> will NOT work
  !> Comments etc in the program are ignored.
  !>                                               1.6.99
  !>
  !> @warning
  !>  - The error handler is not working yet in all cases ...
  !>  - In this version only files 5000 lines long can be read in
  !> @endwarning
  !-------------------------------------------------------------------------------
  subroutine ioinput(key_in, value_out, skiplines, ifile, ierror)

    implicit none

    integer :: ncolio, nlinio
    parameter (ncolio=256, nlinio=5000)

    character (len=*), intent(in) :: key_in
    character (len=:), intent(out), allocatable :: value_out
    integer, intent(in) :: skiplines
    integer, intent(in) :: ifile
    integer, intent(out) :: ierror

    character(len=:), allocatable :: key
    character(len=1) :: testchar

    integer :: nlines, ios, ier, i_start_key, i_check_eq, key_len, iline
    character(len=ncolio) :: linetxt(nlinio)

    character(len=*), parameter :: abc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890-_<>'
    integer, parameter :: nlen = len(abc)

    ierror = 0

    !--------------------------------------------------------------------------------
    ! Read in of inputcard
    !--------------------------------------------------------------------------------
    ier = 0
    open (unit=ifile, status='OLD', file='inputcard', iostat=ios)
    if (ios>0) then
      write (*, *) 'Error in reading the inputcard file'
      stop
    end if

    nlines = 1
    do
      read (ifile, fmt='(A256)', iostat=ios) linetxt(nlines)
      if (ios<0 .or. nlines>=nlinio) exit
      nlines = nlines + 1
    end do
    nlines = nlines - 1
    ! write(6,*) 'LINES :',npt
    if (nlines>=nlinio) write (1337, *) 'Not all lines are read in from inputcard'
    !--------------------------------------------------------------------------------


    !bring keyword into canonical form (i.e. no blanks and all upper case)
    key = convert_to_uppercase(trim(adjustl(key_in)))
    key_len = len(key)

    do iline = 1, nlines

      !look if keyword (without the '=') exists in line
      i_start_key = index(linetxt(iline), key)

      if (i_start_key>0) then
        !the keyword was found in the line
        !try if the keyword with '=' exists in this line as well; note: achar(61): '='
        i_check_eq = index(linetxt(iline), key//achar(61))

        if (i_check_eq>0) then

          !yes, it exists: return the string after the keyword (incl. the '=')
          value_out = linetxt(iline)(i_start_key+key_len+1:)
          close (ifile)
          return

        else!i_check_eq>0

          !the value is expected skiplines lines below the keyword
          if (iline+skiplines>nlines) then
            write (1337, *) "ERROR in IoInput scanning for '", key, "':"
            write (1337, *) "No more lines in file! skiplines = ", skiplines
            ierror = 1
            close (ifile)
            return
          else
            value_out = linetxt(iline+skiplines)(i_start_key:)

            !make sanity test: character just before the keyword should be blank
            if (i_start_key>1) then
              testchar = linetxt(iline+skiplines)(i_start_key-1:i_start_key)
              if (testchar /= ' ') then
                write (1337, *) 'WARNING from IoInput: Possible ERROR!!!'
                write (1337, *) "Value for '", key, "' maybe read in incorrectly."
              end if!testchar
            end if!i_start_key>1

            close (ifile)
            return

          end if!iline+skiplines

        end if!i_check_eq>0

      end if!i_start_key>0

    end do!iline

    ierror = 1
    close (ifile)
    return

  end subroutine ioinput


  !-------------------------------------------------------------------------------
  !> Summary: Return position of first whitespace and letter
  !> Author: 
  !> Category: KKRhost, input-output
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> This sub returns the position of the first space character
  !> in ipos2, and the position of the first letter in the string
  !> STR1
  !-------------------------------------------------------------------------------
  subroutine verify77(nabc, abc, nchar, str1, ipos1, ipos2)

    implicit none
    integer :: nchar, nabc
    character (len=nchar) :: str1
    character (len=nabc) :: abc
    character (len=1) :: char
    integer :: ipos, ipos1, ipos2, i, j

    ipos2 = 0

    ipos1 = index(str1, ' ')
    do j = 1, 10
      char = str1(j:j+1)
      ipos = 0
      do i = 1, 40
        ipos = index(char, abc(i:i))
        if (ipos>0) then
          ipos2 = j
          return
        end if
      end do
    end do
    return
  end subroutine verify77

  !-------------------------------------------------------------------------------
  !> Summary: Convert a sting to upper case
  !> Author: Bernd Zimmermann
  !> Category: KKRhost, input-output
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Convert a sting to upper case.
  !> 
  !> Adjusted from https://en.wikibooks.org/wiki/Fortran/strings#Approaches_to_Case_Conversion 
  !>
  !-------------------------------------------------------------------------------
  function convert_to_uppercase(in) result(out)

    implicit none

    character(*), intent(in)  :: in
    character(len(in))        :: out
    integer                   :: i, j
    character(*), parameter   :: upp = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(*), parameter   :: low = 'abcdefghijklmnopqrstuvwxyz'

    out = in                            !transfer all characters and left-indent
    do i = 1, LEN_TRIM(out)             !all non-blanks
      j = INDEX(low, out(i:i))          !is ith character in low
    if (j > 0) out(i:i) = upp(j:j)      !yes, then subst with upp
    end do

  end function convert_to_uppercase

end module mod_ioinput
