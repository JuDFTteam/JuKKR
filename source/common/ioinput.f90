!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_ioinput

  private
  public :: ioinput, convert_to_uppercase

  logical, save :: inputcard_is_read = .false.
  integer, save :: ncols  = 0!! Number of columns of inputcard
  integer, save :: nlines = 0!! Number of lines   of inputcard
  character(len=:), dimension (:), allocatable, save :: linetxt !! Saved inputcard

contains

  !-------------------------------------------------------------------------------
  !> Summary: Read value of keyword from inputcard
  !> Author: Bernd Zimmermann
  !> Category: KKRhost, input-output
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> This subroutine is responsible for the I/O with the input file.
  !>
  !> Given a keyword `key_in` the subroutine positions the
  !> reading just after the keyword if this
  !> includes a '=', or `skipline` lines after the
  !> occurence of the keyword.
  !> Usage:
  !> To read `LMAX` include in the inputcard
  !>
  !>     `LMAX= 3`                CORRECT!
  !>  or `LMAX=3`                 CORRECT!
  !>  or `LMAX=   3  my comment`  CORRECT!
  !>
  !>     or
  !>
  !>     LMAX                    CORRECT for skipline=1!
  !>      3                      (note the missing '=')
  !>   
  !>     or
  !>
  !>     LMAX
  !>  ---------                  CORRECT for skipline=2!
  !>      3
  !>
  !> Be carefull in the latter case to put the value not too far to the left
  !> keyword. Example:
  !>
  !>    LMAX
  !>   3               WRONG!
  !>
  !>                                               1.6.99
  !>         Refactored by Bernd Zimmermann on 13.11.2018
  !> @todo
  !>   `ifile` should be made a module-local parameter. Ensures better readability on calling routines. (BZ)
  !> @endtodo
  !-------------------------------------------------------------------------------
  subroutine ioinput(key_in, value_out, skiplines, ifile, ierror)

    implicit none

    character (len=*), intent(in) :: key_in !!string of arbitrary length specifying the keyword. Leading and trailing blanks will be ignored.
    character (len=:), intent(out), allocatable :: value_out !!string containing the rest of the line (i.e. after the mark as determined by the keyword).
    integer, intent(in) :: skiplines !!integer specifying which line (after the occurrence of the keyword) is read in (only effective if the keyword is NOT directly wollowed by a '=').
    integer, intent(in) :: ifile     !!unit of the inputcard
    integer, intent(out) :: ierror   !!error state. 0 = no error, 1 = error

    character(len=:), allocatable :: key
    character(len=1) :: testchar

    integer :: i_start_key, i_check_eq, key_len, iline

    ierror = 0

    call read_inputcard(ifile)

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
          return

        else!i_check_eq>0

          !the value is expected skiplines lines below the keyword
          if (iline+skiplines>nlines) then
            write (1337, *) "ERROR in IoInput scanning for '", key, "':"
            write (1337, *) "Not enough lines in inputcard! skiplines = ", skiplines
            ierror = 1
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

            return

          end if!iline+skiplines

        end if!i_check_eq>0

      end if!i_start_key>0

    end do!iline

    ierror = 1
    return

  end subroutine ioinput

  !-------------------------------------------------------------------------------
  !> Summary: Read inputcard and store to array
  !> Author: Bernd Zimmermann
  !> Category: KKRhost, input-output
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Read inputcard once and store into module-local array `linetxt`.
  !> The most efficient size of `linetxt` is automatically determined. All lines of 
  !>   the inputcard are read in, but the length of each line is limited by the 
  !>   constant parameter `maxcol`.
  !> 
  !> @warning
  !>   The length of each line is limited by the constant parameter `maxcol`.
  !> @endwarning
  !-------------------------------------------------------------------------------
  subroutine read_inputcard(ifile)

    implicit none
    integer, intent(in) :: ifile
    integer, parameter :: maxcol = 2048

    integer :: ierror, ios, nlinetmp, ncoltmp, iline, ipos
    character(len=maxcol) :: line
    character(len=:), allocatable :: linetmp

    if (.not.inputcard_is_read) then
      !--------------------------------------------------------------------------------
      ! open inputcard
      !--------------------------------------------------------------------------------
      open (unit=ifile, status='OLD', file='inputcard', iostat=ios)
      if (ios>0) then
        write (*, *) 'Error reading the inputcard file.'
        stop
      end if

      !--------------------------------------------------------------------------------
      ! determine array sizes
      !--------------------------------------------------------------------------------
      nlinetmp = 0
      ncoltmp  = 0
      do
        read (ifile, fmt='(A)', iostat=ios) line
        if (ios<0 ) exit

        nlinetmp = nlinetmp + 1
        ncoltmp = max(ncoltmp,len_trim(line))
      end do

      !--------------------------------------------------------------------------------
      ! modify module-wide variables and allocate 
      !--------------------------------------------------------------------------------
      if(allocated(linetxt)) deallocate(linetxt)
      nlines = nlinetmp
      ncols  = ncoltmp
      allocate(character(len=ncols) :: linetxt(nlines))
      allocate(character(len=ncols) :: linetmp)

      !--------------------------------------------------------------------------------
      ! read in txt from inputcard 
      !--------------------------------------------------------------------------------
      rewind (ifile)
      do iline = 1, nlines
        read (ifile, fmt='(A)') linetmp

        ! now remove comments which follow the '#' sign
        ! to do this we overwrite linetmp starting from the '#' sign with spaces
        ipos = index(linetmp, '#')
        if (ipos>0) then
          write(1337,'(A)') 'WARNING: Found "#" in inputcard and change line for reading'
          write(1337,'(2A)') 'Original line:', linetmp
          linetmp(ipos:) = ' '
          write(1337,'(2A)') 'Changed to:', linetmp
        end if
 
        call capitalize_escaped_keywords(linetmp, linetxt(iline), ierror)
        if (ierror==1) then
          write(1337,'(A,I5,A)') 'Warning capitalizing line ', iline, ' of inputcard.'
          write(*,   '(A,I5,A)') 'Warning capitalizing line ', iline, ' of inputcard.'
        else if (ierror==2) then
          write(1337,'(A,I5,A)') 'Error capitalizing line ', iline, ' of inputcard.'
          write(*,   '(A,I5,A)') 'Error capitalizing line ', iline, ' of inputcard.'
          stop
        end if

      end do ! iline

      !--------------------------------------------------------------------------------
      ! close inputcard and finish up
      !--------------------------------------------------------------------------------
      close (ifile)
      inputcard_is_read = .true.

    end if ! .not.inputcard_is_read
    
  end subroutine read_inputcard

  !-------------------------------------------------------------------------------
  !> Summary: Capitalize all escaped keywords in a string
  !> Author: Bernd Zimmermann
  !> Category: KKRhost, input-output
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Capitalized all escaped keywords in a string, i.e. it converts `<lmax>` to `<LMAX>`, but not `lmax`.
  !>
  !> `error=0` on sucessful return
  !> `error=1` signals warnings (i.e. hatd-coded threshold reached)
  !> `error=2` signals problems with order of escaping characters
  !>
  !> @note
  !>   Can only treat up to `nkeymax` escaped keywords per line and prints a warning if this threshold is reached.
  !> @endnote
  !> @warning
  !>   The routine might get confused if addidional characters `<` or `>`, not used for escaping, are preset in `string`. Some of this is tried to be detected and sugnaled via `error=1`
  !> @endwarning
  !-------------------------------------------------------------------------------
  subroutine capitalize_escaped_keywords(line_in, line_out, ierror)

    implicit none

    character, parameter :: cstart='<', cstop='>'
    integer, parameter   :: nkeymax = 1000

    character(len=*), intent(in) :: line_in
    character(len=len(line_in)), intent(out) :: line_out
    integer,          intent(out) :: ierror

    integer :: nlen, isafe, ipos1, ipos2, nkeys, ikey
    character(len=:), allocatable :: word
    integer, dimension(1:nkeymax) :: istart, istop

    !init
    ierror = 0
    isafe = 0 !lower bound for scan
    nkeys = 0
    nlen  = len_trim(line_in)
    line_out = line_in

    do while ( (isafe>=0) .and. (nkeys<nkeymax) )
      ipos1 = index(line_in(isafe+1:nlen),cstart)
      ipos2 = index(line_in(isafe+1:nlen),cstop)
      if ((ipos1>0) .and. (ipos2>0) .and. ipos1<ipos2) then
        nkeys = nkeys + 1
        istart(nkeys) = isafe + ipos1 + 1 !+1 to remove leading `cstart`
        istop(nkeys)  = isafe + ipos2 - 1 !-1 to remove trailing `cstop`
        isafe         = isafe + ipos2 + 1 !count lower bound up
      else if ((ipos1>0) .and. (ipos2>0)) then
        ierror = 2
        write (*,*) "Problem with order of '<' or '>' in 'capitalize_escaped_keywords'."
        write (*,*) "Are additional '<' or '>' present?"
        return
      else
        isafe  = -1
      end if
    end do

    if (nkeys==nkeymax) then
      ierror = 1
      write (*,*) "WARNING! Hit 'nkeymax' threshold in 'capitalize_escaped_keywords'."
      write (*,*) "         Maybe not all keywords have been capitalized."
    end if

    do ikey = 1, nkeys
      word = line_in(istart(ikey):istop(ikey))
      if (is_keyword(word)) then
        line_out(istart(ikey):istop(ikey)) = convert_to_uppercase(word)
      end if

    end do!ikey

  end subroutine capitalize_escaped_keywords

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

  !-------------------------------------------------------------------------------
  !> Summary: Checks if a string is a valid keyword
  !> Author: Bernd Zimmermann
  !> Category: KKRhost, input-output
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Checks if a string is a valid keyword, meanig that it only consists of 
  !> alpha-numeric characters, as well as the dash (-) and underscore (_).
  !> 
  !-------------------------------------------------------------------------------
  function is_keyword(word) result(l_kw)
    implicit none
    character(len=*), intent(in) :: word
    logical :: l_kw
    integer :: nlen, ii
    character(len=*), parameter :: abcbig = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890-_'
    character :: ch

    l_kw = .true.

    nlen = len(word)
    do ii = 1, nlen
      ch = word(ii:ii)
      if (index(abcbig,ch)==0) then
        l_kw = .false.
        return
      end if
    end do

    return

  end function is_keyword

end module mod_ioinput
