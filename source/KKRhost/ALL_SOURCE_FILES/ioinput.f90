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
  subroutine ioinput(charkey, char, iline, ifile, ierror)

    implicit none
    integer :: nchar, nabc, ncolio, nlinio
    parameter (nchar=16, nabc=40, ncolio=256, nlinio=5000)
    character (len=nchar) :: charkey
    character (len=ncolio) :: char
    integer :: iline, ierror, ifile
    integer :: i, ios, ier, npt, ilen, ipos, ipos1, iklen
    character (len=ncolio) :: string(nlinio)
    character (len=ncolio) :: string1
    character (len=nabc) :: abc
    character :: atest
    data abc/'ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890-_<>'/


    ierror = 0
    ier = 0
    char(1:50) = '                                                   '
    open (unit=ifile, status='OLD', file='inputcard', iostat=ios, err=100)
    if (ios>0) then
      write (6, *) 'Error in reading the inputcard file'
      stop
    end if

    npt = 1
    do
      read (ifile, fmt='(A256)', iostat=ios) string(npt)
      if (ios<0 .or. npt>=nlinio) exit
      npt = npt + 1
    end do
    npt = npt - 1
    ! write(6,*) 'LINES :',npt
    if (npt>=nlinio) write (1337, *) 'Not all lines are read in from inputcard'

    ! 2 lines below where changed
    ! ILEN = VERIFY(CHARKEY,ABC)
    ! IKLEN= VERIFY(CHARKEY,' ')
    ! for linux
    call verify77(nabc, abc, nchar, charkey, ilen, iklen)
    ! for linux
    ! write(6,*) CHARKEY(1:ILEN-1),ILEN,IKLEN
    if (ilen<1) then
      write (1337, *) 'Input ERROR!'
      write (1337, *) 'Cannot evaluate : ', charkey
      write (1337, *) 'IoInput is returning no value! '
      return
    end if

    do i = 1, npt                  ! loop in all line
      string1 = '   ' // string(i) ! shift by 2 characters
      ipos = index(string1, charkey(1:ilen-1))
      ! return the position of occurence
      if (ipos/=0) then
        if (ipos<4) then
          write (6, *) 'CONSISTENCY ERROR IOINPUT!'
          stop
        end if
        ! write(6,*) 'ipos is not zero',CHARKEY//'=','**'
        ipos1 = index(string1, charkey(1:ilen-1)//achar(61))
        if (ipos1/=0) then
          ! return the string after 'CHARKEY='
          char = string1(ipos1+ilen:)
          ! write(6,*) CHARKEY,CHAR ! test
          close (ifile)
          return
        else
          ! return the ILINE line below this CHARKEY
          if (i+iline<=npt) then
            ! write(6,*) IPOS,ILEN

            char = string(i+iline)(ipos-3:)
            if (ipos-4>0) then     ! Changed on 28.01.2000
              atest = string(i+iline)(ipos-4:ipos-3)
              if (atest/=' ') then
                write (1337, *) 'Possible ERROR !!!'
                write (1337, *) 'Parameter ', charkey, ' maybe read in incorrectrly'
              end if
            end if
            ! write(6,*) CHARKEY,CHAR ! test
            close (ifile)
            return
          else
            write (1337, *) 'IoInput : No more lines in file '
          end if
        end if
      end if
    end do                         ! i=1,npt
    ier = 1
    ierror = ierror + ier
    ! ccc       if (CHAR(1:20).eq.'                    ') then
    ! ccc       write(6,*) 'Parameter ........ ',CHARKEY , ' NOT found'
    ! ccc       write(6,*) 'Check your inputcard'
    ! ccc       end if
    close (ifile)
    return
100 write (6, *) ' Error while reading..... ', charkey
    write (6, *) ' Check your  inputcard ! '
    stop
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

end module mod_ioinput
