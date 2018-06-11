    Subroutine ioinput(charkey, char, iline, ifile, ierror)
! *********************************************************
! *  This subroutine is responsible for the I/O
! *  with the input file.
! *
! *  GIVEN a KEYWORD: CHARKEY it positions the
! *  reading just after the CHARKEY if this
! *  includes a '=', or ILINE lines after the
! *  occurence of THE CHARKEY.
! *  USAGE :
! *  To read lmax include in the input card (ifile)
! *
! *      LMAX= 3      CORRECT!
! *
! *      or
! *
! *      LMAX         CORRECT!     (iline=1)
! *       3
! *    (without  the '=' )
! *      LMAX
! *   ---------                    (iline=2),etc
! *       3
! *  be carefull in this case to put the value after the
! *  keyword example:
! *
! *     LMAX
! *   3               WRONG!
! *
! * will NOT work
! * Comments etc in the program are ignored.
! *                                               1.6.99
! *
! * The error handler is not working yet in all cases ....
! * In this version only files 5000 lines long can be read in
! *******************************************************
      Implicit None
      Integer :: nchar, nabc, ncolio, nlinio
      Parameter (nchar=16, nabc=40, ncolio=256, nlinio=5000)
      Character (Len=nchar) :: charkey !*NCHAR

      Character (Len=ncolio) :: char !*NCOLIO

      Integer :: iline, ierror, ifile
      Integer :: i, ios, ier, npt, ilen, ipos, ipos1, iklen
      Character (Len=ncolio) :: string(nlinio) !*NCOLIO

      Character (Len=ncolio) :: string1 !*NCOLIO

      Character (Len=nabc) :: abc !*NABC

      Character :: atest
      Data abc/'ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890-_<>'/

      ierror = 0
      ier = 0
      char(1:50) = '                                                   '
      Open (Unit=ifile, Status='OLD', File='inputcard', Iostat=ios, Err=100)
      If (ios>0) Then
        Write (6, *) 'Error in reading the inputcard file'
        Stop
      End If


      npt = 1
      Do
        Read (ifile, Fmt='(A256)', Iostat=ios) string(npt)
        If (ios<0 .Or. npt>=nlinio) Exit
        npt = npt + 1
      End Do
      npt = npt - 1
!          write(6,*) 'LINES :',npt
      If (npt>=nlinio) Write (1337, *) &
        'Not all lines are read in from inputcard'

! 2 lines below where changed
!       ILEN = VERIFY(CHARKEY,ABC)
!       IKLEN= VERIFY(CHARKEY,' ')
! for linux
      Call verify77(nabc, abc, nchar, charkey, ilen, iklen)
! for linux
!        write(6,*) CHARKEY(1:ILEN-1),ILEN,IKLEN
      If (ilen<1) Then
        Write (1337, *) 'Input ERROR!'
        Write (1337, *) 'Cannot evaluate : ', charkey
        Write (1337, *) 'IoInput is returning no value! '
        Return
      End If

      Do i = 1, npt ! loop in all line
        string1 = '   ' // string(i) ! shift by 2 characters
        ipos = index(string1, charkey(1:ilen-1))
! return the position of occurence
        If (ipos/=0) Then
          If (ipos<4) Then
            Write (6, *) 'CONSISTENCY ERROR IOINPUT!'
            Stop
          End If
!                write(6,*) 'ipos is not zero',CHARKEY//'=','**'
          ipos1 = index(string1, charkey(1:ilen-1)//achar(61))
          If (ipos1/=0) Then
! return the string after 'CHARKEY='
            char = string1(ipos1+ilen:)
!                    write(6,*) CHARKEY,CHAR ! test
            Close (ifile)
            Return
          Else
! return the ILINE line below this CHARKEY
            If (i+iline<=npt) Then
!                    write(6,*) IPOS,ILEN

              char = string(i+iline)(ipos-3:)
              If (ipos-4>0) Then ! Changed on 28.01.2000
                atest = string(i+iline)(ipos-4:ipos-3)
                If (atest/=' ') Then
                  Write (1337, *) 'Possible ERROR !!!'
                  Write (1337, *) 'Parameter ', charkey, &
                    ' maybe read in incorrectrly'
                End If
              End If
!                   write(6,*) CHARKEY,CHAR ! test
              Close (ifile)
              Return
            Else
              Write (1337, *) 'IoInput : No more lines in file '
            End If
          End If
        End If
      End Do ! i=1,npt
      ier = 1
      ierror = ierror + ier
!ccc       if (CHAR(1:20).eq.'                    ') then
!ccc       write(6,*) 'Parameter ........ ',CHARKEY , ' NOT found'
!ccc       write(6,*) 'Check your inputcard'
!ccc       end if
      Close (ifile)
      Return
100   Write (6, *) ' Error while reading..... ', charkey
      Write (6, *) ' Check your  inputcard ! '
      Stop
    End Subroutine

    Subroutine verify77(nabc, abc, nchar, str1, ipos1, ipos2)
! This sub returns the position of the first space character
! in ipos2, and the position of the first letter in the string
! STR1
      Implicit None
      Integer :: nchar, nabc
      Character (Len=nchar) :: str1 !*NCHAR

      Character (Len=nabc) :: abc !*NABC

      Character (Len=1) :: char !*1

      Integer :: ipos, ipos1, ipos2, i, j
!        DATA ABC/'ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890-_<>'/
      ipos2 = 0

      ipos1 = index(str1, ' ')
      Do j = 1, 10
        char = str1(j:j+1)
!           write(6,*) 'char : ',j, char
        ipos = 0
        Do i = 1, 40
          ipos = index(char, abc(i:i))
          If (ipos>0) Then
            ipos2 = j
            Return
          End If
        End Do

      End Do
      Return
    End Subroutine
