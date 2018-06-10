SUBROUTINE ioinput(charkey,CHAR,iline,ifile,ierror)
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
       implicit none
       INTEGER NCHAR,NABC,NCOLIO,NLINIO
       PARAMETER(NCHAR=16,NABC=40,NCOLIO=256,NLINIO=5000)
       CHARACTER (len=nchar) CHARKEY !*NCHAR
       CHARACTER (len=ncolio) CHAR !*NCOLIO
       INTEGER ILINE,IERROR,IFILE
       integer i,ios,ier,npt,ilen,ipos,ipos1,iklen
       CHARACTER (len=ncolio) STRING(NLINIO) !*NCOLIO
       CHARACTER (len=ncolio) STRING1 !*NCOLIO
       CHARACTER (len=nabc) ABC !*NABC
       CHARACTER ATEST
       DATA ABC/'ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890-_<>'/

ierror = 0
ier = 0
CHAR(1:50)='                                                   '
OPEN(UNIT=ifile,STATUS='OLD',FILE='inputcard',IOSTAT=ios, ERR=2000)
IF (ios > 0) THEN
  WRITE(6,*) "Error in reading the inputcard file"
  STOP
END IF


npt = 1
DO
  READ(ifile,FMT='(A256)',IOSTAT=ios) string(npt)
  IF(ios < 0.OR.npt >= nlinio) EXIT
  npt = npt + 1
END DO
npt = npt - 1
!          write(6,*) 'LINES :',npt
IF (npt >= nlinio) WRITE(1337,*)'Not all lines are read in from inputcard'

! 2 lines below where changed
!       ILEN = VERIFY(CHARKEY,ABC)
!       IKLEN= VERIFY(CHARKEY,' ')
! for linux
CALL verify77(nabc,abc,nchar,charkey,ilen,iklen)
! for linux
!        write(6,*) CHARKEY(1:ILEN-1),ILEN,IKLEN
IF(ilen < 1) THEN
  WRITE(1337,*) 'Input ERROR!'
  WRITE(1337,*) 'Cannot evaluate : ',charkey
  WRITE(1337,*) 'IoInput is returning no value! '
  RETURN
END IF

DO i=1,npt       ! loop in all line
  string1 = '   ' // string(i)     ! shift by 2 characters
  ipos = INDEX(string1,charkey(1:ilen-1))
! return the position of occurence
  IF (ipos /= 0) THEN
    IF (ipos < 4) THEN
      WRITE(6,*) 'CONSISTENCY ERROR IOINPUT!'
      STOP
    END IF
!                write(6,*) 'ipos is not zero',CHARKEY//'=','**'
    ipos1= INDEX(string1,charkey(1:ilen-1)//achar(61))
    IF (ipos1 /= 0) THEN
! return the string after 'CHARKEY='
      CHAR = string1(ipos1+ilen:)
!                    write(6,*) CHARKEY,CHAR ! test
      CLOSE(ifile)
      RETURN
    ELSE
! return the ILINE line below this CHARKEY
      IF (i+iline <= npt) THEN
!                    write(6,*) IPOS,ILEN
        
        CHAR = string(i+iline)(ipos-3:)
        IF (ipos-4 > 0) THEN    ! Changed on 28.01.2000
          atest= string(i+iline)(ipos-4:ipos-3)
          IF (atest /= ' ') THEN
            WRITE(1337,*) 'Possible ERROR !!!'
            WRITE(1337,*) 'Parameter ',charkey, ' maybe read in incorrectrly'
          END IF
        END IF
!                   write(6,*) CHARKEY,CHAR ! test
        CLOSE(ifile)
        RETURN
      ELSE
        WRITE(1337,*) 'IoInput : No more lines in file '
      END IF
    END IF
  END IF
END DO  ! i=1,npt
ier = 1
ierror = ierror + ier
!ccc       if (CHAR(1:20).eq.'                    ') then
!ccc       write(6,*) 'Parameter ........ ',CHARKEY , ' NOT found'
!ccc       write(6,*) 'Check your inputcard'
!ccc       end if
CLOSE(ifile)
RETURN
2000 WRITE(6,*) ' Error while reading..... ',charkey
WRITE(6,*) ' Check your  inputcard ! '
STOP
END SUBROUTINE ioinput

SUBROUTINE verify77(nabc,abc,nchar,str1,ipos1,ipos2)
! This sub returns the position of the first space character
! in ipos2, and the position of the first letter in the string
! STR1
        implicit none  
         INTEGER NCHAR,NABC
         CHARACTER (len=nchar) STR1 !*NCHAR
         CHARACTER (len=nabc) ABC !*NABC
         CHARACTER (len=1) CHAR !*1
         integer ipos,ipos1,ipos2,i,j
!        DATA ABC/'ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890-_<>'/
ipos2 =0

ipos1 = INDEX(str1,' ')
DO j=1,10
  CHAR = str1(j:j+1)
!           write(6,*) 'char : ',j, char
  ipos = 0
  DO i=1,40
    ipos = INDEX(CHAR,abc(i:i))
    IF (ipos > 0) THEN
      ipos2 = j
      RETURN
    END IF
  END DO
  
END DO
RETURN
END SUBROUTINE verify77

