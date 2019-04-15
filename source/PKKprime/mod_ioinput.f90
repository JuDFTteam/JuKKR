!-----------------------------------------------------------------------------------------!
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of kk-prime@juKKR and available as free software under the conditions !
! of the MIT license as expressed in the LICENSE file in more detail.                     !
!-----------------------------------------------------------------------------------------!


module mod_ioinput

  implicit none

contains

      !-------------------------------------------------------------------------------
      !> Summary: This subroutine is responsible for the I/O with the input file.
      !> Author: 
      !> Category: PKKprime, input-output
      !> Deprecated: False ! This needs to be set to True for deprecated subroutines
      !>
      !>  GIVEN a KEYWORD: CHARKEY it positions the
      !>  reading just after the CHARKEY if this
      !>  includes a '=', or ILINE lines after the
      !>  occurence of THE CHARKEY.
      !>  USAGE :
      !>  To read lmax include in the input card (ifile)
      !>
      !>      LMAX= 3      CORRECT!
      !>
      !>      or
      !>
      !>      LMAX         CORRECT!     (iline=1)
      !>       3
      !>    (without  the '=' )
      !>      LMAX
      !>   ---------                    (iline=2),etc
      !>       3
      !>  be carefull in this case to put the value after the
      !>  keyword example:
      !>
      !>     LMAX
      !>   3               WRONG!
      !>
      !> will NOT work
      !> Comments etc in the program are ignored.
      !>                                               1.6.99
      !>
      !> @warning 
      !> - The error handler is not working yet in all cases ...
      !> - In this version only files 500 lines long can be read in
      !> @endwarning
      !-------------------------------------------------------------------------------
      SUBROUTINE IOinput(CHARKEY,CHAR,ILINE,IFILE,IERROR)

      use mod_verify77, only : VERIFY77

      implicit none
      CHARACTER(len=10) :: CHARKEY
      CHARACTER CHAR*80
      INTEGER ILINE,IERROR,IFILE
      integer i,ios,ier,npt,ilen,ipos,ipos1,iklen
      CHARACTER STRING(500)*80
      CHARACTER STRING1*80
      CHARACTER ATEST
      character(len=*), parameter :: abc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890-_<>'
      integer, parameter :: NABC = len(abc)
      integer, parameter :: NCHAR = len(CHARKEY)

      IERROR = 0
      CHAR(1:50)='                                                   '
      OPEN(UNIT=ifile,status='OLD',FILE='inputFS',iostat=ios,err=2000)
         if (ios.gt.0) then
            write(6,*) "Error in reading the inputcard file"
            STOP
         end if

         npt = 1
         do
            read(ifile,FMT='(A80)',IOSTAT=ios) STRING(npt)
            IF(ios.lt.0.or.npt.ge.500) EXIT
            npt = npt + 1
         end do
          npt = npt - 1
!          write(6,*) 'LINES :',npt
          if (NPT.GE.500) WRITE(6,*)'Not all lines are read in from inputcard'

! 2 lines below where changed
!       ILEN = VERIFY(CHARKEY,ABC)
!       IKLEN= VERIFY(CHARKEY,' ')
! for linux
        CALL VERIFY77(NABC,ABC,NCHAR,CHARKEY,ILEN,IKLEN)
! for linux
!        write(6,*) CHARKEY(1:ILEN-1),ILEN,IKLEN
          IF(ILEN.LT.1) THEN 
             write(6,*) 'Input ERROR!'
             write(6,*) 'Cannot evaluate : ',CHARKEY
             write(6,*) 'IoInput is returning no value! '
             RETURN
          END IF
!
       DO i=1,NPT       ! loop in all line
         STRING1 = '   ' // STRING(I)     ! shift by 2 characters
          IPOS = INDEX(STRING1,CHARKEY(1:ILEN-1)) ! return the position of occurence
             if (ipos.ne.0) then
                 if (ipos.lt.4) then 
                  write(6,*) 'CONSISTENCY ERROR IOINPUT!'
                  STOP
                 end if
!                write(6,*) 'ipos is not zero',CHARKEY//'=','**'
                ipos1= INDEX(STRING1,CHARKEY(1:ILEN-1)//ACHAR(61))
                if (IPOS1.NE.0) then        ! return the string after 'CHARKEY='
                   CHAR = STRING1(ipos1+ilen:)
!                    write(6,*) CHARKEY,CHAR ! test
                   close(IFILE)
                   return
                else
! return the ILINE line below this CHARKEY
                  if (I+ILINE.LE.NPT) then
!                    write(6,*) IPOS,ILEN

                   CHAR = STRING(I+ILINE)(IPOS-3:)
                   if (ipos-4.gt.0) then    ! Changed on 28.01.2000
                   ATEST= STRING(I+ILINE)(IPOS-4:IPOS-3)
                   IF (ATEST.NE.' ') THEN
                   write(6,*) 'Possible ERROR !!!'
                   write(6,*) 'Parameter ',CHARKEY, ' maybe read in incorrectrly'
                   END IF
                   end if
!                   write(6,*) CHARKEY,CHAR ! test
                  close(ifile)
                  return
                  else
                  write(6,*) 'IoInput : No more lines in file '
                  end if
               end if
            end if
       END DO  ! i=1,npt
       IERROR = 1
       if (CHAR(1:20).eq.'                    ') then
       write(6,*) 'Parameter ........ ',CHARKEY , ' NOT found'
       write(6,*) 'Check your inputcard'
       end if 
       close(IFILE)
      RETURN
 2000 write(6,*) ' Error while reading..... ',CHARKEY
      write(6,*) ' Check your  inputcard ! '
      STOP
 1000 FORMAT(10I4)
 1001 FORMAT(3F12.9)
 1002 FORMAT(80A1)
 1003 FORMAT(10L4)
 1004 FORMAT(I4)
 1005 FORMAT(4I4)
      END SUBROUTINE IOinput

end module mod_ioinput
