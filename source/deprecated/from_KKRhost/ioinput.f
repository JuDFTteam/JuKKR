      SUBROUTINE IOinput(CHARKEY,CHAR,ILINE,IFILE,IERROR)
c *********************************************************
c *  This subroutine is responsible for the I/O
c *  with the input file.
c *
c *  GIVEN a KEYWORD: CHARKEY it positions the
c *  reading just after the CHARKEY if this 
c *  includes a '=', or ILINE lines after the 
c *  occurence of THE CHARKEY.
c *  USAGE :
c *  To read lmax include in the input card (ifile)
c *
c *      LMAX= 3      CORRECT!          
c *
c *      or 
c *  
c *      LMAX         CORRECT!     (iline=1)
c *       3
c *    (without  the '=' )
c *      LMAX
c *   ---------                    (iline=2),etc
c *       3
c *  be carefull in this case to put the value after the 
c *  keyword example:
c *
c *     LMAX
c *   3               WRONG!
c *
c * will NOT work 
c * Comments etc in the program are ignored.
c *                                               1.6.99
c *
c * The error handler is not working yet in all cases ....
c * In this version only files 5000 lines long can be read in
c *******************************************************
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
c
      IERROR = 0
      IER = 0
      CHAR(1:50)='                                                   '
      OPEN(UNIT=ifile,status='OLD',FILE='inputcard',iostat=ios,
     &           err=2000)
         if (ios.gt.0) then
            write(6,*) "Error in reading the inputcard file"
            STOP
         end if
c
c
         npt = 1
         do
            read(ifile,FMT='(A256)',IOSTAT=ios) STRING(npt)
            IF(ios.lt.0.or.npt.ge.NLINIO) EXIT
            npt = npt + 1
         end do
          npt = npt - 1
c          write(6,*) 'LINES :',npt
          if (NPT.GE.NLINIO) 
     &         WRITE(1337,*)'Not all lines are read in from inputcard'

c 2 lines below where changed
c       ILEN = VERIFY(CHARKEY,ABC)
c       IKLEN= VERIFY(CHARKEY,' ')
c for linux
        CALL VERIFY77(NABC,ABC,NCHAR,CHARKEY,ILEN,IKLEN)
c for linux
c        write(6,*) CHARKEY(1:ILEN-1),ILEN,IKLEN
          IF(ILEN.LT.1) THEN 
             write(1337,*) 'Input ERROR!' 
             write(1337,*) 'Cannot evaluate : ',CHARKEY
             write(1337,*) 'IoInput is returning no value! '
             RETURN              
          END IF
c
       DO i=1,NPT       ! loop in all line
         STRING1 = '   ' // STRING(I)     ! shift by 2 characters   
          IPOS = INDEX(STRING1,CHARKEY(1:ILEN-1)) 
             ! return the position of occurence
             if (ipos.ne.0) then
                 if (ipos.lt.4) then 
                  write(6,*) 'CONSISTENCY ERROR IOINPUT!'
                  STOP
                 end if
c                write(6,*) 'ipos is not zero',CHARKEY//'=','**'
                ipos1= INDEX(STRING1,CHARKEY(1:ILEN-1)//ACHAR(61))
                if (IPOS1.NE.0) then
                   ! return the string after 'CHARKEY=' 
                   CHAR = STRING1(ipos1+ilen:)
c                    write(6,*) CHARKEY,CHAR ! test 
                   close(IFILE)
                   return                                      
                else
c return the ILINE line below this CHARKEY
                  if (I+ILINE.LE.NPT) then
c                    write(6,*) IPOS,ILEN
                   
                   CHAR = STRING(I+ILINE)(IPOS-3:)
                   if (ipos-4.gt.0) then    ! Changed on 28.01.2000
                   ATEST= STRING(I+ILINE)(IPOS-4:IPOS-3)
                   IF (ATEST.NE.' ') THEN
                   write(1337,*) 'Possible ERROR !!!'
                   write(1337,*) 'Parameter ',CHARKEY, 
     &                        ' maybe read in incorrectrly'
                   END IF
                   end if
c                   write(6,*) CHARKEY,CHAR ! test
                  close(ifile)
                  return 
                  else
                  write(1337,*) 'IoInput : No more lines in file '  
                  end if
               end if
            end if
       END DO  ! i=1,npt
       IER = 1
       IERROR = IERROR + IER
Cccc       if (CHAR(1:20).eq.'                    ') then
Cccc       write(6,*) 'Parameter ........ ',CHARKEY , ' NOT found'
Cccc       write(6,*) 'Check your inputcard'
Cccc       end if 
       close(IFILE)
      RETURN
 2000 write(6,*) ' Error while reading..... ',CHARKEY
      write(6,*) ' Check your  inputcard ! '
      STOP
      END

        SUBROUTINE VERIFY77(NABC,ABC,NCHAR,STR1,ipos1,ipos2)
        implicit none  
c This sub returns the position of the first space character
c in ipos2, and the position of the first letter in the string
c STR1
         INTEGER NCHAR,NABC
         CHARACTER (len=nchar) STR1 !*NCHAR
         CHARACTER (len=nabc) ABC !*NABC
         CHARACTER (len=1) CHAR !*1
         integer ipos,ipos1,ipos2,i,j
 !        DATA ABC/'ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890-_<>'/
          ipos2 =0
c
         ipos1 = INDEX(STR1,' ')
         do j=1,10
            char = str1(j:j+1)
c           write(6,*) 'char : ',j, char 
            ipos = 0
            do i=1,40
               ipos = INDEX(CHAR,ABC(I:I))
               if (IPOS.GT.0) THEN
                  ipos2 = j
                  RETURN
               end if 
            end do
            
         end do   
         RETURN
         END 

