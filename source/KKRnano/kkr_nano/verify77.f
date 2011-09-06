        SUBROUTINE VERIFY77(STR1,ipos1,ipos2)
        implicit none  
c This sub returns the position of the first space character
c in ipos2, and the position of the first letter in the string
c STR1
c
        CHARACTER STR1*10
        CHARACTER ABC*37
        CHARACTER CHAR*1
        integer ipos,ipos1,ipos2,i,j
        DATA ABC/'ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890-'/
         ipos2 =0
c
         ipos1 = INDEX(STR1,' ')
         do j=1,10
            char = str1(j:j+1)
c           write(6,*) 'char : ',j, char 
            ipos = 0
            do i=1,37
               ipos = INDEX(CHAR,ABC(I:I))
               if (IPOS.GT.0) THEN
                  ipos2 = j
                  RETURN
               end if 
            end do
            
         end do   
         RETURN
         END 

