module mod_log
contains
   !-------------------------------------------------------------------------------
   !> Summary: Logs date and time
   !> Author: 
   !> Category: KKRimp, input-output
   !>
   !-------------------------------------------------------------------------------
 subroutine log_write(string1)
 use mod_types, only: t_inc
 implicit none
        character(len=*)     :: string1
        character(len=22)    :: prefix
        character(len=10)    :: time2, date2
        call DATE_AND_TIME(date=date2,time=time2)
        prefix     =     date2(7:8) // '.' // date2(5:6) // '.' // date2(1:4)// ' ' // & 
                           time2(1:2) // ":" // time2(3:4) // ":" // time2(5:6) // ":   "
!         write(1337,'(2A)') prefix 
        if (t_inc%i_write>0) write(1337,'(2A)') prefix, string1
!         write(1337,'(2A)') prefix

end subroutine log_write



end module mod_log
