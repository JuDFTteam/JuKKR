module mod_logo

!           \_       \_   \_\_   \_\_\_\_  
!           \_       \_  \_ \_  \_ \_   \_ 
!           \_\_   \_\_\_\_ \_\_\_ \_\_\_  
!      \_   \_\_   \_\_  \_ \_  \_ \_  \_  
!       \_\_    \_\_ \_   \_\_   \_\_   \_


#define host     ! for use in host code
! #define impurity ! for use in KKRimpurity code
! #define scatter  ! for use in KKRscatter (zulapi) code
! #define KKprime  ! for use in KKprime (pkkr) code


private
public :: JuKKRlogo, print_logo, logo_line_width, logo_num_lines

implicit none

integer, parameter :: logo_line_width = 56
integer, parameter :: logo_num_lines = 10
character(len=logo_line_width):: JuKKRlogo(logo_num_lines)

JuKKRlogo = ['********************************************************', &
#ifdef host
          &  '*                  --- KKRhost ---                     *', &
#elif impurity
          &  '*                  --- KKRimp ---                      *', &
#elif scatter
          &  '*                 --- KKRscatter ---                   *', &
#elif KKprime
          &  '*                  --- KKprime ---                     *', &
#endif
          &  '*                         @                            *', &
          &  '*                                                      *', &
          &  '*                v        v  vv v  vv vvv              *', &
          &  '*                v        v v   v v   v  v             *', &
          &  '*                v v   v  vvv   vvv   vvv              *', &
          &  '*           v   v  v   v  v  v  v  v  v  v             *', &
          &  '*            vv     vv    v   v v   v v   v            *', &
          &  '********************************************************']

subroutine print_logo(logo, llen, nl)
  use mod_mympi, only: myrank, master
  implicit none
  !input
  integer, intent(in) :: llen, nl
  character(len=llen), intent(in):: logo(nl)
  !local
  integer :: i
  
  do i=1,nl
    ! write to std out and output file
    if(myrank==master) write(*,"(A)") logo(i)
    write(1337,"(A)") logo(i)
  end do
end subroutine print_logo

end module mod_logo