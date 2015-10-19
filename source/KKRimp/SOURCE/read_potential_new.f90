module read_potentialnew
integer                 ::      ifile=23425566






subroutine read_label()


line=this_readline(ifile,ios)
if line=='<'///>

end subroutine read_label


function this_readline(ifile,ios)
!--------------------------------------------------------
!--  reads the next line in file unit ifile            --
!--------------------------------------------------------
!--  files starting with a dash (#) are treated as     --
!--  a comment !!!                                     --
!--  OUTPUT: next line which is not commented out      --
!--          error variable IOS (should be zero)       --
!--------------------------------------------------------
! input variables
  implicit none
integer,intent(in)               :: ifile
integer,intent(out)              :: ios
! local variables
character(len=200)  ::this_readline
do
  read(unit=ifile,fmt='(A)', iostat=ios) this_readline
  if (ios/=0 .or. this_readline(1:1)/='#') exit
end do
end function this_readline

end module read_potentialnew
