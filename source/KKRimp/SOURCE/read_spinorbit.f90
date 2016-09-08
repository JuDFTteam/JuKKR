module mod_read_spinorbit
contains
subroutine read_spinorbit(natom,cellorbit,myrank)
  use nrtype
  use type_cellorbit
  use mod_version_info
  implicit none
  integer                           ::  natom
  integer                           ::  myrank
  type(cell_typeorbit)              ::  cellorbit
! local variables
integer                                  :: iatom,ios,ival

open(unit=1000001, file='kkrflex_spinorbitperatom', status='old', iostat=ios)
call version_check_header(1000001)
if (ios/=0) then
   write(1337,*) '[read_spinorbit] no kkrflex_spinorbitperatom exists'
   if (myrank==0) then
     write(*,*) '[read_spinorbit] no kkrflex_spinorbitperatom exists'
   end if
!set default value
do iatom = 1, natom
   cellorbit%use_spinorbit(iatom)=1
enddo
else 
  if (myrank==0) then
    print *, 'Using kkrflex_spinorbitperatom '
    print *, '---------------------------------'
  end if
  do iatom = 1, natom
    read(unit=1000001, fmt=*) ival 
    if (ival/=1 .and. ival/= 0) then
      print *, '[read_spinorbit] kkrflex_spinorbitperatom inconsistent'
    end if
    cellorbit%use_spinorbit(iatom)=ival
    if (myrank==0) then
      write(*,*) 'Atom ',iatom, ' : ',cellorbit%use_spinorbit(iatom)
    end if
  end do
end if
close(1000001)

end subroutine read_spinorbit


end module mod_read_spinorbit
