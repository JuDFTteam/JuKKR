module mod_read_spinorbit
contains
subroutine read_spinorbit(natom,cellnew,myrank)
  use nrtype
  use type_cellnew
  implicit none
  integer                           ::  natom
  type(cell_typenew)                 ::  cellnew(natom)
  integer                           ::  myrank
! local variables
integer                                  :: iatom,ios,ival


open(unit=1000001, file='kkrflex_spinorbitperatom', status='old', iostat=ios)
if (ios/=0) then
   write(1337,*) '[read_spinorbit] no kkrflex_spinorbitperatom exists'
   write(1337,*) '[read_spinorbit] setting cellnew(1:natom)%use_spinorbit=1 for all atoms'
   cellnew(1:natom)%use_spinorbit=1
   if (myrank==0) then
     write(*,*) '[read_spinorbit] no kkrflex_spinorbitperatom exists'
   end if
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
    cellnew(iatom)%use_spinorbit=ival
    if (myrank==0) then
      write(*,*) 'Atom ',iatom, ' : ',cellnew(iatom)%use_spinorbit
    end if
  end do
end if
close(1000001)

end subroutine read_spinorbit


end module mod_read_spinorbit
