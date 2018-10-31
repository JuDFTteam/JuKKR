module mod_change_nrmin

contains 

subroutine change_nrmin(natom,cell)
use type_cell
implicit none

integer  :: natom
type(cell_type)  :: cell(natom)

integer  :: iatom

do iatom=1,natom
  write(*,*) 'Atom ',iatom
  write(*,*) 'current NRNS is',cell(iatom)%NRNS
  write(*,*) 'current NRMIN_NS is',cell(iatom)%nrmin_ns
  write(*,*) 'current RMIN_NS is',cell(iatom)%rmesh(cell(iatom)%nrmin_ns)
  write(*,*) 'new NRNS ? '
  read(*,*)  cell(iatom)%NRNS
  cell(iatom)%nrmin_ns=cell(iatom)%nrmax-cell(iatom)%nrns
  write(*,*) 'new NRNS is',cell(iatom)%NRNS
  write(*,*) 'new NRMIN_NS is',cell(iatom)%nrmin_ns
  write(*,*) 'new RMIN_NS is',cell(iatom)%rmesh(cell(iatom)%nrmin_ns)
end do

end subroutine change_nrmin











end module mod_change_nrmin