module mod_wavefunctodisc
integer             :: first=1
integer,parameter   :: nlabel=300
integer             :: label(2,nlabel)
integer             :: icounter=1

contains

subroutine wavefunctodisc_write(wavefunction,cellnew,iatom,ispin,my_rank)
use type_wavefunction
use type_cellnew
use nrtype, only: wlength
implicit none

type(wavefunction_type) :: wavefunction
type(cell_typenew) :: cellnew

integer                 :: iatom
integer                 :: ispin
integer                 :: my_rank
!local
integer                 :: ilabel
character(len=20)       :: ctemp
if (first==1) then
  write(ctemp,'(I03.3)') my_rank
  open(unit=2343363,file='temp_rll_'    //trim(ctemp)//'.txt',access='direct',&
       recl=wlength*4*wavefunction%lmsize2*wavefunction%lmsize*cellnew%nrmaxnew,form='unformatted')
  open(unit=2343364,file='temp_sll_'    //trim(ctemp)//'.txt',access='direct',&
       recl=wlength*4*wavefunction%lmsize2*wavefunction%lmsize*cellnew%nrmaxnew,form='unformatted')
  open(unit=2343365,file='temp_rllleft_'//trim(ctemp)//'.txt',access='direct',&
       recl=wlength*4*wavefunction%lmsize2*wavefunction%lmsize*cellnew%nrmaxnew,form='unformatted')
  open(unit=2343366,file='temp_sllleft_'//trim(ctemp)//'.txt',access='direct',&
       recl=wlength*4*wavefunction%lmsize2*wavefunction%lmsize*cellnew%nrmaxnew,form='unformatted')
end if

call setlabel(iatom,ispin,ilabel)

if ( allocated(wavefunction%rll) ) then
  write(2343363,rec=ilabel) wavefunction%rll
  wavefunction%rll_saved=1
  deallocate(wavefunction%rll)
end if

if ( allocated(wavefunction%sll) ) then
  write(2343364,rec=ilabel) wavefunction%sll
  wavefunction%sll_saved=1
  deallocate(wavefunction%sll)
end if

if ( allocated(wavefunction%rllleft) ) then
  write(2343365,rec=ilabel) wavefunction%rllleft
  wavefunction%rllleft_saved=1
  deallocate(wavefunction%rllleft)
end if

if ( allocated(wavefunction%sllleft) ) then
  write(2343366,rec=ilabel) wavefunction%sllleft
  wavefunction%sllleft_saved=1
  deallocate(wavefunction%sllleft)
end if

first=0
end subroutine wavefunctodisc_write

subroutine wavefunctodisc_read(wavefunction,cellnew,iatom,ispin)
use type_wavefunction
use type_cellnew
implicit none

type(wavefunction_type) :: wavefunction
type(cell_typenew) :: cellnew
integer                 :: iatom
integer                 :: ispin
!local
integer                 :: ilabel

call getlabel(iatom,ispin,ilabel)

if ( .not. allocated(wavefunction%rll) .and. wavefunction%rll_saved==1) then
  allocate(wavefunction%rll(wavefunction%lmsize2,wavefunction%lmsize,cellnew%nrmaxnew,1))
  read(2343363,rec=ilabel) wavefunction%rll
end if
if ( .not. allocated(wavefunction%sll) .and. wavefunction%sll_saved==1) then
  allocate(wavefunction%sll(wavefunction%lmsize2,wavefunction%lmsize,cellnew%nrmaxnew,1))
  read(2343364,rec=ilabel) wavefunction%sll
end if
if ( .not. allocated(wavefunction%rllleft) .and. wavefunction%rllleft_saved==1) then
  allocate(wavefunction%rllleft(wavefunction%lmsize2,wavefunction%lmsize,cellnew%nrmaxnew,1))
  read(2343365,rec=ilabel) wavefunction%rllleft
end if
if ( .not. allocated(wavefunction%sllleft) .and. wavefunction%sllleft_saved==1) then
  allocate(wavefunction%sllleft(wavefunction%lmsize2,wavefunction%lmsize,cellnew%nrmaxnew,1))
  read(2343366,rec=ilabel) wavefunction%sllleft
end if

end subroutine wavefunctodisc_read


subroutine setlabel(iatom,ispin,ival)
integer :: iatom,ispin,ival
if (label(ispin,iatom)==0) then
  label(ispin,iatom)= icounter
  ival              = icounter
  icounter          = icounter+1
else
  ival=label(ispin,iatom)
end if

end subroutine setlabel

subroutine getlabel(iatom,ispin,ival)
integer :: iatom,ispin,ival
if (label(ispin,iatom)==0) then
  stop '[wavefunctodisc] getlabel: label error'
else
  ival=label(ispin,iatom)
end if

end subroutine getlabel



end module mod_wavefunctodisc
