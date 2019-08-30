module mod_read_atominfo
contains
subroutine read_atominfo(cmode,cfilename_atom,natom,ntotatom,ratom,zatom,lmaxd,lmaxatom,killatom,vtotatom)
  use nrtype
  use mod_version_info
  use mod_types, only: t_inc
  implicit none
! the routine reads the atom information using 2 modes
! - 'total' reads the information of all atoms, also the once which are beiing killed. 
!           Arrays have dimension NTOTATOM
! - 'imp'  reads information only for impurity atoms
!           Array dimension is NATOM

! interface variables
character(len=*),          intent(in)    :: cmode
character(len=*),          intent(in)    :: cfilename_atom    ! filename of the atominfo file
integer,                   intent(out)   :: natom             ! number of impurity atoms
integer,                   intent(out)   :: ntotatom          ! total number of atoms impurity + surrounding
real(kind=dp),allocatable, intent(out)   :: ratom(:,:)     ! spacial positions
real(kind=dp),allocatable, intent(out)   :: zatom(:)          ! core charge of each atom
integer,                   intent(out)   :: lmaxd             ! maximum number of lmaxd
integer,allocatable,       intent(out)   :: lmaxatom(:)       ! lmax of each individual atoms
! real(kind=dp),allocatable, intent(out)   :: rmt(:)            ! muffin tin radius of each atom
integer,allocatable,       intent(out)   :: killatom(:)     ! 1= is screening atom 0= is not
integer,allocatable                      :: vtotatom(:)     ! 1= is screening atom 0= is not

! local variables
integer                                  :: ifile_atom,ios,numb
integer                                  :: iatom, jatom, icountvatom
real(kind=dp)                            :: temp1,temp2,temp3,temp4
integer                                  :: temp5,temp6,temp7
character(len=200)  ::string1

!--------------------------------------------------------
!--  open the file                        --
!--------------------------------------------------------
ifile_atom=1000000 
icountvatom=0
open(unit=ifile_atom, file=cfilename_atom, status='old', iostat=ios)
call version_check_header(ifile_atom)

  if (ios/=0) then
     write(*,*) '[read_atominfo] config file does not exist'
     stop
  end if !ios/=0


string1=this_readline(ifile_atom,ios)
read(unit=string1,fmt=*) ntotatom,natom
! write(*,*) ntotatom,natom
if (cmode=='total') then
  numb=ntotatom
elseif (cmode=='imp') then
  numb=ntotatom !Phivos 2.6.14: before this was natom, too small in case of "killed" atoms
end if


allocate(lmaxatom(numb), killatom(numb), ratom(3,numb), &
         zatom(numb),vtotatom(numb) )
! Initialize
lmaxatom(:) = 0
killatom(:) = 0
ratom(:,:) = 0.d0
zatom(:) = 0.d0
vtotatom(:) = 0

!--------------------------------------------------------
!--  read in all impurity atoms                        --
!--------------------------------------------------------
jatom=1
do iatom=1,ntotatom
   string1=this_readline(ifile_atom,ios)
   if (ios/=0) stop '[read_atominfo] Error reading atom info2'
   if (ios==-1) stop '[read_atominfo] EOF'
   if (cmode=='total') then
     read(unit=string1,iostat=ios,fmt=*) ratom(1,iatom),ratom(2,iatom),ratom(3,iatom), &
                                         zatom(iatom),vtotatom(iatom),killatom(iatom),lmaxatom(iatom)
     if (killatom(iatom)==1) icountvatom=icountvatom+1
   elseif (cmode=='imp') then
     read(unit=string1,iostat=ios,fmt=*) temp1,temp2,temp3,temp4,temp5,temp6,temp7
     if (temp6/=1) then !if killatom is not 1
       read(unit=string1,iostat=ios,fmt=*) ratom(1,jatom),ratom(2,jatom),ratom(3,jatom), &
                                         zatom(jatom),vtotatom(jatom),killatom(jatom),lmaxatom(jatom)
     jatom=jatom+1
     else
       icountvatom=icountvatom+1
     end if
   end if
!    write(*,*) iatom,ios
   if (ios/=0) stop '[read_atominfo] Error reading atom info'
end do !iatom
! if (ntotatom-icountvatom/=natom) stop '[read_atominfo] error defining virtual atoms'
close(ifile_atom)

!--------------------------------------------------------
!--  calculate the maximum lm cut-off                  --
!--------------------------------------------------------
lmaxd = maxval(lmaxatom)


!--------------------------------------------------------
!--  write out information                             --
!--------------------------------------------------------
if (t_inc%i_write>0) write(1337,*) ''
if (t_inc%i_write>0) write(1337,*) '-------------------------------------------'
if (t_inc%i_write>0) write(1337,*) '-----      read atominfo               ----'
if (t_inc%i_write>0) write(1337,*) '-------------------------------------------'
if (t_inc%i_write>0) write(1337,*) 'NATOM is    ', natom 
if (t_inc%i_write>0) write(1337,*) 'NTOTATOM is ', ntotatom 
if (t_inc%i_write>0) write(1337,*) 'LMAXD is    ', lmaxd
if (t_inc%i_write>0) write(1337,*) '-------------------------------------------'
if (t_inc%i_write>0) write(1337,'(2A)') '  No.   x        y      z',&
                  '            Z    VATOM?  KILLATOM?  LMAX'
if (t_inc%i_write>0) write(1337,*) '-------------------------------------------'
do iatom=1,numb
   if (t_inc%i_write>0) write(1337,'(I4,8g10.3,F10.2,4I6)') iatom, ratom(:,iatom), zatom(iatom), &
                                      vtotatom(iatom),killatom(iatom),lmaxatom(iatom)
end do !iatom
if (t_inc%i_write>0) write(1337,*) '-------------------------------------------'


end subroutine read_atominfo

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


end module mod_read_atominfo
