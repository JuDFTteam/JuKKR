module mod_mpienergy
contains

subroutine mpienergy_distribute(myrank,numbproc,ielast,mpi_iebounds)
! ##############################################################
! This routine tries to distribute the energy points equally
! to all processors used in parallel MPI calculations
! ##############################################################

!interface
 implicit none
integer,intent(in)          :: myrank
integer,intent(in)                     :: numbproc
integer,intent(in)                     :: ielast
integer,allocatable         :: mpi_iebounds(:,:)
!local
integer                  :: mpi_ielast(0:numbproc)
integer                  :: ie,iestart,iestop,irank,ranktemp,ieleft


! kkrsusc
if (.not.allocated(mpi_iebounds)) then
  allocate( mpi_iebounds(2,0:numbproc-1) )
end if
! mpi_ielast=0

do irank=0,numbproc-1
  mpi_ielast(irank)=ielast/numbproc
end do

ieleft=mod(ielast,numbproc)
do irank=0,ieleft-1
  mpi_ielast(irank)=mpi_ielast(irank)+1
end do

!   do irank=0,numbproc-1
! write(*,*) irank,mpi_ielast(irank)
! end do
! stop
! ieleft=ielast;ie=0
! do while(ieleft/=0)
!   ranktemp = mod(ie,numbproc)
!   mpi_ielast(ranktemp)=mpi_ielast(ranktemp)+1
!   ieleft= ieleft-1
!   ie=ie+1
! end do

! do ie=1,ielast
!   ranktemp = mod(ie,numbproc) + 1
!   mpi_ielast(ranktemp) = ie
! end do

! iestart=1-mpi_ielast(0)
iestop = 0
do irank=0,numbproc-1
 iestart = iestop+1
 iestop  = iestop+mpi_ielast(irank)
 mpi_iebounds(1,irank)= iestart
 mpi_iebounds(2,irank)= iestop
end do

!   do irank=0,numbproc-1
! write(*,*) irank,mpi_iebounds(1,irank),mpi_iebounds(2,irank)
! end do
! stop


write(1337,*) '###################################'
write(1337,*) '###  mpi_energy distribution    ###'
write(1337,*) '###################################'
write(1337,*) 'number of threads       :',numbproc
write(1337,*) 'number of energy points :',ielast
write(1337,*) 'Distribution to threads:'
Do irank=0,numbproc-1
  write(1337,'(4(A,I4))') ' Thread',irank,' IE=',mpi_iebounds(1,irank),'..',mpi_iebounds(2,irank),' Total: ',mpi_ielast(irank)
end do
end subroutine mpienergy_distribute

end module mod_mpienergy
