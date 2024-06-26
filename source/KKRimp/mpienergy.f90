module mod_mpienergy
contains

   !-------------------------------------------------------------------------------
   !> Summary: Tries to distribute the energy points equally to all processors
   !> Author:
   !> Category: KKRimp, communication
   !> Deprecated: False ! This needs to be set to True for deprecated subroutines
   !>
   !-------------------------------------------------------------------------------
subroutine mpienergy_distribute(myrank,numbproc,ielast,mpi_iebounds)
! ##############################################################
! This routine tries to distribute the energy points equally
! to all processors used in parallel MPI calculations
! ##############################################################
  use mod_types, only: t_inc

!interface
 implicit none
integer,intent(in)          :: myrank
integer,intent(in)                     :: numbproc
integer,intent(in)                     :: ielast
integer,allocatable         :: mpi_iebounds(:,:)
!local
integer                  :: mpi_ielast(0:numbproc)
integer                  :: iestart,iestop,irank,ieleft

allocate( mpi_iebounds(2,0:numbproc-1) )
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


if (t_inc%i_write>0) write(1337,*) '###################################'
if (t_inc%i_write>0) write(1337,*) '###  mpi_energy distribution    ###'
if (t_inc%i_write>0) write(1337,*) '###################################'
if (t_inc%i_write>0) write(1337,*) 'number of threads       :',numbproc
if (t_inc%i_write>0) write(1337,*) 'number of energy points :',ielast
if (t_inc%i_write>0) write(1337,*) 'Distribution to threads:'
Do irank=0,numbproc-1
  if (t_inc%i_write>0) write(1337,'(4(A,I4))') ' Thread',irank,' IE=',mpi_iebounds(1,irank),'..',mpi_iebounds(2,irank),' Total: ',mpi_ielast(irank)
end do
end subroutine mpienergy_distribute

end module mod_mpienergy
