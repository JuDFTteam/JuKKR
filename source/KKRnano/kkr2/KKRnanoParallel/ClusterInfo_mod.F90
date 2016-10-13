!------------------------------------------------------------------------------
!> Module for collecting information about all reference clusters of atoms
!> in truncation zone.
!>
!> The information collected determines the sparsity structure of the
!> coefficient matrix used in solving the multiple scattering problem.
!>
!> @author Elias Rabel


! Some macros for checked allocation/deallocation
! they need an integer variable named memory_stat declared in each routine
! they are used.

#define CHECKALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define CHECKDEALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Deallocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)
#define DEALLOCATECHECK(X) deallocate(X, stat=memory_stat); CHECKDEALLOC(memory_stat)

module ClusterInfo_mod
#include "macros.h"
  implicit none
  private
  public :: ClusterInfo, create, destroy

  type ClusterInfo
    integer :: naclsd !< maximal number of cluster atoms
    integer :: naez_trc
    integer, allocatable :: nacls_trc(:) !> dim(naez_trc) number of target atoms in local interaction cluster
    integer, allocatable :: numn0_trc(:) !> dim(naez_trc) number of inequivent target atoms
    integer(kind=2), allocatable :: indn0_trc(:,:) !> dim(naclsd,naez_trc) ! can be stored in 16bit since it encodes only cluster indices
    integer(kind=2), allocatable ::  atom_trc(:,:) !> dim(naclsd,naez_trc)
    integer(kind=2), allocatable ::  ezoa_trc(:,:) !> dim(naclsd,naez_trc)
  endtype

  interface create
    module procedure createClusterInfo
  endinterface
  
  interface destroy
    module procedure destroyClusterInfo
  endinterface
  
  integer, parameter, private :: MAGIC = 385306
  
  contains

  !----------------------------------------------------------------------------
  !> Communicates and creates cluster info.
  !>
  !> Note: The cluster info determines the sparsity structure of
  !> the multiple scattering matrix
  !> The cluster information is communicated within each truncation zone
  !>
  !> @param ref_clusters    all the local ref. clusters
  subroutine createClusterInfo(self, ref_clusters, trunc_zone, communicator)
    use RefCluster_mod, only: RefCluster
    use TruncationZone_mod, only: TruncationZone
    use one_sided_commI_mod, only: copyFromI_com

    include 'mpif.h'

    type(ClusterInfo), intent(inout)    :: self
    type(RefCluster), intent(in)        :: ref_clusters(:)
    type(TruncationZone), intent(in)    :: trunc_zone
    integer, intent(in)                 :: communicator

    integer :: ii, jj, cnt, ila
    integer :: OFFSET_INDN, OFFSET_ATOM, OFFSET_EZOA
    integer :: nacls, numn0, naclsd, naez_trc, num_local_atoms, blocksize
    integer :: memory_stat ! needed in allocatecheck
    integer :: ierr
    integer(kind=2) :: ezoa, ind ! local atom index
    integer(kind=4) :: atom_id, atom, indn0 ! global atom ids
    integer(kind=4), allocatable :: send_buf(:,:), recv_buf(:,:) ! global atom ids, require more than 16bit
#ifdef DEBUG
    integer(kind=4), allocatable :: global_target_atom_id(:) !!! DEBUG
#endif
    num_local_atoms = size(ref_clusters)

    nacls = maxval(ref_clusters(:)%nacls) ! find maximum for locally stored clusters
    ! determine global maximal number of cluster atoms
    call MPI_Allreduce(nacls, naclsd, 1, MPI_INTEGER, MPI_MAX, communicator, ierr)
    ! here, additional alignment could be imposed into naclsd, e.g. this: 
    !   naclsd = 4*((naclsd - 1)/4 + 1)
    self%naclsd = naclsd

    naez_trc = trunc_zone%naez_trc
    self%naez_trc = naez_trc

    OFFSET_INDN = 0*naclsd + 3
    OFFSET_ATOM = 1*naclsd + 3
    OFFSET_EZOA = 2*naclsd + 3
    blocksize   = 3*naclsd + 3 + 1
    
    ALLOCATECHECK(send_buf(blocksize,num_local_atoms))
    send_buf(:,:) = -1 ! init

    do ila = 1, num_local_atoms
      nacls = ref_clusters(ila)%nacls
      numn0 = ref_clusters(ila)%numn0
      CHECKASSERT( nacls <= naclsd )
      CHECKASSERT( numn0 <= nacls )
      send_buf(1,ila) = ref_clusters(ila)%atom_id ! global atom index
      send_buf(2,ila) = nacls
      send_buf(3,ila) = numn0
      send_buf(OFFSET_INDN + 1:numn0 + OFFSET_INDN,ila) = ref_clusters(ila)%indn0(:) ! indn0 has dim(numn0) now, however, numn0 <= nacls holds
      send_buf(OFFSET_ATOM + 1:nacls + OFFSET_ATOM,ila) = ref_clusters(ila)%atom(:)
      send_buf(OFFSET_EZOA + 1:nacls + OFFSET_EZOA,ila) = ref_clusters(ila)%ezoa(:)
      send_buf(blocksize,ila) = MAGIC
    enddo ! ila

    ALLOCATECHECK(recv_buf(blocksize,naez_trc))

    call copyFromI_com(recv_buf, send_buf, trunc_zone%global_atom_id, blocksize, num_local_atoms, communicator)

    DEALLOCATECHECK(send_buf)
#ifdef DEBUG
    ALLOCATECHECK(global_target_atom_id(naclsd)) !!! DEBUG
#endif
    ALLOCATECHECK(self%nacls_trc(naez_trc))
    ALLOCATECHECK(self%numn0_trc(naez_trc))
    ALLOCATECHECK(self%indn0_trc(naclsd,naez_trc))
    ALLOCATECHECK(self%atom_trc(naclsd,naez_trc))
    ALLOCATECHECK(self%ezoa_trc(naclsd,naez_trc))
    self%nacls_trc = 0
    self%numn0_trc = 0
    self%indn0_trc = -1
    self%atom_trc = -1
    self%ezoa_trc = -1

    do ii = 1, naez_trc
      atom_id = recv_buf(1,ii)
      CHECKASSERT( atom_id == trunc_zone%global_atom_id(ii) ) ! check if send_buf from right atom was received

      numn0 = recv_buf(3,ii)
      ! indn0 and atom have to be transformed to 'truncation-zone-indices'
      cnt = 0
      do jj = 1, numn0

        indn0 = recv_buf(OFFSET_INDN + jj,ii) ! indn0 received
! ! ! ! write(*,'(9(a,i0))') __FILE__,__LINE__,' ii=',ii,' jj=',jj,' ind=',indn0
        ind = trunc_zone%local_atom_idx(indn0)

        if (ind > 0) then ! ind == -1 means that this atom is outside of truncation zone
          cnt = cnt + 1
          self%indn0_trc(cnt,ii) = ind ! indn0 translated into local indices of the trunc_zone
          ! if trunc_zone%global_atom_id(:) is in ascending order (strictly monotonous), also indn0_trc will get that property
#ifdef DEBUG
          global_target_atom_id(cnt) = indn0 !!! DEBUG
#endif         
        endif ! ind > 0
      enddo ! jj

      CHECKASSERT( 0 < cnt .and. cnt <= numn0 )
      self%numn0_trc(ii) = cnt
#ifdef DEBUG
      write(*,'(999(a,i0))') __FILE__,__LINE__,' atom_id=',atom_id,' numn0_trc=',cnt,' indn0= ',(global_target_atom_id(jj),' ', jj=1,cnt) !!! DEBUG
#endif

      nacls = recv_buf(2,ii)
      cnt = 0
      do jj = 1, nacls

        atom = recv_buf(OFFSET_ATOM + jj,ii) ! atom received
        ezoa = recv_buf(OFFSET_EZOA + jj,ii) ! ezoa received (proposed new code)
! ! ! ! write(*,'(9(a,i0))') __FILE__,__LINE__,' ii=',ii,' jj=',jj,' ind=',atom
        ind = trunc_zone%local_atom_idx(atom)

        if (ind > 0) then ! ind == -1 means that this atom is outside of truncation zone
          cnt = cnt + 1
          self%atom_trc(cnt,ii) = ind ! atom translated into local indices of the trunc_zone
          self%ezoa_trc(cnt,ii) = ezoa ! does not need translation (proposed new code)
#ifdef DEBUG
          global_target_atom_id(cnt) = atom !!! DEBUG
#endif          
        endif ! ind > 0
      enddo ! jj

      CHECKASSERT( 0 < cnt .and. cnt <= nacls )
      self%nacls_trc(ii) = cnt
#ifdef DEBUG
      write(*,'(999(a,i0))') __FILE__,__LINE__,' atom_id=',atom_id,' nacls_trc=',cnt,' atom@ezoa= ',(global_target_atom_id(jj),'@',self%ezoa_trc(jj,ii),' ', jj=1,cnt) !!! DEBUG
#endif
! ! ! write(*,'(9(a,i0))') __FILE__,__LINE__,' naclsd=',naclsd

!!    ! ToDo: discuss if we have to treat ezoa_trc in sync with atom_trc concerning the (ind == -1) case
!!    self%ezoa_trc(1:nacls,ii) = recv_buf(OFFSET_EZOA + 1:nacls + OFFSET_EZOA,ii) ! (old code)
#ifdef DEBUG
      write(*,'(999(a,i0))') __FILE__,__LINE__,' atom_id=',atom_id,' old nacls=',cnt,' atom@ezoa= ',(global_target_atom_id(jj),'@',recv_buf(OFFSET_EZOA+jj,ii),' ', jj=1,cnt) !!! DEBUG
#endif
      CHECKASSERT( recv_buf(blocksize,ii) == MAGIC ) ! check if the end of send_buf seems correct
    enddo ! ii
    
#ifdef DEBUG
    DEALLOCATECHECK(global_target_atom_id) !!! DEBUG
#endif    
    DEALLOCATECHECK(recv_buf)
  endsubroutine ! create
  
  elemental subroutine destroyClusterInfo(self)
    type(ClusterInfo), intent(inout) :: self

    integer :: ist ! ignore status
    deallocate(self%nacls_trc, self%numn0_trc, self%indn0_trc, self%atom_trc, self%ezoa_trc, stat=ist)
    self%naclsd = 0
    self%naez_trc = 0
  endsubroutine ! destroy

endmodule ! ClusterInfo_mod
