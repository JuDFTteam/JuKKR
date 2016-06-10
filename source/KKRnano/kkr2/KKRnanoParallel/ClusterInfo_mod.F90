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

#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; STOP; endif

module ClusterInfo_mod
  implicit none
  private
  public :: ClusterInfo, create, destroy

  type ClusterInfo
    integer :: naclsd !< maximal number of cluster atoms
    integer :: naez_trc
    integer, allocatable :: nacls_trc(:)
    integer, allocatable :: numn0_trc(:)
    integer, allocatable :: indn0_trc(:,:)
    integer, allocatable :: atom_trc(:,:)
    integer, allocatable :: ezoa_trc(:,:)
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

    integer :: ii, jj, cnt, ind
    integer :: nacls, numn0, naclsd, naez_trc, num_local_atoms, blocksize
    integer :: memory_stat ! needed in allocatecheck
    integer :: ierr
    integer, allocatable :: send_buf(:,:), recv_buf(:,:)

    num_local_atoms = size(ref_clusters)

    nacls = maxval(ref_clusters(:)%nacls) ! find maximum for locally stored clusters
    ! determine global maximal number of cluster atoms
    call MPI_Allreduce(nacls, naclsd, 1, MPI_INTEGER, MPI_MAX, communicator, ierr)
    self%naclsd = naclsd

    naez_trc = trunc_zone%naez_trc
    self%naez_trc = naez_trc

    ALLOCATECHECK(self%nacls_trc(naez_trc))
    self%nacls_trc = 0
    ALLOCATECHECK(self%numn0_trc(naez_trc))
    self%numn0_trc = 0
    ALLOCATECHECK(self%indn0_trc(naez_trc,naclsd))
    self%indn0_trc = -1
    ALLOCATECHECK(self%atom_trc(naclsd,naez_trc))
    self%atom_trc = 0
    ALLOCATECHECK(self%ezoa_trc(naclsd,naez_trc))
    self%ezoa_trc = -1

    blocksize = 3*naclsd+4
    ALLOCATECHECK(send_buf(blocksize,num_local_atoms))
    send_buf(:,:) = -1

    do ii = 1, num_local_atoms
      nacls = ref_clusters(ii)%nacls
      numn0 = ref_clusters(ii)%numn0
      CHECKASSERT( nacls <= naclsd )
      send_buf(1,ii) = ref_clusters(ii)%atom_index ! global atom index
      send_buf(2,ii) = ref_clusters(ii)%nacls
      send_buf(3,ii) = ref_clusters(ii)%numn0
      send_buf(0*naclsd+3+1:0*naclsd+3+numn0,ii) = ref_clusters(ii)%indn0(:) ! indn0 is dim(numn0) now, however, numn0 <= nacls holds
      send_buf(1*naclsd+3+1:1*naclsd+3+nacls,ii) = ref_clusters(ii)%atom(:)
      send_buf(2*naclsd+3+1:2*naclsd+3+nacls,ii) = ref_clusters(ii)%ezoa(:)
      send_buf(3*naclsd+3+1,ii) = MAGIC
    enddo ! ii

    ALLOCATECHECK(recv_buf(blocksize,naez_trc))

    call copyFromI_com(recv_buf, send_buf, trunc_zone%trunc2atom_index, blocksize, num_local_atoms, communicator)

    do ii = 1, naez_trc
      CHECKASSERT( recv_buf(1,ii) == trunc_zone%trunc2atom_index(ii) ) ! check if send_buf from right atom was received

      numn0 = recv_buf(3,ii)
      ! indn0 and atom have to be transformed to 'truncation-zone-indices'
      cnt = 0
      do jj = 1, numn0
! ! ! ! write(*,'(9(a,i0))') __FILE__,__LINE__,' ii=',ii,' jj=',jj,' ind=',recv_buf(3+jj,ii)   
!       ind = translateInd(trunc_zone, recv_buf(0*naclsd+3+jj,ii)) ! indn0 received
        ind = trunc_zone%index_map(recv_buf(0*naclsd+3+jj,ii)) ! indn0 received

        ! ind = -1 means that this atom is outside of truncation zone
        if (ind > 0) then
          cnt = cnt + 1
          self%indn0_trc(ii,cnt) = ind
        endif ! ind > 0
      enddo ! jj

      CHECKASSERT( 0 < cnt .and. cnt <= naclsd )
      self%numn0_trc(ii) = cnt

      nacls = recv_buf(2,ii)
      cnt = 0
      do jj = 1, nacls
! ! ! ! write(*,'(9(a,i0))') __FILE__,__LINE__,' ii=',ii,' jj=',jj,' ind=',recv_buf(naclsd+3+jj,ii)
!       ind = translateInd(trunc_zone, recv_buf(1*naclsd+3+jj,ii)) ! atom received
        ind = trunc_zone%index_map(recv_buf(1*naclsd+3+jj,ii)) ! atom received

        if (ind > 0) then
          cnt = cnt + 1
          self%atom_trc(cnt,ii) = ind
        endif ! ind > 0
      enddo ! jj

      CHECKASSERT( 0 < cnt .and. cnt <= naclsd )
      self%nacls_trc(ii) = cnt

! ! ! write(*,'(9(a,i0))') __FILE__,__LINE__,' naclsd=',naclsd

      self%ezoa_trc(:,ii) = recv_buf(2*naclsd+4:3*naclsd+3,ii)

      CHECKASSERT( recv_buf(3*naclsd+4,ii) == MAGIC ) ! check if end of send_buf is correct
    enddo ! ii

    DEALLOCATECHECK(recv_buf)
    DEALLOCATECHECK(send_buf)
  endsubroutine ! create
  
  elemental subroutine destroyClusterInfo(self)
    type(ClusterInfo), intent(inout) :: self

    integer :: ist ! ignore status
    deallocate(self%nacls_trc, self%numn0_trc, self%indn0_trc, self%atom_trc, self%ezoa_trc, stat=ist)
    self%naclsd = 0
    self%naez_trc = 0
  endsubroutine ! destroy

endmodule ! ClusterInfo_mod
