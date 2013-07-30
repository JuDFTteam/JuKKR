! Some macros for checked allocation/deallocation
! they need an integer variable named memory_stat declared in each routine
! they are used.

#define CHECKALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define CHECKDEALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Deallocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)
#define DEALLOCATECHECK(X) deallocate(X, stat=memory_stat); CHECKDEALLOC(memory_stat)

#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; STOP; endif

! TODO: add reference to LatticeVectors and remove reference from RefClusters ???

module ClusterInfo_mod
  implicit none

  type ClusterInfo
    integer :: naclsd !< maximal number of cluster atoms
    integer :: naez_trc
    integer, dimension(:), allocatable :: nacls_trc
    integer, dimension(:), allocatable :: numn0_trc
    integer, dimension(:,:), allocatable :: indn0_trc
    integer, dimension(:,:), allocatable :: atom_trc
    integer, dimension(:,:), allocatable :: ezoa_trc
  end type

  public :: createClusterInfo_com
  private :: constructIndices

  CONTAINS

  !----------------------------------------------------------------------------
  !> Communicates and creates cluster info.
  !>
  !> Note: The cluster info determines the sparsity structure of
  !> the multiple scattering matrix
  !> The cluster information is communicated within each truncation zone
  !>
  !> @param ref_cluster_array    all the local ref. clusters
  subroutine createClusterInfo_com(self, ref_cluster_array, trunc_zone, communicator)
    use RefCluster_mod
    use TruncationZone_mod,  only: TruncationZone, translateInd
    use one_sided_commI_mod, only: copyFromI_com
    implicit none

    include 'mpif.h'

    type (ClusterInfo),    intent(inout) :: self
    type (RefCluster),     intent(in), dimension(:) :: ref_cluster_array
    type (TruncationZone), intent(in)    :: trunc_zone
    integer,               intent(in)    :: communicator

    integer :: ii
    integer :: ierr
    integer :: nacls_loc, nacls
    integer :: naclsd
    integer :: naez_trc
    integer :: num_local_atoms
    integer :: blocksize
    integer :: memory_stat
    integer, allocatable, dimension(:,:) :: buffer
    integer, allocatable, dimension(:,:) :: recv_buf

    integer, parameter :: MAGIC = 385306

    num_local_atoms = size(ref_cluster_array)

    ! find maximum for locally stored clusters
    nacls_loc = 0
    do ii = 1, num_local_atoms
      nacls_loc = max(nacls_loc, ref_cluster_array(ii)%nacls)
    end do

    ! determine maximal number of cluster atoms
    call MPI_Allreduce(nacls_loc, naclsd, 1, MPI_INTEGER, &
                       MPI_MAX, communicator, ierr)
    self%naclsd = naclsd

    naez_trc = trunc_zone%naez_trc
    self%naez_trc = naez_trc

    ALLOCATECHECK(self%nacls_trc(naez_trc))
    self%nacls_trc = 0
    ALLOCATECHECK(self%numn0_trc(naez_trc))
    self%numn0_trc = 0
    ALLOCATECHECK(self%indn0_trc(naez_trc, naclsd))
    self%indn0_trc = -1
    ALLOCATECHECK(self%atom_trc(naclsd, naez_trc))
    self%atom_trc = 0
    ALLOCATECHECK(self%ezoa_trc(naclsd, naez_trc))
    self%ezoa_trc = -1

    blocksize = 3*naclsd + 5
    ALLOCATECHECK(buffer(blocksize, num_local_atoms) )
    buffer = 0
    ALLOCATECHECK(recv_buf(blocksize, naez_trc) )

    do ii = 1, num_local_atoms
      nacls = ref_cluster_array(ii)%nacls
      CHECKASSERT(nacls <= naclsd)
      buffer(1, ii) = ref_cluster_array(ii)%atom_index
      buffer(2, ii) = ref_cluster_array(ii)%nacls
      buffer(3, ii) = ref_cluster_array(ii)%numn0
      !write(*,*) buffer(4             :(3 + nacls)          , ii), nacls
      buffer(4             :(3 + nacls)          , ii) = ref_cluster_array(ii)%indn0
      buffer((naclsd + 5)  :(naclsd + 4 + nacls) , ii) = ref_cluster_array(ii)%atom
      buffer((2*naclsd + 5):(2*naclsd+4 + nacls) , ii) = ref_cluster_array(ii)%ezoa
      buffer((3*naclsd + 5), ii) = MAGIC
    end do

    call copyFromI_com(recv_buf, buffer, trunc_zone%trunc2atom_index, &
                       blocksize, num_local_atoms, communicator)

    call constructIndices(self, trunc_zone, naez_trc, recv_buf, naclsd)

    DEALLOCATECHECK(recv_buf)
    DEALLOCATECHECK(buffer)

  end subroutine

  !----------------------------------------------------------------------------
  ! Helper routine
  subroutine constructIndices(self, trunc_zone, naez_trc, recv_buf, naclsd)
    use TruncationZone_mod,  only: TruncationZone, translateInd
    implicit none

    type (ClusterInfo), intent(inout) :: self
    type (TruncationZone), intent(in) :: trunc_zone
    integer, intent(in) :: naez_trc
    integer, intent(in), dimension(:,:) :: recv_buf
    integer, intent(in) :: naclsd

    integer :: ii
    integer, parameter :: MAGIC = 385306

    do ii = 1, naez_trc
      ! check if buffer from right atom was received
      CHECKASSERT( recv_buf(1, ii) == trunc_zone%trunc2atom_index(ii) )

      self%nacls_trc(ii) = recv_buf(2, ii)
      CHECKASSERT( self%nacls_trc(ii) <= naclsd .and. self%nacls_trc(ii) > 0 )
      self%numn0_trc(ii) = recv_buf(3, ii)
      CHECKASSERT( self%numn0_trc(ii) <= naclsd .and. self%numn0_trc(ii) > 0 )

      ! indn0 and atom have to be transformed to 'truncation-zone-indices'
      self%indn0_trc(ii,:) = translateInd(trunc_zone, &
                             recv_buf(4:(naclsd + 4), ii))

      self%atom_trc(:,ii)  = translateInd(trunc_zone, &
                             recv_buf((naclsd   + 5) : (2*naclsd + 4), ii))

      self%ezoa_trc(:,ii)  = recv_buf((2*naclsd + 5) : (3*naclsd + 4), ii)

      ! check if end of buffer is correct
      CHECKASSERT( recv_buf((3*naclsd + 5), ii) == MAGIC )
    end do
  end subroutine

end module ClusterInfo_mod
