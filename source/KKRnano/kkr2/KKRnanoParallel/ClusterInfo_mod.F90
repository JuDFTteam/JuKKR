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
  public :: createClusterInfo_com, destroyClusterInfo_do_nothing

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
    module procedure createClusterInfo_com
  endinterface
  
  interface destroy
    module procedure destroyClusterInfo_do_nothing
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
  !> @param ref_cluster_array    all the local ref. clusters
  subroutine createClusterInfo_com(self, ref_cluster_array, trunc_zone, communicator)
    use RefCluster_mod, only: RefCluster
    use TruncationZone_mod,  only: TruncationZone
    use one_sided_commI_mod, only: copyFromI_com

    include 'mpif.h'

    type(ClusterInfo),    intent(inout) :: self
    type(RefCluster),     intent(in)    :: ref_cluster_array(:)
    type(TruncationZone), intent(in)    :: trunc_zone
    integer,              intent(in)    :: communicator

    integer :: ii
    integer :: ierr
    integer :: nacls, numn0
    integer :: naclsd
    integer :: naez_trc
    integer :: num_local_atoms
    integer :: blocksize
    integer :: memory_stat
    integer, allocatable :: buffer(:,:)
    integer, allocatable :: recv_buf(:,:)


    num_local_atoms = size(ref_cluster_array)

    ! find maximum for locally stored clusters
    nacls = 0
    do ii = 1, num_local_atoms
      nacls = max(nacls, ref_cluster_array(ii)%nacls)
    enddo ! ii

    ! determine maximal number of cluster atoms
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
    ALLOCATECHECK(buffer(blocksize, num_local_atoms) )
    buffer(:,:) = -1
    ALLOCATECHECK(recv_buf(blocksize, naez_trc) )

    do ii = 1, num_local_atoms
      nacls = ref_cluster_array(ii)%nacls
      numn0 = ref_cluster_array(ii)%numn0
      CHECKASSERT( nacls <= naclsd )
!!!!  buffer(:,:) = -1 ! introducing this line leads to an error -- suspicious
      buffer(1,ii) = ref_cluster_array(ii)%atom_index
      buffer(2,ii) = ref_cluster_array(ii)%nacls
      buffer(3,ii) = ref_cluster_array(ii)%numn0
      buffer(4         :         3+numn0,ii) = ref_cluster_array(ii)%indn0 ! indn0 is dim(numn0) now, however, numn0 <= nacls holds
      buffer(4+naclsd  :  naclsd+3+nacls,ii) = ref_cluster_array(ii)%atom
      buffer(4+2*naclsd:2*naclsd+3+nacls,ii) = ref_cluster_array(ii)%ezoa
      buffer(4+3*naclsd,ii) = MAGIC
    enddo ! ii

    call copyFromI_com(recv_buf, buffer, trunc_zone%trunc2atom_index, blocksize, num_local_atoms, communicator)

    call constructIndices(self, trunc_zone, naez_trc, recv_buf, naclsd)

    DEALLOCATECHECK(recv_buf)
    DEALLOCATECHECK(buffer)

  endsubroutine ! create

  
  subroutine destroyClusterInfo_do_nothing(self)
    type(ClusterInfo), intent(inout) :: self
    ! nothing to be done
  endsubroutine ! destroy

  !----------------------------------------------------------------------------
  !> Helper routine.
  subroutine constructIndices(self, trunc_zone, naez_trc, recv_buf, naclsd)
    use TruncationZone_mod, only: TruncationZone, translateInd

    type(ClusterInfo), intent(inout) :: self
    type(TruncationZone), intent(in) :: trunc_zone
    integer, intent(in) :: naez_trc
    integer, intent(in) :: recv_buf(:,:)
    integer, intent(in) :: naclsd

    integer :: ii, jj, counter, ind

    do ii = 1, naez_trc
      ! check if buffer from right atom was received
      CHECKASSERT( recv_buf(1,ii) == trunc_zone%trunc2atom_index(ii) )

      ! indn0 and atom have to be transformed to 'truncation-zone-indices'
      counter = 0
      do jj = 1, naclsd
! ! ! ! write(*,'(9(a,i0))') __FILE__,__LINE__,' ii=',ii,' jj=',jj,' ind=',recv_buf(3+jj,ii)   
        ind = translateInd(trunc_zone, recv_buf(3+jj,ii))

        ! ind = -1 means that atom is outside of truncation zone
        if (ind > 0) then
          counter = counter + 1
          self%indn0_trc(ii,counter) = ind
        endif ! ind > 0
      enddo ! jj

      self%numn0_trc(ii) = counter

      counter = 0
      do jj = 1, naclsd
! ! ! ! write(*,'(9(a,i0))') __FILE__,__LINE__,' ii=',ii,' jj=',jj,' ind=',recv_buf(naclsd+3+jj,ii)
        ind = translateInd(trunc_zone, recv_buf(naclsd+3+jj,ii))
        
        if (ind > 0) then
          counter = counter + 1
          self%atom_trc(counter,ii) = ind
        endif ! ind > 0
      enddo ! jj

      self%nacls_trc(ii) = counter
      
! ! ! write(*,'(9(a,i0))') __FILE__,__LINE__,' naclsd=',naclsd
      CHECKASSERT( self%nacls_trc(ii) <= naclsd .and. self%nacls_trc(ii) > 0 )
      CHECKASSERT( self%numn0_trc(ii) <= naclsd .and. self%numn0_trc(ii) > 0 )

      self%ezoa_trc(:,ii) = recv_buf(2*naclsd+4:3*naclsd+3,ii)

      ! check if end of buffer is correct
      CHECKASSERT( recv_buf(3*naclsd+4,ii) == MAGIC )
    enddo ! ii
  endsubroutine constructIndices

endmodule ClusterInfo_mod
