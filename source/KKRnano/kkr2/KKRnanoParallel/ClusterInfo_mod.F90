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
    integer :: naclsd !< maximal number of cluster atoms,              maybe padded for memory alignment
    integer :: numn0d !< maximal number of inequivalent cluster atoms, maybe padded for memory alignment
    integer :: naez_trc !< number of atoms in the truncation zone
    integer, allocatable :: nacls(:) !> dim(naez_trc) number of target atoms in local interaction cluster
    integer, allocatable :: numn0(:) !> dim(naez_trc) number of inequivalent target atoms
    !> some quantities can be stored in 16bit since it encodes only cluster indices
    integer(kind=2), allocatable :: indn0(:,:) !> dim(numn0d,naez_trc) ! ordered set of indices of inequivalent atoms in the cluster
    integer(kind=2), allocatable ::  atom(:,:) !> dim(naclsd,naez_trc) ! indices of target atoms (periodic images included)
    integer(kind=2), allocatable ::  ezoa(:,:) !> dim(naclsd,naez_trc) ! index of periodic image into array lattice_vectors%rr
  endtype

  interface create
    module procedure createClusterInfo
  endinterface
  
  interface destroy
    module procedure destroyClusterInfo
  endinterface
  
  integer, parameter, private :: MAGIC_NUMBER = 385306
  logical, parameter, private :: dbg = .false.
#define cDBG if(dbg) 

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

    integer :: isa, nn0, ila, ineqv, ist
    integer :: OFFSET_INDN, OFFSET_ATOM, OFFSET_EZOA, POS_MAGIC
    integer, parameter :: POS_INDEX=1, POS_NACLS=2, POS_NUMN0=3
    integer :: nacls, numn0, naclsd, numn0d, num_local_atoms, blocksize
    integer :: memory_stat ! needed in allocatecheck
    integer :: ierr ! MPI error status
    integer(kind=2) :: ind ! local atom index
    integer(kind=4) :: atom_id, indn0 ! global atom ids
    integer(kind=4), allocatable :: send_buf(:,:), recv_buf(:,:) ! global atom ids, require more than 16bit
    integer(kind=4), allocatable :: global_target_atom_id(:) !!! DEBUG
    integer, parameter :: N_ALIGN = 1 ! 1:no alignment, 4: 8 Byte alignment, 32: 64 Byte alignment 

    num_local_atoms = size(ref_clusters)

    nacls = maxval(ref_clusters(:)%nacls) ! find maximum for locally stored clusters
    numn0 = maxval(ref_clusters(:)%numn0) !
    CHECKASSERT( numn0 <= nacls )

    ! determine global maximal number of cluster atoms
    call MPI_Allreduce(nacls, naclsd, 1, MPI_INTEGER, MPI_MAX, communicator, ierr)
    call MPI_Allreduce(numn0, numn0d, 1, MPI_INTEGER, MPI_MAX, communicator, ierr)
    CHECKASSERT( numn0d <= naclsd )

    self%naez_trc = trunc_zone%naez_trc

    ! start packing a send buffer
    OFFSET_ATOM = POS_NUMN0 ! the first three numbers are [source_atom_index, nacls, numn0]
    OFFSET_EZOA = OFFSET_ATOM + naclsd
    OFFSET_INDN = OFFSET_EZOA + naclsd
    POS_MAGIC   = OFFSET_INDN + numn0d + 1
    blocksize   = POS_MAGIC ! this many numbers are communicated per site

    ALLOCATECHECK(send_buf(blocksize,num_local_atoms))
    do ila = 1, num_local_atoms
      send_buf(:,ila) = -1 ! init
      nacls = ref_clusters(ila)%nacls
      numn0 = ref_clusters(ila)%numn0
      CHECKASSERT( nacls <= naclsd )
      CHECKASSERT( numn0 <= nacls )
      send_buf(POS_INDEX,ila) = ref_clusters(ila)%source_atom_index ! global atom index
      send_buf(POS_NACLS,ila) = nacls
      send_buf(POS_NUMN0,ila) = numn0
      send_buf(OFFSET_ATOM + 1:nacls + OFFSET_ATOM,ila) = ref_clusters(ila)%atom(:)
      send_buf(OFFSET_EZOA + 1:nacls + OFFSET_EZOA,ila) = ref_clusters(ila)%ezoa(:)
      send_buf(OFFSET_INDN + 1:numn0 + OFFSET_INDN,ila) = ref_clusters(ila)%indn0(:) ! indn0 has dim(numn0) now, however, numn0 <= nacls holds
      send_buf(POS_MAGIC,ila) = MAGIC_NUMBER
    enddo ! ila
    ! send buffer is packed

    ALLOCATECHECK(recv_buf(blocksize,self%naez_trc))

    ! communication
    call copyFromI_com(recv_buf, send_buf, trunc_zone%global_atom_id, blocksize, num_local_atoms, communicator)

    DEALLOCATECHECK(send_buf)
cDBG  ALLOCATECHECK(global_target_atom_id(naclsd)) !!! DEBUG
    ALLOCATECHECK(self%nacls(self%naez_trc))
    ALLOCATECHECK(self%numn0(self%naez_trc))
    self%nacls = 0
    self%numn0 = 0

    CHECKASSERT( N_ALIGN >= 1 )
    self%naclsd = ((naclsd - 1)/N_ALIGN + 1)*N_ALIGN ! alignment
    self%numn0d = ((numn0d - 1)/N_ALIGN + 1)*N_ALIGN ! alignment

    ALLOCATECHECK(self%indn0(self%numn0d,self%naez_trc))
    ALLOCATECHECK( self%atom(self%naclsd,self%naez_trc))
    ALLOCATECHECK( self%ezoa(self%naclsd,self%naez_trc))
    self%indn0 = -1
    self%atom = -1
    self%ezoa = -1

    do isa = 1, self%naez_trc
      atom_id = recv_buf(POS_INDEX,isa)
      CHECKASSERT( atom_id == trunc_zone%global_atom_id(isa) ) ! check if send_buf from right atom was received

      numn0 = recv_buf(POS_NUMN0,isa)
      ! the indices in indn0 have to be transformed to truncation zone indices
      nn0 = 0 ! create a new version of numn0
      do ineqv = 1, numn0 ! loop over all inequivalent atoms in the cluster

        indn0 = recv_buf(OFFSET_INDN + ineqv,isa) ! indn0 received
cDBG    write(*,'(9(a,i0))') __FILE__,__LINE__,' isa=',isa,' ineqv=',ineqv,' indn0=',indn0
        ind = trunc_zone%local_atom_idx(indn0) ! translate into a local index of the truncation zone

        if (ind > 0) then ! ind == -1 means that this atom is outside of truncation zone
          nn0 = nn0 + 1
          self%indn0(nn0,isa) = ind ! indn0 translated into local indices of the truncation zone
          ! if trunc_zone%global_atom_id(:) is in ascending order (strictly monotonous), also indn0 will inherit that property
cDBG      global_target_atom_id(nn0) = indn0 !!! DEBUG
        endif ! ind > 0
        
      enddo ! ineqv

      CHECKASSERT( 0 < nn0 .and. nn0 <= numn0 )
      self%numn0(isa) = nn0
cDBG  write(*,'(4(a,i0),999(" ",i0))') __FILE__,__LINE__,' atom_id=',atom_id,' numn0=',nn0,' indn0= ',global_target_atom_id(1:nn0) !!! DEBUG

      nacls = recv_buf(POS_NACLS,isa) ! nacls
      self%nacls(isa) = nacls

      self%atom(1:nacls,isa) = trunc_zone%local_atom_idx(recv_buf(OFFSET_ATOM + 1:nacls + OFFSET_ATOM,isa)) ! translate into truncation zone indices here
cDBG  global_target_atom_id(1:nacls) = recv_buf(OFFSET_ATOM + 1:nacls + OFFSET_ATOM,isa) !!! DEBUG
cDBG  write(*,'(999(a,i0))') __FILE__,__LINE__,' atom_id=',atom_id,' nacls=',nacls,' atom@ezoa= ',(global_target_atom_id(ist),'@',self%ezoa(ist,isa),' ', ist=1,nacls) !!! DEBUG

      ! the index of the periodic image does not need translation
      self%ezoa(1:nacls,isa) = recv_buf(OFFSET_EZOA + 1:nacls + OFFSET_EZOA,isa) ! convert to integer(kind=2)
cDBG  write(*,'(999(a,i0))') __FILE__,__LINE__,' atom_id=',atom_id,' old nacls=',nacls,' atom@ezoa= ',(global_target_atom_id(ist),'@',recv_buf(OFFSET_EZOA+ist,isa),' ', ist=1,nacls) !!! DEBUG

      CHECKASSERT( recv_buf(POS_MAGIC,isa) == MAGIC_NUMBER ) ! check if the end of recv_buf seems correct
    enddo ! isa

cDBG  deallocate(global_target_atom_id, stat=ist) ! ignore status !!! DEBUG
    DEALLOCATECHECK(recv_buf)
  endsubroutine ! create
  
  elemental subroutine destroyClusterInfo(self)
    type(ClusterInfo), intent(inout) :: self

    integer :: ist
    deallocate(self%nacls, self%numn0, self%indn0, self%atom, self%ezoa, stat=ist) ! ignore status
!   deallocate(self%jacls, stat=ist) ! ignore status
    self%naclsd = 0
    self%naez_trc = 0
  endsubroutine ! destroy

endmodule ! ClusterInfo_mod
