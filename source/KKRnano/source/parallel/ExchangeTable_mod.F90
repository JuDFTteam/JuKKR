!------------------------------------------------------------------------------
!> Module for two-sided MPI communication
!>
!> @author Paul Baumeister
!>
!> Only for equally sized chunks of data which have basic numeric datatypes.
!>
!> Data layout:
!>
!> \verbatim
!>
!> A distributed array is created (very similar to Co-array Fortran),
!> which consists of equally sized chunks.
!>
!> Each rank then copies a slice of this array into its local memory.
!>
!> Note: this corresponds to the "partitioned global address space" (PGAS)
!>       programming model - however in our case only READ access is allowed.
!>
!> Each chunk has an "owner"-index (MPI-rank starting at 0)
!> and a local index (starting at 1 !)
!>
!> chunk index = (rank, local index)
!>
!>  ___________|_____ _____|
!> |     |     |     |     |
!> |(0,1)|(0,2)|(1,1)|(1,2)| usw...
!> |_____|_____|_____|_____|
!>    rank 0   |  rank 1   |
!>
!> |atm 1|atm 2|atm 3|atm 4| ...
!>
!> To access the chunks by a continuous index (called 'atom index')
!> one can use the routines getOwner and getLocalInd
!> to convert to a chunk index
!>
!> \endverbatim

module ExchangeTable_mod
! use ChunkIndex_mod, only: getRankAndLocalIndex

#define assert(condition) if(.not. (condition)) stop __LINE__

  implicit none
  private
  
  public :: ExchangeTable, create, destroy
  
  type ExchangeTable
    integer              :: comm !> MPI communicator handle
    integer              :: nPairs
    integer, allocatable :: pair2rank(:) ! dim(nrows)
    integer              :: recv_n
    integer, allocatable :: recv_start(:) ! dim(nrows + 1)
    integer, allocatable :: recv_index(:)  ! dim(recv_n)
    integer              :: send_n
    integer, allocatable :: send_start(:) ! dim(nrows + 1)
    integer, allocatable :: send_index(:)  ! dim(nsend_n)
    ! ToDo:
    ! this could be implemented as we can save some communication volume by 
    ! transferring nacls(a) instead of naclsd == maxval(nacls) elements
!   integer, allocatable :: recv_size(:) ! dim(recv_n) 
!   integer, allocatable :: send_size(:) ! dim(nsend_n)
  endtype

  interface create
    module procedure create_ExchangeTable
  endinterface
  
  interface destroy
    module procedure destroy_ExchangeTable
  endinterface

  contains

  subroutine create_ExchangeTable(self, global_atom_id, comm, max_local_atoms)
    use ChunkIndex_mod, only: getRankAndLocalIndex
    type(ExchangeTable), intent(out) :: self
    integer, intent(in) :: global_atom_id(:) !> dim(num_trunc_atoms) mapping trunc index -> global atom index
    integer, intent(in) :: comm !> MPI communicator handle
    integer, intent(in) :: max_local_atoms !> needed only for the mapping from global atom IDs to (rank, local index)

    ! locals
    integer :: iout, iinp, rank, tag, myrank, nranks, ierr, ist
    integer :: num_trunc_atoms ! number of atoms inside the truncation zone
    integer(kind=4) :: gid ! global atom id
    integer, parameter :: TAGMOD = 2**15
    integer :: ipair, inz, nPairs
    
    integer, allocatable :: nrecv(:), nsend(:), temp_recv_index(:), reqs(:,:), stats(:,:,:)
    integer(kind=4), allocatable :: chunk_inds(:,:) ! dim(2,num_trunc_atoms)
    integer(kind=2), allocatable :: nreq(:), rank2pair(:) ! assume we will not communicate with more than 32000 MPI processes
    
    include 'mpif.h' ! only: MPI_STATUS_SIZE, MPI_INTEGER, MPI_REQUEST_NULL
    
    num_trunc_atoms = size(global_atom_id)

    call MPI_Comm_size(comm, nranks, ierr)
    call MPI_Comm_rank(comm, myrank, ierr)

    ! ===============================================================================
    ! part 1: exchange the Ginp arrays between the MPI processes 
    ! ===============================================================================
    
    allocate(nreq(0:nranks-1), chunk_inds(2,num_trunc_atoms))!, stat=ist) ! ToDo: catch status
    
    ! loop over the sites sending the information
    nreq(:) = 0
    do iout = 1, num_trunc_atoms
      gid = global_atom_id(iout) ! get the global atom id

!     rank = (gid - 1)/max_local_atoms ! block distribution of atoms to ranks
      chunk_inds(:,iout) = getRankAndLocalIndex(gid, max_local_atoms*nranks, nranks)
      rank = chunk_inds(1,iout) ! remote rank
!     iinp = chunk_inds(2,iout) ! local index on remote rank, index not used here
      
      nreq(rank) = nreq(rank) + 1 ! add one request
    enddo ! iout
    
    
    nPairs = count(nreq > 0)
    allocate(self%pair2rank(nPairs))
    allocate(rank2pair(0:nranks-1)) ; rank2pair(:) = -1
    
    assert( sum(nreq) == num_trunc_atoms )
    self%recv_n = num_trunc_atoms

    ! create prefix sum
    allocate(self%recv_start(nPairs + 1), self%recv_index(self%recv_n))
!   allocate(self%recv_size(self%recv_n)) ! ToDo
    allocate(temp_recv_index(self%recv_n))
    ipair = 0
    inz = 0
    do rank = 0, nranks - 1
      if (nreq(rank) > 0) then
        ipair = ipair + 1 ! create a new pair
        self%pair2rank(ipair) = rank
        rank2pair(rank) = ipair
        self%recv_start(ipair) = inz + 1
        inz = inz + nreq(rank)
      endif
    enddo ! rank
    assert( inz == self%recv_n )
    self%recv_start(nPairs + 1) = self%recv_n + 1

    deallocate(nreq, stat=ist) ! ignore status

    allocate(nrecv(nPairs)) ; nrecv(:) = 0
    do iout = 1, num_trunc_atoms
      rank = chunk_inds(1,iout) ! remote rank
      iinp = chunk_inds(2,iout) ! local index on remote rank
      
      ipair = rank2pair(rank)
      assert( ipair > 0 ) ! this pair must exist, otherwise we made a fatal error
      
      inz = self%recv_start(ipair) + nrecv(ipair)
      nrecv(ipair) = nrecv(ipair) + 1 ! create a new request

      temp_recv_index(inz) = iinp ! store the local owner index in remote sender rank which I want to receive
      self%recv_index(inz) = iout

!     self%recv_index(inz) = global_atom_id(iout) ! get the global atom id, needed for DEBUG only
!     self%recv_size(inz)  = nacls_trunc(iout) ! ToDo
    enddo ! atoms

    deallocate(chunk_inds, stat=ist) ! ignore status
    
    allocate(reqs(2,nPairs), stats(MPI_STATUS_SIZE,2,nPairs))!, stat=ist) ! ToDo: catch status
    allocate(nsend(nPairs)) ; nsend(:) = 0
    
    reqs(:,:) = MPI_REQUEST_NULL
    do ipair = 1, nPairs
      rank = self%pair2rank(ipair)
    
      tag = modulo(rank*myrank, TAGMOD)
      call MPI_Isend(nrecv(ipair), 1, MPI_INTEGER, rank, tag, comm, reqs(1,ipair), ierr)
      call MPI_Irecv(nsend(ipair), 1, MPI_INTEGER, rank, tag, comm, reqs(2,ipair), ierr)
      
    enddo ! ipair
    
! #define DEBUG
    
    ! we have to assume that the communication pattern is symmetric
    ! i.e. each pair of two ranks either communicates into both directions or not at all
    ! we assume a symmetric communication pattern, i.e. if 
    ! rank myrank has one or more requests towards rank remote,
    ! rank remote will also have one or more requests rank myrank. 
    ! otherwise we will hang in the following call
#ifdef DEBUG
    if (myrank == 0) write(*,*) "symmetric communication pattern assumed ..."
#endif
    call MPI_Waitall(2*nPairs, reqs, stats, ierr) ! wait until all sends and all receives have finished
    call MPI_Barrier(comm, ierr) ! wait until all sends and all receives have finished on all processes
#ifdef DEBUG
    if (myrank == 0) write(*,*) "symmetric communication pattern found --> ok"
#endif

    self%send_n = sum(nsend)
    allocate(self%send_start(nPairs + 1), self%send_index(self%send_n))
!   allocate(self%send_size(self%send_n)) ! ToDo

    ! prefix sum
    inz = 0
    do ipair = 1, nPairs
      self%send_start(ipair) = inz + 1
      inz = inz + nsend(ipair)
    enddo ! ipair
    assert( inz == self%send_n )
    self%send_start(nPairs + 1) = self%send_n + 1

    ! communicate lists of iinp-indices
    reqs(:,:) = MPI_REQUEST_NULL
    do ipair = 1, nPairs
      rank = self%pair2rank(ipair)

      tag = modulo(rank*myrank, TAGMOD)
      call MPI_Isend(temp_recv_index(self%recv_start(ipair):), nrecv(ipair), MPI_INTEGER, rank, tag, comm, reqs(1,ipair), ierr)
      call MPI_Irecv(self%send_index(self%send_start(ipair):), nsend(ipair), MPI_INTEGER, rank, tag, comm, reqs(2,ipair), ierr)

    enddo ! ipair
#ifdef DEBUG
    if (myrank == 0) write(*,*) "exchange indices  in communication pattern ..."
#endif
    call MPI_Waitall(2*nPairs, reqs, stats, ierr) ! wait until all sends and all receives have finished
    call MPI_Barrier(comm, ierr) ! wait until all sends and all receives have finished on all processes
#ifdef DEBUG
    if (myrank == 0) write(*,*) "indices exchanged in communication pattern --> ok"
#endif

    self%npairs = nPairs ! store and mark the object as usable
    self%comm = comm ! store a copy of the MPI communicator handle

!   self%send_size(inz) = nacls_local(self%send_index(inz)) ! ToDo
    
    deallocate(nrecv, nsend, rank2pair, temp_recv_index, reqs, stats, stat=ist) ! ignore status
  endsubroutine ! create_ExchangeTable


  elemental subroutine destroy_ExchangeTable(self)
    type(ExchangeTable), intent(inout) :: self
    integer :: ist
    self%comm = 0
    self%nPairs = 0
    self%recv_n = 0
    self%send_n = 0
    deallocate(self%pair2rank, self%recv_start, self%recv_index, self%send_start, self%send_index, stat=ist)
!   deallocate(self%send_size, self%recv_size, stat=ist) ! ignore status
  endsubroutine ! destroy

endmodule ! ExchangeTable_mod
