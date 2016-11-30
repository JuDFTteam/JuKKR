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

#define NUMBER_TYPE double complex
#define NUMBERMPI_TYPE MPI_DOUBLE_COMPLEX

module two_sided_comm_TYPE_mod
! use ChunkIndex_mod, only: getRankAndLocalIndex

#define assert(condition) if(.not. (condition)) stop __LINE__

  implicit none
  private
  
  public :: DataExchangeTable, create, destroy
  public :: reference_sys_com
  
  type DataExchangeTable
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
    module procedure create_DataExchangeTable
  endinterface
  
  interface destroy
    module procedure destroy_DataExchangeTable
  endinterface

  contains

  subroutine create_DataExchangeTable(self, max_local_atoms, natoms, global_atom_id, comm)
    use ChunkIndex_mod, only: getRankAndLocalIndex
    type(DataExchangeTable), intent(out) :: self
    integer, intent(in) :: max_local_atoms
    integer, intent(in) :: natoms ! number of atoms inside the truncation zone
    integer, intent(in) :: global_atom_id(:) !> mapping trunc index -> global atom index
    integer, intent(in) :: comm

    ! locals
    integer(kind=4), allocatable :: chunk_inds(:,:) ! dim(2,natoms) 
    integer :: iout, iinp, rank, tag, myrank, nranks, ierr, ist
    integer(kind=4) :: gid ! global atom id
    integer, parameter :: TAGMOD = 2**15
    integer, allocatable :: nrecv(:), nsend(:), temp_recv_index(:)
    integer :: ipair, inz, nPairs
    integer(kind=2), allocatable :: nreq(:), rank2pair(:) ! assume we will not communicate with more than 32000 MPI processes
    integer, allocatable :: reqs(:,:), stats(:,:,:)
    
    include 'mpif.h' ! only: MPI_STATUS_SIZE, MPI_INTEGER, MPI_REQUEST_NULL
    
    assert( natoms == size(global_atom_id) )

    call MPI_Comm_size(comm, nranks, ierr)
    call MPI_Comm_rank(comm, myrank, ierr)

    ! ===============================================================================
    ! part 1: exchange the Ginp arrays between the MPI processes 
    ! ===============================================================================
    
    allocate(nreq(0:nranks-1), chunk_inds(2,natoms))!, stat=ist) ! ToDo: catch status
    
    ! loop over the sites sending the information
    nreq(:) = 0
    do iout = 1, natoms
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
    
    assert( sum(nreq) == natoms )
    self%recv_n = natoms

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
    do iout = 1, natoms
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
    ! we have to assume that the communication pattern is symmetric
    ! i.e. each pair of two ranks either communicates into both directions or not at all
    ! we assume a symmetric communication pattern, i.e. if 
    ! rank myrank has one or more requests towards rank remote,
    ! rank remote will also have one or more requests rank myrank. 
    ! otherwise we will hang in the following call
    write(*,*) "symmetric communication pattern assumed ..."
    call MPI_Waitall(2*nPairs, reqs, stats, ierr) ! wait until all sends and all receives have finished
    write(*,*) "symmetric communication pattern found --> ok"
    
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
    call MPI_Waitall(2*nPairs, reqs, stats, ierr) ! wait until all sends and all receives have finished
    
    self%npairs = nPairs ! store and mark the object as usable
    self%comm = comm ! store a copy of the MPI communicator handle

!   self%send_size(inz) = nacls_local(self%send_index(inz)) ! ToDo
    
    deallocate(nrecv, nsend, rank2pair, temp_recv_index, stat=ist) ! ignore status
  endsubroutine ! create_DataExchangeTable


  elemental subroutine destroy_DataExchangeTable(self)
    type(DataExchangeTable), intent(inout) :: self
    integer :: ist
    self%comm = 0
    self%nPairs = 0
    self%recv_n = 0
    self%send_n = 0
    deallocate(self%pair2rank, self%recv_start, self%recv_index, self%send_start, self%send_index, stat=ist)
!   deallocate(self%send_size, self%recv_size, stat=ist) ! ignore status
  endsubroutine ! destroy

  
  subroutine reference_sys_com(self, Gout, Ginp)
    type(DataExchangeTable), intent(in) :: self
    NUMBER_TYPE, intent(out) :: Gout(:,:,:,:,:) !> dim(lmmaxd,lmmaxd,0:LLy,naclsd,num_trunc_atoms)
    NUMBER_TYPE, intent(in)  :: Ginp(:,:,:,:,:) !< dim(lmmaxd,lmmaxd,0:LLy,naclsd,num_local_atoms)
    
    ! In order to implement point-to-point communication also with max_local_atoms > 1,
    ! we need to set up a sparse matrix that describes 
    !   which references system (col) we need to receive from which rank (pair)
    ! then, we need to communicate these matrices
    !   and create the tables about which local reference system we need to send to which rank (pair)

    integer, allocatable :: rreq(:), sreq(:), rstats(:,:), sstats(:,:)
    integer :: ncount, ipair, rank, tag, ierr, ist, inz, iinp, iout
    include 'mpif.h' ! only: MPI_STATUS_SIZE, MPI_INTEGER, MPI_REQUEST_NULL

    assert( self%comm /= 0 )
    
    ncount = size(Ginp(:,:,:,:,1))
    do ist = 1, 4
      assert( size(Ginp, ist) == size(Gout, ist) )
    enddo ! ist

    allocate(sreq(self%send_n), sstats(MPI_STATUS_SIZE,self%send_n), &
             rreq(self%recv_n), rstats(MPI_STATUS_SIZE,self%recv_n), stat=ist) ! ToDo: catch status
    sreq(:) = MPI_REQUEST_NULL
    rreq(:) = MPI_REQUEST_NULL

    ! now communication works like this
    do ipair = 1, self%npairs
      rank = self%pair2rank(ipair)
      
      do inz = self%send_start(ipair), self%send_start(ipair + 1) - 1
        tag = inz - self%send_start(ipair)
        iinp = self%send_index(inz)
        ! ncount = size(Ginp, 1)*size(Ginp, 2)*size(Ginp, 3)*self%send_size(inz) ! ToDo
        call MPI_Isend(Ginp(:,:,:,:,iinp), ncount, NUMBERMPI_TYPE, rank, tag, self%comm, sreq(inz), ierr)
      enddo ! inz

      do inz = self%recv_start(ipair), self%recv_start(ipair + 1) - 1
        tag = inz - self%recv_start(ipair)
        iout = self%recv_index(inz)
        ! ncount = size(Gout, 1)*size(Gout, 2)*size(Gout, 3)*self%recv_size(inz) ! ToDo
        call MPI_Irecv(Gout(:,:,:,:,iout), ncount, NUMBERMPI_TYPE, rank, tag, self%comm, rreq(inz), ierr)
      enddo ! inz

    enddo ! ipair
    call MPI_Waitall(self%send_n, sreq, sstats, ierr) ! wait until all sends have finished
    call MPI_Waitall(self%recv_n, rreq, rstats, ierr) ! wait until all receives have finished

    deallocate(rreq, rstats, sreq, sstats, stat=ist) ! ignore status
  endsubroutine ! reference_sys_com

endmodule ! two_sided_comm_TYPE_mod

#ifdef TEST_TWO_SIDED_COMM__TYPE__
! a test program - not compiled due to conditional compilation
program test
  use two_sided_comm_TYPE_mod
  implicit none

  include 'mpif.h'

  integer :: myrank
  integer :: num_ranks
  integer :: nchunks_total
  integer :: partner_rank

  integer, parameter :: CHUNKSIZE=6, CHUNKSPERPROC=3, NUMREQUESTED=7

  NUMBER_TYPE :: buffer(CHUNKSIZE,CHUNKSPERPROC)
  NUMBER_TYPE :: dest_buffer(CHUNKSIZE, NUMREQUESTED)
  integer :: chunks_req(NUMREQUESTED)
  integer(kind=4) :: chunk_inds(2,NUMREQUESTED)
  integer :: win, ierr, ii, local_ind
  
#define COMMCHECK(X) if ( (X) /= 0 ) then; write(*,*) "Comm failure", X, __LINE__; STOP; endif
  
  call MPI_Init(ierr)
  COMMCHECK(ierr)

  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
  COMMCHECK(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, num_ranks, ierr)
  COMMCHECK(ierr)

  write(*,*) "Rank ", myrank, " has finished." ! correct if EVERY rank prints this message

  call MPI_Finalize(ierr)
  COMMCHECK(ierr)
  
  !write(*,*) "OK: TEST successful"
endprogram
#endif
