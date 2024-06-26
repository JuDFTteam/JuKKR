!------------------------------------------------------------------------------
!> Module for one-sided MPI communication
!>
!> @author Elias Rabel
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

!> Note by Paul F Baumeister, 2019
!
!> The file extension .F95 indicates that sed will be applied to this source file
!> to replace a missing template feature in Fortran.
!> Do not use more than one name with _TYPE per line!

#define COMMCHECK(X) if ( (X) /= 0 ) then; write(*,*) "Comm failure", X, __LINE__; STOP; endif
#define CHECK(X) if ( .not. (X) ) then; write(*,*) "FAIL: ", __LINE__; STOP; endif

#define NUMBERZ double complex
#define NUMBERMPIZ MPI_DOUBLE_COMPLEX
#define NUMBERC complex
#define NUMBERMPIC MPI_COMPLEX
#define NUMBERD double precision
#define NUMBERMPID MPI_DOUBLE_PRECISION
#define NUMBERI integer
#define NUMBERMPII MPI_INTEGER

#ifndef NUMBER_TYPE
! defaults
#define NUMBER_TYPE integer
#define NUMBERMPI_TYPE MPI_INTEGER
#endif

module one_sided_comm_TYPE_mod
  implicit none
  private
  
  ! these statements must be each one in one line as we use the sed text replacement
  public :: copyFrom_TYPE_com
  public :: exposeBuffer_TYPE
  public :: copyChunks_TYPE
  public :: copyChunksNoSync_TYPE
  public :: fence_TYPE
  public :: hideBuffer_TYPE
  
  include 'mpif.h'

  contains

  !------------------------------------------------------------------------------
  !> Routine for exchanging data within certain atoms only.
  !>
  !> Convenience function
  !>
  !> size of local buffer  :  chunk_size*num_local_atoms
  !> size of receive buffer:  chunk_size*size(atom_ids)
  !> Uses MPI-RMA
  subroutine copyFrom_TYPE_com(receive_buf, local_buf, atom_ids, chunk_size, num_local_atoms, comm)
    use ChunkIndex_mod, only: getRankAndLocalIndex
    NUMBER_TYPE, intent(inout) :: receive_buf(*)    ! receive
    NUMBER_TYPE, intent(inout) :: local_buf(*) ! send
    integer, intent(in) :: atom_ids(:)
    integer, intent(in) :: chunk_size
    integer, intent(in) :: num_local_atoms
    integer, intent(in) :: comm

    integer(kind=4), allocatable :: chunk_inds(:,:)
    integer :: ii, ierr, nranks, atom_requested, naez_trc, win, max_local_atoms

    naez_trc = size(atom_ids) ! number of truncated atoms

    call MPI_Comm_size(comm, nranks, ierr)
    call MPI_Allreduce(num_local_atoms, max_local_atoms, 1, MPI_INTEGER, MPI_MAX, comm, ierr)

    allocate(chunk_inds(2,naez_trc))

    do ii = 1, naez_trc
      ! get 'real' atom index
      atom_requested = atom_ids(ii)

      chunk_inds(:,ii) = getRankAndLocalIndex(atom_requested, max_local_atoms*nranks, nranks)

    enddo ! ii

    call exposeBuffer_TYPE(win, local_buf, chunk_size*num_local_atoms, chunk_size, comm)
    call copyChunks_TYPE(receive_buf, win, chunk_inds, chunk_size)
    call hideBuffer_TYPE(win)

    deallocate(chunk_inds)

  endsubroutine ! copyFrom_TYPE_com



  subroutine exposeBuffer_TYPE(win, buffer, bsize, chunk_size, comm)
    integer, intent(inout) :: win 
    NUMBER_TYPE, intent(in) :: buffer(*)
    integer, intent(in) :: bsize
    integer, intent(in) :: chunk_size 
    integer, intent(in) :: comm

    integer(kind=MPI_ADDRESS_KIND) :: typesize, lowerbound
    integer :: disp_unit, ierr

    call MPI_Type_get_extent(NUMBERMPI_TYPE, lowerbound, typesize, ierr)
    COMMCHECK(ierr)

    disp_unit = typesize * chunk_size ! has to be plain integer!!

    ! Measure in units of chunks here disp_unit = CHUNKSIZE
    call MPI_Win_create(buffer, typesize*bsize, disp_unit, MPI_INFO_NULL, comm, win, ierr)
    COMMCHECK(ierr)
    
  endsubroutine ! exposeBuffer_TYPE

  !------------------------------------------------------------------------------
  !> Copy chunks of size 'chunk_size' located at 
  !> (rank/local index) locations given in 'chunk_inds' into 'dest_buffer'.
  !> On output dest_buffer contains chunks in 
  !> the order as specified in 'chunk_inds' 
  subroutine copyChunks_TYPE(dest_buffer, win, chunk_inds, chunk_size)
    NUMBER_TYPE, intent(out) :: dest_buffer(*)
    integer, intent(inout) :: win
    integer(kind=4), intent(in) :: chunk_inds(:,:) ! (2,*)
    integer, intent(in) :: chunk_size

    call fence_TYPE(win)
    call copyChunksNoSync_TYPE(dest_buffer, win, chunk_inds, chunk_size)
    ! ensure that Get has completed and dest_buffer is valid
    call fence_TYPE(win)

  endsubroutine ! copyChunks_TYPE

  !------------------------------------------------------------------------------
  !> Copy chunks of size 'chunk_size' located at
  !> (rank/local index) locations given in 'chunk_inds' into 'dest_buffer'.
  !> On output dest_buffer contains chunks in
  !> the order as specified in 'chunk_inds' - NO FENCE CALLS
  !>
  !> NOTE: MUST do fence calls before and after one or many calls to this routine.
  subroutine copyChunksNoSync_TYPE(dest_buffer, win, chunk_inds, chunk_size)
    NUMBER_TYPE, intent(out) :: dest_buffer(*)
    integer, intent(inout) :: win
    integer(kind=4), intent(in) :: chunk_inds(:,:) ! (2,*)
    integer, intent(in) :: chunk_size

    integer :: owner_rank, local_ind, ii, ierr
    integer(kind=MPI_ADDRESS_KIND) :: disp

    do ii = 1, size(chunk_inds, 2)
      owner_rank = chunk_inds(1,ii)
      local_ind  = chunk_inds(2,ii)

      disp = local_ind - 1 ! Measure in units of chunks here disp_unit = CHUNKSIZE

      call MPI_Get(dest_buffer((ii-1)*chunk_size+1), chunk_size, NUMBERMPI_TYPE, &
                        owner_rank, disp, chunk_size, NUMBERMPI_TYPE, win, ierr)

    enddo ! ii

  endsubroutine ! copyChunksNoSync_TYPE

  !------------------------------------------------------------------------------
  !> Wrapper for fence call.
  subroutine fence_TYPE(win)
    integer, intent(inout) :: win

    integer :: ierr
    call MPI_Win_fence(0, win, ierr)
    COMMCHECK(ierr)
    
  endsubroutine ! fence_TYPE

  !------------------------------------------------------------------------------
  !> Hide buffer after completing one-sided communication.
  subroutine hideBuffer_TYPE(win)
    integer, intent(inout) :: win

    integer :: ierr
    call MPI_Win_free(win, ierr)
    COMMCHECK(ierr)

  endsubroutine ! hideBuffer_TYPE

endmodule ! one_sided_comm_TYPE_mod

#ifdef TEST_ONE_SIDED_COMM__TYPE__
! a test program - not compiled due to conditional compilation
program test_one_sided_comm_TYPE
  use ChunkIndex_mod, only: getRankAndLocalIndex
  use one_sided_comm_TYPE_mod
  implicit none

  include 'mpif.h'

  integer :: myrank, num_ranks, nchunks_total, partner_rank
  integer, parameter :: CHUNKSIZE=6, CHUNKSPERPROC=3, NUMREQUESTED=7

  NUMBER_TYPE :: buffer(CHUNKSIZE,CHUNKSPERPROC)
  NUMBER_TYPE :: dest_buffer(CHUNKSIZE, NUMREQUESTED)
  integer :: chunks_req(NUMREQUESTED)
  integer(kind=4) :: chunk_inds(2,NUMREQUESTED)
  integer :: win, ierr, ii, local_ind

  call MPI_Init(ierr)
  COMMCHECK(ierr)

  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
  COMMCHECK(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, num_ranks, ierr)
  COMMCHECK(ierr)

  !buffer = dcmplx( dble(myrank+1), dble(myrank+1) )
  do ii = 1, CHUNKSPERPROC
    buffer(:,ii) = (myrank * CHUNKSPERPROC + ii - 1) ! encode rank and local index in 1 number
  enddo ! ii
  !write(*,*) buffer

  call MPI_Comm_size(MPI_COMM_WORLD, num_ranks, ierr)
  COMMCHECK(ierr)

  nchunks_total = CHUNKSPERPROC * num_ranks ! every rank owns CHUNKSPERPROC chunks

  do ii = 1, NUMREQUESTED
    chunks_req(ii) = mod((myrank + 1) * CHUNKSPERPROC + ii - 1, CHUNKSPERPROC*num_ranks) + 1
    chunk_inds(:,ii) = getRankAndLocalIndex(chunks_req(ii), nchunks_total, num_ranks)
  enddo ! ii

  call exposeBuffer_TYPE(win, buffer, size(buffer), CHUNKSIZE, MPI_COMM_WORLD)
  call copyChunks_TYPE(dest_buffer, win, chunk_inds, CHUNKSIZE)
  call hideBuffer_TYPE(win)
 
  !write(*,*) "Rank ", myrank, " wants ", chunks_req 
  !write(*,*) "Rank ", myrank, dest_buffer
  
  ! check if right buffer was received
  do ii = 1, NUMREQUESTED
    partner_rank = getOwner(chunks_req(ii), nchunks_total, num_ranks)
    local_ind = getLocalInd(chunks_req(ii), nchunks_total, num_ranks)
    CHECK( sum( abs( dest_buffer(:, ii) - (partner_rank * CHUNKSPERPROC + local_ind - 1) ) ) < 1e-10 )
  enddo ! ii  

  write(*,*) "Rank ", myrank, " has finished." ! correct if EVERY rank prints this message

  call MPI_Finalize(ierr)
  COMMCHECK(ierr)
  
  !write(*,*) "OK: TEST successful"
endprogram ! test
#endif
