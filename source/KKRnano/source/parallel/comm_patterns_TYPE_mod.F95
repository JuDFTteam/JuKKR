#define COMMCHECK(IERR) if((IERR) /= 0) then; write(*,*) "ERROR: Communication failure: ", __FILE__, __LINE__; STOP; endif

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

!------------------------------------------------------------------------------
!> Module that implements common communication patterns.
!> Purpose: For use when collective communication is not possible
!> Author: Elias Rabel, 2012
!
!> Change only comm_patterns_T Y P E_mod.F90, then run create_comm_patterns.sh
!> This creates the files for the different datatypes needed.
!
module comm_patterns_TYPE_mod
  implicit none
  private
  
!   !!! interfacing these functions does not work because of (*) interfaces
!   public :: comm_gather, comm_redistribute, comm_bcast, send_array
  
  public :: comm_gather_TYPE
  public :: comm_redistribute_TYPE
  public :: comm_redistributeV_TYPE
  public :: comm_bcast_TYPE
  public :: comm_bcast2_TYPE ! is not used in the code so far
  public :: send_array_TYPE
  
  contains

  !--------------------------------------------------------------------------
  !> Communicates distributed array-parts to receiver.
  subroutine comm_gather_TYPE(my_world_rank, array, blocksize, owning_ranks, receiver)
    integer, intent(in) :: my_world_rank
    NUMBER_TYPE, intent(inout) :: array(*)
    integer, intent(in) :: blocksize
    integer, intent(in) :: owning_ranks(:)
    integer, intent(in) :: receiver

    integer :: rank_index
    integer :: start
    integer :: sender

    do rank_index = 1, size(owning_ranks)
      sender = owning_ranks(rank_index)

      start = (rank_index - 1) * blocksize + 1

      call send_array_TYPE(my_world_rank, array(start), blocksize, sender, receiver)

    enddo ! rank_index

  endsubroutine ! comm_gather

  !--------------------------------------------------------------------------
  !> Redistributes array-parts among groups of ranks.
  subroutine comm_redistribute_TYPE(my_world_rank, array, blocksize, old_owners, new_owners)
    integer, intent(in) :: my_world_rank
    NUMBER_TYPE, intent(inout) :: array(*)
    integer, intent(in) :: blocksize
    integer, intent(in) :: old_owners(:)
    integer, intent(in) :: new_owners(:)

    integer :: rank_index
    integer :: start
    integer :: sender
    integer :: receiver

    do rank_index = 1, size(old_owners)
      sender = old_owners(rank_index)
      receiver = new_owners(rank_index)

      start = (rank_index - 1) * blocksize + 1

      call send_array_TYPE(my_world_rank, array(start), blocksize, sender, receiver)

    enddo ! rank_index

  endsubroutine ! comm_redistribute

  !--------------------------------------------------------------------------
  !> Redistributes array-parts of different sizes (given by array 'blocksizes'
  !> among groups of ranks.
  !> Note: each process needs a buffer that is large enough to hold ALL parts.
  subroutine comm_redistributeV_TYPE(my_world_rank, array, blocksizes, old_owners, new_owners)
    integer, intent(in) :: my_world_rank
    NUMBER_TYPE, intent(inout) :: array(*)
    integer, intent(in) :: blocksizes(:)
    integer, intent(in) :: old_owners(:)
    integer, intent(in) :: new_owners(:)

    integer :: rank_index
    integer :: start
    integer :: sender
    integer :: receiver

    start = 1
    do rank_index = 1, size(old_owners)
      sender = old_owners(rank_index)
      receiver = new_owners(rank_index)

      call send_array_TYPE(my_world_rank, array(start), blocksizes(rank_index), sender, receiver)
      start = start + blocksizes(rank_index)

    enddo ! rank_index

  endsubroutine ! comm_redistributeV

  !-------------------------------------------------------------------
  !> communicate 'array' starting from 'owner' to all the ranks
  !> in 'ranks' in a round-robin-fashion
  !> owner --> 1 --> 2 --> 3 --> 4
  !> Also accepts single entry in ranks
  !> It is no problem if the owner is the first entry in 'ranks'
  !> Also accepts duplicate entries - but: unnecessary communication
  subroutine comm_bcast_TYPE(my_world_rank, array, length, ranks, owner)
    integer, intent(in) :: my_world_rank
    NUMBER_TYPE, intent(inout) :: array(*)
    integer, intent(in) :: length
    integer, intent(in) :: ranks(:)
    integer, intent(in) :: owner

    integer :: rank_index
    integer :: sender, receiver

    sender = owner
    receiver = ranks(1)

    call send_array_TYPE(my_world_rank, array, length, sender, receiver)

    do rank_index = 1, size(ranks)-1
      sender = ranks(rank_index)
      receiver = ranks(rank_index+1)

      call send_array_TYPE(my_world_rank, array, length, sender, receiver)

    enddo ! rank_index

  endsubroutine ! comm_bcast

  !-------------------------------------------------------------------
  !> communicate 'array' starting from the first entry in 'ranks'
  !> to all ranks in 'ranks' in a round-robin-fashion
  !> 1 --> 2 --> 3 --> 4
  !> Also accepts single entry in ranks
  !> Also accepts duplicate entries - but: unnecessary communication
  subroutine comm_bcast2_TYPE(my_world_rank, array, length, ranks)
    integer, intent(in) :: my_world_rank
    NUMBER_TYPE, intent(inout) :: array(*)
    integer, intent(in) :: length
    integer, intent(in) :: ranks(:)

    integer :: rank_index
    integer :: sender, receiver

    do rank_index = 1, size(ranks)-1
      sender = ranks(rank_index)
      receiver = ranks(rank_index+1)

      call send_array_TYPE(my_world_rank, array, length, sender, receiver)

    enddo ! rank_index

  endsubroutine ! comm_bcast2

  !-----------------------------------------------------------------
  !> Helper routine. Sends 'length' entries of 'array' from sender to
  !> receiver (ranks in MPI_COMM_WORLD).
  subroutine send_array_TYPE(my_world_rank, array, length, sender, receiver)
    integer, intent(in) :: my_world_rank
    NUMBER_TYPE, intent(inout) :: array(*)
    integer, intent(in) :: length
    integer, intent(in) :: sender
    integer, intent(in) :: receiver

    include 'mpif.h'
    
    integer :: ierr
    integer :: tag
    integer :: status(MPI_STATUS_SIZE)

    tag = sender

    if (sender /= receiver) then

      if (my_world_rank == sender) &
        call MPI_Send(array, length, NUMBERMPI_TYPE, receiver, tag, MPI_COMM_WORLD, ierr)

      if (my_world_rank == receiver) &
        call MPI_Recv(array, length, NUMBERMPI_TYPE, sender, tag, MPI_COMM_WORLD, status, ierr)
        
      COMMCHECK(ierr) ! only MPI_Send or MPI_Recv has been called, never both
    endif ! sender /= receiver

  endsubroutine ! send_array

endmodule ! comm_patterns_mod
