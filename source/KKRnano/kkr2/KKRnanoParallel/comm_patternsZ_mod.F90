#define COMMCHECK(IERR) if((IERR) /= 0) then; write(*,*) "ERROR: Communication failure: ", __FILE__, __LINE__; STOP; endif

#define NUMBERZ double complex
#define NUMBERMPIZ MPI_DOUBLE_COMPLEX
#define NUMBERC complex
#define NUMBERMPIC MPI_COMPLEX
#define NUMBERD double precision
#define NUMBERMPID MPI_DOUBLE_PRECISION
#define NUMBERI integer
#define NUMBERMPII MPI_INTEGER

#ifndef NUMBERZ
! defaults
#define NUMBERZ integer
#define NUMBERMPIZ MPI_INTEGER
#endif

!------------------------------------------------------------------------------
!> Module that implements common communication patterns.
!> Purpose: For use when collective communication is not possible
!> Author: Elias Rabel, 2012
!
!> Change only comm_patterns_T Y P E_mod.F90, then run create_comm_patterns.sh
!> This creates the files for the different datatypes needed.
!
module comm_patternsZ_mod
  implicit none
  private
  public :: comm_gather, comm_redistribute, comm_bcast, send_array
  
  interface comm_gather
    module procedure comm_gatherZ
  endinterface

  interface comm_redistribute
    module procedure comm_redistributeZ, &
                     comm_redistributeVZ
  endinterface
  
  interface comm_bcast
    module procedure comm_bcastZ, &
                     comm_bcast2Z
  endinterface

  interface send_array
    module procedure send_arrayZ
  endinterface
  
  ! deprecated public statements (to be private in the future)
  public :: comm_gatherZ
  public :: comm_redistributeZ
  public :: comm_redistributeVZ
  public :: comm_bcastZ
  public :: comm_bcast2Z
  public :: send_arrayZ
  
  include 'mpif.h'
  
  contains

!--------------------------------------------------------------------------
!> Communicates distributed array-parts to receiver.
  subroutine comm_gatherZ(my_world_rank, array, blocksize, owning_ranks, receiver)
    integer, intent(in) :: my_world_rank
    NUMBERZ, intent(inout) :: array(*)
    integer, intent(in) :: blocksize
    integer, intent(in) :: owning_ranks(:)
    integer, intent(in) :: receiver

    integer :: rank_index
    integer :: start
    integer :: sender

    do rank_index = 1, size(owning_ranks)
      sender = owning_ranks(rank_index)

      start = (rank_index - 1) * blocksize + 1

      call send_arrayZ(my_world_rank, array(start), blocksize, sender, receiver)

    enddo ! rank_index

  endsubroutine comm_gatherZ

!--------------------------------------------------------------------------
!> Redistributes array-parts among groups of ranks.
subroutine comm_redistributeZ(my_world_rank, array, blocksize, old_owners, new_owners)
  integer, intent(in) :: my_world_rank
  NUMBERZ, intent(inout) :: array(*)
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

    call send_arrayZ(my_world_rank, array(start), blocksize, sender, receiver)

  enddo ! rank_index

endsubroutine comm_redistributeZ

!--------------------------------------------------------------------------
!> Redistributes array-parts of different sizes (given by array 'blocksizes'
!> among groups of ranks.
!> Note: each process needs a buffer that is large enough to hold ALL parts.
subroutine comm_redistributeVZ(my_world_rank, array, blocksizes, old_owners, new_owners)
  integer, intent(in) :: my_world_rank
  NUMBERZ, intent(inout) :: array(*)
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

    call send_arrayZ(my_world_rank, array(start), blocksizes(rank_index), sender, receiver)
    start = start + blocksizes(rank_index)

  enddo ! rank_index

endsubroutine comm_redistributeVZ

!-------------------------------------------------------------------
!> communicate 'array' starting from 'owner' to all the ranks
!> in 'ranks' in a round-robin-fashion
!> owner --> 1 --> 2 --> 3 --> 4
!> Also accepts single entry in ranks
!> It is no problem if the owner is the first entry in 'ranks'
!> Also accepts duplicate entries - but: unnecessary communication
subroutine comm_bcastZ(my_world_rank, array, length, ranks, owner)
  integer, intent(in) :: my_world_rank
  NUMBERZ, intent(inout) :: array(*)
  integer, intent(in) :: length
  integer, intent(in) :: ranks(:)
  integer, intent(in) :: owner

  integer :: rank_index
  integer :: sender, receiver

  sender = owner
  receiver = ranks(1)

  call send_arrayZ(my_world_rank, array, length, sender, receiver)

  do rank_index = 1, size(ranks)-1
    sender = ranks(rank_index)
    receiver = ranks(rank_index+1)

    call send_arrayZ(my_world_rank, array, length, sender, receiver)

  enddo ! rank_index

endsubroutine comm_bcastZ

!-------------------------------------------------------------------
!> communicate 'array' starting from the first entry in 'ranks'
!> to all ranks in 'ranks' in a round-robin-fashion
!> 1 --> 2 --> 3 --> 4
!> Also accepts single entry in ranks
!> Also accepts duplicate entries - but: unnecessary communication
subroutine comm_bcast2Z(my_world_rank, array, length, ranks)
  integer, intent(in) :: my_world_rank
  NUMBERZ, intent(inout) :: array(*)
  integer, intent(in) :: length
  integer, intent(in) :: ranks(:)

  integer :: rank_index
  integer :: sender, receiver

  do rank_index = 1, size(ranks)-1
    sender = ranks(rank_index)
    receiver = ranks(rank_index+1)

    call send_arrayZ(my_world_rank, array, length, sender, receiver)

  enddo ! rank_index

endsubroutine comm_bcast2Z

!-----------------------------------------------------------------
!> Helper routine. Sends 'length' entries of 'array' from sender to
!> receiver (ranks in MPI_COMM_WORLD).
subroutine send_arrayZ(my_world_rank, array, length, sender, receiver)
  integer, intent(in) :: my_world_rank
  NUMBERZ, intent(inout) :: array(*)
  integer, intent(in) :: length
  integer, intent(in) :: sender
  integer, intent(in) :: receiver

  integer :: ierr
  integer :: tag
  integer :: status(MPI_STATUS_SIZE)

  tag = sender

  if (sender /= receiver) then

    if (my_world_rank == sender) then
      call MPI_Send(array, length, NUMBERMPIZ, receiver, tag, MPI_COMM_WORLD, ierr)
      COMMCHECK(ierr)
    endif

    if (my_world_rank == receiver) then
      call MPI_Recv(array, length, NUMBERMPIZ, sender, tag, MPI_COMM_WORLD, status, ierr)
      COMMCHECK(ierr)
    endif

  endif

endsubroutine send_arrayZ

endmodule comm_patternsZ_mod
