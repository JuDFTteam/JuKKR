#define COMMCHECK(IERR) if((IERR) /= 0) then; write(*,*) "ERROR: Communication failure: ", __FILE__, __LINE__; STOP; endif
#define NUMBERD integer
#define NUMBERMPID MPI_INTEGER
#define NUMBERZ double complex
#define NUMBERMPIZ MPI_DOUBLE_COMPLEX
#define NUMBERC complex
#define NUMBERMPIC MPI_COMPLEX
#define NUMBERD double precision
#define NUMBERMPID MPI_DOUBLE_PRECISION
#define NUMBERI integer
#define NUMBERMPII MPI_INTEGER

!------------------------------------------------------------------------------
!> Module that implements common communication patterns.
!> Purpose: For use when collective communication is not possible
!> Author: Elias Rabel, 2012
!
!> Change only comm_patterns_T Y P E_mod.F90, then run create_comm_patterns.sh
!> This creates the files for the different datatypes needed.
!
module comm_patternsD_mod
  implicit none
  private
  public :: comm_gatherD
  public :: comm_redistributeD
  public :: comm_redistributeVD
  public :: comm_bcastD
  public :: comm_bcast2D
  public :: send_arrayD
  
  include 'mpif.h'
  
  contains

!--------------------------------------------------------------------------
!> Communicates distributed array-parts to receiver.
  subroutine comm_gatherD(my_world_rank, array, blocksize, owning_ranks, receiver)
    integer, intent(in) :: my_world_rank
    NUMBERD, intent(inout), dimension(*) :: array 
    integer, intent(in) :: blocksize
    integer, intent(in), dimension(:) :: owning_ranks
    integer, intent(in) :: receiver

    integer :: rank_index
    integer :: start
    integer :: sender

    do rank_index = 1, size(owning_ranks)
      sender = owning_ranks(rank_index)

      start = (rank_index - 1) * blocksize + 1

      call send_arrayD(my_world_rank, array(start), blocksize, sender, receiver)

    end do

  end subroutine

!--------------------------------------------------------------------------
!> Redistributes array-parts among groups of ranks.
subroutine comm_redistributeD(my_world_rank, array, blocksize, old_owners, new_owners)
  integer, intent(in) :: my_world_rank
  NUMBERD, intent(inout), dimension(*) :: array 
  integer, intent(in) :: blocksize
  integer, intent(in), dimension(:) :: old_owners
  integer, intent(in), dimension(:) :: new_owners

  integer :: rank_index
  integer :: start
  integer :: sender
  integer :: receiver

  do rank_index = 1, size(old_owners)
    sender = old_owners(rank_index)
    receiver = new_owners(rank_index)

    start = (rank_index - 1) * blocksize + 1

    call send_arrayD(my_world_rank, array(start), blocksize, sender, receiver)

  end do

end subroutine

!--------------------------------------------------------------------------
!> Redistributes array-parts of different sizes (given by array 'blocksizes'
!> among groups of ranks.
!> Note: each process needs a buffer that is large enough to hold ALL parts.
subroutine comm_redistributeVD(my_world_rank, array, blocksizes, old_owners, new_owners)
  integer, intent(in) :: my_world_rank
  NUMBERD, intent(inout), dimension(*) :: array
  integer, intent(in), dimension(:) :: blocksizes
  integer, intent(in), dimension(:) :: old_owners
  integer, intent(in), dimension(:) :: new_owners

  integer :: rank_index
  integer :: start
  integer :: sender
  integer :: receiver

  start = 1
  do rank_index = 1, size(old_owners)
    sender = old_owners(rank_index)
    receiver = new_owners(rank_index)

    call send_arrayD(my_world_rank, array(start), blocksizes(rank_index), sender, receiver)
    start = start + blocksizes(rank_index)

  end do

end subroutine

!-------------------------------------------------------------------
!> communicate 'array' starting from 'owner' to all the ranks
!> in 'ranks' in a round-robin-fashion
!> owner --> 1 --> 2 --> 3 --> 4
!> Also accepts single entry in ranks
!> It is no problem if the owner is the first entry in 'ranks'
!> Also accepts duplicate entries - but: unnecessary communication
subroutine comm_bcastD(my_world_rank, array, length, ranks, owner)
  integer, intent(in) :: my_world_rank
  NUMBERD, intent(inout), dimension(*) :: array
  integer, intent(in) :: length
  integer, intent(in), dimension(:) :: ranks
  integer, intent(in) :: owner

  integer :: rank_index
  integer :: sender, receiver

  sender = owner
  receiver = ranks(1)

  call send_arrayD(my_world_rank, array, length, sender, receiver)

  do rank_index = 1, size(ranks)-1
    sender = ranks(rank_index)
    receiver = ranks(rank_index+1)

    call send_arrayD(my_world_rank, array, length, sender, receiver)

  end do

end subroutine

!-------------------------------------------------------------------
!> communicate 'array' starting from the first entry in 'ranks'
!> to all ranks in 'ranks' in a round-robin-fashion
!> 1 --> 2 --> 3 --> 4
!> Also accepts single entry in ranks
!> Also accepts duplicate entries - but: unnecessary communication
subroutine comm_bcast2D(my_world_rank, array, length, ranks)
  integer, intent(in) :: my_world_rank
  NUMBERD, intent(inout), dimension(*) :: array
  integer, intent(in) :: length
  integer, intent(in), dimension(:) :: ranks

  integer :: rank_index
  integer :: sender, receiver

  do rank_index = 1, size(ranks)-1
    sender = ranks(rank_index)
    receiver = ranks(rank_index+1)

    call send_arrayD(my_world_rank, array, length, sender, receiver)

  end do

end subroutine

!-----------------------------------------------------------------
!> Helper routine. Sends 'length' entries of 'array' from sender to
!> receiver (ranks in MPI_COMM_WORLD).
subroutine send_arrayD(my_world_rank, array, length, sender, receiver)
  integer, intent(in) :: my_world_rank
  NUMBERD, intent(inout), dimension(*) :: array
  integer, intent(in) :: length
  integer, intent(in) :: sender
  integer, intent(in) :: receiver

  integer :: ierr
  integer :: tag
  integer, dimension(MPI_STATUS_SIZE) :: status

  tag = sender

  if (sender /= receiver) then

    if (my_world_rank == sender) then
      call MPI_Send(array, length, NUMBERMPID, receiver, &
                    tag, MPI_COMM_WORLD, ierr)

      COMMCHECK(ierr)

    end if

    if (my_world_rank == receiver) then
      call MPI_Recv(array, length, NUMBERMPID, sender, &
                    tag, MPI_COMM_WORLD, status, ierr)

      COMMCHECK(ierr)

    end if

  end if

end subroutine

end module
