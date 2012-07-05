!> Author: Elias Rabel, 2012

module KKRnanoParallel_mod

  !> Opaque object. Use getter routines to access member variables
  type KKRnanoParallel
    PRIVATE
    
    !> MPI communicator handle
    !> This communicator contains all atom process belonging
    !> to the same (Spin,Energy) - group
    integer :: my_SE_communicator
    
    integer :: my_active_communicator
    integer :: active
    
    integer :: my_world_rank
    integer :: my_atom_rank  ! Note ranks start from 0 !
          
    integer :: my_spin_id    ! Note ids start from 1 !
    integer :: my_energy_id
    integer :: my_atom_id
    integer :: my_SE_id

    integer :: num_comms_
    integer :: procs_per_comm_

    integer :: num_atom_ranks_
    integer :: num_spin_ranks_
    integer :: num_energy_ranks_
  end type

  contains

  !---------------------------------------------------------------------------------------
  subroutine createKKRnanoParallel(my_mpi, num_atom_ranks, num_spin_ranks, num_energy_ranks)
    implicit none
    include 'mpif.h'

    type (KKRnanoParallel), intent(inout) :: my_mpi
    integer, intent(in) :: num_atom_ranks
    integer, intent(in) :: num_spin_ranks
    integer, intent(in) :: num_energy_ranks

    integer :: num_comms
    integer :: procs_per_comm
    integer :: ierr
    integer :: color
    integer :: key
    integer :: num_ranks    
 
    num_comms = num_spin_ranks * num_energy_ranks
    procs_per_comm = num_atom_ranks

    my_mpi%num_comms_ = num_comms
    my_mpi%procs_per_comm_ = procs_per_comm
    my_mpi%num_atom_ranks_ = num_atom_ranks
    my_mpi%num_energy_ranks_ = num_energy_ranks
    my_mpi%num_spin_ranks_ = num_spin_ranks

    call MPI_Init(ierr)

    call MPI_Comm_size(MPI_COMM_WORLD, num_ranks, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, my_mpi%my_world_rank, ierr)

    ! Check if there are enough ranks.
    if (num_ranks < num_comms*procs_per_comm) then
      if (my_mpi%my_world_rank == 0) then
        write(*,*) "Not enough MPI ranks."
      end if
      call MPI_Finalize(ierr)
      stop
    end if
    
    color = my_mpi%my_world_rank / procs_per_comm
    key = my_mpi%my_world_rank - (color * procs_per_comm)
    
    !---------- Define ids ---------------
    ! Note: ids start with 1
    ! TODO
    my_mpi%my_atom_id = key + 1
    my_mpi%my_SE_id = color + 1
    my_mpi%my_spin_id = (color / num_energy_ranks) + 1
    my_mpi%my_energy_id = mod(color, num_energy_ranks) + 1
    
    if (color >= num_comms) then
      color = MPI_UNDEFINED ! inactive
    end if

    call MPI_Comm_split(MPI_COMM_WORLD, color, key, my_mpi%my_SE_communicator, ierr)
    
    if (color /= MPI_UNDEFINED) then
      call MPI_Comm_rank(my_mpi%my_SE_communicator, my_mpi%my_atom_rank, ierr)
    else
      my_mpi%my_atom_rank = -1
    end if

    ! Assertion:
    if ((color /= MPI_UNDEFINED) .and. &
        (my_mpi%my_atom_rank /= key)) then
        
      write(*,*) "Inconsistency in MPI rank ordering."
      stop 
    end if
    
    ! define active and inactive ranks
    my_mpi%active = 1
    key = my_mpi%my_world_rank
    if (my_mpi%my_world_rank >= &
        num_atom_ranks*num_spin_ranks*num_energy_ranks) then
      my_mpi%active = 0    
    end if
    
    call MPI_Comm_split(MPI_COMM_WORLD, my_mpi%active, key, my_mpi%my_active_communicator, ierr)
    
    ! Assertion: check if ids are correct
    if (my_mpi%my_world_rank /= mapToWorldRank(my_mpi, &
                                               my_mpi%my_atom_id, &
                                               my_mpi%my_spin_id, &
                                               my_mpi%my_energy_id)) then
      
      write(*,*) "Inconsistency in KKRnano process ids."
      stop
    end if

  end subroutine

  !--------------------------------------------------------------
  subroutine destroyKKRnanoParallel(my_mpi)
    implicit none
    include 'mpif.h'

    type (KKRnanoParallel), intent(inout) :: my_mpi

    integer :: ierr

    if (my_mpi%my_SE_communicator /= MPI_COMM_NULL) then
      call MPI_Comm_free(my_mpi%my_SE_communicator, ierr)
    end if
    call MPI_Comm_free(my_mpi%my_active_communicator, ierr)
    call MPI_Finalize(ierr)
  end subroutine

  !================== Getter routines =======================

  !--------------------------------------------------------------
  !> Returns rank of process in MPI_COMM_WORLD.
  integer function getMyWorldRank(my_mpi)
    implicit none
    type (KKRnanoParallel), intent(in) :: my_mpi
    getMyWorldRank = my_mpi%my_world_rank
  end function
  
  !--------------------------------------------------------------
  !> Returns rank of process in (Spin,Energy)-communicator.
  !> Equals (getMyAtomId() - 1)
  integer function getMyAtomRank(my_mpi)
    implicit none
    type (KKRnanoParallel), intent(in) :: my_mpi
    getMyAtomRank = my_mpi%my_atom_rank
  end function
  
  ! ----------------- Id getters --------------------------------
  
  !--------------------------------------------------------------
  !> Returns atom id.
  integer function getMyAtomId(my_mpi)
    implicit none
    type (KKRnanoParallel), intent(in) :: my_mpi
    getMyAtomId = my_mpi%my_atom_id
  end function
  
  !--------------------------------------------------------------
  !> Returns spin id.
  integer function getMySpinId(my_mpi)
    implicit none
    type (KKRnanoParallel), intent(in) :: my_mpi
    getMySpinId = my_mpi%my_spin_id
  end function
  
  !--------------------------------------------------------------
  !> Returns energy id.
  integer function getMyEnergyId(my_mpi)
    implicit none
    type (KKRnanoParallel), intent(in) :: my_mpi
    getMyEnergyId = my_mpi%my_energy_id
  end function
  
  !--------------------------------------------------------------
  !> Returns (spin,energy) id.
  integer function getMySEId(my_mpi)
    implicit none
    type (KKRnanoParallel), intent(in) :: my_mpi
    getMySEId = my_mpi%my_SE_id
  end function
  
  ! ----------- Number of ranks getters -------------------------
  
  !--------------------------------------------------------------
  !> Returns number of atom ranks.
  integer function getNumAtomRanks(my_mpi)
    implicit none
    type (KKRnanoParallel), intent(in) :: my_mpi
    getNumAtomRanks = my_mpi%num_atom_ranks_
  end function
  
  !--------------------------------------------------------------
  !> Returns number of spin ranks.
  integer function getNumSpinRanks(my_mpi)
    implicit none
    type (KKRnanoParallel), intent(in) :: my_mpi
    getNumSpinRanks = my_mpi%num_spin_ranks_
  end function
  
  !--------------------------------------------------------------
  !> Returns number of spin ranks.
  integer function getNumEnergyRanks(my_mpi)
    implicit none
    type (KKRnanoParallel), intent(in) :: my_mpi
    getNumEnergyRanks = my_mpi%num_energy_ranks_
  end function
  
  !--------------------------------------------------------------
  !> Returns number of (Spin,Energy)-ranks/groups.
  integer function getNumSERanks(my_mpi)
    implicit none
    type (KKRnanoParallel), intent(in) :: my_mpi
    getNumSERanks = my_mpi%num_energy_ranks_ * &
                    my_mpi%num_spin_ranks_
  end function
  
  !--------------------------------------------------------------
  !> Returns .true. if rank is the master rank.
  logical function isMasterRank(my_mpi)
    implicit none
    type (KKRnanoParallel), intent(in) :: my_mpi
    isMasterRank = (my_mpi%my_world_rank == 0)
  end function
  
  !------------ Other -------------------------------------------
  
  !--------------------------------------------------------------
  !> Returns .true. if rank is in the master group.
  !> means that SE-Id is 1
  logical function isInMasterGroup(my_mpi)
    implicit none
    type (KKRnanoParallel), intent(in) :: my_mpi
    isInMasterGroup = (my_mpi%my_SE_id == 1)
  end function
  
  !--------------------------------------------------------------
  !> Returns (Spin,Energy)-communicator handle of the communicator
  !> the process belongs to
  integer function getMySEcommunicator(my_mpi)
    implicit none
    type (KKRnanoParallel), intent(in) :: my_mpi
    getMySEcommunicator = my_mpi%my_SE_communicator
  end function
  
  !--------------------------------------------------------------
  !> Returns communicator handle of active/inactive ranks,
  !> according to the group the calling process belongs to
  integer function getMyActiveCommunicator(my_mpi)
    implicit none
    type (KKRnanoParallel), intent(in) :: my_mpi
    getMyActiveCommunicator = my_mpi%my_active_communicator
  end function
  
  !--------------------------------------------------------------
  !> Returns .true. if rank is an active rank
  logical function isActiveRank(my_mpi)
    implicit none
    type (KKRnanoParallel), intent(in) :: my_mpi
    isActiveRank = (my_mpi%active == 1)
  end function

  ! START RANKS FROM 0 or FROM 1 ???
  ! rank ... Naming as in MPI, starts from 0
  ! id   ... starts from 1
  !---------------------------------------------------------------
  !> Given atom_id, spin_id and energy_id, return rank of
  !> corresponding process in MPI_COMM_WORLD.
  !> Note id-s start from 1 !
  integer function mapToWorldRank(my_mpi, atom_id, spin_id, energy_id)
    implicit none
    type (KKRnanoParallel), intent(in) :: my_mpi
    integer, intent(in) :: spin_id
    integer, intent(in) :: atom_id  ! 1..N
    integer, intent(in) :: energy_id

    mapToWorldRank = mapToWorldRankSE(my_mpi, atom_id, &
                     (spin_id - 1)*my_mpi%num_energy_ranks_ + energy_id)
    
    !(((spin_id - 1) * my_mpi%num_energy_ranks_ + &
    !                    energy_id) - 1) * my_mpi%num_atom_ranks_ + atom_id - 1

  end function
  
  !---------------------------------------------------------------
  !> Given atom_id, SE_id (spin,energy)-id, return rank of
  !> corresponding process in MPI_COMM_WORLD.
  !> This is an alternative to mapToWorldRank using the
  !> combined (spin,energy)-id
  !> Note id-s start from 1 !
  integer function mapToWorldRankSE(my_mpi, atom_id, SE_id)
    implicit none
    type (KKRnanoParallel), intent(in) :: my_mpi
    integer, intent(in) :: SE_id
    integer, intent(in) :: atom_id  ! 1..N

    mapToWorldRankSE = (SE_id - 1) * my_mpi%num_atom_ranks_ + atom_id - 1

  end function


end module

! Example of usage.

!program xy
!  use KKRnano_mpi_new_mod
!  implicit none
!  type (KKRnanoParallel) :: my_mpi
!
!  call createKKRnanoParallel(my_mpi, 2, 2, 3)
!
!  write (*,*) "Worldrank", getMyWorldRank(my_mpi), &
!                           isMasterRank(my_mpi), &
!                           isActiveRank(my_mpi), &
!                           getMyAtomId(my_mpi), &
!                           getMySpinId(my_mpi), &
!                           getMyEnergyId(my_mpi), &
!                           getMySEId(my_mpi), &
!                           isInMasterGroup(my_mpi)
!
!  if (isMasterRank(my_mpi)) then
!    write(*,*) "Num Ranks: ", getNumAtomRanks(my_mpi), &
!                              getNumSpinRanks(my_mpi), &
!                              getNumEnergyRanks(my_mpi), &
!                              getNumSERanks(my_mpi)
!  end if
!
!  call destroyKKRnanoParallel(my_mpi);
!end program
