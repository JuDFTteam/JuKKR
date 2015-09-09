!------------------------------------------------------------------------------
!> This module is used to set up the rank informations+communicators
!> (world-rank, atom-rank, spin-rank, energy-rank) for KKRnano.
!>
!> It replaces the complicated setup in the original KKRnano.
!> Note that also this code could be simplified using
!> standard MPI routines (MPI_Cart_create etc.)
!>
!> For setup one creates an object of type 'KKRnanoParallel'
!> using the 'createKKRnanoParallel' routine.
!> This object is immutable to ensure that the rank information
!> is constant during the program run.
!> Use exclusively the getXY routines to retrieve rank information
!>
!> After use call 'destroyKKRnanoParallel'
!>
!> Author: Elias Rabel, 2012

! Example of usage.

!program xy
!  use KKRnanoParallel_mod
!  implicit none
!  type(KKRnanoParallel) :: my_mpi
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
!  endif
!
!  call destroyKKRnanoParallel(my_mpi);
!endprogram

module KKRnanoParallel_mod
implicit none
  private
  public :: KKRnanoParallel, create, destroy
  public :: createKKRnanoParallel, destroyKKRnanoParallel ! deprecated
  public :: getMyWorldRank, getMyAtomRank, getMasterRank 
  public :: getMySEId, getMyAtomId, getMySpinId, getMyEnergyId, getResponsibleSpinId    
  public :: getNumAtomRanks, getNumSpinRanks, getNumEnergyRanks, getNumSERanks, getNumWorldRanks
  public :: getMySEcommunicator, getMyActiveCommunicator        
  public :: isMasterRank, isInMasterGroup, isActiveRank, isWorkingSpinRank        
  public :: mapToWorldRank, mapToWorldRankSE       
  
  !> Opaque object. Use getter routines to access member variables
  type KKRnanoParallel
    PRIVATE
    
    !> MPI communicator handle
    !> This communicator contains all atom process belonging
    !> to the same (Spin,Energy) - group
    integer :: my_SE_communicator
    
    integer :: my_active_communicator
    integer :: active
    integer :: active_rank
    
    integer :: my_world_rank
    integer :: my_atom_rank  ! Note ranks start from 0 !
          
    integer :: my_spin_id    ! Note ids start from 1 !
    integer :: my_energy_id
    integer :: my_atom_id
    integer :: my_SE_id

    integer :: num_comms_
    integer :: procs_per_comm_

    integer :: num_ranks_      !< Total number of ranks in MPI_COMM_WORLD
    integer :: num_atom_ranks_
    integer :: num_spin_ranks_
    integer :: num_energy_ranks_
  endtype

  interface create
    module procedure createKKRnanoParallel
  endinterface
  
  interface destroy
    module procedure destroyKKRnanoParallel
  endinterface
  
  contains

  !---------------------------------------------------------------------------------------
  subroutine createKKRnanoParallel(my_mpi, num_atom_ranks, num_spin_ranks, num_energy_ranks)
    include 'mpif.h'

    type(KKRnanoParallel), intent(inout) :: my_mpi
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
    my_mpi%num_ranks_ = num_ranks
    call MPI_Comm_rank(MPI_COMM_WORLD, my_mpi%my_world_rank, ierr)

    ! Check if there are enough ranks.
    if (num_ranks < num_comms*procs_per_comm) then
      if (my_mpi%my_world_rank == 0) then
        write(*,*) "Not enough MPI ranks."
      endif
      call MPI_Finalize(ierr)
      stop
    endif
    
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
    endif

    call MPI_Comm_split(MPI_COMM_WORLD, color, key, my_mpi%my_SE_communicator, ierr)
    
    if (color /= MPI_UNDEFINED) then
      call MPI_Comm_rank(my_mpi%my_SE_communicator, my_mpi%my_atom_rank, ierr)
    else
      my_mpi%my_atom_rank = -1
    endif

    ! Assertion:
    if ((color /= MPI_UNDEFINED) .and. &
        (my_mpi%my_atom_rank /= key)) then
        
      write(*,*) "Inconsistency in MPI rank ordering."
      stop 
    endif
    
    ! define active and inactive ranks
    my_mpi%active = 1
    key = my_mpi%my_world_rank
    if (my_mpi%my_world_rank >= &
        num_atom_ranks*num_spin_ranks*num_energy_ranks) then
      my_mpi%active = 0    
    endif
    
    call MPI_Comm_split(MPI_COMM_WORLD, my_mpi%active, key, my_mpi%my_active_communicator, ierr)
    
    call MPI_Comm_rank(my_mpi%my_active_communicator, my_mpi%active_rank, ierr)

    ! Assertion: check if ids are correct
    if (my_mpi%my_world_rank /= mapToWorldRank(my_mpi, &
                                               my_mpi%my_atom_id, &
                                               my_mpi%my_spin_id, &
                                               my_mpi%my_energy_id)) then
      
      write(*,*) "Inconsistency in KKRnano process ids."
      stop
    endif

    ! Assertion: check if atom-rank is correct
    if (my_mpi%active == 1) then
      if (my_mpi%my_atom_rank /=  my_mpi%my_atom_id - 1) then
        write(*,*) "Inconsistency in KKRnano atom rank."
        stop
      endif
    endif

    if (my_mpi%active == 1 .and. my_mpi%active_rank /= my_mpi%my_world_rank) then
      write(*,*) "ERROR: active rank not equal to world rank."
      stop
    endif

  endsubroutine

  !--------------------------------------------------------------
  subroutine destroyKKRnanoParallel(my_mpi)
    include 'mpif.h'

    type(KKRnanoParallel), intent(inout) :: my_mpi

    integer :: ierr

    if (my_mpi%my_SE_communicator /= MPI_COMM_NULL) then
      call MPI_Comm_free(my_mpi%my_SE_communicator, ierr)
    endif
    call MPI_Comm_free(my_mpi%my_active_communicator, ierr)
    call MPI_Finalize(ierr)
  endsubroutine

  !================== Getter routines =======================

  !--------------------------------------------------------------
  !> Returns rank of process in MPI_COMM_WORLD.
  integer function getMyWorldRank(my_mpi)
    type(KKRnanoParallel), intent(in) :: my_mpi
    getMyWorldRank = my_mpi%my_world_rank
  endfunction
  
  !--------------------------------------------------------------
  !> Returns rank of process in (Spin,Energy)-communicator.
  !> Equals (getMyAtomId() - 1)
  integer function getMyAtomRank(my_mpi)
    type(KKRnanoParallel), intent(in) :: my_mpi
    getMyAtomRank = my_mpi%my_atom_rank
  endfunction
  
  ! ----------------- Id getters --------------------------------
  
  !--------------------------------------------------------------
  !> Returns atom id.
  integer function getMyAtomId(my_mpi)
    type(KKRnanoParallel), intent(in) :: my_mpi
    getMyAtomId = my_mpi%my_atom_id
  endfunction
  
  !--------------------------------------------------------------
  !> Returns spin id.
  integer function getMySpinId(my_mpi)
    type(KKRnanoParallel), intent(in) :: my_mpi
    getMySpinId = my_mpi%my_spin_id
  endfunction
  
  !--------------------------------------------------------------
  !> Returns energy id.
  integer function getMyEnergyId(my_mpi)
    type(KKRnanoParallel), intent(in) :: my_mpi
    getMyEnergyId = my_mpi%my_energy_id
  endfunction
  
  !--------------------------------------------------------------
  !> Returns (spin,energy) id.
  integer function getMySEId(my_mpi)
    type(KKRnanoParallel), intent(in) :: my_mpi
    getMySEId = my_mpi%my_SE_id
  endfunction
  
  ! ----------- Number of ranks getters -------------------------
  
  !--------------------------------------------------------------
  !> Returns number of atom ranks.
  integer function getNumAtomRanks(my_mpi)
    type(KKRnanoParallel), intent(in) :: my_mpi
    getNumAtomRanks = my_mpi%num_atom_ranks_
  endfunction
  
  !--------------------------------------------------------------
  !> Returns number of spin ranks.
  integer function getNumSpinRanks(my_mpi)
    type(KKRnanoParallel), intent(in) :: my_mpi
    getNumSpinRanks = my_mpi%num_spin_ranks_
  endfunction
  
  !--------------------------------------------------------------
  !> Returns number of energy ranks.
  integer function getNumEnergyRanks(my_mpi)
    type(KKRnanoParallel), intent(in) :: my_mpi
    getNumEnergyRanks = my_mpi%num_energy_ranks_
  endfunction
  
  !--------------------------------------------------------------
  !> Returns number of (Spin,Energy)-ranks/groups.
  integer function getNumSERanks(my_mpi)
    type(KKRnanoParallel), intent(in) :: my_mpi
    getNumSERanks = my_mpi%num_energy_ranks_ * my_mpi%num_spin_ranks_
  endfunction
  
  !--------------------------------------------------------------
  !> Returns number ranks in MPI_COMM_WORLD.
  integer function getNumWorldRanks(my_mpi)
    type(KKRnanoParallel), intent(in) :: my_mpi
    getNumWorldRanks = my_mpi%num_ranks_
  endfunction

  !--------------------------------------------------------------
  !> Returns rank number of MasterRank in MPI_COMM_WORLD.
  integer function getMasterRank(my_mpi) ! independent of my_mpi
    type(KKRnanoParallel), intent(in) :: my_mpi
    getMasterRank = 0
  endfunction

  !--------------------------------------------------------------
  !> Returns .true. if rank is the master rank.
  logical function isMasterRank(my_mpi)
    type(KKRnanoParallel), intent(in) :: my_mpi
    isMasterRank = (my_mpi%my_world_rank == 0)
  endfunction
  
  !------------ Other -------------------------------------------
  
  !--------------------------------------------------------------
  !> Returns .true. if rank is in the master group.
  !> means that SE-Id is 1
  logical function isInMasterGroup(my_mpi)
    type(KKRnanoParallel), intent(in) :: my_mpi
    isInMasterGroup = (my_mpi%my_SE_id == 1)
  endfunction
  
  !--------------------------------------------------------------
  !> Returns (Spin,Energy)-communicator handle of the communicator
  !> the process belongs to
  integer function getMySEcommunicator(my_mpi)
    type(KKRnanoParallel), intent(in) :: my_mpi
    getMySEcommunicator = my_mpi%my_SE_communicator
  endfunction
  
  !--------------------------------------------------------------
  !> Returns communicator handle of active/inactive ranks,
  !> according to the group the calling process belongs to
  integer function getMyActiveCommunicator(my_mpi)
    type(KKRnanoParallel), intent(in) :: my_mpi
    getMyActiveCommunicator = my_mpi%my_active_communicator
  endfunction
  
  !--------------------------------------------------------------
  !> Returns .true. if rank is an active rank
  logical function isActiveRank(my_mpi)
    type(KKRnanoParallel), intent(in) :: my_mpi
    isActiveRank = (my_mpi%active == 1)
  endfunction

  !--------------------------------------------------------------
  !> Given the spin index, determine which Spin-Id is responsible
  !> for this spin index.
  integer function getResponsibleSpinId(my_mpi, ispin)
    type(KKRnanoParallel), intent(in) :: my_mpi
    integer, intent(in) :: ispin

    !-------
    integer :: mapspin

    if (my_mpi%num_spin_ranks_ == 1) then
      mapspin = 1
    else
      mapspin = ispin
    endif

    getResponsibleSpinId = mapspin
  endfunction

  !--------------------------------------------------------------
  !> Given the spin index, determine whether process has to
  !> work or not.
  logical function isWorkingSpinRank(my_mpi, ispin)
    type(KKRnanoParallel), intent(in) :: my_mpi
    integer, intent(in) :: ispin

    !-------
    integer :: mapspin

    if (my_mpi%num_spin_ranks_ == 1) then
      mapspin = 1
    else
      mapspin = ispin
    endif

    isWorkingSpinRank = (my_mpi%my_spin_id == mapspin)
  endfunction

  ! START RANKS FROM 0 or FROM 1 ???
  ! rank ... Naming as in MPI, starts from 0
  ! id   ... starts from 1
  !---------------------------------------------------------------
  !> Given atom_id, spin_id and energy_id, return rank of
  !> corresponding process in MPI_COMM_WORLD.
  !> Note id-s start from 1 !
  integer function mapToWorldRank(my_mpi, atom_id, spin_id, energy_id)
    type(KKRnanoParallel), intent(in) :: my_mpi
    integer, intent(in) :: spin_id
    integer, intent(in) :: atom_id  ! 1..N
    integer, intent(in) :: energy_id

    mapToWorldRank = mapToWorldRankSE(my_mpi, atom_id, (spin_id - 1)*my_mpi%num_energy_ranks_ + energy_id)
    
    !(((spin_id - 1) * my_mpi%num_energy_ranks_ + &
    !                    energy_id) - 1) * my_mpi%num_atom_ranks_ + atom_id - 1
  endfunction
  
  !---------------------------------------------------------------
  !> Given atom_id, SE_id (spin,energy)-id, return rank of
  !> corresponding process in MPI_COMM_WORLD.
  !> This is an alternative to mapToWorldRank using the
  !> combined (spin,energy)-id
  !> Note id-s start from 1 !
  integer function mapToWorldRankSE(my_mpi, atom_id, SE_id)
    type(KKRnanoParallel), intent(in) :: my_mpi
    integer, intent(in) :: SE_id
    integer, intent(in) :: atom_id  ! 1..N

    mapToWorldRankSE = (SE_id - 1) * my_mpi%num_atom_ranks_ + atom_id - 1
  endfunction


endmodule

