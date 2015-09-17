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
!  type(KKRnanoParallel) :: self
!
!  call createKKRnanoParallel(self, 2, 2, 3)
!
!  write (*,*) "Worldrank", getMyWorldRank(self), &
!                           isMasterRank(self), &
!                           isActiveRank(self), &
!                           getMyAtomId(self), &
!                           getMySpinId(self), &
!                           getMyEnergyId(self), &
!                           getMySEId(self), &
!                           isInMasterGroup(self)
!
!  if (isMasterRank(self)) then
!    write(*,*) "Num Ranks: ", getNumAtomRanks(self), &
!                              getNumSpinRanks(self), &
!                              getNumEnergyRanks(self), &
!                              getNumSERanks(self)
!  endif
!
!  call destroyKKRnanoParallel(self);
!endprogram

module KKRnanoParallel_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
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
    private
    
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

    integer :: num_comms
    integer :: procs_per_comm

    integer :: num_all_ranks      !< Total number of ranks in MPI_COMM_WORLD
    integer :: num_atom_ranks
    integer :: num_spin_ranks
    integer :: num_energy_ranks
  endtype

  interface create
    module procedure createKKRnanoParallel
  endinterface
  
  interface destroy
    module procedure destroyKKRnanoParallel
  endinterface
  
  contains

  !---------------------------------------------------------------------------------------
  subroutine createKKRnanoParallel(self, num_atom_ranks, num_spin_ranks, num_energy_ranks)
    include 'mpif.h'

    type(KKRnanoParallel), intent(inout) :: self
    integer, intent(in) :: num_atom_ranks
    integer, intent(in) :: num_spin_ranks
    integer, intent(in) :: num_energy_ranks

    integer :: n, ierr, color, key
 
!   call MPI_Init(ierr) ! must have been call earlier
    call MPI_Comm_size(MPI_COMM_WORLD, self%num_all_ranks, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, self%my_world_rank,  ierr)
 
    self%num_spin_ranks   = max(1, num_spin_ranks)
    self%num_energy_ranks = max(1, num_energy_ranks)
    self%num_comms = self%num_spin_ranks * self%num_energy_ranks
     
    if (num_atom_ranks < 1) then
      ! the number of atomic ranks is chosen atomatically
      self%num_atom_ranks = self%num_all_ranks / self%num_comms
      ! todo: warn if the numbers are not divisible
      if (self%num_atom_ranks * self%num_comms /= self%num_all_ranks) &
        warn(6, "Number of all MPI ranks ="+self%num_all_ranks+"is not divisible by"+self%num_comms)
    else
      self%num_atom_ranks = num_atom_ranks ! controlled via the input file
    endif
    self%procs_per_comm = self%num_atom_ranks

    n = self%num_comms * self%procs_per_comm ! abbrev.
    if (self%num_all_ranks < n) & ! Check if there are enough ranks
      die_here("found only "+self%num_all_ranks+"MPI ranks, but requires"+self%num_comms+"*"+self%procs_per_comm+"="+n+"ranks!")
    
    color = self%my_world_rank / self%procs_per_comm
    key = self%my_world_rank - (color * self%procs_per_comm)
    
    !---------- Define ids ---------------
    ! Note: ids start with 1
    ! TODO
    self%my_atom_id = key + 1
    self%my_SE_id = color + 1
    self%my_spin_id = (color / self%num_energy_ranks) + 1
    self%my_energy_id = mod(color, self%num_energy_ranks) + 1
    
    if (color >= self%num_comms) color = MPI_UNDEFINED ! inactive process

    call MPI_Comm_split(MPI_COMM_WORLD, color, key, self%my_SE_communicator, ierr)
    
    self%my_atom_rank = -1
    if (color /= MPI_UNDEFINED) then
      call MPI_Comm_rank(self%my_SE_communicator, self%my_atom_rank, ierr)
      if (self%my_atom_rank /= key) die_here("Inconsistency in MPI rank ordering!") ! Assertion:
    endif ! color is defined
    
    self%active = 1 ! define active and inactive ranks
    key = self%my_world_rank
    if (self%my_world_rank >= self%num_atom_ranks*self%num_spin_ranks*self%num_energy_ranks) self%active = 0    
    
    call MPI_Comm_split(MPI_COMM_WORLD, self%active, key, self%my_active_communicator, ierr)
    
    call MPI_Comm_rank(self%my_active_communicator, self%active_rank, ierr)

    ! Assertion: check if ids are correct
    if (self%my_world_rank /= mapToWorldRank(self, self%my_atom_id, self%my_spin_id, self%my_energy_id)) &
      die_here("Inconsistency in KKRnano process ids!")

    ! Assertion: check if atom-rank is correct
    if (self%active == 1 .and. self%my_atom_rank /= self%my_atom_id-1) die_here("Inconsistency in KKRnano atom rank!")
    if (self%active == 1 .and. self%active_rank /= self%my_world_rank) die_here("Active rank not equal to world rank!")

  endsubroutine

  !--------------------------------------------------------------
  subroutine destroyKKRnanoParallel(self)
    include 'mpif.h'

    type(KKRnanoParallel), intent(inout) :: self
    integer :: ierr

    if (self%my_SE_communicator /= MPI_COMM_NULL) &
      call MPI_Comm_free(self%my_SE_communicator, ierr)
    call MPI_Comm_free(self%my_active_communicator, ierr)
    call MPI_Finalize(ierr)
  endsubroutine

  !================== Getter routines =======================

  !--------------------------------------------------------------
  !> Returns rank of process in MPI_COMM_WORLD.
  integer function getMyWorldRank(self)
    type(KKRnanoParallel), intent(in) :: self
    getMyWorldRank = self%my_world_rank
  endfunction
  
  !--------------------------------------------------------------
  !> Returns rank of process in (Spin,Energy)-communicator.
  !> Equals (getMyAtomId() - 1)
  integer function getMyAtomRank(self)
    type(KKRnanoParallel), intent(in) :: self
    getMyAtomRank = self%my_atom_rank
  endfunction
  
  ! ----------------- Id getters --------------------------------
  
  !--------------------------------------------------------------
  !> Returns atom id.
  integer function getMyAtomId(self)
    type(KKRnanoParallel), intent(in) :: self
    getMyAtomId = self%my_atom_id
  endfunction
  
  !--------------------------------------------------------------
  !> Returns spin id.
  integer function getMySpinId(self)
    type(KKRnanoParallel), intent(in) :: self
    getMySpinId = self%my_spin_id
  endfunction
  
  !--------------------------------------------------------------
  !> Returns energy id.
  integer function getMyEnergyId(self)
    type(KKRnanoParallel), intent(in) :: self
    getMyEnergyId = self%my_energy_id
  endfunction
  
  !--------------------------------------------------------------
  !> Returns (spin,energy) id.
  integer function getMySEId(self)
    type(KKRnanoParallel), intent(in) :: self
    getMySEId = self%my_SE_id
  endfunction
  
  ! ----------- Number of ranks getters -------------------------
  
  !--------------------------------------------------------------
  !> Returns number of atom ranks.
  integer function getNumAtomRanks(self)
    type(KKRnanoParallel), intent(in) :: self
    getNumAtomRanks = self%num_atom_ranks
  endfunction
  
  !--------------------------------------------------------------
  !> Returns number of spin ranks.
  integer function getNumSpinRanks(self)
    type(KKRnanoParallel), intent(in) :: self
    getNumSpinRanks = self%num_spin_ranks
  endfunction
  
  !--------------------------------------------------------------
  !> Returns number of energy ranks.
  integer function getNumEnergyRanks(self)
    type(KKRnanoParallel), intent(in) :: self
    getNumEnergyRanks = self%num_energy_ranks
  endfunction
  
  !--------------------------------------------------------------
  !> Returns number of (Spin,Energy)-ranks/groups.
  integer function getNumSERanks(self)
    type(KKRnanoParallel), intent(in) :: self
    getNumSERanks = self%num_energy_ranks * self%num_spin_ranks
  endfunction
  
  !--------------------------------------------------------------
  !> Returns number ranks in MPI_COMM_WORLD.
  integer function getNumWorldRanks(self)
    type(KKRnanoParallel), intent(in) :: self
    getNumWorldRanks = self%num_all_ranks
  endfunction

  !--------------------------------------------------------------
  !> Returns rank number of MasterRank in MPI_COMM_WORLD.
  integer function getMasterRank(self) ! independent of self
    type(KKRnanoParallel), intent(in) :: self
    getMasterRank = 0
  endfunction

  !--------------------------------------------------------------
  !> Returns .true. if rank is the master rank.
  logical function isMasterRank(self)
    type(KKRnanoParallel), intent(in) :: self
    isMasterRank = (self%my_world_rank == 0)
  endfunction
  
  !------------ Other -------------------------------------------
  
  !--------------------------------------------------------------
  !> Returns .true. if rank is in the master group.
  !> means that SE-Id is 1
  logical function isInMasterGroup(self)
    type(KKRnanoParallel), intent(in) :: self
    isInMasterGroup = (self%my_SE_id == 1)
  endfunction
  
  !--------------------------------------------------------------
  !> Returns (Spin,Energy)-communicator handle of the communicator
  !> the process belongs to
  integer function getMySEcommunicator(self)
    type(KKRnanoParallel), intent(in) :: self
    getMySEcommunicator = self%my_SE_communicator
  endfunction
  
  !--------------------------------------------------------------
  !> Returns communicator handle of active/inactive ranks,
  !> according to the group the calling process belongs to
  integer function getMyActiveCommunicator(self)
    type(KKRnanoParallel), intent(in) :: self
    getMyActiveCommunicator = self%my_active_communicator
  endfunction
  
  !--------------------------------------------------------------
  !> Returns .true. if rank is an active rank
  logical function isActiveRank(self)
    type(KKRnanoParallel), intent(in) :: self
    isActiveRank = (self%active == 1)
  endfunction

  !--------------------------------------------------------------
  !> Given the spin index, determine which Spin-Id is responsible
  !> for this spin index.
  integer function getResponsibleSpinId(self, ispin)
    type(KKRnanoParallel), intent(in) :: self
    integer, intent(in) :: ispin
    getResponsibleSpinId = ispin; if (self%num_spin_ranks == 1) getResponsibleSpinId = 1
  endfunction

  !--------------------------------------------------------------
  !> Given the spin index, determine whether process has to
  !> work or not.
  logical function isWorkingSpinRank(self, ispin)
    type(KKRnanoParallel), intent(in) :: self
    integer, intent(in) :: ispin
    isWorkingSpinRank = (self%my_spin_id == getResponsibleSpinId(self, ispin))
  endfunction

  ! START RANKS FROM 0 or FROM 1 ???
  ! rank ... Naming as in MPI, starts from 0
  ! id   ... starts from 1
  !---------------------------------------------------------------
  !> Given atom_id, spin_id and energy_id, return rank of
  !> corresponding process in MPI_COMM_WORLD.
  !> Note id-s start from 1 !
  integer function mapToWorldRank(self, atom_id, spin_id, energy_id)
    type(KKRnanoParallel), intent(in) :: self
    integer, intent(in) :: spin_id, energy_id, atom_id  ! 1..N
    mapToWorldRank = mapToWorldRankSE(self, atom_id, (spin_id - 1)*self%num_energy_ranks + energy_id)
    !(((spin_id - 1) * self%num_energy_ranks + energy_id) - 1) * self%num_atom_ranks + atom_id - 1
  endfunction
  
  !---------------------------------------------------------------
  !> Given atom_id, SE_id (spin,energy)-id, return rank of
  !> corresponding process in MPI_COMM_WORLD.
  !> This is an alternative to mapToWorldRank using the
  !> combined (spin,energy)-id
  !> Note id-s start from 1 !
  integer function mapToWorldRankSE(self, atom_id, SE_id)
    type(KKRnanoParallel), intent(in) :: self
    integer, intent(in) :: SE_id, atom_id  ! 1..N
    mapToWorldRankSE = (SE_id - 1) * self%num_atom_ranks + atom_id - 1
  endfunction

endmodule

