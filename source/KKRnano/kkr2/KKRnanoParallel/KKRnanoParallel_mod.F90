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
!  use KKRnanoParallel_mod, only: KKRnanoParallel
!  implicit none
!  type(KKRnanoParallel) :: mp
!
!  call create(mp, 2, 2, 3)
!  write (*,*) "Worldrank", mp%myWorldRank, mp%isMasterRank, mp%isActiveRank, mp%myAtomId, mp%mySpinId, mp%myEnergyId, mp%mySEId, mp%isInMasterGroup 
!
!  if (mp%isMasterRank) &
!    write(*,*) "Num Ranks: ", mp%numAtomRanks, mp%numSpinRanks, mp%numEnergyRanks, mp%getNumSERanks()
!  call destroy(mp);
!endprogram

module KKRnanoParallel_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
implicit none
  private
  public :: KKRnanoParallel, create, destroy
  public :: getResponsibleSpinId, getNumSERanks, isWorkingSpinRank, mapToWorldRank, mapToWorldRankSE

  type KKRnanoParallel
    !> MPI communicator handle
    !> This communicator contains all atom process belonging
    !> to the same (Spin,Energy) - group
    integer :: mySEComm
    
    integer :: myWorldRank
    integer :: myAtomRank  ! Note ranks start from 0 !
          
    integer :: mySpinId    ! Note Ids start from 1 !
    integer :: myEnergyId
    integer :: myAtomId
    integer :: mySEId

    integer :: procsPerComm
    integer :: numComms
    integer :: numWorldRanks      !< Total number of ranks in MPI_COMM_WORLD
    integer :: numAtomRanks
    integer :: numSpinRanks
    integer :: numEnergyRanks
    
    integer :: myActiveComm
    
    logical :: isActiveRank, isMasterRank, isInMasterGroup
    
  endtype

  interface create
    module procedure createKKRnanoParallel
  endinterface
  
  interface destroy
    module procedure destroyKKRnanoParallel
  endinterface
  
  contains

  !---------------------------------------------------------------------------------------
  subroutine createKKRnanoParallel(self, numAtomRanks, numSpinRanks, numEnergyRanks)
    include 'mpif.h'

    type(KKRnanoParallel), intent(inout) :: self
    integer, intent(in) :: numAtomRanks
    integer, intent(in) :: numSpinRanks
    integer, intent(in) :: numEnergyRanks

    integer :: n, ierr, color, key, active_rank, active
 
!   call MPI_Init(ierr) ! must have been call earlier
    call MPI_Comm_size(MPI_COMM_WORLD, self%numWorldRanks, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, self%myWorldRank,  ierr)
    self%isMasterRank = (self%myWorldRank == 0)
 
    self%numSpinRanks   = max(1, numSpinRanks)
    self%numEnergyRanks = max(1, numEnergyRanks)
    self%numComms = self%numSpinRanks * self%numEnergyRanks
     
    if (numAtomRanks < 1) then
      ! the number of atomic ranks is chosen atomatically
      self%numAtomRanks = self%numWorldRanks / self%numComms
      ! todo: warn if the numbers are not divisible
      if (self%numAtomRanks * self%numComms /= self%numWorldRanks) &
        warn(6, "Number of all MPI ranks ="+self%numWorldRanks+"is not divisible by"+self%numComms)
    else
      self%numAtomRanks = numAtomRanks ! controlled via the input file
    endif
    self%procsPerComm = self%numAtomRanks

    n = self%numComms * self%procsPerComm ! abbrev.
    if (self%numWorldRanks < n) & ! Check if there are enough ranks
      die_here("found only "+self%numWorldRanks+"MPI ranks, but requires"+self%numComms+"*"+self%procsPerComm+"="+n+"ranks!")
    
    color = self%myWorldRank / self%procsPerComm
    key = self%myWorldRank - (color * self%procsPerComm)
    
    !---------- Define ids ---------------
    ! Note: ids start with 1
    ! TODO
    self%myAtomId = key + 1
    self%mySEId = color + 1
    self%isInMasterGroup = (self%mySEId == 1)

    self%mySpinId = (color / self%numEnergyRanks) + 1
    self%myEnergyId = mod(color, self%numEnergyRanks) + 1
    
    if (color >= self%numComms) color = MPI_UNDEFINED ! inactive process

    call MPI_Comm_split(MPI_COMM_WORLD, color, key, self%mySEComm, ierr)
    
    self%myAtomRank = -1
    if (color /= MPI_UNDEFINED) then
      call MPI_Comm_rank(self%mySEComm, self%myAtomRank, ierr)
      if (self%myAtomRank /= key) die_here("Inconsistency in MPI rank ordering!") ! Assertion:
    endif ! color is defined
    
    key = self%myWorldRank
    self%isActiveRank = (self%myWorldRank < self%numAtomRanks*self%numSpinRanks*self%numEnergyRanks) ! define active and inactive ranks
    active = 0 ; if( self%isActiveRank ) active = 1 ! convert logical to integer
    
    call MPI_Comm_split(MPI_COMM_WORLD, active, key, self%myActiveComm, ierr)
    
    call MPI_Comm_rank(self%myActiveComm, active_rank, ierr)

    ! Assertion: check if ids are correct
    if (self%myWorldRank /= mapToWorldRank(self, self%myAtomId, self%mySpinId, self%myEnergyId)) &
      die_here("Inconsistency in KKRnano process ids!")

    ! Assertion: check if atom-rank is correct
    if (self%isActiveRank) then
      if (self%myAtomRank /= self%myAtomId - 1) die_here("Inconsistency in KKRnano atom rank!")
      if (active_rank /= self%myWorldRank)      die_here("Active rank not equal to world rank!")
    endif ! active
    
  endsubroutine ! create

  !--------------------------------------------------------------
  subroutine destroyKKRnanoParallel(self)
    include 'mpif.h'

    type(KKRnanoParallel), intent(inout) :: self
    integer :: ierr

    if (self%mySEComm /= MPI_COMM_NULL) call MPI_Comm_free(self%mySEComm, ierr)
    call MPI_Comm_free(self%myActiveComm, ierr)
    call MPI_Finalize(ierr)
  endsubroutine ! destroy

  !================== Getter routines =======================

  !--------------------------------------------------------------
  !> Returns rank of process in MPI_COMM_WORLD.
!   integer function getMyWorldRank(self)
!     type(KKRnanoParallel), intent(in) :: self
!     getMyWorldRank = self%myWorldRank
!   endfunction
  
  !--------------------------------------------------------------
  !> Returns rank of process in (Spin,Energy)-communicator.
  !> Equals (getMyAtomId() - 1)
!   integer function getMyAtomRank(self)
!     type(KKRnanoParallel), intent(in) :: self
!     getMyAtomRank = self%myAtomRank
!   endfunction
  
  ! ----------------- Id getters --------------------------------
  
  !--------------------------------------------------------------
  !> Returns atom id.
!   integer function getMyAtomId(self)
!     type(KKRnanoParallel), intent(in) :: self
!     getMyAtomId = self%myAtomId
!   endfunction
  
  !--------------------------------------------------------------
  !> Returns spin id.
!   integer function getMySpinId(self)
!     type(KKRnanoParallel), intent(in) :: self
!     getMySpinId = self%mySpinId
!   endfunction
  
  !--------------------------------------------------------------
  !> Returns energy id.
!   integer function getMyEnergyId(self)
!     type(KKRnanoParallel), intent(in) :: self
!     getMyEnergyId = self%myEnergyId
!   endfunction
  
  !--------------------------------------------------------------
  !> Returns (spin,energy) id.
!   integer function getMySEId(self)
!     type(KKRnanoParallel), intent(in) :: self
!     getMySEId = self%mySEId
!   endfunction
  
  ! ----------- Number of ranks getters -------------------------
  
  !--------------------------------------------------------------
  !> Returns number of atom ranks.
!   integer function getNumAtomRanks(self)
!     type(KKRnanoParallel), intent(in) :: self
!     getNumAtomRanks = self%numAtomRanks
!   endfunction
  
  !--------------------------------------------------------------
  !> Returns number of spin ranks.
!   integer function getNumSpinRanks(self)
!     type(KKRnanoParallel), intent(in) :: self
!     getNumSpinRanks = self%numSpinRanks
!   endfunction
  
  !--------------------------------------------------------------
  !> Returns number of energy ranks.
!   integer function getNumEnergyRanks(self)
!     type(KKRnanoParallel), intent(in) :: self
!     getNumEnergyRanks = self%numEnergyRanks
!   endfunction
  
  !--------------------------------------------------------------
  !> Returns number of (Spin,Energy)-ranks/groups.
  integer function getNumSERanks(self)
    type(KKRnanoParallel), intent(in) :: self
    getNumSERanks = self%numEnergyRanks * self%numSpinRanks
  endfunction
  
  !--------------------------------------------------------------
  !> Returns number ranks in MPI_COMM_WORLD.
!   integer function getNumWorldRanks(self)
!     type(KKRnanoParallel), intent(in) :: self
!     getNumWorldRanks = self%numWorldRanks
!   endfunction

  !--------------------------------------------------------------
  !> Returns .true. if rank is the master rank.
!   logical function isMasterRank(self)
!     type(KKRnanoParallel), intent(in) :: self
!     isMasterRank = (self%myWorldRank == 0)
!   endfunction
  
  !------------ Other -------------------------------------------
  
  !--------------------------------------------------------------
  !> Returns .true. if rank is in the master group.
  !> means that SE-Id is 1
!   logical function isInMasterGroup(self)
!     type(KKRnanoParallel), intent(in) :: self
!     isInMasterGroup = (self%mySEId == 1)
!   endfunction
  
  !--------------------------------------------------------------
  !> Returns (Spin,Energy)-communicator handle of the communicator
  !> the process belongs to
!   integer function getMySEcommunicator(self)
!     type(KKRnanoParallel), intent(in) :: self
!     getMySEcommunicator = self%mySEComm
!   endfunction
  
  !--------------------------------------------------------------
  !> Returns communicator handle of active/inactive ranks,
  !> according to the group the calling process belongs to
!   integer function getMyActiveCommunicator(self)
!     type(KKRnanoParallel), intent(in) :: self
!     getMyActiveCommunicator = self%myActiveComm
!   endfunction
  
  !--------------------------------------------------------------
  !> Returns .true. if rank is an active rank
!   logical function isActiveRank(self)
!     type(KKRnanoParallel), intent(in) :: self
!     isActiveRank = (self%active == 1)
!   endfunction

  !--------------------------------------------------------------
  !> Given the spin index, determine which Spin-Id is responsible
  !> for this spin index.
  integer function getResponsibleSpinId(self, ispin)
    type(KKRnanoParallel), intent(in) :: self
    integer, intent(in) :: ispin
    getResponsibleSpinId = ispin; if (self%numSpinRanks == 1) getResponsibleSpinId = 1
  endfunction

  !--------------------------------------------------------------
  !> Given the spin index, determine whether process has to
  !> work or not.
  logical function isWorkingSpinRank(self, ispin)
    type(KKRnanoParallel), intent(in) :: self
    integer, intent(in) :: ispin
    isWorkingSpinRank = (self%mySpinId == getResponsibleSpinId(self, ispin))
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
    mapToWorldRank = mapToWorldRankSE(self, atom_id, (spin_id - 1)*self%numEnergyRanks + energy_id)
    !(((spin_id - 1) * self%numEnergyRanks + energy_id) - 1) * self%numAtomRanks + atom_id - 1
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
    mapToWorldRankSE = (SE_id - 1) * self%numAtomRanks + atom_id - 1
  endfunction

endmodule ! KKRnanoParallel_mod

