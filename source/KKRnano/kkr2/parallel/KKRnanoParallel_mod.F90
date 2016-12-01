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
 
!   call MPI_Init(ierr) ! must have been called earlier
    call MPI_Comm_size(MPI_COMM_WORLD, self%numWorldRanks, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, self%myWorldRank, ierr)
    self%isMasterRank = (self%myWorldRank == 0)
 
    self%numSpinRanks = max(1, numSpinRanks)
    self%numEnergyRanks = max(1, numEnergyRanks)
    self%numComms = self%numSpinRanks * self%numEnergyRanks
     
    if (numAtomRanks < 1) then
      ! the number of atomic ranks is chosen automatically
      self%numAtomRanks = self%numWorldRanks / self%numComms
      if (self%numAtomRanks * self%numComms /= self%numWorldRanks) &
        warn(6, "Number of all MPI ranks ="+self%numWorldRanks+"is not divisible by"+self%numComms)
    else
      self%numAtomRanks = numAtomRanks ! controlled via the input file
    endif

    n = self%numComms * self%numAtomRanks ! abbrev.
    if (self%numWorldRanks < n) & ! Check if there are enough ranks
      die_here("found only "+self%numWorldRanks+"MPI ranks, but requires"+self%numComms+"*"+self%numAtomRanks+"="+n+"ranks!")
    
    color = self%myWorldRank / self%numAtomRanks
    key = self%myWorldRank - (color * self%numAtomRanks)
    
    !---------- Define Ids ---------------
    ! Note: Ids start with 1
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
    active = 0 ; if (self%isActiveRank) active = 1 ! convert logical to integer
 
    call MPI_Comm_split(MPI_COMM_WORLD, active, key, self%myActiveComm, ierr)
 
    call MPI_Comm_rank(self%myActiveComm, active_rank, ierr)

    ! Assertion: check if Ids are correct
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

  
  !--------------------------------------------------------------
  !> Returns number of (Spin,Energy)-ranks/groups.
  integer function getNumSERanks(self)
    type(KKRnanoParallel), intent(in) :: self
    getNumSERanks = self%numEnergyRanks * self%numSpinRanks
  endfunction

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

