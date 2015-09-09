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
!> ChunkIndex = (rank, local index)
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
!> to convert to a 'ChunkIndex'
!>
!> \endverbatim


module ChunkIndex_mod
  implicit none
  private
  
  public :: ChunkIndex, getOwner, getLocalInd, getChunkIndex
  type ChunkIndex
    integer :: owner
    integer :: local_ind
  endtype

  contains

!------------------------------------------------------------------------------
!> Returns number of rank that owns atom/matrix/chunk with index 'ind'.
!> @param num Total number of chunks/atoms/matrices
!> @param nranks total number of ranks
integer function getOwner(ind, num, nranks)
  integer, intent(in) :: ind, num, nranks

  integer :: atoms_per_proc  

  atoms_per_proc = num / nranks  
  getOwner = (ind - 1) / atoms_per_proc
  ! 0 ... nranks-1  
endfunction getOwner

!------------------------------------------------------------------------------
!> Returns local index (on owning rank) of atom/matrix/chunk with index 'ind'.
!> @param num Total number of chunks/atoms/matrices
integer function getLocalInd(ind, num, nranks)
  integer, intent(in) :: ind, num, nranks

  integer :: atoms_per_proc  

  atoms_per_proc = num / nranks
  getLocalInd = mod((ind - 1), atoms_per_proc) + 1

  ! 1 ... atoms_per_proc
  
endfunction getLocalInd

!------------------------------------------------------------------------------
!> Returns chunk index of atom/matrix/chunk with index 'ind'.
!> @param ind    "atom"-index
!> @param num    Total number of chunks/atoms/matrices
!> @param nranks number of ranks
function getChunkIndex(ind, num, nranks)
  type (ChunkIndex) :: getChunkIndex
  integer, intent(in) :: ind, num, nranks

  getChunkIndex%owner = getOwner(ind, num, nranks)
  getChunkIndex%local_ind = getLocalInd(ind, num, nranks)

endfunction getChunkIndex

endmodule ! ChunkIndex_mod
 
