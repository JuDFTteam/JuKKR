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
  
  public :: getRankAndLocalIndex
  
  contains

  function getRankAndLocalIndex(ind, num, nranks) result(rank_loc)
    integer, intent(in) :: ind, num, nranks
    integer(kind=4) :: rank_loc(2) ! result

    integer :: atoms_per_proc ! == max_local_atoms
    atoms_per_proc = (num - 1) / nranks + 1 ! == ceiling(num/dble(nranks))
    
    rank_loc(1) = (ind - 1) / atoms_per_proc ! owner rank
    rank_loc(2) = ind - rank_loc(1)*atoms_per_proc ! local index
!   rank_loc(2) = mod(ind - 1, atoms_per_proc) + 1 ! same with mod
    
#ifndef NDEBUG    
    if (ind > num) stop "ChunkIndex: fatal error: ind > num"
    if (ind < 1)   stop "ChunkIndex: fatal error: ind < 1"
    if (rank_loc(1) >= nranks) stop "ChunkIndex: fatal error: computed rank >= nranks"
    if (rank_loc(2) > atoms_per_proc) stop "ChunkIndex: fatal error: computed local index > max_local_atoms"
#endif
  endfunction ! get

endmodule ! ChunkIndex_mod

