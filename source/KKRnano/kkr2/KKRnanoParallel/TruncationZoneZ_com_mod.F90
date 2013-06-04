!------------------------------------------------------------------------------
!> Module for communication within truncation zone.
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
!> \endverbatim

#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; STOP; endif

! the following macro is redefined later
#define NUMBERZ integer

#define NUMBERZ double complex
#define NUMBERC complex
#define NUMBERD double precision
#define NUMBERI integer

module TruncationZoneZ_com_mod

  implicit none

  CONTAINS

!------------------------------------------------------------------------------
!> Routine for exchanging data within truncation zone only
!> size of local buffer  :  chunk_size*num_local_atoms
!> size of receive buffer:  chunk_size*naez_trc
!> Uses MPI-RMA
subroutine exchangeZ_com(receive_buf, local_buf, trunc_zone, chunk_size, num_local_atoms, communicator)
  use CalculationData_mod
  use TruncationZone_mod
  use one_sided_commZ_mod
  implicit none

  include 'mpif.h'

  NUMBERZ, dimension(*), intent(inout) :: receive_buf     ! receive
  NUMBERZ, dimension(*), intent(inout) :: local_buf ! send
  integer, intent(in) :: communicator
  type (TruncationZone), intent(in) :: trunc_zone
  integer, intent(in) :: chunk_size
  integer, intent(in) :: num_local_atoms

  type (ChunkIndex), dimension(:), allocatable :: chunk_inds
  integer :: ii
  integer :: ierr
  integer :: naez_trc ! number of atoms in trunc. zone
  integer :: naez
  integer :: nranks
  integer :: atom_requested
  integer :: win

  naez_trc = trunc_zone%naez_trc

  call MPI_Comm_size(communicator, nranks, ierr)

  naez = num_local_atoms * nranks
  CHECKASSERT( naez == size(trunc_zone%index_map) ) !TODO: remove

  allocate(chunk_inds(naez_trc))

  do ii = 1, naez_trc
    ! get 'real' atom index, not truncation zone atom index!
    atom_requested = trunc_zone%trunc2atom_index(ii)
    chunk_inds(ii)%owner = getOwner(atom_requested, naez, nranks)
    chunk_inds(ii)%local_ind = getLocalInd(atom_requested, naez, nranks)
  end do

  !call exposeBufferZ(win, local_buf, size(local_buf), chunk_size, communicator)
  call exposeBufferZ(win, local_buf, chunk_size*num_local_atoms, chunk_size, communicator)
  call copyChunksZ(receive_buf, win, chunk_inds, chunk_size)
  call hideBufferZ(win)

  deallocate(chunk_inds)

end subroutine


end module TruncationZoneZ_com_mod
