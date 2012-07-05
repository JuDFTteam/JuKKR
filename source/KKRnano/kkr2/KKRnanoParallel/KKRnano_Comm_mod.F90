#include "../DebugHelpers/test_macros.h"

module KKRnano_Comm_mod
  implicit none

  CONTAINS

  !----------------------------------------------------------------------------
  !> Collect the results from the multiple scattering part at
  !> the corresponding atom process of the master group.
  !> Master Group: (Spin, Energy)-id = 1
  subroutine collectMultScatteringResults(my_mpi, GMATN_ALL, LLY_GRDT_ALL, EPROC)
    use KKRnanoParallel_mod
    use comm_patternsZ_mod
    implicit none

    type (KKRnanoParallel), intent(in) :: my_mpi

    double complex, dimension(:,:,:,:), intent(inout) ::  GMATN_ALL
    double complex, dimension(:,:), intent(inout) ::  LLY_GRDT_ALL
    integer, dimension(:), intent(in) :: EPROC

    integer :: iemxd
    integer :: nspind
    integer :: lmmaxd

    integer :: ie
    integer :: ispin

    integer :: my_atom_id
    integer :: my_world_rank

    integer, dimension(:), allocatable :: owning_ranks
    integer :: receiver

    lmmaxd = size(GMATN_ALL,1)
    iemxd = size(GMATN_ALL,3)
    nspind = size(GMATN_ALL,4)

    ASSERT(lmmaxd == size(GMATN_ALL, 2))
    ASSERT(iemxd == size(LLY_GRDT_ALL, 1))
    ASSERT(nspind == size(LLY_GRDT_ALL, 2))

    ! TODO: check allocate
    allocate(owning_ranks(iemxd * nspind))

    my_atom_id = getMyAtomId(my_mpi)

    do ispin = 1, nspind
      do ie = 1, iemxd

        owning_ranks( (ispin-1)*iemxd + ie ) = mapToWorldRank(my_mpi, my_atom_id, ispin, EPROC(ie))

      end do
    end do

    receiver = mapToWorldRankSE(my_mpi, my_atom_id, 1)

    my_world_rank = getMyWorldRank(my_mpi)

    call comm_gatherZ(my_world_rank, GMATN_ALL, lmmaxd*lmmaxd, owning_ranks, receiver)
    call comm_gatherZ(my_world_rank, LLY_GRDT_ALL, 1, owning_ranks, receiver)

    ! TODO: check dimensions

    deallocate(owning_ranks)

  end subroutine

  !----------------------------------------------------------------------------
  subroutine redistributeInitialGuess(my_mpi, PRSC, EPROC, EPROCO, KMESH, NofKs)
    use KKRnanoParallel_mod
    use comm_patternsC_mod
    implicit none

    type (KKRnanoParallel), intent(in) :: my_mpi
    complex, dimension(:,:), intent(inout) :: PRSC

    integer, dimension(:), intent(in) :: EPROC
    integer, dimension(:), intent(in) :: EPROCO
    integer, dimension(:), intent(in) :: KMESH
    integer, dimension(:), intent(in) :: NofKs

    !--------
    integer, dimension(:), allocatable :: blocksizes
    integer, dimension(:), allocatable :: old_owners
    integer, dimension(:), allocatable :: new_owners

    integer :: my_world_rank
    integer :: my_atom_id
    integer :: my_spin_id
    integer :: single_block_size
    integer :: num_energy_iemxd
    integer :: ie

    my_atom_id = getMyAtomId(my_mpi)
    my_world_rank = getMyWorldRank(my_mpi)
    my_spin_id  = getMySpinId(my_mpi)

    num_energy_iemxd = size(EPROC)
    single_block_size = size(PRSC, 1)

    ASSERT(num_energy_iemxd == size(EPROCO))
    ASSERT(num_energy_iemxd == size(KMESH))

    allocate(blocksizes(num_energy_iemxd))
    allocate(old_owners(num_energy_iemxd))
    allocate(new_owners(num_energy_iemxd))

    do ie = 1, num_energy_iemxd
      old_owners(ie) = mapToWorldRank(my_mpi, my_atom_id, my_spin_id, EPROCO(ie))
      new_owners(ie) = mapToWorldRank(my_mpi, my_atom_id, my_spin_id, EPROC(ie))
      ASSERT(KMESH(ie) < size(NofKs))
      blocksizes(ie) = single_block_size * NofKs(KMESH(IE))
    end do

    call comm_redistributeVC(my_world_rank, PRSC, blocksizes, old_owners, new_owners)

    deallocate(blocksizes)
    deallocate(old_owners)
    deallocate(new_owners)

  end subroutine

end module KKRnano_Comm_mod
