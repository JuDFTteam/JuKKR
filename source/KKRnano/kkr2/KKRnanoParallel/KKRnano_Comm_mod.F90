#include "../DebugHelpers/test_macros.h"

module KKRnano_Comm_mod
  implicit none

  CONTAINS

  !----------------------------------------------------------------------------
  !> Print informative message about KKRnano parallelisation.
  subroutine printKKRnanoInfo(my_mpi, nthrds)
    use KKRnanoParallel_mod
    implicit none

    type (KKRnanoParallel), intent(in) :: my_mpi
    integer, intent(in) :: nthrds

    !$ integer::omp_get_max_threads
    integer :: natoms, nspin, nenergy

    if (isMasterRank(my_mpi)) then

      natoms = getNumAtomRanks(my_mpi)
      nspin  = getNumSpinRanks(my_mpi)
      nenergy = getNumEnergyRanks(my_mpi)

      write(*,'(79("="))')
      write(*,*) '  total no. of MPI ranks  = ', getNumWorldRanks(my_mpi)
      !$omp parallel
      !$omp single
      !$  write(*,*) '  OMP max. threads        = ',omp_get_max_threads()
      !$omp end single
      !$omp end parallel
      write(*,*) '  groups of processes created'
      write(*,*) '  NMPI                    = ',natoms
      write(*,*) '  SMPI                    = ',nspin
      write(*,*) '  EMPI                    = ',nenergy
      write(*,*) '  NTHRDS                  = ',NTHRDS
      write(*,*) '  MPI-processes (active)  = ',natoms*nspin*nenergy
      write(*,*) '  total no. of tasks      = ',natoms*nspin*nenergy*nthrds
      write(*,'(79("="))')

    end if

  end subroutine

  !----------------------------------------------------------------------------
  !> Collect the results from the multiple scattering part at
  !> the corresponding atom process of the master group.
  !> Master Group: (Spin, Energy)-id = 1
  subroutine collectMSResults_com(my_mpi, GMATN_ALL, LLY_GRDT_ALL, EPROC)
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
  subroutine redistributeInitialGuess_com(my_mpi, PRSC, EPROC, EPROCO, KMESH, NofKs)
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

  !----------------------------------------------------------------------------
  !> Communicate Potential from Master-Group to all other (Spin,Energy)-Groups.
  subroutine communicatePotential(my_mpi, VISP, VINS, ECORE)
    use KKRnanoParallel_mod
    use comm_patternsD_mod
    implicit none

    type (KKRnanoParallel), intent(in) :: my_mpi

    double precision, dimension(:,:,:), intent(inout) :: VINS        ! .. input potential
    double precision, dimension(:,:), intent(inout) :: VISP
    double precision, intent(inout) :: ECORE(20,2)

    !----------
    integer, dimension(:), allocatable :: ranks
    integer :: my_world_rank
    integer :: my_atom_id
    integer :: owner
    integer :: num_SE_ranks
    integer :: ind

    num_SE_ranks = getNumSERAnks(my_mpi)
    my_world_rank = getMyWorldRank(my_mpi)
    my_atom_id = getMyAtomId(my_mpi)

    allocate(ranks(num_SE_ranks))

    owner = mapToWorldRankSE(my_mpi, my_atom_id, 1)

    do ind = 1, num_SE_ranks
      ranks(ind) = mapToWorldRankSE(my_mpi, my_atom_id, ind)
    end do

    call comm_bcastD(my_world_rank, VINS, size(VINS), ranks, owner)
    call comm_bcastD(my_world_rank, VISP, size(VISP), ranks, owner)
    call comm_bcastD(my_world_rank, ECORE, size(ECORE), ranks, owner)

    deallocate(ranks)

  end subroutine

! ---------------------- Jij --------------------------------------------------

  !----------------------------------------------------------------------------
  !> Communicates the off-diagonal Green's function elements to ranks with
  !> Spin_id = 1.
  subroutine jijSpinCommunication_com(my_mpi, GMATXIJ, nspind)
    use KKRnanoParallel_mod
    use comm_patternsZ_mod
    implicit none

    type (KKRnanoParallel), intent(in) :: my_mpi
    double complex, dimension(:,:,:,:), intent(inout) :: GMATXIJ
    integer, intent(in) :: nspind

    integer :: blocksize
    integer :: my_world_rank
    integer :: my_atom_id
    integer :: my_energy_id
    integer :: ranks(nspind)
    integer :: receiver
    integer :: ispin

    my_world_rank = getMyWorldRank(my_mpi)
    my_atom_id = getMyAtomId(my_mpi)
    my_energy_id = getMyEnergyId(my_mpi)

    blocksize = size(GMATXIJ, 1) * size(GMATXIJ, 2) * size(GMATXIJ, 3)

    ASSERT(nspind == size(GMATXIJ, 4))

    do ispin = 1, nspind
      ranks(ispin) = mapToWorldRank(my_mpi, my_atom_id, ispin, my_energy_id)
    end do

    receiver = ranks(1)   ! S=1 processes receive

    call comm_gatherZ(my_world_rank, GMATXIJ, blocksize, ranks, receiver)

  end subroutine

  !----------------------------------------------------------------------------
  !> Performs the energy integration over energy points that are locally known
  !> to rank.
  !> Only Spin_Id == 1 ranks work.
  !> Start with JXCIJINT set to zero then call for each energy point.
  !> Wrapper for XCCPLJIJ_START
  subroutine jijLocalEnergyIntegration(my_mpi, energy_weight, GMATXIJ, DTIXIJ, RXIJ, NXIJ, IXCP, RXCCLS, JXCIJINT)
    use KKRnanoParallel_mod
    implicit none

    type (KKRnanoParallel), intent(in) :: my_mpi
    double complex, intent(in) :: energy_weight
    double complex, dimension(:,:,:,:), intent(inout) :: GMATXIJ
    double complex, dimension(:,:) :: DTIXIJ
    double precision, dimension(:) :: RXIJ
    integer, intent(in) :: NXIJ
    double precision, dimension(:,:) :: RXCCLS
    integer, dimension(:) :: IXCP
    double complex, dimension(:) :: JXCIJINT


    logical :: ERESJIJ
    integer :: I1
    integer :: IER
    integer :: communicator
    integer :: naez, lmmaxd, nxijd, nspind

    ERESJIJ = .false.  ! not supported
    I1 = 0             ! just a dummy value
    IER = 2            ! dummy value != 1
    communicator = getMySEcommunicator(my_mpi)
    naez = getNumAtomRanks(my_mpi)

    lmmaxd = size(GMATXIJ, 1)
    ASSERT(lmmaxd == size(GMATXIJ, 2))
    nxijd = size(GMATXIJ, 3)
    nspind = size(GMATXIJ, 4)
    ASSERT(nspind == 2)

    ASSERT(size(DTIXIJ, 1) == lmmaxd)
    ASSERT(size(DTIXIJ, 2) == lmmaxd)
    ASSERT(size(RXIJ) == nxijd)
    ASSERT(size(RXCCLS) == nxijd*3)
    ASSERT(size(IXCP) == nxijd)
    ASSERT(size(JXCIJINT) == nxijd)
    ASSERT(NXIJ <= nxijd)

    if (getMySpinId(my_mpi) == 1) then
      call XCCPLJIJ_START(I1,IER,energy_weight,RXIJ,NXIJ,IXCP,RXCCLS,GMATXIJ,DTIXIJ, &
                          communicator, JXCIJINT,ERESJIJ, naez, lmmaxd, nxijd, nspind)
    end if

  end subroutine

  !----------------------------------------------------------------------------
  !> Communicate and sum results from all energy processes
  subroutine jijReduceIntResults_com(my_mpi, JXCIJINT)
    use KKRnanoParallel_mod
    use comm_patternsZ_mod
    implicit none

    type (KKRnanoParallel), intent(in) :: my_mpi
    double complex, dimension(:), intent(inout) :: JXCIJINT

    !-----------
    double complex, parameter :: CZERO = (0.0d0, 0.0d0)
    double complex, dimension(:), allocatable :: recv
    double complex, dimension(:), allocatable :: summed
    integer :: length

    integer :: my_world_rank
    integer :: my_atom_id
    integer :: my_spin_id
    integer :: my_energy_id
    integer :: num_eranks
    integer :: ind

    integer :: receiver
    integer :: sender

    my_world_rank = getMyWorldRank(my_mpi)
    my_atom_id = getMyAtomId(my_mpi)
    my_spin_id = getMySpinId(my_mpi)
    my_energy_id = getMyEnergyId(my_mpi)

    length = size(JXCIJINT)

    receiver = mapToWorldRank(my_mpi, my_atom_id, my_spin_id, 1)
    num_eranks = getNumEnergyRanks(my_mpi)

    if (my_spin_id == 1) then
      allocate(recv(length))
      allocate(summed(length))

      recv = CZERO

      summed = JXCIJINT

      do ind = 2, num_eranks
        sender = mapToWorldRank(my_mpi, my_atom_id, my_spin_id, ind)
        call send_arrayZ(my_world_rank, recv, length, sender, receiver)
        summed = summed + recv
      end do

      JXCIJINT = summed

      deallocate(recv)
      deallocate(summed)
    end if

  end subroutine

end module KKRnano_Comm_mod
