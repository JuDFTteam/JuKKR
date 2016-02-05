!------------------------------------------------------------------------------
!> This module contains the routines that communicate intermediate results
!> between different processors.
!> Author: Elias Rabel, 2012
!------------------------------------------------------------------------------

#include "../DebugHelpers/test_macros.h"

! Some macros for checked allocation/deallocation
! they need an integer variable named memory_stat declared in each routine
! they are used.

#define CHECKALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define CHECKDEALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Deallocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)
#define DEALLOCATECHECK(X) deallocate(X, stat=memory_stat); CHECKDEALLOC(memory_stat)


module KKRnano_Comm_mod
  implicit none
  private
  public :: jijSpinCommunication_com, jijLocalEnergyIntegration, jijReduceIntResults_com
  public :: collectMSResults_com, redistributeInitialGuess_com
  public :: setKKRnanoNumThreads, printKKRnanoInfo, communicatePotential
  
  contains

  !----------------------------------------------------------------------------
  !> Print informative message about KKRnano parallelisation.
  subroutine printKKRnanoInfo(my_mpi, nthrds)
    use KKRnanoParallel_mod, only: KKRnanoParallel, getNumWorldRanks, getNumAtomRanks, getNumSpinRanks, getNumEnergyRanks, isMasterRank

    type(KKRnanoParallel), intent(in) :: my_mpi
    integer, intent(in) :: nthrds ! requested number of threads

    integer :: natoms, nspin, nenergy, nthrds_requested, num_threads ! actual number of threads
    !$ integer, external :: omp_get_num_threads
    !$ integer, external :: omp_get_max_threads

    if (isMasterRank(my_mpi)) then

      natoms  = getNumAtomRanks(my_mpi)
      nspin   = getNumSpinRanks(my_mpi)
      nenergy = getNumEnergyRanks(my_mpi)

      num_threads = 1
      nthrds_requested = num_threads
      if (nthrds > 0) nthrds_requested = nthrds
      
      write(*,'(79("="))')
      write(*,*) '  total no. of MPI ranks  = ', getNumWorldRanks(my_mpi)
      !$omp parallel
      !$omp single
      !$  num_threads = omp_get_num_threads()
      !$  write(*,*) '  OMP max. threads        = ',omp_get_max_threads()
      !$  write(*,*) '  OMP num. threads        = ',num_threads
      !$omp endsingle
      !$omp endparallel
      write(*,*) '  groups of processes created'
      write(*,*) '  NMPI                    = ',natoms
      write(*,*) '  SMPI                    = ',nspin
      write(*,*) '  EMPI                    = ',nenergy
      write(*,*) '  NTHRDS                  = ',nthrds_requested
      write(*,*) '  MPI-processes (active)  = ',natoms*nspin*nenergy
      write(*,*) '  total no. of tasks      = ',natoms*nspin*nenergy*num_threads
      write(*,'(79("="))')

    endif ! is master

  endsubroutine ! print

  !----------------------------------------------------------------------------
  !> Set the number of OpenMP threads to nthrds.
  subroutine setKKRnanoNumThreads(nthrds)
    integer, intent(in) :: nthrds
    if (nthrds > 0) then
      !$ call OMP_SET_NUM_THREADS(NTHRDS)
    endif ! nthrds > 0
  endsubroutine ! set

  !----------------------------------------------------------------------------
  !> Collect the results from the multiple scattering part at
  !> the corresponding atom process of the master group.
  !> Master Group: (Spin, Energy)-id = 1
  subroutine collectMSResults_com(my_mpi, GMATN_ALL, LLY_GRDT_ALL, EPROC)
    use KKRnanoParallel_mod, only: KKRnanoParallel, getNumSpinRanks, getResponsibleSpinId, mapToWorldRank, getMyAtomId, mapToWorldRankSE, getMyWorldRank
    use comm_patternsZ_mod, only: comm_gatherZ

    type(KKRnanoParallel), intent(in) :: my_mpi
    double complex, intent(inout) :: GMATN_ALL(:,:,:,:), LLY_GRDT_ALL(:,:)
    integer, intent(in) :: EPROC(:)

    integer :: memory_stat, iemxd, num_spin_ranks, lmmaxd, ie, spin_id, ispin, nspind, my_atom_id, my_world_rank, receiver
    integer, allocatable :: owning_ranks(:)

    lmmaxd = size(GMATN_ALL,1)
    iemxd = size(GMATN_ALL,3)
    nspind = size(GMATN_ALL, 4)

    num_spin_ranks = getNumSpinRanks(my_mpi)

    ASSERT(lmmaxd == size(GMATN_ALL, 2))
    ASSERT(iemxd == size(LLY_GRDT_ALL, 1))
    ASSERT(nspind == size(LLY_GRDT_ALL, 2))
    ASSERT(num_spin_ranks >= 1)
    ASSERT(num_spin_ranks <= 2)
    ASSERT(nspind >= 1)
    ASSERT(nspind <= 2)

    ! TODO: check allocate
    ALLOCATECHECK(owning_ranks(iemxd * nspind))

    my_atom_id = getMyAtomId(my_mpi)

    do ispin = 1, nspind
      spin_id = getResponsibleSpinId(my_mpi, ispin)
      do ie = 1, iemxd
        owning_ranks((ispin-1)*iemxd + ie) = mapToWorldRank(my_mpi, my_atom_id, spin_id, EPROC(ie))
      enddo ! ie
    enddo ! ispin

    receiver = mapToWorldRankSE(my_mpi, my_atom_id, 1)

    my_world_rank = getMyWorldRank(my_mpi)

    call comm_gatherZ(my_world_rank, GMATN_ALL, lmmaxd*lmmaxd, owning_ranks, receiver)
    call comm_gatherZ(my_world_rank, LLY_GRDT_ALL, 1, owning_ranks, receiver)

    ! TODO: check dimensions

    DEALLOCATECHECK(owning_ranks)

  endsubroutine ! collect

  !----------------------------------------------------------------------------
  subroutine redistributeInitialGuess_com(my_mpi, PRSC, EPROC, EPROCO, KMESH, NofKs)
    use KKRnanoParallel_mod, only: KKRnanoParallel, getMyAtomId, getMyWorldRank, getMySpinId, mapToWorldRank
    use comm_patternsC_mod, only: comm_redistributeVC

    type(KKRnanoParallel), intent(in) :: my_mpi
    complex, intent(inout) :: PRSC(:,:)
    integer, intent(in) :: EPROC(:), EPROCO(:), KMESH(:), NofKs(:)

    integer, allocatable :: blocksizes(:), old_owners(:), new_owners(:)

    integer :: memory_stat, my_world_rank, my_atom_id, my_spin_id
    integer :: single_block_size, num_energy_iemxd, ie

    my_atom_id = getMyAtomId(my_mpi)
    my_world_rank = getMyWorldRank(my_mpi)
    my_spin_id  = getMySpinId(my_mpi)

    num_energy_iemxd = size(EPROC)
    single_block_size = size(PRSC, 1)

    ASSERT(num_energy_iemxd == size(EPROCO))
    ASSERT(num_energy_iemxd == size(KMESH))

    ALLOCATECHECK(blocksizes(num_energy_iemxd))
    ALLOCATECHECK(old_owners(num_energy_iemxd))
    ALLOCATECHECK(new_owners(num_energy_iemxd))

    do ie = 1, num_energy_iemxd
      old_owners(ie) = mapToWorldRank(my_mpi, my_atom_id, my_spin_id, EPROCO(ie))
      new_owners(ie) = mapToWorldRank(my_mpi, my_atom_id, my_spin_id, EPROC(ie))
      ASSERT(KMESH(ie) < size(NofKs))
      blocksizes(ie) = single_block_size*NofKs(KMESH(IE))
    enddo ! ie

    call comm_redistributeVC(my_world_rank, PRSC, blocksizes, old_owners, new_owners)

    DEALLOCATECHECK(blocksizes)
    DEALLOCATECHECK(old_owners)
    DEALLOCATECHECK(new_owners)

  endsubroutine ! redistribute

  !----------------------------------------------------------------------------
  !> Communicate Potential from Master-Group to all other (Spin,Energy)-Groups.
  subroutine communicatePotential(my_mpi, VISP, VINS, ECORE)
    use KKRnanoParallel_mod, only: KKRnanoParallel, getNumSERanks, getMyWorldRank, getMyAtomId, mapToWorldRankSE
    use comm_patternsD_mod, only: comm_bcastD

    type(KKRnanoParallel), intent(in) :: my_mpi
    double precision, intent(inout) :: VINS(:,:,:) ! .. input potential
    double precision, intent(inout) :: VISP(:,:)
    double precision, intent(inout) :: ECORE(20,2)

    integer, allocatable :: ranks(:)
    integer :: memory_stat, my_world_rank, my_atom_id, owner, num_SE_ranks, ind

    num_SE_ranks = getNumSERanks(my_mpi)
    my_world_rank = getMyWorldRank(my_mpi)
    my_atom_id = getMyAtomId(my_mpi)

    ALLOCATECHECK(ranks(num_SE_ranks))

    owner = mapToWorldRankSE(my_mpi, my_atom_id, 1)

    do ind = 1, num_SE_ranks
      ranks(ind) = mapToWorldRankSE(my_mpi, my_atom_id, ind)
    enddo ! ind

    call comm_bcastD(my_world_rank, VINS, size(VINS), ranks, owner)
    call comm_bcastD(my_world_rank, VISP, size(VISP), ranks, owner)
    call comm_bcastD(my_world_rank, ECORE, size(ECORE), ranks, owner)

    DEALLOCATECHECK(ranks)

  endsubroutine ! communicate

! ---------------------- Jij --------------------------------------------------

  !----------------------------------------------------------------------------
  !> Communicates the off-diagonal Green's function elements and t-matrices
  !> to ranks with Spin_id = 1.
  subroutine jijSpinCommunication_com(my_mpi, GMATXIJ, DTIXIJ)
    use KKRnanoParallel_mod, only: KKRnanoParallel, getMyWorldRank, getMyAtomId, getMyEnergyId, getNumSpinRanks, mapToWorldRank, getMySpinId
    use comm_patternsZ_mod, only: comm_gatherZ

    type(KKRnanoParallel), intent(in) :: my_mpi
    double complex, intent(inout) :: GMATXIJ(:,:,:,:)
    double complex, intent(inout) :: DTIXIJ(:,:,:)

    integer :: blocksize, my_world_rank, my_atom_id, my_energy_id
    integer :: ranks(2), receiver, spin_id, num_spin_ranks

    my_world_rank = getMyWorldRank(my_mpi)
    my_atom_id = getMyAtomId(my_mpi)
    my_energy_id = getMyEnergyId(my_mpi)
    num_spin_ranks = getNumSpinRanks(my_mpi)

    ASSERT(size(GMATXIJ, 4) == 2)
    ASSERT(num_spin_ranks >= 1)
    ASSERT(num_spin_ranks <= 2)

    do spin_id = 1, num_spin_ranks
      ranks(spin_id) = mapToWorldRank(my_mpi, my_atom_id, spin_id, my_energy_id)
    enddo ! spin_id

    receiver = ranks(1)   ! S=1 processes receive

    !--------------------------- GMATXIJ --------------------------------------
    if (num_spin_ranks == 1) then
      ! "communicate" both spin-channels
      ! Note: no communication necessary in this case - comm_gatherZ detects this
      blocksize = size(GMATXIJ, 1) * size(GMATXIJ, 2) * size(GMATXIJ, 3) * 2
    else
      ! communicate only one spin channel
      blocksize = size(GMATXIJ, 1) * size(GMATXIJ, 2) * size(GMATXIJ, 3)
    endif

    call comm_gatherZ(my_world_rank, GMATXIJ, blocksize, ranks(1:num_spin_ranks), receiver)

    !--------------------------- DTIXIJ ---------------------------------------

    if (num_spin_ranks == 1) then
      ! "communicate" both spin-channels
      ! Note: no communication necessary in this case - comm_gatherZ detects this
      blocksize = size(DTIXIJ, 1) * size(DTIXIJ, 2) * 2
    else
      ! communicate only one spin channel
      blocksize = size(DTIXIJ, 1) * size(DTIXIJ, 2)
    endif

    call comm_gatherZ(my_world_rank, DTIXIJ, blocksize, ranks(1:num_spin_ranks), receiver)

    if (getMySpinId(my_mpi) /= 1) then
      ! invalidate results for other ranks to be able to detect errors
      GMATXIJ = dcmplx(1d9, 1d9)
      DTIXIJ  = dcmplx(1d9, 1d9)
    endif

  endsubroutine ! jij

  !----------------------------------------------------------------------------
  !> Performs the energy integration over energy points that are locally known
  !> to rank.
  !> Only Spin_Id == 1 ranks work.
  !> Start with JXCIJINT set to zero then call for each energy point.
  !> Wrapper for XCCPLJIJ_START
  subroutine jijLocalEnergyIntegration(my_mpi, energy_weight, GMATXIJ, DTIXIJ, RXIJ, NXIJ, IXCP, RXCCLS, JXCIJINT)
    use KKRnanoParallel_mod, only: KKRnanoParallel, getMySEcommunicator, getNumAtomRanks, getMySpinId
    use jij_calc_mod, only: XCCPLJIJ_START

    type(KKRnanoParallel), intent(in) :: my_mpi
    double complex, intent(in) :: energy_weight
    double complex, intent(inout) :: GMATXIJ(:,:,:,:)
    double complex :: DTIXIJ(:,:)
    double precision :: RXIJ(:)
    integer, intent(in) :: NXIJ
    double precision :: RXCCLS(:,:)
    integer :: IXCP(:)
    double complex :: JXCIJINT(:)

    logical :: ERESJIJ
    integer :: I1, IER, naez, lmmaxd, nxijd, nspind
    integer :: communicator

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

    ! Only ranks with Spin-Id=1 work!!!
    if (getMySpinId(my_mpi) == 1) then
      call XCCPLJIJ_START(I1, IER, energy_weight, RXIJ, NXIJ, IXCP, RXCCLS, GMATXIJ, DTIXIJ, &
                          communicator, JXCIJINT, ERESJIJ, naez, lmmaxd, nxijd, nspind)
    else
      JXCIJINT = dcmplx(1d9, 1d9) ! invalidate results for other ranks to be able to detect errors
    endif

  endsubroutine ! jij

  !----------------------------------------------------------------------------
  !> Communicate and sum results from all energy processes to ranks with
  !> (S-Id, E-Id) = (1,1)
  !> Only those ranks hold the correct result!!!
  subroutine jijReduceIntResults_com(my_mpi, JXCIJINT)
    use KKRnanoParallel_mod, only: KKRnanoParallel, getMyWorldRank, getMyAtomId, getMySpinId, getMyEnergyId, mapToWorldRank, getNumEnergyRanks
    use comm_patternsZ_mod, only: send_arrayZ

    type(KKRnanoParallel), intent(in) :: my_mpi
    double complex, intent(inout) :: JXCIJINT(:)

    double complex, parameter :: CZERO = (0.0d0, 0.0d0)
    double complex, allocatable :: sendrecv(:), summed(:)
    integer :: memory_stat, length
    integer :: my_atom_id, my_spin_id, my_energy_id
    integer :: num_eranks, my_world_rank
    integer :: ind, receiver, sender

    my_world_rank = getMyWorldRank(my_mpi)
    my_atom_id = getMyAtomId(my_mpi)
    my_spin_id = getMySpinId(my_mpi)
    my_energy_id = getMyEnergyId(my_mpi)

    length = size(JXCIJINT)

    receiver = mapToWorldRank(my_mpi, my_atom_id, my_spin_id, 1)
    num_eranks = getNumEnergyRanks(my_mpi)

    ! Only ranks with Spin-Id=1 work!!!

    if (my_spin_id == 1) then
      ALLOCATECHECK(sendrecv(length))
      ALLOCATECHECK(summed(length))

      ! JXCIJINT is the array to send
      sendrecv = JXCIJINT

      ! start with 2nd energy group, because 1st already knows its contribution
      ! therefore initialise 'summed' with JXCIJINT (contribution of 1st group)
      summed = JXCIJINT

      do ind = 2, num_eranks
        sender = mapToWorldRank(my_mpi, my_atom_id, my_spin_id, ind)
        call send_arrayZ(my_world_rank, sendrecv, length, sender, receiver)

        if (my_world_rank == receiver) then  ! only sum when (S-id, E-id) = (1,1) !
          summed = summed + sendrecv
        endif

      enddo ! ind

      JXCIJINT = summed

      if (my_energy_id /= 1) &
        JXCIJINT = dcmplx(1d9, 1d9) ! invalidate results for other ranks to be able to detect errors

      DEALLOCATECHECK(sendrecv)
      DEALLOCATECHECK(summed)
    endif  ! my_spin_id == 1

  endsubroutine ! jij

endmodule KKRnano_Comm_mod
