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
#include "../macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  public :: jijSpinCommunication_com, jijLocalEnergyIntegration, jijReduceIntResults_com
  public :: collectMSResults_com, redistributeInitialGuess
  public :: setKKRnanoNumThreads, printKKRnanoInfo, communicatePotential, communicateNoncoBfields

  interface redistributeInitialGuess
    module procedure redistributeInitialGuess_c, redistributeInitialGuess_z
  endinterface
  
  contains

  !----------------------------------------------------------------------------
  !> Print informative message about KKRnano parallelisation.
  subroutine printKKRnanoInfo(mp, nthrds)
    use KKRnanoParallel_mod, only: KKRnanoParallel

    type(KKRnanoParallel), intent(in) :: mp
    integer, intent(in) :: nthrds ! requested number of threads

    integer :: natoms, nspin, nenergy, nthrds_requested, num_threads ! actual number of threads
    !$ integer, external :: omp_get_num_threads
    !$ integer, external :: omp_get_max_threads

    if (mp%isMasterRank) then

      natoms  = mp%numAtomRanks
      nspin   = mp%numSpinRanks
      nenergy = mp%numEnergyRanks

      num_threads = 1
      nthrds_requested = num_threads
      if (nthrds > 0) nthrds_requested = nthrds
      
      write(*,'(79("="))')
      write(*,*) '  total no. of MPI ranks  = ', mp%numWorldRanks
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
!$  external :: OMP_SET_NUM_THREADS
    if (nthrds > 0) then
      !$ call OMP_SET_NUM_THREADS(NTHRDS)
    endif ! nthrds > 0
  endsubroutine ! set

  !----------------------------------------------------------------------------
  !> Collect the results from the multiple scattering part at
  !> the corresponding atom process of the master group.
  !> Master Group: (Spin, Energy)-id = 1
  subroutine collectMSResults_com(mp, GMATN_ALL, LLY_GRDT_ALL, EPROC)
    use KKRnanoParallel_mod, only: KKRnanoParallel, getResponsibleSpinId, mapToWorldRank, mapToWorldRankSE
    use comm_patternsZ_mod, only: comm_gatherZ

    type(KKRnanoParallel), intent(in) :: mp
    double complex, intent(inout) :: GMATN_ALL(:,:,:,:), LLY_GRDT_ALL(:,:)
    integer, intent(in) :: EPROC(:)

    integer :: memory_stat, iemxd, lmmaxd, ie, spin_id, ispin, nspind, receiver
    integer, allocatable :: owning_ranks(:)

    lmmaxd = size(GMATN_ALL, 1)
    iemxd  = size(GMATN_ALL, 3)
    nspind = size(GMATN_ALL, 4)

    ASSERT(lmmaxd == size(GMATN_ALL, 2))
    ASSERT(iemxd  == size(LLY_GRDT_ALL, 1))
    ASSERT(nspind == size(LLY_GRDT_ALL, 2))
    ASSERT(mp%numSpinRanks >= 1)
    ASSERT(mp%numSpinRanks <= 2)
    ASSERT(nspind >= 1)
    ASSERT(nspind <= 2)

    ! TODO: check allocate
    ALLOCATECHECK(owning_ranks(iemxd*nspind))

    do ispin = 1, nspind
      spin_id = getResponsibleSpinId(mp, ispin)
      do ie = 1, iemxd
        owning_ranks((ispin-1)*iemxd + ie) = mapToWorldRank(mp, mp%myAtomId, spin_id, EPROC(ie))
      enddo ! ie
    enddo ! ispin

    receiver = mapToWorldRankSE(mp, mp%myAtomId, 1)

    call comm_gatherZ(mp%myWorldRank, GMATN_ALL, lmmaxd*lmmaxd, owning_ranks, receiver)
    call comm_gatherZ(mp%myWorldRank, LLY_GRDT_ALL, 1, owning_ranks, receiver)

    ! TODO: check dimensions
    
    DEALLOCATECHECK(owning_ranks)

  endsubroutine ! collect

  !----------------------------------------------------------------------------
  subroutine redistributeInitialGuess_c(mp, prs, EPROC, EPROCO, KMESH, NofKs)
    use KKRnanoParallel_mod, only: KKRnanoParallel, mapToWorldRank
    use comm_patternsC_mod, only: comm_redistributeVC

    type(KKRnanoParallel), intent(in) :: mp
    complex, intent(inout) :: prs(:,:) ! single precision complex
    integer, intent(in) :: EPROC(:), EPROCO(:), KMESH(:), NofKs(:)

    integer, allocatable :: blocksizes(:), old_owners(:), new_owners(:)
    integer :: memory_stat
    integer :: single_block_size, num_energy_iemxd, ie

    num_energy_iemxd = size(EPROC)
    single_block_size = size(prs, 1)

    ASSERT(num_energy_iemxd == size(EPROCO))
    ASSERT(num_energy_iemxd == size(KMESH))

    ALLOCATECHECK(blocksizes(num_energy_iemxd))
    ALLOCATECHECK(old_owners(num_energy_iemxd))
    ALLOCATECHECK(new_owners(num_energy_iemxd))

    do ie = 1, num_energy_iemxd
      old_owners(ie) = mapToWorldRank(mp, mp%myAtomId, mp%mySpinId, EPROCO(ie))
      new_owners(ie) = mapToWorldRank(mp, mp%myAtomId, mp%mySpinId, EPROC(ie))
      ASSERT(KMESH(ie) < size(NofKs))
      blocksizes(ie) = single_block_size*NofKs(KMESH(IE))
    enddo ! ie

    call comm_redistributeVC(mp%myWorldRank, prs, blocksizes, old_owners, new_owners)

    DEALLOCATECHECK(blocksizes)
    DEALLOCATECHECK(old_owners)
    DEALLOCATECHECK(new_owners)

  endsubroutine ! redistribute

  !----------------------------------------------------------------------------
  subroutine redistributeInitialGuess_z(mp, prs, EPROC, EPROCO, KMESH, NofKs)
    use KKRnanoParallel_mod, only: KKRnanoParallel, mapToWorldRank
    use comm_patternsZ_mod, only: comm_redistributeVZ

    type(KKRnanoParallel), intent(in) :: mp
    double complex, intent(inout) :: prs(:,:) ! double precision complex
    integer, intent(in) :: EPROC(:), EPROCO(:), KMESH(:), NofKs(:)

    integer, allocatable :: blocksizes(:), old_owners(:), new_owners(:)
    integer :: memory_stat
    integer :: single_block_size, num_energy_iemxd, ie

    num_energy_iemxd = size(EPROC)
    single_block_size = size(prs, 1)

    ASSERT(num_energy_iemxd == size(EPROCO))
    ASSERT(num_energy_iemxd == size(KMESH))

    ALLOCATECHECK(blocksizes(num_energy_iemxd))
    ALLOCATECHECK(old_owners(num_energy_iemxd))
    ALLOCATECHECK(new_owners(num_energy_iemxd))

    do ie = 1, num_energy_iemxd
      old_owners(ie) = mapToWorldRank(mp, mp%myAtomId, mp%mySpinId, EPROCO(ie))
      new_owners(ie) = mapToWorldRank(mp, mp%myAtomId, mp%mySpinId, EPROC(ie))
      ASSERT(KMESH(ie) < size(NofKs))
      blocksizes(ie) = single_block_size*NofKs(KMESH(IE))
    enddo ! ie

    call comm_redistributeVZ(mp%myWorldRank, prs, blocksizes, old_owners, new_owners)

    DEALLOCATECHECK(blocksizes)
    DEALLOCATECHECK(old_owners)
    DEALLOCATECHECK(new_owners)

  endsubroutine ! redistribute
  
  
  !----------------------------------------------------------------------------
  !> Communicate Potential from Master-Group to all other (Spin,Energy)-Groups.
  subroutine communicatePotential(mp, VISP, VINS, ECORE)
    use KKRnanoParallel_mod, only: KKRnanoParallel, getNumSERanks, mapToWorldRankSE
    use comm_patternsD_mod, only: comm_bcastD

    type(KKRnanoParallel), intent(in) :: mp
    double precision, intent(inout) :: VINS(:,:,:) ! .. input potential
    double precision, intent(inout) :: VISP(:,:)
    double precision, intent(inout) :: ECORE(20,2)

    integer, allocatable :: ranks(:)
    integer :: memory_stat, owner, numSERanks, ind

    numSERanks = getNumSERanks(mp)

    ALLOCATECHECK(ranks(numSERanks))

    owner = mapToWorldRankSE(mp, mp%myAtomId, 1)

    do ind = 1, numSERanks
      ranks(ind) = mapToWorldRankSE(mp, mp%myAtomId, ind)
    enddo ! ind

    call comm_bcastD(mp%myWorldRank, VINS, size(VINS), ranks, owner)
    call comm_bcastD(mp%myWorldRank, VISP, size(VISP), ranks, owner)
    call comm_bcastD(mp%myWorldRank, ECORE, size(ECORE), ranks, owner)

    DEALLOCATECHECK(ranks)

  endsubroutine ! communicate

  subroutine communicateNoncoBfields(mp, bfields)
    use KKRnanoParallel_mod, only: KKRnanoParallel, getNumSERanks, mapToWorldRankSE
    use comm_patternsD_mod, only: comm_bcastD
    use mod_bfield, only: bfield_data

    type(KKRnanoParallel), intent(in) :: mp
    type(bfield_data), intent(inout)  :: bfields

    integer, allocatable :: ranks(:)
    integer :: memory_stat, owner, numSERanks, ind

    ! Get number of MPI ranks for this atom
    numSERanks = getNumSERanks(mp)

    ! Get owner rank for this atom
    owner = mapToWorldRankSE(mp, mp%myAtomId, 1)

    ! Get all other ranks for this atom
    ALLOCATECHECK(ranks(numSERanks))
    do ind = 1, numSERanks
      ranks(ind) = mapToWorldRankSE(mp, mp%myAtomId, ind)
    enddo

    ! Brodcast data from owner of this atom to all other ranks dealing with this atom
    call comm_bcastD(mp%myWorldRank, bfields%bfield_constr, 3, ranks, owner)

    DEALLOCATECHECK(ranks)

  end subroutine

! ---------------------- Jij --------------------------------------------------

  !----------------------------------------------------------------------------
  !> Communicates the off-diagonal Green's function elements and t-matrices
  !> to ranks with Spin_id = 1.
  subroutine jijSpinCommunication_com(mp, GMATXIJ, DTIXIJ)
    use KKRnanoParallel_mod, only: KKRnanoParallel, mapToWorldRank
    use comm_patternsZ_mod, only: comm_gatherZ

    type(KKRnanoParallel), intent(in) :: mp
    double complex, intent(inout) :: GMATXIJ(:,:,:,:)
    double complex, intent(inout) :: DTIXIJ(:,:,:)

    integer :: blocksize, ranks(2), receiver, spin_id

    ASSERT(size(GMATXIJ, 4) == 2)
    ASSERT(mp%numSpinRanks >= 1)
    ASSERT(mp%numSpinRanks <= 2)

    do spin_id = 1, mp%numSpinRanks
      ranks(spin_id) = mapToWorldRank(mp, mp%myAtomId, spin_id, mp%myEnergyId)
    enddo ! spin_id

    receiver = ranks(1) ! S=1 processes receive

    !--------------------------- GMATXIJ --------------------------------------
    if (mp%numSpinRanks == 1) then
      ! "communicate" both spin-channels
      ! Note: no communication necessary in this case - comm_gatherZ detects this
      blocksize = size(GMATXIJ, 1) * size(GMATXIJ, 2) * size(GMATXIJ, 3) * 2
    else
      ! communicate only one spin channel
      blocksize = size(GMATXIJ, 1) * size(GMATXIJ, 2) * size(GMATXIJ, 3)
    endif

    call comm_gatherZ(mp%myWorldRank, GMATXIJ, blocksize, ranks(1:mp%numSpinRanks), receiver)

    !--------------------------- DTIXIJ ---------------------------------------

    if (mp%numSpinRanks == 1) then
      ! "communicate" both spin-channels
      ! Note: no communication necessary in this case - comm_gatherZ detects this
      blocksize = size(DTIXIJ, 1) * size(DTIXIJ, 2) * 2
    else
      ! communicate only one spin channel
      blocksize = size(DTIXIJ, 1) * size(DTIXIJ, 2)
    endif

    call comm_gatherZ(mp%myWorldRank, DTIXIJ, blocksize, ranks(1:mp%numSpinRanks), receiver)

    if (mp%mySpinId /= 1) then
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
  subroutine jijLocalEnergyIntegration(mp, energy_weight, GMATXIJ, DTIXIJ, RXIJ, NXIJ, IXCP, RXCCLS, JXCIJINT)
    use KKRnanoParallel_mod, only: KKRnanoParallel
    use jij_calc_mod, only: XCCPLJIJ_START

    type(KKRnanoParallel), intent(in) :: mp
    double complex, intent(in) :: energy_weight
    double complex, intent(inout) :: GMATXIJ(:,:,:,:)
    double complex :: DTIXIJ(:,:)
    double precision :: RXIJ(:)
    integer, intent(in) :: NXIJ
    double precision :: RXCCLS(:,:)
    integer :: IXCP(:)
    double complex :: JXCIJINT(:)

    logical :: ERESJIJ
    integer :: I1, IER, lmmaxd, nxijd, nspind

    ERESJIJ = .false.  ! not supported
    I1 = 0             ! just a dummy value
    IER = 2            ! dummy value != 1

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
    if (mp%mySpinId == 1) then
      call XCCPLJIJ_START(IER, energy_weight, NXIJ, IXCP, GMATXIJ, DTIXIJ, &
                          mp%mySEComm, JXCIJINT, ERESJIJ, mp%numAtomRanks, lmmaxd, nxijd, nspind)
    else
      JXCIJINT = dcmplx(1d9, 1d9) ! invalidate results for other ranks to be able to detect errors
    endif

  endsubroutine ! jij

  !----------------------------------------------------------------------------
  !> Communicate and sum results from all energy processes to ranks with
  !> (S-Id, E-Id) = (1,1)
  !> Only those ranks hold the correct result!!!
  subroutine jijReduceIntResults_com(mp, JXCIJINT)
    use KKRnanoParallel_mod, only: KKRnanoParallel, mapToWorldRank
    use comm_patternsZ_mod, only: send_arrayZ

    type(KKRnanoParallel), intent(in) :: mp
    double complex, intent(inout) :: JXCIJINT(:)

    double complex, parameter :: CZERO = (0.0d0, 0.0d0)
    double complex, allocatable :: sendrecv(:), summed(:)
    integer :: memory_stat
    integer :: ind, receiver, sender, length

    length = size(JXCIJINT)

    receiver = mapToWorldRank(mp, mp%myAtomId, mp%mySpinId, 1)

    ! Only ranks with Spin-Id=1 work!!!

    if (mp%mySpinId == 1) then
      ALLOCATECHECK(sendrecv(length))
      ALLOCATECHECK(summed(length))

      ! JXCIJINT is the array to send
      sendrecv = JXCIJINT

      ! start with 2nd energy group, because 1st already knows its contribution
      ! therefore initialise 'summed' with JXCIJINT (contribution of 1st group)
      summed = JXCIJINT

      do ind = 2, mp%numEnergyRanks
        sender = mapToWorldRank(mp, mp%myAtomId, mp%mySpinId, ind)
        call send_arrayZ(mp%myWorldRank, sendrecv, length, sender, receiver)

        if (mp%myWorldRank == receiver) summed = summed + sendrecv  ! only sum when (S-id, E-id) = (1,1) !
      enddo ! ind

      JXCIJINT = summed

      if (mp%myEnergyId /= 1) JXCIJINT = dcmplx(1d9, 1d9) ! invalidate results for other ranks to be able to detect errors

      DEALLOCATECHECK(sendrecv)
      DEALLOCATECHECK(summed)
    endif  ! mp%mySpinId == 1

  endsubroutine ! jij

endmodule KKRnano_Comm_mod
