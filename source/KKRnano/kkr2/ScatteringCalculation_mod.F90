#include "DebugHelpers/logging_macros.h"
#include "DebugHelpers/test_array_log.h"
#include "DebugHelpers/test_macros.h"

!> @author Modularisation: Elias Rabel
module ScatteringCalculation_mod
  use Logging_mod, only:    !import no name here, just mention it for the module dependency 
  use arraytest2_mod, only: !import no name here, just mention it for the module dependency 
implicit none
  private
  public  :: energyLoop

  contains

  !> Input:  *) ebalance_handler (properly initialised)
  !>         *) jij_data, ldau_data (properly initialised !!!)
  !>
  !> Output: ebalance_handler (changed, updated timings and process distribution)
  !>         KKRResults
  !>         jij_data results of jij-calculation
  !>         ldau_data LDA+U results
  !>
  !>         FILES WRITTEN:
  !>         *) Logfiles (if requested)
  !>         *) JIJ-Files (if requested)
  !>         *) matrix dump (if requested)
  subroutine energyLoop(iter, calc, emesh, params, dims, ebalance_handler, my_mpi, arrays)

    USE_LOGGING_MOD
    USE_ARRAYLOG_MOD

    use TimerMpi_mod, only: TimerMpi, getElapsedTime, resetTimer, stopTimer, resumeTimer, outtime
    use EBalanceHandler_mod, only: EBalanceHandler
    use TEST_lcutoff_mod, only: DEBUG_dump_matrix

    use DimParams_mod, only: DimParams
    use InputParams_mod, only: InputParams
    use Main2Arrays_mod, only: Main2Arrays

    use CalculationData_mod, only: CalculationData, getKKR, getAtomData, getLDAUData, getNumLocalAtoms, getAtomIndexOfLocal

    use KKRresults_mod, only: KKRresults
    use BasisAtom_mod, only: BasisAtom
    use JijData_mod, only: JijData
    use LDAUData_mod, only: LDAUData
    use EnergyMesh_mod, only: EnergyMesh
    
    use KKRnanoParallel_mod, only: KKRnanoParallel, isMasterRank, getMyEnergyId, getMySEcommunicator, isWorkingSpinRank, getMyAtomRank, isInMasterGroup
    use KKRnano_Comm_mod, only: jijSpinCommunication_com, jijLocalEnergyIntegration, jijReduceIntResults_com, collectMSResults_com, redistributeInitialGuess_com

    use EBalanceHandler_mod, only: EBalanceHandler, startEBalanceTiming, stopEBalanceTiming, updateEBalance_com

    use kloopz1_mod, only: kloopz1_new
    use InitialGuess_mod, only: InitialGuess, iguess_set_energy_ind, iguess_set_spin_ind

    use wrappers_mod, only: calctmat_wrapper, calcdtmat_wrapper
    use jij_calc_mod, only: clsjij, writejijs, jij_data => global_jij_data

    use TFQMRSolver_mod, only: TFQMRSolver
    use BCPOperator_mod, only: BCPOperator
    use KKROperator_mod, only: KKROperator

    use SingleSiteRef_mod, only: TREF, GREF
    
    integer, intent(in) :: iter
    type(CalculationData), intent(inout) :: calc
    type(KKRnanoParallel), intent(in)    :: my_mpi
    type(EBalanceHandler), intent(inout) :: ebalance_handler
    type(EnergyMesh), intent(in)         :: emesh
    type(Main2Arrays), intent(in)        :: arrays
    type(DimParams), intent(in)          :: dims
    type(InputParams), intent(in)        :: params

    ! locals

    type(BasisAtom), pointer             :: atomdata  ! referenced data does not change
    type(KKRresults), pointer            :: kkr       ! changes
    type(LDAUData), pointer              :: ldau_data ! changes

    type(TFQMRSolver), target :: solv
    type(KKROperator), target :: kkr_op
    type(BCPOperator), target :: precond

    double complex, parameter :: ZERO = (0.d0, 0.d0)
    type(TimerMpi) :: mult_scattering_timer, single_site_timer
    double complex :: JSCAL ! scaling factor for Jij calculation
    integer, allocatable :: atom_indices(:)
    integer :: ie, ispin, prspin, nmesh
    integer :: I1, ilocal, num_local_atoms
    integer :: lmmaxd
    logical :: xccpl

    double complex, allocatable :: Tref_local(:,:,:)  !< local tref-matrices
    double complex, allocatable :: dTref_local(:,:,:) !< local deriv. tref-matrices
    double complex, allocatable :: tmatll(:,:,:) !< all t-matrices
    double complex, allocatable :: GmatN_buffer(:,:,:) !< GmatN for all local atoms
    double complex, allocatable :: GrefN_buffer(:,:,:,:) !< GrefN for all local atoms

    lmmaxd = (dims%lmaxd+1)**2

    atomdata  => getAtomData(calc, 1)
    I1 = atomdata%atom_index
    kkr       => null()
    ldau_data => getLDAUData(calc, 1)
#define clusters calc%clusters
#define lattice_vectors calc%lattice_vectors
#define trunc_zone calc%trunc_zone

    jij_data => calc%jij_data_a(1) ! global, jij works only with max. 1 local atom

    num_local_atoms = getNumLocalAtoms(calc)

    
    allocate(tmatll(lmmaxd,lmmaxd,trunc_zone%naez_trc)) ! allocate buffer for t-matrices
    allocate(Tref_local(lmmaxd,lmmaxd,num_local_atoms)) ! allocate buffers for reference t-matrices
    allocate(dTref_local(lmmaxd,lmmaxd,num_local_atoms))
    allocate(GmatN_buffer(lmmaxd,lmmaxd,num_local_atoms))
    allocate(GrefN_buffer(lmmaxd,lmmaxd,clusters%naclsd,num_local_atoms))

    allocate(atom_indices(num_local_atoms))

    if (params%jij  .and. num_local_atoms > 1) stop "Jij and num_local_atoms > 1 not supported."
    if (params%ldau .and. num_local_atoms > 1) stop "LDA+U and num_local_atoms > 1 not supported."

    xccpl = .false.

    call resetTimer(mult_scattering_timer)
    call stopTimer(mult_scattering_timer)
    
    call resetTimer(single_site_timer)

    prspin = 1

    ! calculate exchange couplings only at last self-consistency step and when Jij=true
    if (ITER == params%SCFSTEPS .and. params%JIJ) XCCPL = .true. ! activate in last SCF iteration

    if (XCCPL) then
      jij_data%do_jij_calculation = .true. ! Trigger jij-calculation

      call CLSJIJ(I1, dims%NAEZ, lattice_vectors%RR, lattice_vectors%nrd, arrays%RBASIS, &
                  jij_data%RCUTJIJ, arrays%NSYMAT, arrays%ISYMINDEX, &
                  jij_data%IXCP, jij_data%NXCP, jij_data%NXIJ, jij_data%RXIJ, &
                  jij_data%RXCCLS, jij_data%ZKRXIJ, &
                  lattice_vectors%nrd, jij_data%nxijd)

      jij_data%JXCIJINT = ZERO
      jij_data%GMATXIJ = ZERO
    endif

    ! get the indices of atoms that shall be treated at once by the process
    ! = truncation zone indices of local atoms
    do ilocal = 1, num_local_atoms
      atom_indices(ilocal) = trunc_zone%index_map(getAtomIndexOfLocal(calc, ilocal))
      CHECKASSERT(atom_indices(ilocal) > 0)
    enddo ! ilocal

    ! setup the solver + bcp preconditioner, allocates a lot of memory
    ! it is good to do these allocations outside of energy loop
    call setup_solver(solv, kkr_op, precond, dims, clusters, lmmaxd, params%qmrbound, atom_indices)

    
  ! IE ====================================================================
  !     BEGIN do loop over energies (EMPID-parallel)
  ! IE ====================================================================
    do IE = 1, emesh%ielast
  ! IE ====================================================================
      if (getMyEnergyId(my_mpi) == ebalance_handler%EPROC(IE)) then
  ! IE ====================================================================
        call startEBalanceTiming(ebalance_handler, IE)

        WRITELOG(2, *) "Working on energy point ", IE

         Tref_local = ZERO
        dTref_local = ZERO
  !------------------------------------------------------------------------------
        !$omp parallel do private(ilocal, atomdata)
        do ilocal = 1, num_local_atoms
          atomdata => getAtomData(calc, ilocal)

          call TREF(emesh%EZ(IE), params%vref, dims%LMAXD, atomdata%rMTref, &
                    Tref_local(:,:,ilocal), dTref_local(:,:,ilocal), derive=(dims%LLY > 0))

        enddo  ! ilocal
        !$omp endparallel do
  !------------------------------------------------------------------------------

        ! Note: ref. system has to be recalculated at each iteration since energy mesh changes
        ! Note: TrefLL is diagonal - however full matrix is stored
        ! Note: Gref is calculated in real space - usually only a few shells

        ! Exchange the reference t-matrices within reference clusters
        ! ToDo: discuss if we can compute them once we know rMTref of all atoms in the reference cluster
        do ilocal = 1, num_local_atoms
          kkr => getKKR(calc, ilocal)

          call gatherTrefMatrices_com( Tref_local,  kkr%TrefLL, calc%ref_cluster_a(ilocal), getMySEcommunicator(my_mpi))
          call gatherTrefMatrices_com(dTref_local, kkr%dTrefLL, calc%ref_cluster_a(ilocal), getMySEcommunicator(my_mpi))
        enddo ! ilocal

  !------------------------------------------------------------------------------
        !$omp parallel do private(ilocal, kkr)
        do ilocal = 1, num_local_atoms
          kkr => getKKR(calc, ilocal)
          
          call GREF(emesh%EZ(IE), params%ALAT, calc%gaunts%IEND, &
                    calc%gaunts%CLEB, calc%ref_cluster_a(ilocal)%RCLS, calc%gaunts%ICLEB, &
                    calc%gaunts%LOFLM, calc%ref_cluster_a(ilocal)%NACLS, &
                    kkr%TrefLL, kkr%dTrefLL, GrefN_buffer(:,:,:,ilocal), &
                    kkr%dGrefN, kkr%LLY_G0TR(IE), &
                    dims%lmaxd, kkr%naclsd, calc%gaunts%ncleb, &
                    dims%LLY)

        enddo  ! ilocal
        !$omp endparallel do
  !------------------------------------------------------------------------------

  ! SPIN ==================================================================
  !     BEGIN do loop over spins
  ! SPIN===================================================================
  !------------------------------------------------------------------------------
  !     beginning of SMPID-parallel section
  !------------------------------------------------------------------------------
        spinloop: do ISPIN = 1, dims%NSPIND
          if (isWorkingSpinRank(my_mpi, ispin)) then

            PRSPIN = 1; if (dims%SMPID == 1) PRSPIN = ISPIN

  !------------------------------------------------------------------------------
            !$omp parallel do private(ilocal, kkr, atomdata, ldau_data, I1)
            do ilocal = 1, num_local_atoms
              kkr => getKKR(calc, ilocal)
              atomdata => getAtomData(calc, ilocal)
              ldau_data => getLDAUData(calc, ilocal)
              I1 = getAtomIndexOfLocal(calc, ilocal)

              call CALCTMAT_wrapper(atomdata, emesh, ie, ispin, params%ICST, &
                              params%NSRA, calc%gaunts, kkr%TmatN, kkr%TR_ALPH, ldau_data, params%Volterra)

              jij_data%DTIXIJ(:,:,ISPIN) = kkr%TmatN(:,:,ISPIN)  ! save t-matrix for Jij-calc.

              if (dims%LLY == 1) &  ! calculate derivative of t-matrix for Lloyd's formula
              call CALCDTMAT_wrapper(atomdata, emesh, ie, ispin, params%ICST, params%NSRA, calc%gaunts, kkr%dTdE, kkr%TR_ALPH, ldau_data, params%Volterra)

              ! t_ref-matrix of central cluster atom has index 1
              call substractReferenceTmatrix(kkr%TmatN(:,:,ISPIN), kkr%TrefLL(:,:,1), kkr%lmmaxd)

              ! do the same for derivative of T-matrix
              call substractReferenceTmatrix(kkr%dTdE(:,:,ISPIN), kkr%dTrefLL(:,:,1), kkr%lmmaxd)

              ! TmatN now contains Delta t = t - t_ref !!!
              ! dTdE now contains Delta dt !!!

              ! renormalize TR_ALPH
              kkr%TR_ALPH(ISPIN) = kkr%TR_ALPH(ISPIN) - kkr%LLY_G0TR(IE)

              call rescaleTmatrix(kkr%TmatN(:,:,ISPIN), kkr%lmmaxd, params%alat)

            enddo ! ilocal
            !$omp endparallel do
  !------------------------------------------------------------------------------

            nmesh = emesh%kmesh(IE)

            call stopTimer(single_site_timer)
            call resumeTimer(mult_scattering_timer)

  ! <<>> Multiple scattering part

            ! gather t-matrices from own truncation zone
            call gatherTmatrices_com(calc, tmatll, ispin, getMySEcommunicator(my_mpi))

            TESTARRAYLOG(3, tmatll)

            call iguess_set_energy_ind(calc%iguess_data, ie)
            call iguess_set_spin_ind(calc%iguess_data, PRSPIN)

            jij_data%active_spin = ispin

  !          WRITE(*,'(14i5)') getMyWorldRank(my_mpi),getMyAtomRank(my_mpi),getMyAtomId(my_mpi),getMySpinId(my_mpi),
  !            getMyEnergyId(my_mpi),getMySEId(my_mpi),getNumAtomRanks(my_mpi),getNumSpinRanks(my_mpi),getNumEnergyRanks(my_mpi),
  !            getNumSERanks(my_mpi),getNumWorldRanks(my_mpi),0,isMasterRank(my_mpi),isInMasterGroup(my_mpi)

  !------------------------------------------------------------------------------
            call kloopz1_new(GmatN_buffer, solv, kkr_op, precond, params%ALAT, &
                    arrays%NOFKS(nmesh), arrays%VOLBZ(nmesh), &
                    arrays%BZKP(:,:,nmesh), arrays%VOLCUB(:,nmesh), &
                    lattice_vectors%RR, &
                    GrefN_buffer, arrays%NSYMAT,arrays%DSYMLL, &
                    tmatll, arrays%lmmaxd, &
                    trunc_zone%trunc2atom_index, getMySEcommunicator(my_mpi), &
                    calc%iguess_data)
  !------------------------------------------------------------------------------
  
            if (getMyAtomRank(my_mpi) == 0 .and. params%KTE >= 0) &
              call printEnergyPoint(emesh%EZ(IE), IE, ISPIN, arrays%NOFKS(NMESH), solv%represent_stats())

            ! copy results from buffer: G_LL'^NN (E, spin) = GmatN_buffer_LL'^N(ilocal) N(ilocal)
            do ilocal = 1, num_local_atoms
              kkr => getKKR(calc, ilocal)
              kkr%GmatN(:,:,ie,ispin) = GmatN_buffer(:,:,ilocal)
            enddo ! ilocal

            call stopTimer(mult_scattering_timer)
            call resumeTimer(single_site_timer)

  ! SPIN ==================================================================
          endif ! isWorkingSpinRank
  ! SPIN ==================================================================
        enddo spinloop ! ISPIN = 1, NSPIN
  ! SPIN ==================================================================
  !     enddo loop over spins (SMPID-parallel)
  ! SPIN===================================================================

  ! =====================================================================
  ! Calculate Jij for the in CLSJIJ predefined atom pairs i,j
  ! xccpl

        if (XCCPL) then
          call jijSpinCommunication_com(my_mpi, jij_data%GMATXIJ, jij_data%DTIXIJ)

          ! calculate DTIXIJ = T_down - T_up
          call calcDeltaTupTdown(jij_data%DTIXIJ)

          JSCAL = emesh%WEZ(IE)/DBLE(jij_data%NSPIND)

          call jijLocalEnergyIntegration(my_mpi, JSCAL, jij_data%GMATXIJ, &
                                          jij_data%DTIXIJ(:,:,1), jij_data%RXIJ,&
                                          jij_data%NXIJ, jij_data%IXCP, &
                                          jij_data%RXCCLS, jij_data%JXCIJINT)
        endif ! xccpl

  ! xccpl
  ! endof Jij calculation
  ! =====================================================================

        call stopEBalanceTiming(ebalance_handler, ie)

  ! IE ====================================================================
      endif ! getMyEnergyId(my_mpi) == ebalance_handler%EPROC(IE)
  ! IE ====================================================================
    enddo ! IE = 1,IELAST
  ! IE ====================================================================
  !     enddo loop over energies (EMPID-parallel)
  ! IE ====================================================================

    call stopTimer(single_site_timer)

!     if (isMasterRank(my_mpi)) &
!       write(6, fmt='(A,I4,9A)') 'iter:',ITER,'  solver stats: ',trim(solv%represent_total_stats())
    
    
  !=======================================================================
    ! communicate information of 1..EMPID and 1..SMPID processors to MASTERGROUP
    do ilocal = 1, num_local_atoms
      kkr => getKKR(calc, ilocal)
      call collectMSResults_com(my_mpi, kkr%GmatN, kkr%LLY_GRDT, ebalance_handler%EPROC)
    enddo ! ilocal
  !=======================================================================

  ! TIME
    call OUTTIME(isMasterRank(my_mpi), 'Single Site took.....', getElapsedTime(single_site_timer), ITER)
    call OUTTIME(isMasterRank(my_mpi), 'Mult. Scat. took.....', getElapsedTime(mult_scattering_timer), ITER)

  !=======================================================================
  !     output of Jij's
  !=======================================================================
    if (XCCPL) then
      call jijReduceIntResults_com(my_mpi, jij_data%JXCIJINT)

      if (isInMasterGroup(my_mpi)) &
        call writeJiJs(I1, jij_data%RXIJ, jij_data%NXIJ, jij_data%IXCP, jij_data%RXCCLS, jij_data%JXCIJINT, jij_data%nxijd)
    endif

  !=======================================================================
  !     on the basis of new timings determine now new distribution of
  !     work to 1 .. EMPID processors - all processes SYNCED
  !=======================================================================
    call updateEBalance_com(ebalance_handler, my_mpi)

  !=======================================================================
  !     in case of IGUESS and EMPID > 1 initial guess arrays might
  !     have to be adjusted to new distributions
  !=======================================================================
    if (dims%IGUESSD == 1 .and. dims%EMPID > 1) then

      do ISPIN = 1,dims%NSPIND
        if (isWorkingSpinRank(my_mpi, ispin)) then

          PRSPIN = 1; if (dims%SMPID == 1) PRSPIN = ISPIN

          WRITELOG(3, *) "EPROC:     ", ebalance_handler%EPROC
          WRITELOG(3, *) "EPROC_old: ", ebalance_handler%EPROC_old

          call redistributeInitialGuess_com(my_mpi, calc%iguess_data%PRSC(:,:,PRSPIN), &
              ebalance_handler%EPROC, ebalance_handler%EPROC_old, &
              emesh%kmesh, arrays%NofKs)

        endif ! isWorkingSpinRank
      enddo ! ISPIN

    endif ! IGUESS == 1 .and. EMPID > 1

    call cleanup_solver(solv, kkr_op, precond)

    deallocate(tmatll)
    deallocate(atom_indices)
    deallocate(GrefN_buffer)
    deallocate(GmatN_buffer)
    deallocate(dTref_local)
    deallocate(Tref_local)
#undef clusters
#undef lattice_vectors
#undef trunc_zone

  endsubroutine ! energyLoop

  ! =============================================================================
  ! Helper routines
  ! =============================================================================

  !------------------------------------------------------------------------------
  !> Don't forget to clean up!!!
  !> Sets up TFQMR and preconditioner - matrix not setup yet!
  !> preconditioner not calculated yet
  !> Matrix setup happens later in kkrmat
  subroutine setup_solver(solv, kkr_op, precond, dims, cluster_info, lmmaxd, qmrbound, atom_indices)
    use TFQMRSolver_mod, only: TFQMRSolver
    use KKROperator_mod, only: KKROperator, create
    use BCPOperator_mod, only: BCPOperator, create
    use DimParams_mod, only: DimParams
    use ClusterInfo_mod, only: ClusterInfo
    use MultScatData_mod, only: MultScatData, create

    type(TFQMRSolver), intent(inout) :: solv
    type(KKROperator), intent(inout) :: kkr_op
    type(BCPOperator), intent(inout) :: precond
    type(DimParams), intent(in) :: dims
    type(ClusterInfo), intent(in) :: cluster_info
    integer, intent(in) :: lmmaxd
    double precision, intent(in) :: qmrbound
    integer, intent(in) :: atom_indices(:) !< indices of atoms treated at once

    call create(kkr_op)

    call solv%init(kkr_op) ! register sparse matrix and preconditioner at solver

    if (dims%bcpd == 1) then
      ! set the solver options for BCP preconditioner
      call create(precond, dims%natbld, [dims%xdim, dims%ydim, dims%zdim], cluster_info, lmmaxd)
      call solv%init_precond(precond)
    endif

    call solv%set_qmrbound(qmrbound)

    call create(kkr_op%ms, cluster_info, lmmaxd, atom_indices)

  endsubroutine ! setup_solver

  !------------------------------------------------------------------------------
  subroutine cleanup_solver(solv, kkr_op, precond)
    use TFQMRSolver_mod, only: TFQMRSolver, destroy
    use KKROperator_mod, only: KKROperator, destroy
    use BCPOperator_mod, only: BCPOperator, destroy

    type(TFQMRSolver), intent(inout) :: solv
    type(KKROperator), intent(inout) :: kkr_op
    type(BCPOperator), intent(inout) :: precond

    call destroy(solv)
    call destroy(precond)
    call destroy(kkr_op)

  endsubroutine ! cleanup_solver

  !----------------------------------------------------------------------------
  !> Print info about Energy-Point currently treated.
  subroutine printEnergyPoint(ez_point, ie, ispin, nmesh, solver_stats)
    double complex, intent(in) :: ez_point
    integer, intent(in) :: ie, ispin, nmesh
    character(len=*), intent(in) :: solver_stats
    write(6, fmt='(A,I4,A,2F11.6,A,I0,A,I0,9A)') ' ** IE =',ie,' ENERGY =',ez_point,' ispin = ',ispin,' NofKs = ',nmesh,'  ',trim(solver_stats)
  endsubroutine ! print

  !----------------------------------------------------------------------------
  !> Calculate \Delta T_up - T_down for exchange couplings calculation.
  !> The result is stored in dtixij(:,:,1)
  subroutine calcDeltaTupTdown(dtixij)
    double complex, intent(inout) :: dtixij(:,:,:)
    
    integer :: lmmaxd

    lmmaxd = size(dtixij, 1)
    dtixij(1:lmmaxd,1:lmmaxd,1) = dtixij(1:lmmaxd,1:lmmaxd,2) - dtixij(1:lmmaxd,1:lmmaxd,1)
  endsubroutine ! calc

  !----------------------------------------------------------------------------
  !> Substract diagonal reference T matrix of certain spin channel
  !> from real system's T matrix.
  subroutine substractReferenceTmatrix(TmatN, TrefLL, lmmaxd)
    integer, intent(in) :: lmmaxd
    double complex, intent(inout) :: TmatN(:,:)
    double complex, intent(in) :: TrefLL(:,:)

    integer :: lm1
    ! note: TrefLL is diagonal due to a spherical reference potential
    do lm1 = 1, lmmaxd
      TmatN(lm1,lm1) = TmatN(lm1,lm1) - TrefLL(lm1,lm1)
    enddo ! lm1

  endsubroutine ! subtract

  !------------------------------------------------------------------------------
  !> Rescale and symmetrise T-matrix.
  subroutine rescaleTmatrix(tsst_local, lmmaxd, alat)
    use Constants_mod, only: pi
    double complex, intent(inout) :: tsst_local(lmmaxd,lmmaxd)
    integer, intent(in) :: lmmaxd
    double precision, intent(in) :: alat

    integer :: lm1, lm2
    double precision :: rfctori

    rfctori = pi/alat ! = 0.5*(alat/(2*pi))^(-1)

    ! convert inverted delta_t-matrices to p.u.
    ! also a symmetrisation of the matrix is performed

    do lm2 = 1, lmmaxd
      do lm1 = 1, lm2
        tsst_local(lm1,lm2) = (tsst_local(lm1,lm2) + tsst_local(lm2,lm1))*rfctori
        tsst_local(lm2,lm1) = tsst_local(lm1,lm2) ! symmetric under exchange lm1 <--> lm2
      enddo ! lm1
    enddo ! lm2
    
  endsubroutine ! rescale

  !------------------------------------------------------------------------------
  !> Gather all tref-matrices of reference cluster.
  !> @param Tref_local   all locally calculated tref-matrices
  !> @param TrefLL       on exit all tref-matrices in ref_cluster
  subroutine gatherTrefMatrices_com(Tref_local, TrefLL, ref_cluster, communicator)
    use RefCluster_mod, only: RefCluster
    use one_sided_commZ_mod, only: copyFromZ_com

    double complex, intent(inout) :: Tref_local(:,:,:)
    double complex, intent(inout) :: TrefLL(:,:,:)
    type(RefCluster), intent(in) :: ref_cluster
    integer, intent(in) :: communicator

    integer :: chunk_size, num_local_atoms

    chunk_size = size(Tref_local, 1) * size(Tref_local, 2)
    num_local_atoms = size(Tref_local, 3)

    ASSERT (size(Tref_local, 1) == size(TrefLL, 1))
    ASSERT (size(Tref_local, 2) == size(TrefLL, 2))
    ASSERT (size(TrefLL, 3) >= ref_cluster%nacls)

    TrefLL = dcmplx(0.d0, 0.d0)

    call copyFromZ_com(TrefLL, Tref_local, ref_cluster%atom, chunk_size, num_local_atoms, communicator)

  endsubroutine ! gather

  !------------------------------------------------------------------------------
  !> Gather all t-matrices for 'ispin'-channel (from truncation zone only).
  !>
  !> Uses MPI-RMA
  subroutine gatherTmatrices_com(calc, tmatll, ispin, communicator)
    use CalculationData_mod, only: CalculationData, getNumLocalAtoms, getKKR
    use KKRresults_mod, only: KKRresults
    use one_sided_commZ_mod, only: copyFromZ_com

    type(CalculationData), intent(in) :: calc
    double complex, intent(inout) :: tmatll(:,:,:)
    integer, intent(in) :: ispin
    integer, intent(in) :: communicator

    type(KKRresults), pointer :: kkr

    integer :: ilocal, num_local_atoms, lmmaxd, chunk_size
    double complex, allocatable :: tsst_local(:,:,:)

    num_local_atoms = getNumLocalAtoms(calc)
    lmmaxd = size(tmatll, 1)

    allocate(tsst_local(lmmaxd,lmmaxd,num_local_atoms))

    chunk_size = size(tsst_local, 1)*size(tsst_local, 2)

    do ilocal = 1, num_local_atoms
      kkr => getKKR(calc, ilocal)
      tsst_local(:,:,ilocal) = kkr%TmatN(:,:,ispin)
    enddo ! ilocal

    call copyFromZ_com(tmatll, tsst_local, calc%trunc_zone%trunc2atom_index, chunk_size, num_local_atoms, communicator)

    deallocate(tsst_local)
  endsubroutine ! gather

endmodule ! ScatteringCalculation_mod
