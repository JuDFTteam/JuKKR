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
  subroutine energyLoop(iter, calc, emesh, params, dims, ebalance_handler, mp, arrays)

    USE_LOGGING_MOD
    USE_ARRAYLOG_MOD

    use TimerMpi_mod, only: TimerMpi, getElapsedTime, resetTimer, stopTimer, resumeTimer, outtime
    use EBalanceHandler_mod, only: EBalanceHandler
    use TEST_lcutoff_mod, only: DEBUG_dump_matrix

    use DimParams_mod, only: DimParams
    use InputParams_mod, only: InputParams
    use Main2Arrays_mod, only: Main2Arrays

    use CalculationData_mod, only: CalculationData, getAtomData, getLDAUData

    use KKRresults_mod, only: KKRresults
    use BasisAtom_mod, only: BasisAtom
    use JijData_mod, only: JijData
    use LDAUData_mod, only: LDAUData
    use EnergyMesh_mod, only: EnergyMesh
    
    use KKRnanoParallel_mod, only: KKRnanoParallel, isWorkingSpinRank
    use KKRnano_Comm_mod, only: jijSpinCommunication_com, jijLocalEnergyIntegration, jijReduceIntResults_com, collectMSResults_com, redistributeInitialGuess_com

    use EBalanceHandler_mod, only: EBalanceHandler, startEBalanceTiming, stopEBalanceTiming, update

    use kloopz1_mod, only: kloopz1
    use InitialGuess_mod, only: InitialGuess, iguess_set_energy_ind, iguess_set_spin_ind

    use wrappers_mod, only: calctmat_wrapper, calcdtmat_wrapper
    use jij_calc_mod, only: clsjij, writejijs, jij_data => global_jij_data

    use IterativeSolver_mod, only: IterativeSolver
    use SolverStats_mod, only: represent
    use BCPOperator_mod, only: BCPOperator
    use KKROperator_mod, only: KKROperator

    use SingleSiteRef_mod, only: TREF, GREF
    
    integer, intent(in) :: iter
    type(CalculationData), intent(inout) :: calc
    type(KKRnanoParallel), intent(in)    :: mp
    type(EBalanceHandler), intent(inout) :: ebalance_handler
    type(EnergyMesh), intent(in)         :: emesh
    type(Main2Arrays), intent(in)        :: arrays
    type(DimParams), intent(in)          :: dims
    type(InputParams), intent(in)        :: params

    ! locals

    type(BasisAtom), pointer             :: atomdata  ! referenced data does not change
    type(LDAUData), pointer              :: ldau_data ! changes

    type(IterativeSolver) :: solv
    type(KKROperator), target :: kkr_op
    type(BCPOperator), target :: precond

    double complex, parameter :: ZERO = (0.d0, 0.d0)
    type(TimerMpi) :: mult_scattering_timer, single_site_timer
    double complex :: JSCAL ! scaling factor for Jij calculation
    integer(kind=2), allocatable :: atom_indices(:)
    integer :: ie, ispin, prspin, nmesh, ist
    integer :: i1, ila, num_local_atoms, iacls
    integer :: lmmaxd
    logical :: xccpl

    double complex, allocatable :: Tref_local(:,:,:)  !< local tref-matrices
    double complex, allocatable :: dTref_local(:,:,:) !< local deriv. tref-matrices
    double complex, allocatable :: tmatLL(:,:,:) !< all t-matrices inside the truncation zone
    double complex, allocatable :: dtmatLL(:,:,:) !< all t-matrices inside the truncation zone
    double complex, allocatable :: GmatN_buffer(:,:,:) !< GmatN for all local atoms
    double complex, allocatable :: GrefN_buffer(:,:,:,:) !< GrefN for all local atoms
    double complex, allocatable :: DGrefN_buffer(:,:,:,:) !< DGrefN for all local atoms, LLY

    lmmaxd = (dims%lmaxd+1)**2

    atomdata  => getAtomData(calc, 1)
    i1 = atomdata%atom_index
    ldau_data => getLDAUData(calc, 1)
    jij_data => calc%jij_data_a(1) ! global name jij_data, jij works only with max. 1 local atom
    
#define kkr(ila) calc%kkr_a(ila)

    num_local_atoms = calc%num_local_atoms

    
    allocate(tmatLL(lmmaxd,lmmaxd,calc%trunc_zone%naez_trc)) ! allocate buffer for t-matrices
    allocate(dtmatLL(lmmaxd,lmmaxd,calc%trunc_zone%naez_trc)) ! allocate buffer for derivative of t-matrices, LLY
    allocate(Tref_local(lmmaxd,lmmaxd,num_local_atoms)) ! allocate buffers for reference t-matrices
    allocate(dTref_local(lmmaxd,lmmaxd,num_local_atoms))
    allocate(GmatN_buffer(lmmaxd,lmmaxd,num_local_atoms))
    allocate(GrefN_buffer(lmmaxd,lmmaxd,calc%clusters%naclsd,num_local_atoms))
    allocate(DGrefN_buffer(lmmaxd,lmmaxd,calc%clusters%naclsd,num_local_atoms)) ! LLY
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

      call CLSJIJ(i1, dims%NAEZ, calc%lattice_vectors%RR, calc%lattice_vectors%nrd, arrays%RBASIS, &
                  jij_data%RCUTJIJ, arrays%NSYMAT, arrays%ISYMINDEX, &
                  jij_data%IXCP, jij_data%NXCP, jij_data%NXIJ, jij_data%RXIJ, &
                  jij_data%RXCCLS, jij_data%ZKRXIJ, &
                  calc%lattice_vectors%nrd, jij_data%nxijd)

      jij_data%JXCIJINT = ZERO
      jij_data%GMATXIJ = ZERO
    endif

    ! get the indices of atoms that shall be treated at once by the process
    ! = truncation zone indices of local atoms
    do ila = 1, num_local_atoms
      atom_indices(ila) = calc%trunc_zone%local_atom_idx(calc%atom_ids(ila)) ! get global local truncation zone indices from global atom_ids 
      CHECKASSERT(atom_indices(ila) > 0)
    enddo ! ila

    ! setup the solver + bcp preconditioner, allocates a lot of memory
    ! it is good to do these allocations outside of energy loop: ToDo
    call setup_solver(solv, kkr_op, precond, dims, calc%clusters, lmmaxd, params%qmrbound, atom_indices)

    
  ! IE ====================================================================
  !     BEGIN do loop over energies (EMPID-parallel)
  ! IE ====================================================================
    do IE = 1, emesh%ielast
  ! IE ====================================================================
      if (mp%myEnergyId == ebalance_handler%EPROC(IE)) then
  ! IE ====================================================================
        call startEBalanceTiming(ebalance_handler, IE)

        WRITELOG(2, *) "Working on energy point ", IE

        ! if we had rMTref given for all atoms inside the reference cluster radius
        !   we could compute the Tref on the fly
! #define COMPUTE_tref_LOCALLY
#ifndef COMPUTE_tref_LOCALLY       
         Tref_local = ZERO
        dTref_local = ZERO
  !------------------------------------------------------------------------------
        !$omp parallel do private(ila, atomdata)
        do ila = 1, num_local_atoms
          atomdata => getAtomData(calc, ila)

          call TREF(emesh%EZ(IE), params%vref, dims%lmaxd, atomdata%rMTref, &
                    Tref_local(:,:,ila), dTref_local(:,:,ila), derive=(dims%Lly > 0)) ! ToDo: combine Tref and dTref in one array

        enddo  ! ila
        !$omp endparallel do
  !------------------------------------------------------------------------------

        ! Note: ref. system has to be recalculated at each iteration since energy mesh changes
        ! Note: TrefLL is diagonal - however, a full block matrix is stored
        ! Note: Gref is calculated in real space - usually only a few shells

        ! Exchange the reference t-matrices within the reference clusters
        ! ToDo: discuss if we can compute them once we know rMTref of all atoms in the reference cluster
        do ila = 1, num_local_atoms
          call gatherTrefMatrices_com( Tref_local,  kkr(ila)%TrefLL, calc%ref_cluster_a(ila), mp%mySEComm) ! ToDo: combine Tref and dTref in one array
          call gatherTrefMatrices_com(dTref_local, kkr(ila)%dTrefLL, calc%ref_cluster_a(ila), mp%mySEComm)
        enddo ! ila
#else
        !$omp parallel do private(ila, iacls)
        do ila = 1, num_local_atoms
          do iacls = 1, kkr(ila)%naclsd
            ! this calls tref several times with the same parameters if the local atoms are close to each other
            call TREF(emesh%EZ(IE), params%vref, dims%lmaxd, kkr(ila)%rMTref(iacls), &
                      kkr(ila)%TrefLL(:,:,iacls), kkr(ila)%dTrefLL(:,:,iacls), derive=(dims%Lly > 0))
          enddo ! iacls
        enddo ! ila
#endif
  !------------------------------------------------------------------------------
        !$omp parallel do private(ila)
        do ila = 1, num_local_atoms
          
          call GREF(emesh%EZ(IE), params%ALAT, calc%gaunts%IEND, &
                    calc%gaunts%CLEB, calc%ref_cluster_a(ila)%RCLS, calc%gaunts%ICLEB, &
                    calc%gaunts%LOFLM, calc%ref_cluster_a(ila)%NACLS, &
                    kkr(ila)%TrefLL, kkr(ila)%dTrefLL, GrefN_buffer(:,:,:,ila), &
                    DGrefN_buffer(:,:,:,ila), kkr(ila)%Lly_G0Tr(IE), &
                    dims%lmaxd, kkr(ila)%naclsd, calc%gaunts%ncleb, &
                    dims%Lly)

        enddo  ! ila
        !$omp endparallel do
  !------------------------------------------------------------------------------

  ! SPIN ==================================================================
  !     BEGIN do loop over spins
  ! SPIN===================================================================
  !------------------------------------------------------------------------------
  !     beginning of SMPID-parallel section
  !------------------------------------------------------------------------------
        spinloop: do ISPIN = 1, dims%NSPIND
          if (isWorkingSpinRank(mp, ispin)) then

            PRSPIN = 1; if (dims%SMPID == 1) PRSPIN = ISPIN

  !------------------------------------------------------------------------------
            !$omp parallel do private(ila, atomdata, ldau_data, i1)
            do ila = 1, num_local_atoms
              atomdata => getAtomData(calc, ila)
              ldau_data => getLDAUData(calc, ila)
              i1 = calc%atom_ids(ila) ! get global atom_id from local index

              call CALCTMAT_wrapper(atomdata, emesh, ie, ispin, params%ICST, params%NSRA, calc%gaunts, kkr(ila)%TmatN, kkr(ila)%Tr_alph, ldau_data, params%Volterra)

              jij_data%DTIXIJ(:,:,ISPIN) = kkr(ila)%TmatN(:,:,ISPIN) ! save t-matrix for Jij-calc.

              if (dims%Lly == 1) &  ! calculate derivative of t-matrix for Lloyd's formula
              call CALCDTMAT_wrapper(atomdata, emesh, ie, ispin, params%ICST, params%NSRA, calc%gaunts, kkr(ila)%dTdE, kkr(ila)%Tr_alph, ldau_data, params%Volterra)

              ! t_ref-matrix of central cluster atom has index 1
              call substractReferenceTmatrix(kkr(ila)%lmmaxd, arrays%NSYMAT, arrays%DSYMLL, kkr(ila)%TrefLL(:,:,1), kkr(ila)%TmatN(:,:,ISPIN))
              
              ! do the same for derivative of T-matrix
              call substractReferenceTmatrix(kkr(ila)%lmmaxd, arrays%NSYMAT, arrays%DSYMLL, kkr(ila)%dTrefLL(:,:,1), kkr(ila)%dTdE(:,:,ISPIN))

              ! TmatN now contains Delta t = t - t_ref !!!
              ! dTdE now contains Delta dt !!!

              ! renormalize Tr_alph
              kkr(ila)%Tr_alph(ISPIN) = kkr(ila)%Tr_alph(ISPIN) - kkr(ila)%Lly_G0Tr(IE)

              call rescaleTmatrix(kkr(ila)%TmatN(:,:,ISPIN), kkr(ila)%lmmaxd, params%alat)

            enddo ! ila
            !$omp endparallel do
  !------------------------------------------------------------------------------

            call stopTimer(single_site_timer)
            call resumeTimer(mult_scattering_timer)

  ! <<>> Multiple scattering part

            ! gather t-matrices from own truncation zone
            call gatherTmatrices_com(calc, tmatLL, dtmatLL, ispin, mp%mySEComm)

            TESTARRAYLOG(3, tmatLL)

            call iguess_set_energy_ind(calc%iguess_data, ie)
            call iguess_set_spin_ind(calc%iguess_data, PRSPIN)

            jij_data%active_spin = ispin

  !          WRITE(*,'(14i5)') getMyWorldRank(mp),mp%myAtomRank,getMyAtomId(mp),getMySpinId(mp),
  !            mp%myEnergyId,getMySEId(mp),getNumAtomRanks(mp),getNumSpinRanks(mp),getNumEnergyRanks(mp),
  !            getNumSERanks(mp),getNumWorldRanks(mp),0,mp%isMasterRank,mp%isInMasterGroup

            nmesh = emesh%kmesh(IE)

  !------------------------------------------------------------------------------
            call kloopz1(GmatN_buffer, solv, kkr_op, precond, params%ALAT, &
                    arrays%NOFKS(nmesh), arrays%VOLBZ(nmesh), arrays%BZKP(:,:,nmesh), arrays%VOLCUB(:,nmesh), &
                    calc%lattice_vectors%RR, & ! periodic images
                    GrefN_buffer, arrays%NSYMAT, arrays%DSYMLL, &
                    tmatLL, arrays%lmmaxd, &
                    calc%trunc_zone%global_atom_id, mp%mySEComm, &
                    calc%iguess_data, &
                    DGrefn_buffer, dtmatLL, kkr(1)%tr_alph, kkr(1)%lly_grdt(ie,ispin), calc%atom_ids(1), dims%lly) ! LLY, note: num_local_atoms must be equal to 1 
                   
  !------------------------------------------------------------------------------

            if (mp%myAtomRank == 0 .and. params%KTE >= 0) &
              call printEnergyPoint(emesh%EZ(IE), IE, ISPIN, arrays%NOFKS(NMESH), represent(solv%stats))

            ! copy results from buffer: G_LL'^NN (E, spin) = GmatN_buffer_LL'^N(ila) N(ila)
            do ila = 1, num_local_atoms
              kkr(ila)%GmatN(:,:,ie,ispin) = GmatN_buffer(:,:,ila)
            enddo ! ila

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
          call jijSpinCommunication_com(mp, jij_data%GMATXIJ, jij_data%DTIXIJ)

          ! calculate DTIXIJ = T_down - T_up
          call calcDeltaTupTdown(jij_data%DTIXIJ)

          JSCAL = emesh%WEZ(IE)/DBLE(jij_data%NSPIND)

          call jijLocalEnergyIntegration(mp, JSCAL, jij_data%GMATXIJ, &
                                          jij_data%DTIXIJ(:,:,1), jij_data%RXIJ,&
                                          jij_data%NXIJ, jij_data%IXCP, &
                                          jij_data%RXCCLS, jij_data%JXCIJINT)
        endif ! xccpl

  ! xccpl
  ! endof Jij calculation
  ! =====================================================================

        call stopEBalanceTiming(ebalance_handler, ie)

  ! IE ====================================================================
      endif ! mp%myEnergyId == ebalance_handler%EPROC(IE)
  ! IE ====================================================================
    enddo ! IE = 1,IELAST
  ! IE ====================================================================
  !     enddo loop over energies (EMPID-parallel)
  ! IE ====================================================================

    call stopTimer(single_site_timer)

!     if (mp%isMasterRank) &
!       write(6, fmt='(A,I4,9A)') 'iter:',ITER,'  solver stats: ',trim(solv%represent_total_stats())
    
    
  !=======================================================================
    ! communicate information of 1..EMPID and 1..SMPID processors to MASTERGROUP
    do ila = 1, num_local_atoms
      call collectMSResults_com(mp, kkr(ila)%GmatN, kkr(ila)%LLY_GRDT, ebalance_handler%EPROC)
    enddo ! ila
  !=======================================================================

  ! TIME
    call OUTTIME(mp%isMasterRank, 'Single Site took.....', getElapsedTime(single_site_timer), ITER)
    call OUTTIME(mp%isMasterRank, 'Mult. Scat. took.....', getElapsedTime(mult_scattering_timer), ITER)

  !=======================================================================
  !     output of Jij's
  !=======================================================================
    if (XCCPL) then
      call jijReduceIntResults_com(mp, jij_data%JXCIJINT)

      if (mp%isInMasterGroup) &
        call writeJiJs(i1, jij_data%RXIJ, jij_data%NXIJ, jij_data%IXCP, jij_data%RXCCLS, jij_data%JXCIJINT, jij_data%nxijd)
    endif

  !=======================================================================
  !     on the basis of new timings determine now new distribution of
  !     work to 1 .. EMPID processors - all processes SYNCED
  !=======================================================================
    call update(ebalance_handler, mp)

  !=======================================================================
  !     in case of IGUESS and EMPID > 1 initial guess arrays might
  !     have to be adjusted to new distributions
  !=======================================================================
    if (dims%IGUESSD == 1 .and. dims%EMPID > 1) then

      do ISPIN = 1, dims%NSPIND
        if (isWorkingSpinRank(mp, ispin)) then

          PRSPIN = 1; if (dims%SMPID == 1) PRSPIN = ISPIN

          WRITELOG(3, *) "EPROC:     ", ebalance_handler%EPROC
          WRITELOG(3, *) "EPROC_old: ", ebalance_handler%EPROC_old

          call redistributeInitialGuess_com(mp, calc%iguess_data%PRSC(:,:,PRSPIN), &
              ebalance_handler%EPROC, ebalance_handler%EPROC_old, emesh%kmesh, arrays%NofKs)

        endif ! isWorkingSpinRank
      enddo ! ISPIN

    endif ! IGUESS == 1 .and. EMPID > 1

    call cleanup_solver(solv, kkr_op, precond)

    deallocate(tmatLL, atom_indices, DGrefN_buffer, GrefN_buffer, GmatN_buffer, dTref_local, Tref_local, stat=ist)

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
    use IterativeSolver_mod, only: IterativeSolver, create
    use KKROperator_mod, only: KKROperator, create
    use BCPOperator_mod, only: BCPOperator, create
    use DimParams_mod, only: DimParams
    use ClusterInfo_mod, only: ClusterInfo

    type(IterativeSolver), intent(inout) :: solv
    type(KKROperator), intent(inout) :: kkr_op
    type(BCPOperator), intent(inout) :: precond
    type(DimParams), intent(in) :: dims
    type(ClusterInfo), intent(in) :: cluster_info
    integer, intent(in) :: lmmaxd
    double precision, intent(in) :: qmrbound
    integer(kind=2), intent(in) :: atom_indices(:) !< indices of atoms treated at once

    call create(kkr_op, cluster_info, lmmaxd, atom_indices)
    
    if (dims%bcpd == 1) then
      ! set the solver options for BCP preconditioner
      call create(precond, dims%natbld, [dims%xdim, dims%ydim, dims%zdim], cluster_info, lmmaxd)
      call create(solv, qmrbound, kkr_op, precond)
    else
      call create(solv, qmrbound, kkr_op) ! register sparse matrix and preconditioner at solver
    endif

    
  endsubroutine ! setup_solver

  !------------------------------------------------------------------------------
  subroutine cleanup_solver(solv, kkr_op, precond)
    use IterativeSolver_mod, only: IterativeSolver, destroy
    use KKROperator_mod, only: KKROperator, destroy
    use BCPOperator_mod, only: BCPOperator, destroy

    type(IterativeSolver), intent(inout) :: solv
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
  subroutine substractReferenceTmatrix(lmmaxd, nsymat, dsymll, TrefLL, TmatN)
    integer, intent(in) :: lmmaxd
    integer, intent(in) :: nsymat
    double complex, intent(in) :: dsymll(:,:,:) !> dim(lmmaxd,lmmaxd,nsymat)
    double complex, intent(in) :: TrefLL(:,:) !> dim(lmmaxd,lmmaxd)
    double complex, intent(inout) :: TmatN(:,:) !> dim(lmmaxd,lmmaxd)
    
    double complex :: uTu_sum(lmmaxd,lmmaxd), uT(lmmaxd,lmmaxd)
    double precision :: denom
    integer :: isym, lm
   
    double complex, parameter :: cone=(1.d0, 0.d0), zero=(0.d0, 0.d0)
 
    ! note: TrefLL is diagonal due to a spherical reference potential, therefore, we only subtract the diagonal elements
    do lm = 1, lmmaxd
      TmatN(lm,lm) = TmatN(lm,lm) - TrefLL(lm,lm)
    enddo ! lm

    !------------------------------------------------- SYMMETRISE TMATN
    uTu_sum(:,:) = TmatN(:,:) ! copy, the 1st entry is the unity operation
    do isym = 2, nsymat
      call zgemm('n', 'n', lmmaxd, lmmaxd, lmmaxd, cone, dsymll(1,1,isym), lmmaxd, TmatN, lmmaxd, zero, uT, lmmaxd)
      call zgemm('n', 'c', lmmaxd, lmmaxd, lmmaxd, cone, uT, lmmaxd, dsymll(1,1,isym), lmmaxd, cone, uTu_sum, lmmaxd)
    enddo ! isym

    denom = 1.d0/dble(nsymat)
    TmatN(:,:) = uTu_sum(:,:)*denom ! average
    !------------------------------------------------- SYMMETRISE TMATN

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
      do lm1 = 1, lm2 ! triangular loop including the diagonal
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

    chunk_size = size(Tref_local, 1)*size(Tref_local, 2)
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
  subroutine gatherTmatrices_com(calc, tmatLL, dtde, ispin, communicator)
    use CalculationData_mod, only: CalculationData
    use KKRresults_mod, only: KKRresults
    use one_sided_commZ_mod, only: copyFromZ_com

    type(CalculationData), intent(in) :: calc
    double complex, intent(inout) :: tmatLL(:,:,:)
    double complex, intent(inout) :: dtde(:,:,:) ! LLY
    integer, intent(in) :: ispin
    integer, intent(in) :: communicator

    integer :: ila, num_local_atoms, lmmaxd, chunk_size
    double complex, allocatable :: tsst_local(:,:,:)
    double complex, allocatable :: dtsst_local(:,:,:) ! LLY

    num_local_atoms = calc%num_local_atoms
    lmmaxd = size(tmatLL, 1)

    allocate(tsst_local(lmmaxd,lmmaxd,num_local_atoms))
    allocate(dtsst_local(lmmaxd,lmmaxd,num_local_atoms)) ! LLY

    chunk_size = size(tsst_local, 1)*size(tsst_local, 2)

    do ila = 1, num_local_atoms
      tsst_local(:,:,ila) = kkr(ila)%TmatN(:,:,ispin)
      dtsst_local(:,:,ila) = kkr(ila)%dtde(:,:,ispin) ! LLY
    enddo ! ila

    call copyFromZ_com(tmatLL, tsst_local, calc%trunc_zone%global_atom_id, chunk_size, num_local_atoms, communicator)
    call copyFromZ_com(dtde, dtsst_local, calc%trunc_zone%global_atom_id, chunk_size, num_local_atoms, communicator)

    deallocate(tsst_local)
    deallocate(dtsst_local)
  endsubroutine ! gather
#undef kkr
endmodule ! ScatteringCalculation_mod
