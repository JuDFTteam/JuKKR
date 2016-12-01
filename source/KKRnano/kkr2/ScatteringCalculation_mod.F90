#include "DebugHelpers/logging_macros.h"
#include "DebugHelpers/test_array_log.h"
#include "DebugHelpers/test_macros.h"

!> @author Modularisation: Elias Rabel
module ScatteringCalculation_mod
  use Logging_mod, only:    !import no name here, just mention it for the module dependency 
  use arraytest2_mod, only: !import no name here, just mention it for the module dependency 
implicit none
  private
  public  :: energyLoop, gatherrMTref_com

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

    use TimerMpi_mod, only: TimerMpi, getTime, createTimer, stopTimer, startTimer, outTime, outTimeStats
    use EBalanceHandler_mod, only: EBalanceHandler

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
    use KKRnano_Comm_mod, only: jijSpinCommunication_com, jijLocalEnergyIntegration, jijReduceIntResults_com, collectMSResults_com, redistributeInitialGuess

    use EBalanceHandler_mod, only: EBalanceHandler, startEBalanceTiming, stopEBalanceTiming, update

    use kloopz1_mod, only: kloopz1
    use InitialGuess_mod, only: InitialGuess
    
    use wrappers_mod, only: calctmat_wrapper, calcdtmat_wrapper
    use jij_calc_mod, only: clsjij, writejijs, jij_data => global_jij_data

    use IterativeSolver_mod, only: IterativeSolver
    use SolverStats_mod, only: represent, GiFlops
    use BCPOperator_mod, only: BCPOperator
    use KKROperator_mod, only: KKROperator

    use SingleSiteRef_mod, only: tref, gref
    
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
    type(TimerMpi) :: mult_scattering_timer, single_site_timer, reference_green_timer, kpoint_timer, kernel_timer
    double complex :: JSCAL ! scaling factor for Jij calculation
    integer(kind=2), allocatable :: atom_indices(:)
    integer :: ie, ispin, prspin, nmesh, ist
    integer :: i1, ila, num_local_atoms, iacls
    integer :: lmmaxd
    logical :: xccpl

    double complex, allocatable :: tmatLL(:,:,:,:) !< all t-matrices inside the truncation zone
    double complex, allocatable :: GmatN_buffer(:,:,:) !< GmatN for all local atoms
    double complex, allocatable :: GrefN_buffer(:,:,:,:,:) !< GrefN for all local atoms

    lmmaxd = (dims%lmaxd+1)**2

    atomdata  => getAtomData(calc, 1)
    i1 = atomdata%atom_index
    ldau_data => getLDAUData(calc, 1)
    jij_data => calc%jij_data_a(1) ! global name jij_data, jij works only with max. 1 local atom
    
#define kkr(ila) calc%kkr_a(ila)

    num_local_atoms = calc%num_local_atoms

    
    allocate(tmatLL(lmmaxd,lmmaxd,calc%trunc_zone%naez_trc,0:dims%Lly)) ! allocate buffer for t-matrices
    allocate(GmatN_buffer(lmmaxd,lmmaxd,num_local_atoms))
    allocate(GrefN_buffer(lmmaxd,lmmaxd,0:dims%Lly,calc%clusters%naclsd,num_local_atoms))
    allocate(atom_indices(num_local_atoms))

    if (params%jij  .and. num_local_atoms > 1) stop "Jij and num_local_atoms > 1 not supported."
    if (params%ldau .and. num_local_atoms > 1) stop "LDA+U and num_local_atoms > 1 not supported."

    xccpl = .false.

    call createTimer(mult_scattering_timer)
    call createTimer(single_site_timer)
    call createTimer(reference_green_timer)
    call createTimer(kpoint_timer)
    call createTimer(kernel_timer)

    prspin = 1

    ! calculate exchange couplings only at last self-consistency step and when Jij=true
    if (iter == params%SCFSTEPS .and. params%JIJ) XCCPL = .true. ! activate in last SCF iteration

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
      atom_indices(ila) = calc%trunc_zone%trunc_atom_idx(calc%atom_ids(ila)) ! get global local truncation zone indices from global atom_ids 
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

  !------------------------------------------------------------------------------
        ! if we have rMTref given for all atoms inside the reference cluster radius we can compute the Tref on the fly

        call startTimer(reference_green_timer)          
        
        !$omp parallel do private(ila, iacls)
        do ila = 1, num_local_atoms

          do iacls = 1, calc%ref_cluster_a(ila)%nacls
            ! this calls tref several times with the same parameters if the local atoms are close to each other
            call tref(emesh%EZ(IE), params%vref, dims%lmaxd, kkr(ila)%rMTref(iacls), &
                      kkr(ila)%Tref_ell(:,iacls), kkr(ila)%dTref_ell(:,iacls), derive=(dims%Lly > 0))
          enddo ! iacls
        
          call gref(emesh%EZ(IE), params%ALAT, calc%gaunts%IEND, &
                    calc%gaunts%CLEB, calc%ref_cluster_a(ila)%RCLS, calc%gaunts%ICLEB, &
                    calc%gaunts%LOFLM, calc%ref_cluster_a(ila)%nacls, &
                    kkr(ila)%Tref_ell, kkr(ila)%dTref_ell, GrefN_buffer(:,:,:,:,ila), &
                    kkr(ila)%Lly_G0Tr(IE), &
                    dims%lmaxd, calc%gaunts%ncleb, dims%Lly)

        enddo  ! ila
        !$omp endparallel do
  !------------------------------------------------------------------------------

        call stopTimer(reference_green_timer)          
  
  ! SPIN ==================================================================
  !     BEGIN do loop over spins
  ! SPIN===================================================================
  !------------------------------------------------------------------------------
  !     beginning of SMPID-parallel section
  !------------------------------------------------------------------------------
        spinloop: do ISPIN = 1, dims%NSPIND
          if (isWorkingSpinRank(mp, ispin)) then

            PRSPIN = 1; if (dims%SMPID == 1) PRSPIN = ISPIN
            
            call startTimer(single_site_timer)          

  !------------------------------------------------------------------------------
            !$omp parallel do private(ila, atomdata, ldau_data, i1)
            do ila = 1, num_local_atoms
              atomdata => getAtomData(calc, ila)
              ldau_data => getLDAUData(calc, ila)
              i1 = calc%atom_ids(ila) ! get global atom_id from local index

              call CALCTMAT_wrapper(atomdata, emesh, ie, ispin, params%ICST, params%NSRA, calc%gaunts, kkr(ila)%TmatN, kkr(ila)%Tr_alph, ldau_data, params%Volterra)

              jij_data%DTIXIJ(:,:,ISPIN) = kkr(ila)%TmatN(:,:,ISPIN) ! save t-matrix for Jij-calc.

              ! t_ref-matrix of central cluster atom has index 1
              call substractReferenceTmatrix(arrays%dsymLL(:,:,:arrays%NSYMAT), kkr(ila)%Tref_ell(:,1), kkr(ila)%TmatN(:,:,ISPIN))

              if (dims%Lly == 1) then  ! calculate derivative of t-matrix for Lloyd's formula
                call CALCDTMAT_wrapper(atomdata, emesh, ie, ispin, params%ICST, params%NSRA, calc%gaunts, kkr(ila)%dTmatN, kkr(ila)%Tr_alph, ldau_data, params%Volterra)
               
                ! do the same for derivative of T-matrix
                call substractReferenceTmatrix(arrays%dsymLL(:,:,:arrays%NSYMAT), kkr(ila)%dTref_ell(:,1), kkr(ila)%dTmatN(:,:,ISPIN))
              endif ! Lly

              ! TmatN now contains Delta t = t - t_ref !!!
              ! dTdE now contains Delta dt !!!

              kkr(ila)%Tr_alph(ISPIN) = kkr(ila)%Tr_alph(ISPIN) - kkr(ila)%Lly_G0Tr(IE) ! renormalize Tr_alph

              call rescaleTmatrix(kkr(ila)%TmatN(:,:,ISPIN), params%alat)

            enddo ! ila
            !$omp endparallel do
  !------------------------------------------------------------------------------

            call stopTimer(single_site_timer)
            call startTimer(mult_scattering_timer)

  ! <<>> Multiple scattering part

            ! gather t-matrices from own truncation zone
            call gatherTmatrices_com(calc, tmatLL, ispin, mp%mySEComm)

            TESTARRAYLOG(3, tmatLL)

            jij_data%active_spin = ispin

            nmesh = emesh%kmesh(IE)

  !------------------------------------------------------------------------------
            call kloopz1(GmatN_buffer, solv, kkr_op, precond, params%ALAT, &
                    arrays%NOFKS(nmesh), arrays%VOLBZ(nmesh), arrays%BZKP(:,:,nmesh), arrays%VOLCUB(:,nmesh), &
                    calc%lattice_vectors%RR, & ! periodic images
                    GrefN_buffer, arrays%dsymLL(:,:,:arrays%NSYMAT), &
                    tmatLL, &
                    calc%trunc_zone%global_atom_id, mp%mySEComm, &
                    calc%iguess_data, IE, PRSPIN, &
!                     dtmatLL, &
                    kkr(1)%tr_alph, kkr(1)%Lly_grdt(ie,ispin), calc%atom_ids(1), dims%Lly, & ! LLY, note: num_local_atoms must be equal to 1 
                    params%solver, kpoint_timer, kernel_timer)
  !------------------------------------------------------------------------------

            if (mp%myAtomRank == 0 .and. params%KTE >= 0) &
              call printEnergyPoint(emesh%EZ(IE), IE, ISPIN, arrays%NOFKS(NMESH), represent(solv%stats))

            ! copy results from buffer: G_LL'^NN (E, spin) = GmatN_buffer_LL'^N(ila) N(ila)
            do ila = 1, num_local_atoms
              kkr(ila)%GmatN(:,:,ie,ispin) = GmatN_buffer(:,:,ila)
            enddo ! ila

            call stopTimer(mult_scattering_timer)

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

          jij_data%DTIXIJ(:,:,1) = jij_data%DTIXIJ(:,:,2) - jij_data%DTIXIJ(:,:,1) ! calculate DTIXIJ = T_down - T_up

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

    call outTime(mp%isMasterRank, 'Reference system Green  took', getTime(reference_green_timer), iter)
    call outTime(mp%isMasterRank, 'Single site scattering  took', getTime(single_site_timer), iter)
    call outTime(mp%isMasterRank, 'Multi. site scattering  took', getTime(kpoint_timer), iter)
    if (mp%isMasterRank) call outTimeStats(reference_green_timer, 'Reference G stats:')
    if (mp%isMasterRank) call outTimeStats(single_site_timer,     'Single site stats:')
    if (mp%isMasterRank) call outTimeStats(kpoint_timer,          'Mult. scat. stats:') ! per k-point
    if (mp%isMasterRank) call outTimeStats(kernel_timer,          'KKRoperator stats:') ! per invocation of the KKR_operator
!   if (mp%isMasterRank) call outTimeStats(mult_scattering_timer, 'Multi. site stats:') ! this timer is ...
!   !         ... only energy point resolved, high variance expected due to different k-point mesh sizes
!   call outTime(mp%isMasterRank, 'Multi. site scattering  took', getTime(mult_scattering_timer), iter)

!   if (mp%isMasterRank) write(6, fmt='(A,I4,9A)') 'iter:',iter,'  solver stats: ',trim(solv%represent_total_stats())
    if (mp%isMasterRank) write(6, fmt='(a,i4,9(a,f0.6))') 'iter:',iter,'  aggregate ',GiFlops/1024.d0,' TiFlop on master process' ! useful flops in the iterative solver part

  !=======================================================================
    ! communicate information of 1..EMPID and 1..SMPID processors to MASTERGROUP
    do ila = 1, num_local_atoms
      call collectMSResults_com(mp, kkr(ila)%GmatN, kkr(ila)%LLY_GRDT, ebalance_handler%EPROC)
    enddo ! ila
  !=======================================================================

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

          if (calc%iguess_data%prec == 1) then
            call redistributeInitialGuess(mp, calc%iguess_data%prsc(:,:,PRSPIN), ebalance_handler%EPROC, ebalance_handler%EPROC_old, emesh%kmesh, arrays%NofKs)
          elseif (calc%iguess_data%prec == 2) then
            call redistributeInitialGuess(mp, calc%iguess_data%prsz(:,:,PRSPIN), ebalance_handler%EPROC, ebalance_handler%EPROC_old, emesh%kmesh, arrays%NofKs)
          endif

        endif ! isWorkingSpinRank
      enddo ! ISPIN

    endif ! IGUESS == 1 .and. EMPID > 1

    call cleanup_solver(solv, kkr_op, precond)

    deallocate(tmatLL, atom_indices, GrefN_buffer, GmatN_buffer, stat=ist) ! ignore status

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

    call create(kkr_op, cluster_info, lmmaxd, atom_indices, dims%Lly)
    
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
  !> Substract diagonal reference T matrix of certain spin channel
  !> from real system's T matrix.
  subroutine substractReferenceTmatrix(dsymLL, TrefLL, TmatN)
    double complex, intent(in) :: dsymLL(:,:,:) !> dim(lmmaxd,lmmaxd,nsymat)
    double complex, intent(in) :: TrefLL(0:) !> dim(0:lmax) ! emm-degenerate and diagonal
    double complex, intent(inout) :: TmatN(:,:) !> dim(lmmaxd,lmmaxd)
    
    integer :: lmax, lmmaxd, nsymat, isym, lm, ell, emm, ist
    double complex, allocatable :: uTu_sum(:,:), uT(:,:)
    double precision :: denom
    double complex, parameter :: cone=(1.d0, 0.d0), zero=(0.d0, 0.d0)

    lmax = size(TrefLL, 1) - 1
    lmmaxd = (lmax + 1)**2
    nsymat = size(dsymLL, 3)

    ! note: TrefLL is diagonal due to a spherical reference potential, therefore, we only subtract the diagonal elements
    ASSERT( all(shape(TmatN) == [lmmaxd, lmmaxd]) )
    ASSERT( all(shape(dsymLL) == [lmmaxd, lmmaxd, nsymat]) )

    do ell = 0, lmax
      do emm = -ell, ell
        lm = ell*ell + ell + emm + 1
        TmatN(lm,lm) = TmatN(lm,lm) - TrefLL(ell) ! Tref is stored emm-degenerate and diagonal
      enddo ! emm
    enddo ! ell

    allocate(uTu_sum(lmmaxd,lmmaxd), uT(lmmaxd,lmmaxd))
    !------------------------------------------------- SYMMETRISE TMATN
    uTu_sum(:,:) = TmatN(:,:) ! copy, since the 1st entry is the unity operation, start loop from 2
    do isym = 2, nsymat
      call zgemm('n', 'n', lmmaxd, lmmaxd, lmmaxd, cone, dsymLL(1,1,isym), lmmaxd, TmatN, lmmaxd, zero, uT, lmmaxd)
      call zgemm('n', 'c', lmmaxd, lmmaxd, lmmaxd, cone, uT, lmmaxd, dsymLL(1,1,isym), lmmaxd, cone, uTu_sum, lmmaxd)
    enddo ! isym

    denom = 1.d0/dble(nsymat)
    TmatN(:,:) = uTu_sum(:,:)*denom ! average
    !------------------------------------------------- SYMMETRISE TMATN
    deallocate(uTu_sum, uT, stat=ist) ! ignore status
  endsubroutine ! subtract

  !------------------------------------------------------------------------------
  !> Rescale and symmetrise T-matrix.
  subroutine rescaleTmatrix(tsst_local, alat)
    use Constants_mod, only: pi
    double complex, intent(inout) :: tsst_local(:,:) ! dim(lmmaxd,lmmaxd)
    double precision, intent(in) :: alat

    integer :: lm1, lm2, lmmaxd
    double precision :: rfctori

    rfctori = pi/alat ! = 0.5*(alat/(2*pi))^(-1)
    lmmaxd = size(tsst_local, 2)
    
    ! convert inverted delta_t-matrices to p.u.
    ! also a symmetrisation of the matrix is performed

    do lm2 = 1, lmmaxd
      do lm1 = 1, lm2 ! triangular loop including the diagonal
        tsst_local(lm1,lm2) = (tsst_local(lm1,lm2) + tsst_local(lm2,lm1))*rfctori
        tsst_local(lm2,lm1) =  tsst_local(lm1,lm2) ! symmetric under exchange lm1 <--> lm2
      enddo ! lm1
    enddo ! lm2
    
  endsubroutine ! rescale

  !------------------------------------------------------------------------------
  !> Gather all t-matrices for 'ispin'-channel (from truncation zone only).
  !>
  !> Uses MPI-RMA
  subroutine gatherTmatrices_com(calc, tmatLL, ispin, communicator)
    use CalculationData_mod, only: CalculationData
    use KKRresults_mod, only: KKRresults
!   use one_sided_commZ_mod, only: copyFromZ_com
    
    use two_sided_commZ_mod, only: reference_sys_com
    use ExchangeTable_mod, only: ExchangeTable
    use ExchangeTable_mod, only: create, destroy
    type(ExchangeTable) :: xTable

    type(CalculationData), intent(in) :: calc
    double complex, intent(inout) :: tmatLL(:,:,:,0:) !> dim(lmmaxd,lmmaxd,num_local_atoms,0:Lly)
    integer, intent(in) :: ispin
    integer, intent(in) :: communicator

    integer :: ila, num_local_atoms, lmmaxd, ist, chunk_size, iLly, nLly
    double complex, allocatable :: tsst_local(:,:,:,:)

    num_local_atoms = calc%num_local_atoms
    lmmaxd = size(tmatLL, 1)
    ASSERT(  size(tmatLL, 2) == lmmaxd )
    ASSERT(  size(tmatLL, 3) == size(calc%trunc_zone%global_atom_id) ) ! num_trunc_atoms
    nLly   = size(tmatLL, 4)

    allocate(tsst_local(lmmaxd,lmmaxd,num_local_atoms,0:nLly-1))

    do ila = 1, num_local_atoms
      tsst_local(:,:,ila,0) = kkr(ila)%TmatN(:,:,ispin)
      if (nLly > 1) &
      tsst_local(:,:,ila,1) = kkr(ila)%dTmatN(:,:,ispin) ! LLY
    enddo ! ila

    chunk_size = size(tmatLL, 1)*size(tmatLL, 2)
    
    call create(xTable, calc%trunc_zone%global_atom_id, communicator, max_local_atoms=num_local_atoms)
    do iLly = 0, nLly - 1
!     call copyFromZ_com(tmatLL(:,:,:,iLly), tsst_local(:,:,:,iLly), calc%trunc_zone%global_atom_id, chunk_size, num_local_atoms, communicator)
      call reference_sys_com(xTable, chunk_size, tsst_local(:,:,:,iLly), tmatLL(:,:,:,iLly))
    enddo ! iLly
    call destroy(xTable)

    deallocate(tsst_local, stat=ist) ! ignore status
  endsubroutine ! gather
#undef kkr

  !------------------------------------------------------------------------------
  !> Gather all rMTref values of the reference cluster.
  !> @param rMTref_local all locaLly determined rMTref value
  !> @param rMTref       on exit all rMTref value in ref_cluster
  subroutine gatherrMTref_com(rMTref_local, rMTref, ref_cluster, communicator)
    use RefCluster_mod, only: RefCluster
    use one_sided_commD_mod, only: copyFromD_com

    double precision, intent(in) :: rMTref_local(:) ! (num_local_atoms)
    double precision, intent(out) :: rMTref(:) ! (ref_cluster%nacls)
    type(RefCluster), intent(in) :: ref_cluster
    integer, intent(in) :: communicator
    
    double precision, allocatable :: rMTref_all(:,:,:), rMTref_loc(:,:,:)

    allocate(rMTref_all(1,1,size(rMTref, 1)), rMTref_loc(1,1,size(rMTref_local, 1)))
    
    rMTref_all = 0.d0
    
    rMTref_loc(1,1,:) = rMTref_local(:) ! in
    
    call copyFromD_com(rMTref_all, rMTref_loc, ref_cluster%atom, 1, size(rMTref_local, 1), communicator)
    
    rMTref(:) = rMTref_all(1,1,:) ! out

    deallocate(rMTref_all, rMTref_loc)
  endsubroutine ! gather

endmodule ! ScatteringCalculation_mod
