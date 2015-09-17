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

  CONTAINS

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
subroutine energyLoop(iter, calc_data, emesh, params, dims, ebalance_handler, my_mpi, arrays)

  USE_LOGGING_MOD
  USE_ARRAYLOG_MOD

  use TimerMpi_mod, only: TimerMpi, getElapsedTime, resetTimer, stopTimer, resumeTimer, outtime
  use EBalanceHandler_mod, only: EBalanceHandler
  use TEST_lcutoff_mod, only: DEBUG_dump_matrix

  use DimParams_mod, only: DimParams
  use InputParams_mod, only: InputParams
  use Main2Arrays_mod, only: Main2Arrays

  use CalculationData_mod, only: CalculationData, &
    getTruncationZone, getGaunts, getKKR, getRefCluster, getJijData, getInitialGuessData, getAtomData, &
    getLDAUData, getClusterInfo, getLatticeVectors, getNumLocalAtoms, getAtomIndexOfLocal
    
  use KKRresults_mod, only: KKRresults
  use GauntCoefficients_mod, only: GauntCoefficients
  use BasisAtom_mod, only: BasisAtom
  use JijData_mod, only: JijData
  use LDAUData_mod, only: LDAUData
  use EnergyMesh_mod, only: EnergyMesh
  
  use KKRnanoParallel_mod, only: KKRnanoParallel, isMasterRank, getMyEnergyId, getMySEcommunicator, isWorkingSpinRank, getMyAtomRank, isInMasterGroup
  use KKRnano_Comm_mod, only: jijSpinCommunication_com, jijLocalEnergyIntegration, jijReduceIntResults_com, collectMSResults_com, redistributeInitialGuess_com
  
  use EBalanceHandler_mod, only: EBalanceHandler, startEBalanceTiming, stopEBalanceTiming, updateEBalance_com
  
  use kloopz1_mod, only: kloopz1_new
  use TruncationZone_mod, only: TruncationZone
  use ClusterInfo_mod, only: ClusterInfo
  use RefCluster_mod, only: RefCluster, LatticeVectors
  use InitialGuess_mod, only: InitialGuess, iguess_set_energy_ind, iguess_set_spin_ind

  use wrappers_mod, only: calctmat_wrapper, calcdtmat_wrapper
  use jij_calc_mod, only: clsjij, writejijs, global_jij_data

  use TFQMRSolver_mod, only: TFQMRSolver
  use BCPOperator_mod, only: BCPOperator
  use KKROperator_mod, only: KKROperator

  use SingleSiteRef_mod, only: TREF, GREF
  
  integer, intent(in) :: iter
  type(CalculationData), intent(inout) :: calc_data
  type(KKRnanoParallel), intent(in)    :: my_mpi
  type(EBalanceHandler), intent(inout) :: ebalance_handler
  type(EnergyMesh), intent(in)         :: emesh
  type(Main2Arrays), intent(in)        :: arrays
  type(DimParams), intent(in)          :: dims
  type(InputParams), intent(in)        :: params

  !---------- locals ----------------------------

  type(BasisAtom), pointer             :: atomdata  ! referenced data does not change
  type(KKRresults), pointer            :: kkr       ! changes
  type(GauntCoefficients), pointer     :: gaunts    ! never changes
  type(LDAUData), pointer              :: ldau_data ! changes
  type(JijData), pointer               :: jij_data  ! changes
  type(TruncationZone), pointer        :: trunc_zone  ! never changes
  type(ClusterInfo), pointer           :: clusters    ! never changes
  type(RefCluster), pointer            :: ref_cluster ! never changes
  type(LatticeVectors), pointer        :: lattice_vectors ! never changes
  type(InitialGuess), pointer          :: iguess_data ! changes

  type(TFQMRSolver), target :: solv
  type(KKROperator), target :: kkr_op
  type(BCPOperator), target :: precond

  double complex, parameter :: CZERO = (0.d0, 0.d0)
  type(TimerMpi) :: mult_scattering_timer, single_site_timer
  double complex :: JSCAL ! scaling factor for Jij calculation
  integer, allocatable :: atom_indices(:)
  integer :: ie, ispin, prspin, nmesh
  integer :: I1, ilocal, num_local_atoms
  integer :: lmmaxd
  logical :: xccpl

  double complex, allocatable :: Tref_local(:,:,:)  !< local tref-matrices
  double complex, allocatable :: DTref_local(:,:,:) !< local deriv. tref-matrices
  double complex, allocatable :: tmatll(:,:,:) !< all t-matrices
  double complex, allocatable :: GmatN_buffer(:,:,:) !< GmatN for all local atoms
  double complex, allocatable :: GrefN_buffer(:,:,:,:) !< GrefN for all local atoms

  lmmaxd = (dims%lmaxd+1)**2

  trunc_zone => getTruncationZone(calc_data)
  gaunts    => getGaunts(calc_data)
  atomdata  => getAtomData(calc_data, 1)
  I1 = atomdata%atom_index
  kkr       => null()
  ldau_data => getLDAUData(calc_data, 1)
  jij_data  => getJijData(calc_data, 1)
  clusters  => getClusterInfo(calc_data)
  lattice_vectors  => getLatticeVectors(calc_data)
  iguess_data => getInitialGuessData(calc_data)

  num_local_atoms = getNumLocalAtoms(calc_data)

  global_jij_data => getJijData(calc_data, 1)

  ! allocate buffer for t-matrices
  allocate(tmatll(lmmaxd, lmmaxd, trunc_zone%naez_trc))
  ! allocate buffers for reference t-matrices
  allocate( Tref_local(lmmaxd, lmmaxd, num_local_atoms))
  allocate(DTref_local(lmmaxd, lmmaxd, num_local_atoms))
  allocate(GmatN_buffer(lmmaxd,lmmaxd,num_local_atoms))
  allocate(GrefN_buffer(lmmaxd,lmmaxd, clusters%naclsd, num_local_atoms))

  allocate(atom_indices(num_local_atoms))

  if (params%jij .and. num_local_atoms > 1) then
    if (isMasterRank(my_mpi)) write(*,*) "Jij and num_local_atoms > 1 not supported."
    STOP
  endif

  if (params%ldau .and. num_local_atoms > 1) then
    if (isMasterRank(my_mpi)) write(*,*) "LDA+U and num_local_atoms > 1 not supported."
    STOP
  endif

  xccpl = .false.

  call resetTimer(mult_scattering_timer)
  call stopTimer(mult_scattering_timer)

  call resetTimer(single_site_timer)

  prspin = 1

  ! calculate exchange couplings only at last self-consistency step and when Jij=true
  if ((ITER==params%SCFSTEPS).and.params%JIJ) XCCPL = .true.

  if (XCCPL) then

    ! Trigger jij-calculation
    jij_data%do_jij_calculation = .true.

    call CLSJIJ(I1,dims%NAEZ,lattice_vectors%RR,lattice_vectors%nrd,arrays%RBASIS, &
                jij_data%RCUTJIJ,arrays%NSYMAT,arrays%ISYMINDEX, &
                jij_data%IXCP,jij_data%NXCP,jij_data%NXIJ,jij_data%RXIJ, &
                jij_data%RXCCLS,jij_data%ZKRXIJ, &
                lattice_vectors%nrd, jij_data%nxijd)

    jij_data%JXCIJINT = CZERO
    jij_data%GMATXIJ = CZERO

  endif

  ! get the indices of atoms that shall be treated at once by the process
  ! = truncation zone indices of local atoms
  do ilocal = 1, num_local_atoms
    atom_indices(ilocal) = getAtomIndexOfLocal(calc_data, ilocal)
    atom_indices(ilocal) = trunc_zone%index_map(atom_indices(ilocal))
    CHECKASSERT(atom_indices(ilocal) > 0)
  enddo

  ! setup the solver + bcp preconditioner, allocates a lot of memory
  ! it is good to do these allocations outside of energy loop
  call setup_solver(solv, kkr_op, precond, dims, clusters, lmmaxd, params%qmrbound, atom_indices)

! IE ====================================================================
!     BEGIN do loop over energies (EMPID-parallel)
! IE ====================================================================
  do IE = 1, emesh%ielast
! IE ====================================================================
    if (getMyEnergyId(my_mpi)==ebalance_handler%EPROC(IE)) then
! IE ====================================================================
      call startEBalanceTiming(ebalance_handler, IE)

      WRITELOG(2, *) "Working on energy point ", IE

      Tref_local = CZERO
      DTref_local = CZERO
!------------------------------------------------------------------------------
      !$omp parallel do private(ilocal, kkr, atomdata)
      do ilocal = 1, num_local_atoms
        kkr => getKKR(calc_data, ilocal)
        atomdata  => getAtomData(calc_data, ilocal)
!------------------------------------------------------------------------------

        call TREF(emesh%EZ(IE),arrays%VREF,dims%LMAXD,atomdata%RMTref, &
                  Tref_local(:,:,ilocal), DTref_local(:,:,ilocal), dims%LLY)

!------------------------------------------------------------------------------
      enddo  ! ilocal
      !$omp endparallel do
!------------------------------------------------------------------------------

      ! Note: ref. system has to be recalculated at each iteration
      ! since energy mesh changes
      ! Note: TREFLL is diagonal - however full matrix is stored
      ! Note: Gref is calculated in real space - usually only a few shells

      ! Exchange the reference t-matrices within reference clusters
      do ilocal = 1, num_local_atoms
        kkr => getKKR(calc_data, ilocal)
        ref_cluster => getRefCluster(calc_data, ilocal)

        call gatherTrefMatrices_com( Tref_local,  kkr%TrefLL, ref_cluster, getMySEcommunicator(my_mpi))
        call gatherTrefMatrices_com(DTref_local, kkr%DTrefLL, ref_cluster, getMySEcommunicator(my_mpi))
      enddo

!------------------------------------------------------------------------------
      !$omp parallel do private(ilocal, kkr, ref_cluster)
      do ilocal = 1, num_local_atoms
        kkr => getKKR(calc_data, ilocal)
        ref_cluster => getRefCluster(calc_data, ilocal)
!------------------------------------------------------------------------------

        call GREF(emesh%EZ(IE),params%ALAT,gaunts%IEND, &
                      gaunts%CLEB,ref_cluster%RCLS,gaunts%ICLEB, &
                      gaunts%LOFLM,ref_cluster%NACLS, &
                      kkr%TREFLL,kkr%DTREFLL, GrefN_buffer(:,:,:,ilocal), &
                      kkr%DGREFN, kkr%LLY_G0TR(IE), &
                      dims%lmaxd, kkr%naclsd, gaunts%ncleb, &
                      dims%LLY)

!------------------------------------------------------------------------------
      enddo  ! ilocal
      !$omp endparallel do
!------------------------------------------------------------------------------

! SPIN ==================================================================
!     BEGIN do loop over spins
! SPIN===================================================================
!------------------------------------------------------------------------------
!     beginning of SMPID-parallel section
!------------------------------------------------------------------------------
      spinloop: do ISPIN = 1,dims%NSPIND
        if (isWorkingSpinRank(my_mpi, ispin)) then

          PRSPIN = 1; if (dims%SMPID == 1) PRSPIN = ISPIN

!------------------------------------------------------------------------------
          !$omp parallel do private(ilocal, kkr, atomdata, ldau_data, jij_data, I1)
          do ilocal = 1, num_local_atoms
            kkr => getKKR(calc_data, ilocal)
            atomdata => getAtomData(calc_data, ilocal)
            ldau_data => getLDAUData(calc_data, ilocal)
            jij_data  => getJijData(calc_data, ilocal)
            I1 = getAtomIndexOfLocal(calc_data, ilocal)
!------------------------------------------------------------------------------

            call CALCTMAT_wrapper(atomdata, emesh, ie, ispin, params%ICST, &
                            params%NSRA, gaunts, kkr%TMATN, kkr%TR_ALPH, ldau_data)

            jij_data%DTIXIJ(:,:,ISPIN) = kkr%TMATN(:,:,ISPIN)  ! save t-matrix for Jij-calc.

            if (dims%LLY==1) then  ! calculate derivative of t-matrix for Lloyd's formula
              call CALCDTMAT_wrapper(atomdata, emesh, ie, ispin, params%ICST, &
                            params%NSRA, gaunts, kkr%DTDE, kkr%TR_ALPH, ldau_data)
            endif

            ! t_ref-matrix of central cluster atom has index 1
            call substractReferenceTmatrix(kkr%TMATN(:,:,ISPIN), kkr%TREFLL(:,:,1), kkr%LMMAXD)

            ! do the same for derivative of T-matrix
            call substractReferenceTmatrix(kkr%DTDE(:,:,ISPIN), kkr%DTREFLL(:,:,1), kkr%LMMAXD)

            ! TMATN now contains Delta t = t - t_ref !!!
            ! DTDE now contains Delta dt !!!

            ! renormalize TR_ALPH
            kkr%TR_ALPH(ISPIN) = kkr%TR_ALPH(ISPIN) - kkr%LLY_G0TR(IE)

            call rescaleTmatrix(kkr%TMATN(:,:,ISPIN), kkr%lmmaxd, params%alat)

!------------------------------------------------------------------------------
          enddo ! ilocal
          !$omp endparallel do
!------------------------------------------------------------------------------

          NMESH = arrays%KMESH(IE)

          if ( getMyAtomRank(my_mpi)==0 ) then
            if (params%KTE >= 0) call printEnergyPoint(emesh%EZ(IE), IE, ISPIN, NMESH)
          endif

          call stopTimer(single_site_timer)
          call resumeTimer(mult_scattering_timer)

! <<>> Multiple scattering part

          ! gather t-matrices from own truncation zone
          call gatherTmatrices_com(calc_data, tmatll, ispin, getMySEcommunicator(my_mpi))

          TESTARRAYLOG(3, tmatll)

          call iguess_set_energy_ind(iguess_data, ie)
          call iguess_set_spin_ind(iguess_data, PRSPIN)

          jij_data%active_spin = ispin

!          WRITE(*,'(14i5)') getMyWorldRank(my_mpi),getMyAtomRank(my_mpi),getMyAtomId(my_mpi),getMySpinId(my_mpi),
!            getMyEnergyId(my_mpi),getMySEId(my_mpi),getNumAtomRanks(my_mpi),getNumSpinRanks(my_mpi),getNumEnergyRanks(my_mpi),
!            getNumSERanks(my_mpi),getNumWorldRanks(my_mpi),getMasterRank(my_mpi),isMasterRank(my_mpi),isInMasterGroup(my_mpi)

!------------------------------------------------------------------------------
          call KLOOPZ1_new(GmatN_buffer, solv, kkr_op, precond, params%ALAT, &
          arrays%NOFKS(NMESH),arrays%VOLBZ(NMESH), &
          arrays%BZKP(:,:,NMESH),arrays%VOLCUB(:,NMESH), &
          lattice_vectors%RR, &
          GrefN_buffer, arrays%NSYMAT,arrays%DSYMLL, &
          tmatll, arrays%lmmaxd, lattice_vectors%nrd, &
          trunc_zone%trunc2atom_index, getMySEcommunicator(my_mpi), &
          iguess_data)
!------------------------------------------------------------------------------

          ! copy results from buffer: G_LL'^NN (E, spin) =
          !                           GmatN_buffer_LL'^N(ilocal) N(ilocal)
          do ilocal = 1, num_local_atoms
            kkr => getKKR(calc_data, ilocal)
            kkr%GMATN(:,:,ie,ispin) = GmatN_buffer(:,:,ilocal)
          enddo

          call stopTimer(mult_scattering_timer)
          call resumeTimer(single_site_timer)

        endif
      enddo spinloop                          ! ISPIN = 1,NSPIN
!------------------------------------------------------------------------------
!        endof SMPID-parallel section
!------------------------------------------------------------------------------
! SPIN ==================================================================
!     enddo loop over spins
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
      endif

! xccpl
! endof Jij calculation
! =====================================================================

      call stopEBalanceTiming(ebalance_handler, ie)

! IE ====================================================================
    endif
! IE ====================================================================

  enddo                   ! IE = 1,IELAST

! IE ====================================================================
!     enddo loop over energies (EMPID-parallel)
! IE ====================================================================

  call stopTimer(single_site_timer)

!=======================================================================
!communicate information of 1..EMPID and 1..SMPID processors to MASTERGROUP
  do ilocal = 1, num_local_atoms
    kkr => getKKR(calc_data, ilocal)
    call collectMSResults_com(my_mpi, kkr%GMATN, kkr%LLY_GRDT, ebalance_handler%EPROC)
  enddo
!=======================================================================

! TIME
  call OUTTIME(isMasterRank(my_mpi), 'Single Site took.....', getElapsedTime(single_site_timer), ITER)
  call OUTTIME(isMasterRank(my_mpi), 'Mult. Scat. took.....', getElapsedTime(mult_scattering_timer), ITER)

!=======================================================================
!     output of Jij's
!=======================================================================
  if (XCCPL) then

    call jijReduceIntResults_com(my_mpi, jij_data%JXCIJINT)

    if (isInMasterGroup(my_mpi)) then
      call writeJiJs(I1,jij_data%RXIJ,jij_data%NXIJ,jij_data%IXCP, &
                     jij_data%RXCCLS,jij_data%JXCIJINT, jij_data%nxijd)
    endif
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
  if ((dims%IGUESSD==1) .and. (dims%EMPID>1)) then

    do ISPIN = 1,dims%NSPIND
      if (isWorkingSpinRank(my_mpi, ispin)) then

        if (dims%SMPID==1) then
          PRSPIN = ISPIN
        else
          PRSPIN = 1
        endif

        WRITELOG(3, *) "EPROC:     ", ebalance_handler%EPROC
        WRITELOG(3, *) "EPROC_old: ", ebalance_handler%EPROC_old

        call redistributeInitialGuess_com(my_mpi, iguess_data%PRSC(:,:,PRSPIN), &
             ebalance_handler%EPROC, ebalance_handler%EPROC_old, &
             arrays%KMESH, arrays%NofKs)

      endif
    enddo

  endif  ! IGUESS == 1 .and. EMPID > 1

  call cleanup_solver(solv, kkr_op, precond)

  deallocate(tmatll)
  deallocate(atom_indices)
  deallocate(GrefN_buffer)
  deallocate(GmatN_buffer)
  deallocate(DTref_local)
  deallocate(Tref_local)

endsubroutine energyLoop

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
  use KKROperator_mod, only: KKROperator
  use BCPOperator_mod, only: BCPOperator
  use DimParams_mod, only: DimParams
  use ClusterInfo_mod, only: ClusterInfo
  use SolverOptions_mod, only: SolverOptions
  use MultScatData_mod, only: MultScatData, createMultScatData

  type(TFQMRSolver), intent(inout) :: solv
  type(KKROperator), intent(inout) :: kkr_op
  type(BCPOperator), intent(inout) :: precond
  type(DimParams), intent(in) :: dims
  type(ClusterInfo), intent(in) :: cluster_info
  integer, intent(in) :: lmmaxd
  double precision, intent(in) :: qmrbound
  integer, dimension(:), intent(in) :: atom_indices !< indices of atoms treated at once

  type(SolverOptions) :: solver_opts
  type(MultScatData), pointer :: ms

  ! set the solver options for bcp preconditioner
  solver_opts%bcp = dims%bcpd
  solver_opts%xdim = dims%xdim
  solver_opts%ydim = dims%ydim
  solver_opts%zdim = dims%zdim
  solver_opts%natbld = dims%natbld

  call kkr_op%create()
  ms => kkr_op%get_ms_workspace()

  ! register sparse matrix and preconditioner at solver
  call solv%init(kkr_op)

  if (solver_opts%bcp == 1) then
    call precond%create(solver_opts, cluster_info, lmmaxd)
    call solv%init_precond(precond)
  endif

  call solv%set_qmrbound(qmrbound)

  call createMultScatData(ms, cluster_info, lmmaxd, atom_indices)

endsubroutine

!------------------------------------------------------------------------------
subroutine cleanup_solver(solv, kkr_op, precond)
  use TFQMRSolver_mod, only: TFQMRSolver
  use KKROperator_mod, only: KKROperator
  use BCPOperator_mod, only: BCPOperator
  use MultScatData_mod, only: MultScatData, destroyMultScatData

  type(TFQMRSolver), intent(inout) :: solv
  type(KKROperator), intent(inout) :: kkr_op
  type(BCPOperator), intent(inout) :: precond

  type(MultScatData), pointer :: ms

  ms => kkr_op%get_ms_workspace()

  call solv%destroy()
  call precond%destroy()

  call destroyMultScatData(ms)
  call kkr_op%destroy()

endsubroutine

!----------------------------------------------------------------------------
!> Print info about Energy-Point currently treated.
!>
subroutine printEnergyPoint(ez_point, ie, ispin, nmesh)
  double complex, intent(in) :: ez_point
  integer, intent(in) :: ie, ispin, nmesh
  write (6,'(A,I3,A,2(1X,F10.6),A,I3,A,I3)') ' ** IE = ',ie,' ENERGY =',ez_point,' KMESH = ',nmesh,' ISPIN = ',ispin
endsubroutine

!----------------------------------------------------------------------------
!> Calculate \Delta T_up - T_down for exchange couplings calculation.
!> The result is stored in dtixij(:,:,1)
subroutine calcDeltaTupTdown(dtixij)
  double complex, intent(inout) :: dtixij(:,:,:)
  
  integer :: lmmaxd

  lmmaxd = size(dtixij,1)
  dtixij(1:lmmaxd,1:lmmaxd,1) = dtixij(1:lmmaxd,1:lmmaxd,2) - dtixij(1:lmmaxd,1:lmmaxd,1)
endsubroutine

!----------------------------------------------------------------------------
!> Substract diagonal reference T matrix of certain spin channel
!> from real system's T matrix.
subroutine substractReferenceTmatrix(tmatn, trefll, lmmaxd)
  integer, intent(in) :: lmmaxd
  double complex, intent(inout) :: tmatn(:,:)
  double complex, intent(in) :: trefll(:,:)

  integer :: lm1
  ! note: trefll is diagonal! - spherical reference potential
  do lm1 = 1,lmmaxd
    tmatn(lm1,lm1) =  tmatn(lm1,lm1) - trefll(lm1,lm1)
  enddo ! lm1

endsubroutine

  !------------------------------------------------------------------------------
  !> Rescale and symmetrise T-matrix.
  subroutine rescaleTmatrix(tsst_local, lmmaxd, alat)
    double complex, intent(inout) :: tsst_local(lmmaxd,lmmaxd)
    integer, intent(in) :: lmmaxd
    double precision, intent(in) :: alat

    integer :: lm1, lm2
    double precision :: rfctor

  !     rfctor=a/(2*pi) conversion factor to p.u.
    rfctor = alat/(8.d0*atan(1.d0))           ! = alat/(2*pi)

! --> convert inverted delta_t-matrices to p.u.
!     also a symmetrisation of the matrix is performed

    do lm2 = 1, lmmaxd
      do lm1 = 1, lm2
        tsst_local(lm1,lm2) = 0.5d0/rfctor * (tsst_local(lm1,lm2) + tsst_local(lm2,lm1))
        tsst_local(lm2,lm1) = tsst_local(lm1,lm2) ! symmtric
      enddo ! lm1
    enddo ! lm2
    
  endsubroutine

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

  !-------------------
  integer chunk_size
  integer num_local_atoms

  chunk_size = size(Tref_local, 1) * size(Tref_local, 2)
  num_local_atoms = size(Tref_local, 3)

  ASSERT (size(Tref_local, 1) == size(TrefLL, 1))
  ASSERT (size(Tref_local, 2) == size(TrefLL, 2))
  ASSERT (size(TrefLL, 3) >= ref_cluster%nacls)

  TrefLL = (0.d0, 0.d0)

  call copyFromZ_com(TrefLL, Tref_local, ref_cluster%atom, chunk_size, num_local_atoms, communicator)

endsubroutine

!------------------------------------------------------------------------------
!> Gather all t-matrices for 'ispin'-channel (from truncation zone only).
!>
!> Uses MPI-RMA
subroutine gatherTmatrices_com(calc_data, tmatll, ispin, communicator)
  use CalculationData_mod, only: CalculationData, getNumLocalAtoms, getTruncationZone, getKKR
  use KKRresults_mod, only: KKRresults
  use TruncationZone_mod, only: TruncationZone
  use one_sided_commZ_mod, only: copyFromZ_com

  type(CalculationData), intent(in) :: calc_data
  double complex, intent(inout) :: tmatll(:,:,:)
  integer, intent(in) :: ispin
  integer, intent(in) :: communicator

  type(KKRresults), pointer :: kkr
  type(TruncationZone), pointer :: trunc_zone

  integer :: ilocal
  integer :: num_local_atoms
  integer :: lmmaxd

  integer :: chunk_size
  double complex, allocatable :: tsst_local(:,:,:)

  num_local_atoms = getNumLocalAtoms(calc_data)
  trunc_zone => getTruncationZone(calc_data)
  lmmaxd = size(tmatll, 1)

  allocate(tsst_local(lmmaxd, lmmaxd, num_local_atoms))

  chunk_size = size(tsst_local, 1) * size(tsst_local, 2)

  do ilocal = 1, num_local_atoms
    kkr => getKKR(calc_data, ilocal)
    tsst_local(:,:,ilocal) = kkr%TMATN(:,:,ispin)
  enddo

  call copyFromZ_com(tmatll, tsst_local, trunc_zone%trunc2atom_index, chunk_size, num_local_atoms, communicator)

  deallocate(tsst_local)

endsubroutine

endmodule
