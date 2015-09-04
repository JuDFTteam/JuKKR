! KKRnano
! massive parallel KKR for nanoscaled systems

#include "DebugHelpers/test_macros.h"
#include "DebugHelpers/logging_macros.h"

program KKRnano

  use Logging_mod, only:    !import no name here, just mention it for the module dependency 
  USE_LOGGING_MOD

  use KKRnanoParallel_mod, only: KKRnanoParallel, isMasterRank, isActiveRank, isInMasterGroup, &
    getMyWorldRank, getMyAtomRank, getMyActiveCommunicator, getMySEcommunicator, createKKRnanoParallel, destroyKKRnanoParallel
  use KKRnano_Comm_mod, only: setKKRnanoNumThreads, printKKRnanoInfo, communicatePotential

  use main2_aux_mod, only: printDoubleLineSep, is_abort_by_rank0, writeIterationTimings
  use EnergyMesh_mod, only: EnergyMesh, createEnergyMesh, readEnergyMesh, readEnergyMeshSemi, &
    broadcastEnergyMesh_com, updateEnergyMesh, updateEnergyMeshSemi, destroyEnergyMesh, writeEnergyMesh, writeEnergyMeshSemi

  use RadialMeshData_mod, only: RadialMeshData
  use BasisAtom_mod, only: BasisAtom
  use LDAUData_mod, only: LDAUData

  use TimerMpi_mod, only: TimerMpi, getElapsedTime, resetTimer, outtime
  use EBalanceHandler_mod, only: EBalanceHandler, createEBalanceHandler, initEBalanceHandler, setEqualDistribution, destroyEBalanceHandler

  use wrappers_mod, only: rhocore_wrapper

  use DimParams_mod, only: DimParams, createDimParams, destroyDimParams
  use InputParams_mod, only: InputParams, readInputParamsFromFile
  use Main2Arrays_mod, only: Main2Arrays, createMain2Arrays, readMain2Arrays, destroyMain2Arrays

  use ScatteringCalculation_mod, only: energyloop
  use ProcessKKRresults_mod, only: processKKRresults, output_forces

  use CalculationData_mod, only: CalculationData, createCalculationData, prepareMadelung, &
    getNumLocalAtoms, getAtomData, getLDAUData, getAtomIndexOfLocal, destroyCalculationData
  
  use kkr0_mod, only: main0
  implicit none

  type (CalculationData) :: calc_data

  !     .. Local Scalars ..

  type (TimerMpi) :: program_timer
  type (TimerMpi) :: iteration_timer
  type (EBalanceHandler) :: ebalance_handler

  integer :: ITER
  integer :: I1
  integer :: ilocal
  integer :: num_local_atoms

  type(KKRnanoParallel) :: my_mpi

  integer :: flag

  type (EnergyMesh), target     :: emesh
  type (Main2Arrays), target    :: arrays
  type (DimParams), target      :: dims
  type (InputParams)            :: params

  type (RadialMeshData), pointer :: mesh
  type (BasisAtom), pointer      :: atomdata
  type (LDAUData), pointer       :: ldau_data
  !----------------------------------------------------------------------------
  character(len=16)              :: arg
  integer                        :: ios, ilen
  
  call get_command_argument(1, arg, ilen, ios)
  selectcase (arg)
  case ('--prepare')
    call main0(checkmode=0) ! call former kkr0.exe
    stop
  case ('--check')
    call main0(checkmode=1) ! former kkr0.exe without overwriting the binary files *.unf
    stop
  case ('--help')
    call get_command_argument(0, arg, ilen, ios)
    write(*,'(9A)') 'usage: ',trim(arg),' [options]'
    write(*,'(A)') '  options:', &
    '    --prepare             Former kkr0.exe functionality', &
    '    --help                This command line help function', &
    '    --check               Check input files for errors', &
    '  options to be implemented in the future:', &
    '    --convert             Converter to and from ASCII files', &
    ''
    stop
  case default
    
  endselect ! arg
  !----------------------------------------------------------------------------

  call createDimParams(dims, 'inp0.unf') ! read dim. parameters from file 'inp0.unf'

! -----------------------------------------------------------------------------
  call createKKRnanoParallel(my_mpi, dims%num_atom_procs, dims%SMPID, dims%EMPID)
  call setKKRnanoNumThreads(dims%nthrds)
  call printKKRnanoInfo(my_mpi, dims%nthrds)

  if (isMasterRank(my_mpi)) then
#ifdef NO_LOCKS_MPI
    write(*,*) "NO_LOCKS_MPI defined: Not using MPI RMA locks. Does not scale well."
#endif
#ifdef IDENTICAL_REF
    write(*,*) "IDENTICAL_REF defined: assuming identical reference clusters."
#endif
#ifdef DEBUG_NO_ELECTROSTATICS
    write(*,*) "DEBUG_NO_ELECTROSTATICS: no electrostatics - results are wrong."
#endif
  endif

!------------------------------------------------------------------------------

  if (getMyWorldRank(my_mpi) < 128) then ! max. 128 logfiles
    OPENLOG(getMyWorldRank(my_mpi), 3)
  else
    OPENLOG(getMyWorldRank(my_mpi), 0)
  endif

! ========= TIMING =========================================================
    call resetTimer(program_timer)
    if (isMasterRank(my_mpi)) then
      open (2,file='time-info',form='formatted')
    endif
!========= TIMING END ======================================================

!-----------------------------------------------------------------------------
! Array allocations BEGIN
!-----------------------------------------------------------------------------
  call createMain2Arrays(arrays, dims)
  call readMain2Arrays(arrays, 'arrays.unf') ! every process does this!
!-----------------------------------------------------------------------------
! Array allocations END
!-----------------------------------------------------------------------------

  flag = readInputParamsFromFile(params, 'input.unf')

 ! ======================================================================
 ! =                     End read in variables                          =
 ! ======================================================================

  if (params%JIJ .and. (dims%nspind /= 2)) then
    write(*,*) "ERROR: Jij calculation not possible for spin-unpolarized calc."
    stop
  endif

  if (dims%LLY /= 0) then
    write(*,*) "WARNING: Lloyds formula not supported in this version. Set LLY=0"
    stop
  endif

!=====================================================================
!     processors not fitting in NAEZ*LMPID*SMPID*EMPID do nothing ...
! ... and wait after SC-ITER loop
!=====================================================================

  ! This if closes much later!
  if (isActiveRank(my_mpi)) then

!+++++++++++ pre self-consistency preparation

    !--------------------------------------------------------------------------
    call createCalculationData(calc_data, dims, params, arrays, my_mpi)
    num_local_atoms = getNumLocalAtoms(calc_data)
    !--------------------------------------------------------------------------

    call createEnergyMesh(emesh, dims%iemxd) !!!!

    ! TO DO: Getting rid of the many if-clauses used for semicore contour!
    if(params%use_semicore==1) then
      call readEnergyMeshSemi(emesh) !every process does this!  !!!!
    else
      call readEnergyMesh(emesh)  !every process does this!  !!!!
    endif

    call OUTTIME(isMasterRank(my_mpi),'input files read.....', getElapsedTime(program_timer), 0)

#ifndef DEBUG_NO_ELECTROSTATICS
    call prepareMadelung(calc_data, arrays)
#endif

    call OUTTIME(isMasterRank(my_mpi),'Madelung sums calc...', getElapsedTime(program_timer), 0)

    call createEBalanceHandler(ebalance_handler, emesh%ielast)
    call initEBalanceHandler(ebalance_handler, my_mpi)
    call setEqualDistribution(ebalance_handler, (emesh%NPNT1 == 0))

! ######################################################################
! ######################################################################
    do ITER = 1, params%SCFSTEPS
! ######################################################################
! ######################################################################

      call resetTimer(iteration_timer)

      if (isMasterRank(my_mpi)) then
        call printDoubleLineSep(unit_number = 2)
        call OUTTIME(isMasterRank(my_mpi),'started at ..........', getElapsedTime(program_timer),ITER)
        call printDoubleLineSep(unit_number = 2)
      endif

      WRITELOG(2, *) "Iteration atom-rank ", ITER, getMyAtomRank(my_mpi)

      ! New: instead of reading potential every time, communicate it
      ! between energy and spin processes of same atom
      do ilocal = 1, num_local_atoms
        atomdata => getAtomdata(calc_data, ilocal)
        call communicatePotential(my_mpi, atomdata%potential%VISP, &
                                  atomdata%potential%VINS, atomdata%core%ECORE)
      enddo ! ilocal

      ! Core relaxation - only mastergroup needs results
      if (isInMasterGroup(my_mpi)) then
        ! Not threadsafe: intcor, intin, intout have a save statement
        !!!$omp parallel do private(ilocal, atomdata)
        do ilocal = 1, num_local_atoms
          atomdata => getAtomdata(calc_data, ilocal)
          if (params%use_semicore==1) then
            call RHOCORE_wrapper(emesh%EBOTSEMI, params%NSRA, atomdata)
          else
            call RHOCORE_wrapper(emesh%E1, params%NSRA, atomdata)
          endif
        enddo ! ilocal
        !!!$omp end parallel do
      endif

! LDA+U ! TODO: doesn't work for num_local_atoms > 1
      if (params%LDAU) then
        ! For now only 1 atom per process is supported (1 local atom)
        CHECKASSERT(num_local_atoms == 1)

        atomdata      => getAtomData(calc_data, 1)
        ldau_data     => getLDAUData(calc_data, 1)
        mesh          => atomdata%mesh_ptr
        I1 = getAtomIndexOfLocal(calc_data, 1)

        ldau_data%EREFLDAU = emesh%EFERMI
        ldau_data%EREFLDAU = 0.48    ! ???

        call LDAUINIT(I1,ITER,params%NSRA,ldau_data%NLDAU,ldau_data%LLDAU, &
                      ldau_data%ULDAU,ldau_data%JLDAU,ldau_data%EREFLDAU, &
                      atomdata%potential%VISP,ldau_data%NSPIND,mesh%R,mesh%DRDI, &
                      atomdata%Z_nuclear,mesh%IPAN,mesh%IRCUT, &
                      ldau_data%PHILDAU,ldau_data%UMLDAU,ldau_data%WMLDAU, &
                      ldau_data%lmaxd, mesh%irmd, mesh%ipand)

      endif
! LDA+U

      call OUTTIME(isMasterRank(my_mpi),'initialized .........', getElapsedTime(program_timer),ITER)

      ! Scattering calculations - that is what KKR is all about
      ! output: (some contained as references in calc_data)
      ! ebalance_handler, kkr (!), jij_data, ldau_data
      call energyLoop(iter, calc_data, emesh, params, dims, ebalance_handler, my_mpi, arrays)

      call OUTTIME(isMasterRank(my_mpi),'G obtained ..........', getElapsedTime(program_timer),ITER)

!----------------------------------------------------------------------
! BEGIN only processes in master-group are working
!----------------------------------------------------------------------
      flag = 0
      if (isInMasterGroup(my_mpi)) then
        ! output: (some contained as references in calc_data)
        ! atomdata, densities, broyden, ldau_data,
        ! emesh (only correct for master)
        flag = processKKRresults(iter, calc_data, my_mpi, emesh, dims, params, arrays, program_timer)

      endif
!----------------------------------------------------------------------
! END only processes in master-group are working
!----------------------------------------------------------------------

      if (isMasterRank(my_mpi)) then
        ! only MASTERRANK updates, other ranks get it broadcasted later
        ! (although other processes could update themselves)

        if (params%use_semicore == 1) then
          call updateEnergyMeshSemi(emesh)
        else
          call updateEnergyMesh(emesh)
        endif

        ! write file 'energy_mesh'
        if (emesh%NPOL /= 0) emesh%EFERMI = emesh%E2  ! if not a DOS-calculation E2 coincides with Fermi-Energy

        if (params%use_semicore == 1) then
          call writeEnergyMeshSemi(emesh)
        else
          call writeEnergyMesh(emesh)
        endif


        call printDoubleLineSep()
        call writeIterationTimings(ITER, getElapsedTime(program_timer), getElapsedTime(iteration_timer))

      endif

      if (is_abort_by_rank0(flag, getMyActiveCommunicator(my_mpi))) exit

      call broadcastEnergyMesh_com(my_mpi, emesh)

! ######################################################################
! ######################################################################
    enddo ! ITER ! SC ITERATION LOOP ITER=1, SCFSTEPS
! ######################################################################
! ######################################################################

    ! write forces if requested, master-group only
    if (params%kforce == 1 .and. isInMasterGroup(my_mpi)) then
      call output_forces(calc_data, 0, getMyAtomRank(my_mpi), getMySEcommunicator(my_mpi))
    endif

    if (isMasterRank(my_mpi)) close(2) ! time-info

    call destroyEBalanceHandler(ebalance_handler)
    call destroyEnergyMesh(emesh)
    call destroyCalculationData(calc_data)

  endif ! active Ranks

  CLOSELOG

!------------------------------------------------------------------------------
  call destroyMain2Arrays(arrays)
  call destroyDimParams(dims)
!------------------------------------------------------------------------------

!=====================================================================
!     processors not fitting in NAEZ*LMPID do nothing ...
! ... and wait here
!=====================================================================
! Free KKRnano mpi resources

  call destroyKKRnanoParallel(my_mpi)
  
endprogram KKRnano
