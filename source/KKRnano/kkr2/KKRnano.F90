! KKRnano
! massive parallel KKR for nanoscaled systems

#include "DebugHelpers/test_macros.h"
#include "DebugHelpers/logging_macros.h"

program KKRnano
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)

  use Logging_mod, only:    !import no name here, just mention it for the module dependency 
  USE_LOGGING_MOD

  use KKRnanoParallel_mod, only: KKRnanoParallel, isMasterRank, isActiveRank, isInMasterGroup
  use KKRnanoParallel_mod, only: getMyWorldRank, getMyAtomRank, getMyActiveCommunicator, getMySEcommunicator 
  use KKRnanoParallel_mod, only: createKKRnanoParallel, destroy

  use KKRnano_Comm_mod, only: setKKRnanoNumThreads, printKKRnanoInfo, communicatePotential

  use main2_aux_mod, only: printDoubleLineSep, is_abort_by_rank0, writeIterationTimings
  use EnergyMesh_mod, only: EnergyMesh, create, destroy, load, store, update, broadcast

  use BasisAtom_mod, only: BasisAtom
  use LDAUData_mod, only: LDAUData

  use TimerMpi_mod, only: TimerMpi, getElapsedTime, resetTimer, outTime
  use EBalanceHandler_mod, only: EBalanceHandler, createEBalanceHandler, initEBalanceHandler, setEqualDistribution, destroy

  use wrappers_mod, only: rhocore_wrapper

  use DimParams_mod, only: DimParams, load, destroy
  use InputParams_mod, only: InputParams, readInputParamsFromFile
  use Main2Arrays_mod, only: Main2Arrays, createMain2Arrays, readMain2Arrays, destroy

  use ScatteringCalculation_mod, only: energyloop
  use ProcessKKRresults_mod, only: processKKRresults, output_forces

  use CalculationData_mod, only: CalculationData, create, prepareMadelung, destroy
  use CalculationData_mod, only: getNumLocalAtoms, getAtomData, getLDAUData, getAtomIndexOfLocal
  use BrillouinZoneMesh_mod, only: BrillouinZoneMesh
  
  use KKRzero_mod, only: main0
  use PotentialConverter_mod, only: kkrvform
  
  implicit none

  type(CalculationData) :: calc_data
  type(TimerMpi) :: program_timer
  type(TimerMpi) :: iteration_timer
  type(EBalanceHandler) :: ebalance_handler

  integer :: ITER, I1, ilocal, num_local_atoms, flag, ios, ilen

  type(KKRnanoParallel) :: my_mpi

  type(EnergyMesh)  :: emesh
  type(Main2Arrays) :: arrays
  type(DimParams)   :: dims
  type(InputParams) :: params

  type(BasisAtom), pointer :: atomdata
  type(LDAUData), pointer  :: ldau_data
  
  type(BrillouinZoneMesh) :: kmesh(8)
  
  external :: MPI_Init
  character(len=16)              :: arg
  
  call MPI_Init(ios) ! --> needs to be called here, otherwise MPI_Abort and MPI_Wtime cannot be used during toolbox functionalities
  
  call get_command_argument(1, arg, ilen, ios)
  selectcase (arg)
  case ('--prepare')
    call main0(checkmode=0) ! call former kkr0.exe
    stop
  case ('--check')
    call main0(checkmode=1) ! former kkr0.exe without overwriting the binary files *.unf
    stop
  case ('--convert')
    call kkrvform()
    stop
  case ('--help')
    call get_command_argument(0, arg, ilen, ios)
    write(*,'(9A)') 'usage: ',trim(arg),' [options]'
    write(*,'(A)') '  options:', &
    '    --prepare             Former kkr0.exe functionality', &
    '    --help                This command line help function', &
    '    --check               Check input files for errors', &
    '    --convert             Converter to and from ASCII files', &
    ''
    stop
  case default
    ! start the former kkr2.exe    
  endselect ! arg

  ! from here the former kkr2.exe actions are performed
  
  call load(dims, 'inp0.unf') ! read dimension parameters from file 'inp0.unf'

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

#ifdef DEBUG_NO_VINS
      write(*,*) 'WARNING: DEBUG_NO_VINS enabled; non-spherical part of potential is set to zero!!!'
#endif

#ifdef PRINT_MTRADII
      write(*,*) 'PRINT_MTRADII enabled; rmtref (0.995 * max_muffin_tin by default) and ', &
      'radius_muffin_tin (0.98 * max_muffin_tin by default) are printed consecutively to files mtradii_out'
#endif

#ifdef USE_MTRADII
      write(*,*) 'USE_MTRADII enabled; rmtref (0.995 * max_muffin_tin by default) and ', &
        'radius_muffin_tin (0.98 * max_muffin_tin by default) are taken from files mtradii_in'
#endif
  endif ! is master
  

  if (getMyWorldRank(my_mpi) < 128) then ! max. 128 logfiles
    OPENLOG(getMyWorldRank(my_mpi), 3)
  else
    OPENLOG(getMyWorldRank(my_mpi), 0)
  endif

  call resetTimer(program_timer)
  if (isMasterRank(my_mpi)) open(2, file='time-info', form='formatted', action='write')

  call createMain2Arrays(arrays, dims)
  call readMain2Arrays(arrays, 'arrays.unf') ! every process does this!

  flag = readInputParamsFromFile(params, 'input.unf')
  ! done reading variables


  if (params%JIJ .and. dims%nspind /= 2) die_here("Jij calculation not possible for spin-unpolarized calc.")

  if (dims%LLY /= 0) die_here("Lloyds formula not supported in this version. Set lly=0")

  !=====================================================================
  ! processors not fitting in NAEZ*LMPID*SMPID*EMPID do nothing ... and wait after SC-ITER loop
  !=====================================================================

  if (isActiveRank(my_mpi)) then

    ! pre self-consistency preparations

    assert(dims%naez > 0) 
    call create(calc_data, dims, params, arrays, my_mpi) ! createCalculationData
    num_local_atoms = getNumLocalAtoms(calc_data)

    call create(emesh, dims%iemxd) ! createEnergyMesh
    call load(emesh) ! every process does this!!!

    call outTime(isMasterRank(my_mpi),'input files read.....', getElapsedTime(program_timer), 0)

#ifdef DEBUG_NO_ELECTROSTATICS
    warn(6, "preprocessor define has switched off electrostatics for debugging")
#else
    call prepareMadelung(calc_data, arrays)
#endif

    call outTime(isMasterRank(my_mpi),'Madelung sums calc...', getElapsedTime(program_timer), 0)

    call createEBalanceHandler(ebalance_handler, emesh%ielast)
    call initEBalanceHandler(ebalance_handler, my_mpi)
    call setEqualDistribution(ebalance_handler, (emesh%NPNT1 == 0))

#ifdef DEBUG_NO_VINS
    do ilocal = 1, num_local_atoms
      atomdata => getAtomdata(calc_data, ilocal)
      atomdata%potential%VINS = 0.0   
    enddo ! ilocal
#endif

#ifdef PRINT_MTRADII
   do ilocal = 1, num_local_atoms
      atomdata => getAtomdata(calc_data, ilocal)
      I1 = getAtomIndexOfLocal(calc_data, ilocal)
      write(num, '(A,I7.7)') "mtradii_out.",I1
      open(20, file=num, form='formatted')
      write(20,*) atomdata%rmtref
      write(20,*) atomdata%radius_muffin_tin
      endfile (20)
      close (20)
    enddo ! ilocal
#endif

#ifdef USE_MTRADII
   do ilocal = 1, num_local_atoms
      atomdata => getAtomdata(calc_data, ilocal)
      I1 = getAtomIndexOfLocal(calc_data, ilocal)
      write(num, '(A,I7.7)') "mtradii_in.",I1
      open(20, file=num, form='formatted')
      read(20,*) atomdata%rmtref
      read(20,*) atomdata%radius_muffin_tin
      endfile (20)
      close (20)
    enddo ! ilocal
#endif

    do ITER = 1, params%SCFSTEPS

      call resetTimer(iteration_timer)

      if (isMasterRank(my_mpi)) then
        call printDoubleLineSep(unit_number = 2)
        call outTime(isMasterRank(my_mpi),'started at ..........', getElapsedTime(program_timer),ITER)
        call printDoubleLineSep(unit_number = 2)
      endif ! master

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
          if (params%use_semicore == 1) then
            call RHOCORE_wrapper(emesh%EBOTSEMI, params%NSRA, atomdata)
          else
            call RHOCORE_wrapper(emesh%E1, params%NSRA, atomdata)
          endif
        enddo ! ilocal
        !!!$omp end parallel do
      endif ! in master group

! LDA+U ! TODO: doesn't work for num_local_atoms > 1
      if (params%LDAU) then
        ! For now only 1 atom per process is supported (1 local atom)
        CHECKASSERT(num_local_atoms == 1)

        atomdata      => getAtomData(calc_data, 1)
        ldau_data     => getLDAUData(calc_data, 1)
        I1 = getAtomIndexOfLocal(calc_data, 1)

        ldau_data%EREFLDAU = emesh%EFERMI
        ldau_data%EREFLDAU = 0.48    ! ???

#define mesh atomdata%mesh_ptr
        call LDAUINIT(I1,ITER,params%NSRA,ldau_data%NLDAU,ldau_data%LLDAU, &
                      ldau_data%ULDAU,ldau_data%JLDAU,ldau_data%EREFLDAU, &
                      atomdata%potential%VISP,ldau_data%NSPIND,mesh%R,mesh%DRDI, &
                      atomdata%Z_nuclear,mesh%IPAN,mesh%IRCUT, &
                      ldau_data%PHILDAU,ldau_data%UMLDAU,ldau_data%WMLDAU, &
                      ldau_data%lmaxd, mesh%irmd, mesh%ipand)
#undef mesh                      

      endif ! ldau
! LDA+U

      call outTime(isMasterRank(my_mpi),'initialized .........', getElapsedTime(program_timer),ITER)

      ! Scattering calculations - that is what KKR is all about
      ! output: (some contained as references in calc_data)
      ! ebalance_handler, kkr (!), jij_data, ldau_data
      call energyLoop(iter, calc_data, emesh, params, dims, ebalance_handler, my_mpi, arrays)

      call outTime(isMasterRank(my_mpi),'G obtained ..........', getElapsedTime(program_timer),ITER)

      flag = 0
      if (isInMasterGroup(my_mpi)) then
        ! output: (some contained as references in calc_data)
        ! atomdata, densities, broyden, ldau_data,
        ! emesh (only correct for master)
        flag = processKKRresults(iter, calc_data, my_mpi, emesh, dims, params, arrays, program_timer)

      endif ! in master group

      if (isMasterRank(my_mpi)) then
        ! only MASTERRANK updates, other ranks get it broadcasted later
        ! (although other processes could update themselves)

        call update(emesh)

        ! write file 'energy_mesh'
        if (emesh%NPOL /= 0) emesh%EFERMI = emesh%E2  ! if not a DOS-calculation E2 coincides with Fermi-Energy

        call store(emesh)

        call printDoubleLineSep()
        call writeIterationTimings(ITER, getElapsedTime(program_timer), getElapsedTime(iteration_timer))

      endif ! is master

      if (is_abort_by_rank0(flag, getMyActiveCommunicator(my_mpi))) exit

      call broadcast(emesh, getMyActiveCommunicator(my_mpi))

    enddo ! ITER ! SC ITERATION LOOP ITER=1, SCFSTEPS

    if (params%kforce == 1 .and. isInMasterGroup(my_mpi)) & ! write forces if requested, master-group only
      call output_forces(calc_data, 0, getMyAtomRank(my_mpi), getMySEcommunicator(my_mpi))

    if (isMasterRank(my_mpi)) close(2) ! time-info

    call destroy(ebalance_handler)
    call destroy(calc_data)
    call destroy(emesh)

  else  ! active Ranks
    !=====================================================================
    ! processors not fitting in NAEZ*LMPID do nothing ... and wait here
    !=====================================================================
  endif ! active Ranks

  CLOSELOG

  call destroy(arrays)
  call destroy(dims)
  call destroy(my_mpi) ! Free KKRnano mpi resources
  
endprogram ! KKRnano
