! KKRnano
! massive parallel KKR for nanoscaled systems

#include "DebugHelpers/test_macros.h"
#include "DebugHelpers/logging_macros.h"

program KKRnano
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)

  use Logging_mod, only:    !import no name here, just mention it for the module dependency 
  USE_LOGGING_MOD

  use KKRnanoParallel_mod, only: KKRnanoParallel, create, destroy

  use KKRnano_Comm_mod, only: setKKRnanoNumThreads, printKKRnanoInfo, communicatePotential

  use main2_aux_mod, only: printDoubleLineSep, is_abort_by_rank0, writeIterationTimings
  use EnergyMesh_mod, only: EnergyMesh, create, destroy, load, store, update, broadcast

  use BasisAtom_mod, only: BasisAtom
  use LDAUData_mod, only: LDAUData

  use TimerMpi_mod, only: TimerMpi, getElapsedTime, resetTimer, outTime
  use EBalanceHandler_mod, only: EBalanceHandler, createEBalanceHandler, initEBalanceHandler, setEqualDistribution, destroy

  use AtomicCore_mod, only: rhocore

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
  
  use BrillouinZoneMesh_mod, only: BrillouinZoneMesh, create, load, store, destroy
  implicit none

  type(CalculationData) :: calc_data
  type(TimerMpi) :: program_timer
  type(TimerMpi) :: iteration_timer
  type(EBalanceHandler) :: ebalance_handler

  integer :: ITER, I1, ila, num_local_atoms, flag, ios, ilen
  double precision  :: ebot 

  type(KKRnanoParallel) :: mp

  type(EnergyMesh)  :: emesh
  type(Main2Arrays) :: arrays
  type(DimParams)   :: dims
  type(InputParams) :: params

  type(BasisAtom), pointer :: atomdata
  type(LDAUData), pointer  :: ldau_data
  
! type(BrillouinZoneMesh) :: kmesh(8)

  external :: MPI_Init
  character(len=16)              :: arg
  
  call MPI_Init(ios) ! --> needs to be called here, otherwise MPI_Abort and MPI_Wtime cannot be used during toolbox functionalities
  
  call get_command_argument(1, arg, ilen, ios)
  selectcase (arg)
  case ('--prepare', '-p')
    call main0(checkmode=0,voronano=0) ! call former kkr0.exe
    stop
  case ('--check', '-c')
    call main0(checkmode=1,voronano=0) ! former kkr0.exe without overwriting the binary files *.unf
    stop
  case ('--voronano')
    call main0(checkmode=0,voronano=1) ! former kkr0.exe without reading of potential and shapefunctions
  case ('--convert')
    call kkrvform()
    stop
  case ('--help', '-h')
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

  call create(mp, dims%num_atom_procs, dims%SMPID, dims%EMPID)
  call setKKRnanoNumThreads(dims%nthrds)
  call printKKRnanoInfo(mp, dims%nthrds)

  if (mp%isMasterRank) then
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
  

  if (mp%myWorldRank < 128) then ! max. 128 logfiles
    OPENLOG(mp%myWorldRank, 3)
  else
    OPENLOG(mp%myWorldRank, 0)
  endif

  call resetTimer(program_timer)
  if (mp%isMasterRank) open(2, file='time-info', form='formatted', action='write')

  call createMain2Arrays(arrays, dims)
  call readMain2Arrays(arrays, 'arrays.unf') ! every process does this!

  flag = readInputParamsFromFile(params, 'input.unf')
  ! done reading variables


  if (params%JIJ .and. dims%nspind /= 2) die_here("Jij calculation not possible for spin-unpolarized calc.")

  if (dims%LLY /= 0) die_here("Lloyds formula not supported in this version. Set lly=0")

  !=====================================================================
  ! processors not fitting in NAEZ*LMPID*SMPID*EMPID do nothing ... and wait after SC-ITER loop
  !=====================================================================

  if (mp%isActiveRank) then

    ! pre self-consistency preparations

    assert(dims%naez > 0) 
    selectcase (arg)
    case ('--voronano')
      params%voronano = 1
      call create(calc_data, dims, params, arrays, my_mpi)! calls 'createCalculationData'
      stop ! Voronoi work is done in 'create'
    case default
      ! do not do Voronoi work
      call create(calc_data, dims, params, arrays, my_mpi) ! calls 'createCalculationData'
    endselect ! arg
    
    num_local_atoms = getNumLocalAtoms(calc_data)

    call create(emesh, dims%iemxd) ! createEnergyMesh
    call load(emesh, filename='energy_mesh.0') ! every process does this!!!

    call outTime(mp%isMasterRank,'input files read.....', getElapsedTime(program_timer), 0)

#ifdef DEBUG_NO_ELECTROSTATICS
    warn(6, "preprocessor define has switched off electrostatics for debugging")
#else
    call prepareMadelung(calc_data, arrays)
#endif

    call outTime(mp%isMasterRank,'Madelung sums calc...', getElapsedTime(program_timer), 0)

    call createEBalanceHandler(ebalance_handler, emesh%ielast)
    call initEBalanceHandler(ebalance_handler, mp)
    call setEqualDistribution(ebalance_handler, (emesh%npnt123(1) == 0))

    ! here, we comunicate the rMTref values early, so we can compute tref and dtrefdE locally for each atom inside the reference cluster 
    do ila = 1, num_local_atoms
      atomdata => getAtomdata(calc_data, ila)
#ifdef DEBUG_NO_VINS
      atomdata%potential%VINS = 0.0   
#endif
      I1 = getAtomIndexOfLocal(calc_data, ila)
#ifdef PRINT_MTRADII
      write(num, '(A,I7.7)') "mtradii_out.",I1
      open(20, file=num, form='formatted')
      write(20,*) atomdata%rmtref
      write(20,*) atomdata%radius_muffin_tin
      endfile (20)
      close (20)
#endif
#ifdef USE_MTRADII
      write(num, '(A,I7.7)') "mtradii_in.",I1
      open(20, file=num, form='formatted')
      read(20,*) atomdata%rmtref
      read(20,*) atomdata%radius_muffin_tin
      endfile (20)
      close (20)
#endif

      call gatherrMTref_com(rMTref_local=calc_data%atomdata_a(:)%rMTref, rMTref=calc_data%kkr_a(ila)%rMTref(:), &
                            ref_cluster=calc_data%ref_cluster_a(ila), communicator=mp%mySEComm)
    enddo ! ila

    do ITER = 1, params%SCFSTEPS

      call resetTimer(iteration_timer)

      if (mp%isMasterRank) then
        call printDoubleLineSep(unit_number=2)
        call outTime(mp%isMasterRank,'started at ..........', getElapsedTime(program_timer),ITER)
        call printDoubleLineSep(unit_number=2)
      endif ! master

      WRITELOG(2, *) "Iteration atom-rank ", ITER, mp%myAtomRank

      ! New: instead of reading potential every time, communicate it
      ! between energy and spin processes of same atom
      do ila = 1, num_local_atoms
        atomdata => getAtomdata(calc_data, ila)
        call communicatePotential(mp, atomdata%potential%VISP, atomdata%potential%VINS, atomdata%core%ECORE)
      enddo ! ila

      ! Core relaxation - only mastergroup needs results
      if (mp%isInMasterGroup) then
        ! Not threadsafe: intcor, intin, intout have a save statement
        ebot = emesh%E1; if (any(params%npntsemi > 0)) ebot = emesh%EBOTSEMI
        !!!$omp parallel do private(ila, atomdata)
        do ila = 1, num_local_atoms
          atomdata => getAtomdata(calc_data, ila)
!         call RHOCORE_wrapper(ebot, params%NSRA, atomdata)
#define mesh atomdata%mesh_ptr
          atomdata%core%qc_corecharge = rhocore(ebot, params%nsra, atomdata%nspin, atomdata%atom_index, &  ! atom_index is used only for debugging output
                        mesh%drdi, mesh%r, atomdata%potential%visp(:,:), &
                        mesh%a, mesh%b, atomdata%z_nuclear, &
                        mesh%ircut, atomdata%core%rhocat, &
                        atomdata%core%ecore(:,:), atomdata%core%ncore(:), atomdata%core%lcore(:,:), &
                        mesh%irmd)
        enddo ! ila
        !!!$omp end parallel do
      endif ! in master group

! LDA+U ! TODO: does not work for num_local_atoms > 1
      if (params%LDAU) then
        ! For now only 1 atom per process is supported (1 local atom)
        CHECKASSERT(num_local_atoms == 1)

        atomdata      => getAtomData(calc_data, 1)
        ldau_data     => getLDAUData(calc_data, 1)
        I1 = getAtomIndexOfLocal(calc_data, 1)

        ldau_data%EREFLDAU = emesh%EFERMI
        ldau_data%EREFLDAU = 0.48    ! ???

        call LDAUINIT(I1,ITER,params%NSRA,ldau_data%NLDAU,ldau_data%LLDAU, &
                      ldau_data%ULDAU,ldau_data%JLDAU,ldau_data%EREFLDAU, &
                      atomdata%potential%VISP,ldau_data%NSPIND,mesh%R,mesh%DRDI, &
                      atomdata%Z_nuclear,mesh%IPAN,mesh%IRCUT, &
                      ldau_data%PHILDAU,ldau_data%UMLDAU,ldau_data%WMLDAU, &
                      ldau_data%lmaxd, mesh%irmd, mesh%ipand)
#undef mesh

      endif ! ldau
! LDA+U

      call outTime(mp%isMasterRank,'initialized .........', getElapsedTime(program_timer),ITER)

      ! Scattering calculations - that is what KKR is all about
      ! output: (some contained as references in calc_data)
      ! ebalance_handler, kkr (!), jij_data, ldau_data
      call energyLoop(iter, calc_data, emesh, params, dims, ebalance_handler, mp, arrays)

      call outTime(mp%isMasterRank,'G obtained ..........', getElapsedTime(program_timer),ITER)

      flag = 0
      if (mp%isInMasterGroup) then
        ! output: (some contained as references in calc_data)
        ! atomdata, densities, broyden, ldau_data,
        ! emesh (only correct for master)
        flag = processKKRresults(iter, calc_data, mp, emesh, dims, params, arrays, program_timer)

      endif ! in master group

      if (mp%isMasterRank) then
        ! only MASTERRANK updates, other ranks get it broadcasted later
        ! (although other processes could update themselves)

        call update(emesh)

        ! write file 'energy_mesh'
        if (emesh%NPOL /= 0) emesh%EFERMI = emesh%E2  ! if not a DOS-calculation E2 coincides with Fermi-Energy

        call store(emesh, filename='energy_mesh')

        call printDoubleLineSep()
        call writeIterationTimings(ITER, getElapsedTime(program_timer), getElapsedTime(iteration_timer))

      endif ! is master

      if (is_abort_by_rank0(flag, mp%myActiveComm)) exit

      call broadcast(emesh, mp%myActiveComm)

    enddo ! ITER ! SC ITERATION LOOP ITER=1, SCFSTEPS

    if (params%kforce == 1 .and. mp%isInMasterGroup) & ! write forces if requested, master-group only
      call output_forces(calc_data, 0, mp%myAtomRank, mp%mySEComm)

    if (mp%isMasterRank) close(2) ! time-info

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
  call destroy(mp) ! Free KKRnano mpi resources

  contains
  
    !------------------------------------------------------------------------------
    !> Gather all rMTref values of the reference cluster.
    !> @param rMTref_local all locally determined rMTref value
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

endprogram ! KKRnano
