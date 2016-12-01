! KKRnano
! massive parallel KKR for nanoscaled systems

! #include "DebugHelpers/test_macros.h"
#include "DebugHelpers/logging_macros.h"

program KKRnano
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  use Warnings_mod, only: show_warning_lines

  use Logging_mod, only:    !import no name here, just mention it for the module dependency 
  USE_LOGGING_MOD

  use KKRnanoParallel_mod, only: KKRnanoParallel, create, destroy

  use KKRnano_Comm_mod, only: setKKRnanoNumThreads, printKKRnanoInfo, communicatePotential

  use main2_aux_mod, only: is_abort_by_rank0
  use EnergyMesh_mod, only: EnergyMesh, create, destroy, load, store, update, broadcast

  use BasisAtom_mod, only: BasisAtom
  use LDAUData_mod, only: LDAUData

  use TimerMpi_mod, only: TimerMpi, getTime, createTimer, outTime, startTimer, stopTimer, outTimeStats
  use EBalanceHandler_mod, only: EBalanceHandler, create, init, setEqualDistribution, destroy

  use AtomicCore_mod, only: rhocore

  use DimParams_mod, only: DimParams, load, destroy
  use InputParams_mod, only: InputParams, load
  use Main2Arrays_mod, only: Main2Arrays, create, load, destroy

  use ScatteringCalculation_mod, only: energyLoop
  use ScatteringCalculation_mod, only: gatherrMTref_com
  use ProcessKKRresults_mod, only: processKKRresults, output_forces

  use CalculationData_mod, only: CalculationData, create, prepareMadelung, destroy
  use CalculationData_mod, only: getAtomData, getLDAUData
  
  use KKRzero_mod, only: main0
  use PotentialConverter_mod, only: kkrvform
  
  implicit none

  type(CalculationData) :: calc_data
  type(TimerMpi) :: program_timer, iteration_timer
  type(EBalanceHandler) :: ebalance_handler

  integer :: iter, global_atom_id, ila, num_local_atoms, ios, ilen, voronano
  double precision  :: ebot 

  type(KKRnanoParallel) :: mp

  type(EnergyMesh)  :: emesh
  type(Main2Arrays) :: arrays
  type(DimParams)   :: dims
  type(InputParams) :: params

  type(BasisAtom), pointer :: atomdata
  type(LDAUData), pointer  :: ldau_data

  external :: MPI_Init
  character(len=16) :: arg
#ifdef PRINT_MTRADII
  character(len=40) :: num 
#else  
#ifdef USE_MTRADII
  character(len=40) :: num 
#endif
#endif
 
  call MPI_Init(ios) ! --> needs to be called here, otherwise MPI_Abort and MPI_Wtime cannot be used during toolbox functionalities
  
  voronano = 0
  call get_command_argument(1, arg, ilen, ios)
  selectcase (arg)
  case ('--prepare', '-p')
    call main0(checkmode=0, voronano=voronano) ! call former kkr0.exe
    stop
  case ('--check', '-c')
    call main0(checkmode=1, voronano=voronano) ! former kkr0.exe without overwriting the binary files bin.*
    stop
  case ('--voronano')
    voronano = 1
    call main0(checkmode=0, voronano=voronano) ! former kkr0.exe without reading of potential and shapefunctions
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
    '    --voronano            Produces Voronoi output necessary to start a calculation with KKRnano', &
    ''
    stop
  case default
    ! start the former kkr2.exe    
  endselect ! arg

  ! from here the former kkr2.exe actions are performed
  
  call load(dims, 'bin.dims') ! read dimension parameters from file 'bin.dims'

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
      write(*,*) 'PRINT_MTRADII enabled; rMTref (0.995 * max_muffin_tin by default) and ', &
      'radius_muffin_tin (0.98 * max_muffin_tin by default) are printed consecutively to files mtradii_out'
#endif

#ifdef USE_MTRADII
      write(*,*) 'USE_MTRADII enabled; rMTref (0.995 * max_muffin_tin by default) and ', &
        'radius_muffin_tin (0.98 * max_muffin_tin by default) are taken from files mtradii_in'
#endif
  endif ! is master
  
  OPENLOG(mp%myWorldRank, 3*theta(mp%myWorldRank < 128)) ! max. 128 logfiles

  call createTimer(program_timer)
  call startTimer(program_timer)
  if (mp%isMasterRank) then
    open(2, file='time-info', action='write', iostat=ios)
    if (ios /= 0) warn(6, "unable to create time-info file!") ! but the output will be redirected to fort.2
    write(2,'(79("="))') ! double separator line in time-info file
  endif ! master

  call create(arrays, dims%lmmaxd, dims%naez, dims%kpoibz, dims%maxmshd)
  call load(arrays, 'bin.arrays') ! every process does this!

  call load(params, 'bin.input', ios=ios)
  ! done reading variables

  if (params%JIJ .and. dims%nspind /= 2) die_here("Jij calculation not possible for spin-unpolarized calc.")

  !=====================================================================
  ! processors not fitting in NAEZ*LMPID*SMPID*EMPID do nothing ... and wait after SC-iteration loop
  !=====================================================================

  if (mp%isActiveRank) then

    ! pre self-consistency preparations

    assert(dims%naez > 0) 
    
    call create(emesh, dims%iemxd) ! createEnergyMesh
    call load(emesh, filename='bin.energy_mesh.0') ! every process does this!!!

    call create(calc_data, dims, params, arrays, mp, emesh%kmesh, voronano)
    
    if (voronano == 1) then
      ios = show_warning_lines(unit=6)
      stop ! Voronoi work is done in 'create'
    endif

    num_local_atoms = calc_data%num_local_atoms

    call outTime(mp%isMasterRank, 'input files read ...............', getTime(program_timer), 0)

#ifdef DEBUG_NO_ELECTROSTATICS
    warn(6, "preprocessor define has switched off electrostatics for debugging")
#else
    call prepareMadelung(calc_data, arrays)
#endif

    call outTime(mp%isMasterRank, 'Madelung sums calc .............', getTime(program_timer), 0)

    call create(ebalance_handler, emesh%ielast)
    call init(ebalance_handler, mp)
    call setEqualDistribution(ebalance_handler, (emesh%npnt123(1) == 0))

    ! here, we communicate the rMTref values early, so we could compute tref and dtrefdE locally for each atom inside the reference cluster 
    do ila = 1, num_local_atoms
      atomdata => getAtomdata(calc_data, ila)
#ifdef DEBUG_NO_VINS
      atomdata%potential%VINS = 0.0   
#endif
      global_atom_id = calc_data%atom_ids(ila) ! get global atom_id from local index
#ifdef PRINT_MTRADII
      write(unit=num, '(A,I7.7)') "mtradii_out.",global_atom_id
      open(20, file=num)
      write(20,*) atomdata%rMTref
      write(20,*) atomdata%radius_muffin_tin
      endfile(20)
      close(20)
#endif
#ifdef USE_MTRADII
      write(unit=num, '(A,I7.7)') "mtradii_in.",global_atom_id
      open(20, file=num)
      read(20,*) atomdata%rMTref
      read(20,*) atomdata%radius_muffin_tin
      endfile(20)
      close(20)
#endif
      ! this only works if all process enter the following subroutine the same number of times, so num_local_atoms may not differ among the processes
      call gatherrMTref_com(rMTref_local=calc_data%atomdata_a(:)%rMTref, rMTref=calc_data%kkr_a(ila)%rMTref(:), &
                            ref_cluster=calc_data%ref_cluster_a(ila), communicator=mp%mySEComm)
    enddo ! ila
    
    call createTimer(iteration_timer)
    if (mp%isMasterRank) write(2,'(79("="))') ! double separator line in time-info file
    
    ! start self-consistency loop   
    do iter = 1, params%SCFSTEPS

      call startTimer(iteration_timer)

      if (mp%isMasterRank) call outTime(mp%isMasterRank, 'start iteration ................', getTime(program_timer), iter)

      WRITELOG(2, *) "Iteration atom-rank ", iter, mp%myAtomRank ! write logg message into time-info file

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

        atomdata  => getAtomData(calc_data, 1)
        ldau_data => getLDAUData(calc_data, 1)
        global_atom_id = calc_data%atom_ids(1) ! get global atom_id from local index

        ldau_data%EREFLDAU = emesh%EFERMI
        ldau_data%EREFLDAU = 0.48 ! ???

        call LDAUINIT(global_atom_id,iter,params%NSRA,ldau_data%NLDAU,ldau_data%LLDAU, &
                      ldau_data%ULDAU,ldau_data%JLDAU,ldau_data%EREFLDAU, &
                      atomdata%potential%VISP,ldau_data%NSPIND,mesh%R,mesh%DRDI, &
                      atomdata%Z_nuclear,mesh%IPAN,mesh%IRCUT, &
                      ldau_data%PHILDAU,ldau_data%UMLDAU,ldau_data%WMLDAU, &
                      ldau_data%lmaxd, mesh%irmd, mesh%ipand)
#undef mesh

      endif ! ldau
! LDA+U

      call outTime(mp%isMasterRank, 'initialized ....................', getTime(program_timer), iter)

      ! Scattering calculations - that is what KKR is all about
      ! output: (some contained as references in calc_data)
      ! ebalance_handler, kkr (!), jij_data, ldau_data
      call energyLoop(iter, calc_data, emesh, params, dims, ebalance_handler, mp, arrays)

      call outTime(mp%isMasterRank, 'G obtained .....................', getTime(program_timer), iter)
      ios = 0
      if (mp%isInMasterGroup) then
        ! output: (some contained as references in calc_data)
        ! atomdata, densities, broyden, ldau_data,
        ! emesh (only correct for master)
        ios = processKKRresults(iter, calc_data, mp, emesh, dims, params, arrays, program_timer)

      endif ! in master group

      if (mp%isMasterRank) then
        ! only MASTERRANK updates, other ranks get it broadcasted later
        ! (although other processes could update themselves)

        call update(emesh)

        ! write file 'energy_mesh'
        if (emesh%NPOL /= 0) emesh%EFERMI = emesh%E2 ! if not a DOS-calculation E2 coincides with Fermi-Energy

        call store(emesh, filename='bin.energy_mesh')

        write(*,'(79("="))') ! double separator line

      endif ! master

      if (is_abort_by_rank0(ios, mp%myActiveComm)) exit

      call broadcast(emesh, mp%myActiveComm)

      if (mp%isMasterRank) &
      call outTime(.true. , 'SCF iteration           took', getTime(iteration_timer), iter)
      if (mp%isMasterRank) write(2,'(79("="))') ! double separator line in time-info file
      call stopTimer(iteration_timer)

    enddo ! iter ! SC ITERATION LOOP iter=1, SCFSTEPS
    
    if (params%kforce == 1 .and. mp%isInMasterGroup) & ! write forces if requested, master-group only
      call output_forces(calc_data, 0, mp%myAtomRank, mp%mySEComm)
      
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
  
  if (mp%isMasterRank) ios = show_warning_lines(unit=6)

  call stopTimer(program_timer)
  if (mp%isMasterRank) then
    call outTimeStats(iteration_timer, 'SCF stats:')
    write(2,'(79("="))') ! double separator line in time-info file
    call outTime(.true., 'KKRnano run             took', getTime(program_timer), iter - 1)
    write(2,'(79("="))') ! double separator line in time-info file
    close(2) ! time-info
  endif ! master
  
  call destroy(mp) ! Free KKRnano mpi resources
  
  contains

    integer function theta(condition) result(oneiftrue)
      logical, intent(in) :: condition
      oneiftrue = 0 ; if (condition) oneiftrue = 1
    endfunction ! theta
    
endprogram ! KKRnano
