!--------------------------------------------------------------------------------
! Copyright (c) 2018 Forschungszentrum Juelich GmbH, Juelich, Germany
! This file is part of KKRnano and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!!! main program KKRnano
!-------------------------------------------------------------------------------
!> Summary: massively parallel density functional theory code KKRnano
!>          for nanoscaled systems based on Korringa-Kohn-Rostoker 
!>          multiple scattering theory
!> Author: Marcel Bornemann, Paul F Baumeister, Elias Rabel, Alexander Thiess, 
!>         Rudolf Zeller, Roman Kovacik, et al.
!> Category: KKRnano
!-------------------------------------------------------------------------------
program KKRnano
#include "DebugHelpers/logging_macros.h"
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  use Warnings_mod, only: show_warning_lines

  use Logging_mod, only:    ! import no name, just mention it for the module dependency 
  USE_LOGGING_MOD

  use KKRnanoParallel_mod, only: KKRnanoParallel, create, destroy

  use KKRnano_Comm_mod, only: setKKRnanoNumThreads, printKKRnanoInfo, communicatePotential, communicateNoncoBfields

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
  use ProcessKKRresults_mod, only: processKKRresults, output_forces

  use CalculationData_mod, only: CalculationData, create, prepareMadelung, destroy
  use CalculationData_mod, only: getAtomData, getLDAUData

  use KKRzero_mod, only: main0
  use PotentialConverter_mod, only: kkrvform
#ifdef BENCHMARK_tfQMR
  use tfQMR_mod, only: benchmark_tfQMR
#endif

  implicit none

  type(CalculationData) :: calc_data
  type(TimerMpi) :: program_timer, iteration_timer
  type(EBalanceHandler) :: ebalance_handler

  integer :: iter, global_atom_id, ila, ios, ilen, voronano
  double precision  :: ebot
  double precision :: program_start_time, program_stop_time

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

  call createTimer(program_timer)

  call MPI_Init(ios) ! --> needs to be called now, otherwise MPI_Abort and MPI_Wtime cannot be used during toolbox functionalities

  call startTimer(program_timer, program_start_time)

#ifdef  SUPERCELL_ELECTROSTATICS
  warn(6, "Fix symmetry with SUPERCELL_ELECTROSTATICS ="+SUPERCELL_ELECTROSTATICS)
#endif

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
#ifdef BENCHMARK_tfQMR
    '    --benchmark           Benchmarking of the tfQMR solver', &
#endif
    '    --version             Version and license information', &
    ''
    stop
  case ('--benchmark')
#ifdef BENCHMARK_tfQMR
    call benchmark_tfQMR()
#else
    die_here('For benchmarking of the tfQMR solver, please compile with -D BENCHMARK_tfQMR')
#endif
  case ('--version')
    write(*,'(A)') 'KKRnano', &
    'Copyright (C) 2018 Forschungszentrum Juelich, Juelich, Germany', &
    'This is free software (MIT license). Please see the LICENSE file for more detail.'
    stop
  case default
    ! start the former kkr2.exe    
  endselect ! arg

  ! from now on the former kkr2.exe actions are performed

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

  if (mp%isMasterRank) then
    open(2, file='time-info', action='write', iostat=ios)
    if (ios /= 0) warn(6, "unable to create time-info file!") ! but the output will be redirected to fort.2
    write(2,'(79("="))') ! double separator line in time-info file
  endif ! master
  call outTime(mp%isMasterRank, 'initialization', getTime(program_timer), 0)

  call create(arrays, dims%lmmaxd, dims%lmmaxd_noco, dims%naez, dims%kpoibz, dims%maxmshd)
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

    call outTime(mp%isMasterRank, 'input files read', getTime(program_timer), 0)

#ifdef DEBUG_NO_ELECTROSTATICS
    warn(6, "preprocessor define has switched off electrostatics for debugging")
#else
    call prepareMadelung(calc_data, arrays%rbasis)
#endif

    call outTime(mp%isMasterRank, 'Madelung sums calc', getTime(program_timer), 0)

    call create(ebalance_handler, emesh%ielast)
    call init(ebalance_handler, mp)
    call setEqualDistribution(ebalance_handler, (emesh%npnt123(1) == 0))

    call outTime(mp%isMasterRank, 'Energy balancer initialized', getTime(program_timer), 0)

    do ila = 1, calc_data%num_local_atoms
      atomdata => getAtomdata(calc_data, ila)
#ifdef DEBUG_NO_VINS
      atomdata%potential%VINS = 0.0   
#endif
      global_atom_id = calc_data%atom_ids(ila) ! get global atom_id from local index
#ifdef PRINT_MTRADII
      write(unit=num, fmt="(A,I7.7)") "mtradii_out.",global_atom_id
      open(20, file=num, action='write')
      write(20,*) atomdata%rMTref
      write(20,*) atomdata%radius_muffin_tin
      endfile(20)
      close(20)
#endif
#ifdef USE_MTRADII
      write(unit=num, fmt="(A,I7.7)") "mtradii_in.",global_atom_id
      open(20, file=num, action='read', status='old')
      read(20,*) atomdata%rMTref
      read(20,*) atomdata%radius_muffin_tin
      close(20)
#endif
    enddo ! ila

    call outTime(mp%isMasterRank, 'Muffin-Tin radii scattered', getTime(program_timer), 0)

    call createTimer(iteration_timer)
    if (mp%isMasterRank) write(2,'(79("="))') ! double separator line in time-info file

    ! start self-consistency loop   
    do iter = 1, params%SCFSTEPS
      call startTimer(iteration_timer)

      WRITELOG(2, *) "Iteration atom-rank ", iter, mp%myAtomRank ! write logg message into time-info file

      ! New: instead of reading potential every time, communicate it
      ! between energy and spin processes of same atom
      do ila = 1, calc_data%num_local_atoms
        atomdata => getAtomdata(calc_data, ila)
        call communicatePotential(mp, atomdata%potential%VISP, atomdata%potential%VINS, atomdata%core%ECORE)
      enddo ! ila
      ! If noncollinear magnetic fields are used, communicate them analogously to the potentials
      if (params%noncobfield) then
        do ila = 1, calc_data%num_local_atoms
          call communicateNoncoBfields(mp, calc_data%bfields(ila))
        end do
      end if

      ! Core relaxation - only mastergroup needs results
      if (mp%isInMasterGroup.and.params%npol /= 0) then
        ! Not threadsafe: intcor, intin, intout have a save statement
        ebot = emesh%E1; if (any(params%npntsemi > 0)) ebot = emesh%EBOTSEMI
        !!!$omp parallel do private(ila, atomdata)
        do ila = 1, calc_data%num_local_atoms
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
        CHECKASSERT(calc_data%num_local_atoms == 1)

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

      call outTime(mp%isMasterRank, 'iteration initialized', getTime(program_timer), iter)

      ! Scattering calculations - that is what KKR is all about
      ! output: (some contained as references in calc_data)
      ! ebalance_handler, kkr (!), jij_data, ldau_data
      call energyLoop(iter, calc_data, emesh, params, dims, ebalance_handler, mp, arrays)

      call outTime(mp%isMasterRank, 'energyLoop', getTime(program_timer), iter)
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
      call outTime(.true. , 'SCF iteration', getTime(iteration_timer), iter)
      if (mp%isMasterRank) write(2,'(79("="))') ! double separator line in time-info file
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

#ifdef  SUPERCELL_ELECTROSTATICS
  warn(6, "Fixed symmetry with SUPERCELL_ELECTROSTATICS ="+SUPERCELL_ELECTROSTATICS)
#endif

  if (mp%isMasterRank) ios = show_warning_lines(unit=6)

  call stopTimer(program_timer, program_stop_time)
  if (mp%isMasterRank) then
    call outTimeStats(iteration_timer, 'SCF stats:')
    write(2,'(79("="))') ! double separator line in time-info file
    call outTime(.true., 'KKRnano', program_stop_time - program_start_time, iter - 1)
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
