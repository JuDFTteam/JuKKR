! -------------------------------------------------------------------------------
! PROGRAM: kkrcode
!> @brief Main program for the JM-KKR
!> @details The JM-KKR code is a Density Functional Theory software package,
!> based on the Green function Korringa-Kohn-Rostocker approach.
!> The package allows the calculation of 3D and 2D systems, as well as the
!> determination of the needed parameters for the calculation of impurities and
!> nanoclusters (KKRImp code). The code has also been modified to produce the
!> needed information for the treatment of TD-DFT calculations based in the linear
!> response approach (KKRSusc code).
!> @author Philipp Rüssmann, Bernd Zimmermann, Phivos Mavropoulos, R. Zeller, and many others ...
!! @note
!> - Jonathan Chico Jan. 2018: Removed inc.p dependencies and rewrote to Fortran90
! -------------------------------------------------------------------------------
program kkrcode

#ifdef CPP_MPI
  use :: mpi
  use :: mod_mympi, only: mympi_init, myrank, nranks, master, find_dims_2d, distribute_linear_on_tasks, create_newcomms_group_ie, mpiatom, mpiadapt, check_communication_pattern
  use :: mod_save_wavefun, only: t_wavefunctions, bcast_params_savewf
  use :: godfrin, only: t_godfrin, bcast_params_godfrin ! GODFRIN Flaviano
  use :: mod_wunfiles, only: bcast_t_params_scalars, bcast_t_params_arrays
  use :: mod_types, only: bcast_t_lly_1, bcast_t_inc_tgmat, save_t_mpi_c_grid
  use :: mod_md5sums, only: mympi_bcast_md5sums
#else
  use :: mod_mympi, only: mympi_init, myrank, nranks, master
  use :: mod_save_wavefun, only: t_wavefunctions
#endif
  use :: mod_constants, only: czero, nsymaxd
  use :: mod_profiling, only: memocc
  use :: mod_types, only: t_inc, t_lloyd, t_cpa, t_mpi_c_grid, t_tgmat
  use :: mod_timing, only: timing_start, timing_stop, timing_init, timings_1a, timings_1b, load_imbalance, print_time_and_date
  use :: memoryhandling, only: allocate_cell, allocate_cpa, allocate_soc, allocate_ldau, allocate_magnetization, allocate_potential, &
    allocate_energies, allocate_relativistic, allocate_clusters, allocate_expansion, allocate_mesh, allocate_pannels, allocate_misc, &
    allocate_green, allocate_ldau_potential, allocate_rel_transformations, allocate_semi_inf_host
  use :: mod_version_info, only: version_print_header, construct_serialnr
  use :: mod_wunfiles, only: t_params, init_t_params
  use :: mod_main1a, only: main1a
  use :: mod_main1b, only: main1b
  use :: mod_main1c, only: main1c
  use :: mod_main2, only: main2
  ! array dimensions
  use :: global_variables, only: iemxd, ipand, irid, irmind, irmd, krel, lassld, lm2d, lmaxd, lmmaxd, lmpotd, lmxspd, mmaxd, naclsd, &
    naezd, natomimpd, natypd, ncelld, nchebd, ncleb, nclsd, nembd, ipand, irid, irmd, irmind, krel, nembd1, nfund, ngshd, nofgij, nrd, &
    nprincd, nrefd, nsheld, nspotd, nspind, nspindd, npotd, ntotd
  ! stuff defined in main0 already
  use :: mod_main0, only: main0, a, atom, atomimp, b, btrel, cleb, cls, cmomhost, conc, crel, cscl, dez, drdi, drdirel, dror, drotq, &
    dsymll, dsymll1, ecore, erefldau, ez, ezoa, fpradius, gsh, hostimp, icheck, icleb, icpa, ifunm, ifunm1, ijtabcalc, ijtabcalc_i, &
    ijtabsh, ijtabsym, ilm_map, imaxsh, imt, inipol, iofgij, ipan, ipan_intervall, iqat, iqcalc, irc, ircut, irm, irmin, irns, irrel, &
    irshift, irws, ish, jsh, ititle, itldau, ixipol, jeff, jend, jofgij, jwsrel, kaoez, kfg, kmesh, lcore, lefttinvll, llmsp, lmax, &
    lmpot, lmsp, lmsp1, lmxc, loflm, lopt, mtfac, nacls, naez, natyp, ncheb, ncore, nemb, nfu, noq, npan_eq_at, npan_log_at, nref, &
    nshell, ntcell, nsh1, nsh2, nrrel, npan_tot, phildau, qmgam, qmgamtab, qmphi, qmphitab, qmtet, qmtettab, r2drdirel, ratom, rbasis, &
    rc, rcls, rclsimp, refpot, righttinvll, rmesh, rmtnew, rmtrefat, rnew, rpan_intervall, rr, rrel, rrot, rs, rws, s, rmt, rmtref, &
    socscale, socscl, srrel, thesme, thetas, thetasnew, tleft, tright, rmrel, uldau, vins, visp, vref, vtrel, wez, wg, wldau, &
    yrg, zat, zrel, ueff


  implicit none

  integer :: i_stat, i_all
#ifdef CPP_MPI
  integer :: mympi_comm_at, mympi_comm_ie, nranks_at, myrank_at
  integer :: ierr, myrank_ie, nranks_ie, nranks_atcomm, myrank_atcomm
  integer, dimension (2) :: dims
#endif
  character (len=3) :: ctemp                               ! name for output file

  ! needed to use test('xxxxxxxx'):
  logical :: test, opt
  external :: test, opt

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! initialize MPI >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef CPP_MPI
  ! initialize MPI
  call mpi_init(ierr)
#endif
  ! set variables master, myrank and nranks for serial (myrank=master=0, nranks=1) as well as parallel execution
  call mympi_init()
  ! save myrank in ctemp, needed to open output unit 1337
  write (ctemp, '(I03.3)') myrank
  ! find serial number that is printed to files
  call construct_serialnr()
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< initialize MPI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! start KKR with main0 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! do all this only on master (read in parameters etc. in main0)
  if (myrank==master) then

    ! initialize timing
    call timing_init(myrank)

    ! start KKR program, first do main0, where everything is read in and initialized
    call timing_start('main0')

    ! open output files
    write (*, *) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write (*, *) '!!! Most output written to output.myrank.txt files !!!'
    write (*, *) '!!! please check these files as well               !!!'
    write (*, *) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    open (1337, file='output.'//trim(ctemp)//'.txt')
    ! here version_print_header needs the print_always flag because test options are not read in yet
    call version_print_header(1337, print_always=.true.)

    ! default value on master (needed for writeout in main0)
    t_inc%i_write = 1

    ! now memocc can be started (needs t_inc%i_write to be set)
    call memocc(0, 0, 'count', 'start')

    ! run main0 only serially, then communicate everything that is in here
    call main0()
    call timing_stop('main0')

  end if                           ! myrank ==master
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< start KKR with main0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! full Dirac works only in serial at the moment
  if (krel>0 .and. nranks>1) then
    if (myrank==master) write (*, *) 'Dirac solver does not work with MPI. Please run serially'
#ifdef CPP_MPI
    call mpi_barrier(mpi_comm_world, ierr)
    call mpi_finalize(ierr)
#endif
    stop
  end if


  ! without MPI (serial or openMP) something goes wrong if if files are not written out
  ! this seems to be only the case with the old solver
#ifndef CPP_MPI
  ! if(.not.t_inc%NEWSOSOL) then
  ! t_tgmat%tmat_to_file = .true.
  ! t_tgmat%gmat_to_file = .true.
  ! t_tgmat%gref_to_file = .true.
  ! end if
#endif

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! distribute stuff from main0 >>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef CPP_MPI
  ! now communicate type t_inc and t_tgmat switches (this has an implicit barrier, so that all other processes wait for master to finish with main0)
  if (myrank==master) call timing_start('MPI 1')
  call bcast_t_inc_tgmat(t_inc, t_tgmat, t_cpa, master)

  ! also communicate logicals from t_lloyd
  call bcast_t_lly_1(t_lloyd, master)

  ! communicate parameters that were written in wunfiles into t_params
  call bcast_t_params_scalars(t_params)
  if (myrank/=master) call init_t_params(t_params)
  call bcast_t_params_arrays(t_params)

  ! communicate global variables (previously in inc.p)
  call bcast_global_variables()
#endif


  ! set verbosity levels for each rank set i_write and i_time to 0 or 1,2, depending on verbosity level specified in inputcard and rank
  ! take care of different ranks here for verbose0 and timings1,2
  ! convention: value of 0 mean do not write, value 1 write and reset file, value 2 write everthing
  if (t_inc%i_write==0 .and. myrank==master) t_inc%i_write = 1 ! only written to master, reset file after each iteration, output.init.txt for main0 output
  if (t_inc%i_time==0 .and. myrank==master) t_inc%i_time = 1 ! only master writes timings of current iteration
  if (t_inc%i_time==1) then
    t_inc%i_time = 0               ! all ranks do not write timings but only the master
    if (myrank==master) t_inc%i_time = 2 ! master write all timings of all iterations
  end if                           ! t_inc%i_time==1


#ifdef CPP_MPI
  if (myrank/=master) then

    ! initialize memocc also for the other ranks after t_inc%i_write has been set
    call memocc(0, 0, 'count', 'start')

    call allocate_cell(1, naezd, nembd, natypd, cls, imt, irws, irns, ntcell, refpot, kfg, kaoez, rmt, zat, rws, mtfac, rmtref, rmtrefat, rmtnew, rbasis, lmxc)
    call allocate_semi_inf_host(1, nembd, tleft, tright)
    call allocate_cpa(1, naezd, natypd, noq, icpa, iqat, hostimp, conc)
    call allocate_soc(1, krel, natypd, lmaxd, socscale, cscl, socscl)
    call allocate_ldau(1, natypd, lopt, ueff, jeff, erefldau)
    call allocate_magnetization(1, naezd, natypd, lmmaxd, inipol, ixipol, qmtet, qmphi, drotq)
    call allocate_potential(1, irmd, natypd, npotd, ipand, nfund, lmxspd, lmpotd, irmind, nspotd, nfu, irc, ncore, irmin, lmsp, lmsp1, ircut, lcore, llmsp, ititle, fpradius, visp, &
      ecore, vins)
    call allocate_ldau_potential(1, irmd, natypd, mmaxd, nspind, itldau, wldau, uldau, phildau)
    call allocate_energies(1, iemxd, ez, dez, wez)
    call allocate_relativistic(1, krel, irmd, naezd, natypd, zrel, jwsrel, irshift, vtrel, btrel, rmrel, drdirel, r2drdirel, qmgam, qmgamtab, qmphitab, qmtettab)
    call allocate_rel_transformations(1, lmmaxd, nrrel, irrel, rc, crel, rrel, srrel)
    call allocate_clusters(1, naezd, lmaxd, ncleb, nclsd, nembd1, nsheld, naclsd, lmpotd, natomimpd, nsh1, nsh2, nacls, nshell, atomimp, atom, ezoa, icleb, jend, ratom, rclsimp, &
      cmomhost, rcls)
    call allocate_expansion(1, lm2d, irid, nfund, ntotd, ncleb, lassld, ncelld, nchebd, loflm, wg, cleb, yrg, thetas, thetasnew)
    call allocate_mesh(1, irmd, natypd, a, b, rmesh, drdi)
    call allocate_pannels(1, natypd, ntotd, ipan, npan_tot, npan_eq_at, npan_log_at, ipan_intervall, rpan_intervall)
    call allocate_misc(1, nrd, irmd, irid, lmaxd, naezd, natypd, nfund, nrefd, iemxd, ntotd, nsheld, lmmaxd, nembd1, nchebd, ncelld, lmxspd, nspindd, nsymaxd, nprincd, ifunm, &
      ifunm1, icheck, vref, s, rr, dror, rnew, rs, rrot, thesme, dsymll, dsymll1, lefttinvll, righttinvll)
    call allocate_green(1, naezd, iemxd, ngshd, nsheld, lmpotd, nofgij, ish, jsh, kmesh, imaxsh, iqcalc, iofgij, jofgij, ijtabsh, ijtabsym, ijtabcalc, ijtabcalc_i, ilm_map, gsh)

  end if                           ! ( myrank/=master )

  ! broadcast md5 sums (written to some output files (kkrflex_* etc.)
  call mympi_bcast_md5sums(t_params%ins, myrank, master)

  ! in case of deci-out run everything only serially to write decifile correctly. This will be fixed in a later version when we get rid of the decifile
  if (t_inc%deci_out .and. t_mpi_c_grid%nranks_at>1) then
    stop 'deci-out option chosen. Please run code serially in energy dimension!'
  end if

  ! call myMPI_distribute_ranks(MPIatom, MPIadapt, t_inc, timings_1a, timings_1b, load_imbalance, nranks, myrank, initial=1)
  ! create 2d matrix for processors so that more processors than energy points can be used.
  ! strategy here is to first parallelize the energy points ideally and then give all processors that are left to atom or k-point loop parallelization

  ! communicate logical MPIatom and MPIadapt which determine how the ranks are devided into groups
  call mpi_bcast(mpiatom, 1, mpi_logical, master, mpi_comm_world, ierr)
  if (ierr/=mpi_success) stop 'error broadcasting MPIatom in main_all'
  call mpi_bcast(mpiadapt, 1, mpi_integer, master, mpi_comm_world, ierr)
  if (ierr/=mpi_success) stop 'error broadcasting MPIadapt in main_all'

  ! allocate timing arrays
  if (mpiadapt>0) then
    allocate (timings_1a(t_inc%ielast,t_inc%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(timings_1a))*kind(timings_1a), 'timings_1a', 'main_all')
    allocate (timings_1b(t_inc%ielast), stat=i_stat)
    call memocc(i_stat, product(shape(timings_1b))*kind(timings_1b), 'timings_1b', 'main_all')
    allocate (load_imbalance(t_inc%nkmesh), stat=i_stat)
    call memocc(i_stat, product(shape(load_imbalance))*kind(load_imbalance), 'load_imbalance', 'main_all')
    timings_1a(:, :) = 0.0d0
    timings_1b(:) = 0.0d0
    load_imbalance(:) = 0
  end if                           ! MPIadapt>0

  ! create_subcomms_2d: first find maximal dimensions
  call find_dims_2d(nranks, t_inc%natyp, t_inc%ielast, dims, mpiatom)
  ! save in dims
  t_mpi_c_grid%dims = dims

  ! create communicator for atom/energy matrix
  call create_newcomms_group_ie(master, nranks, myrank, dims(1), dims(2), t_inc%nkmesh, t_inc%kmesh, mympi_comm_ie, myrank_ie, nranks_ie, mympi_comm_at, myrank_at, nranks_at, &
    myrank_atcomm, nranks_atcomm)
  ! save grid info in type 't_mpi_c_grid'
  call save_t_mpi_c_grid(t_mpi_c_grid, dims, mympi_comm_ie, mympi_comm_at, myrank_ie, myrank_at, myrank_atcomm, nranks_ie, nranks_at, nranks_atcomm)
  if (myrank==master) call timing_stop('MPI 1')


  ! communicate parameters for save_wavefunctions
  call bcast_params_savewf(t_wavefunctions)
  ! communicate parameters of godfrin inversion scheme ! GODFRIN Flaviano
  call bcast_params_godfrin(t_godfrin) ! GODFRIN Flaviano
#endif
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< distribute stuff from main0 !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! initialize timing and output files
  ! set if(t_inc%i_write>0) in front of every write(1337) and in mod_timing for timing_stop writeout (if(t_inc%i_time>0))
  ! for i_write (or i_time) =2 do not reset files > here for output.*.txt, after main2, copy writeout after main0 to different file
  if (myrank/=master) call timing_init(myrank)
  if (t_inc%i_write<2) then
#ifdef CPP_OLDCOMP
    if (myrank==master) call system('cp output.000.txt output.0.txt')
#else
    if (myrank==master) call execute_command_line('cp output.000.txt output.0.txt')
#endif
    if (myrank==master) close (1337, status='delete')
    if (t_inc%i_write>0) then
      open (1337, file='output.'//trim(ctemp)//'.txt')
      call version_print_header(1337)
    end if
  end if
  if (t_inc%i_write>0 .and. myrank/=master) then
    open (1337, file='output.'//trim(ctemp)//'.txt')
    call version_print_header(1337)
  end if

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SCF-ITERATION >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Now start scf iterations and do all steps of the KKR formalism until convergence
  if (myrank==master) then
    write (*, *) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write (*, *) '+++            SCF ITERATIONS START                +++'
    write (*, *) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    call print_time_and_date('started')
  end if
  do while ((t_inc%i_iteration<t_inc%n_iteration) .and. (t_inc%n_iteration/=-1))

    ! reset files for t_inc%i_write<2
    ! first copy lat output to output.2.txt so that all information of the precious iteration can be accessed while the next iteration runs
    if (t_inc%i_write<2 .and. t_inc%i_write>0 .and. myrank==master .and. t_inc%i_iteration>1) then
    ! the old intel compiler does not know 'execute_command_line' and has to use the 'system' call
#ifdef CPP_OLDCOMP
      call system('cp output.000.txt output.2.txt')
#else
      call execute_command_line('cp output.000.txt output.2.txt')
#endif
    end if
    ! rewind output.xxx.txt
    if (t_inc%i_write<2 .and. t_inc%i_write>0) then
      rewind (1337)
      if (.not. test('noserial')) read (1337, *)    ! skip first line to keep serial number
    end if
    ! rewind timing files if t_inc%i_time<2 (see mod_timing)
    if (t_inc%i_time<2 .and. t_inc%i_time>0) then
      rewind (43234059)
      if (.not. test('noserial')) read (43234059, *) ! skip first line to keep serial number
    end if

    call timing_start('Time in Iteration')

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculate tmat and gref
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call timing_start('main1a')
    call main1a()
    call timing_stop('main1a')
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculate gmat
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call timing_start('main1b')
    call main1b()
    call timing_stop('main1b')
    if (test('STOP1B  ')) then
      if (.not. opt('WRTGREEN') .and. myrank==master) write (*, *) 'done with WRTGREEN step'
      if (myrank==master) write (*, *) 'Stop after main1b'
#ifdef CPP_MPI
      call mpi_barrier(mpi_comm_world, ierr)
      call mpi_finalize(ierr)
#endif
      stop
    end if                         ! test

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculate density
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call timing_start('main1c')
    call main1c()
    call timing_stop('main1c')
    if (test('STOP1C  ')) then
      if (.not. opt('qdos    ') .and. myrank==master) write (*, *) 'done with qdos steps'
      if (myrank==master) write (*, *) 'Stop after main1c'
#ifdef CPP_MPI
      call mpi_barrier(mpi_comm_world, ierr)
      call mpi_finalize(ierr)
#endif
      stop
    end if                         ! test

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculate DFT stuff (potential from density, exc-potential, calculate total energy, ...)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call timing_start('main2')
    if (myrank==master) then
      call main2()
    end if
    call timing_stop('main2')

    ! reset arrays for next iteration
    if (t_params%lly/=0) then
      if (.not. t_lloyd%dtmat_to_file) t_lloyd%dtmat = czero
      if (.not. t_lloyd%tralpha_to_file) t_lloyd%tralpha = czero
      if (.not. t_lloyd%g0tr_to_file) t_lloyd%g0tr = czero
    end if

#ifdef CPP_MPI
    if (myrank==master) call timing_start('MPI 2')

    ! call myMPI_update_iteration()
    ! update i_iteration after check for convergence in main2:
    call mpi_bcast(t_inc%i_iteration, 1, mpi_integer, master, mpi_comm_world, ierr)
    if (ierr/=mpi_success) stop 'error broadcasting i_iteration in main_all'
    ! broadcast parameter arrays from master (e.g. update nonco angles etc.)
    call bcast_t_params_scalars(t_params)
    call bcast_t_params_arrays(t_params)

    ! find out if MPI_communication pattern should be modified: (with test option 'MPIadapt' the program will be forced to change the communication grid after the first iteration and then compares the timings
    ! call myMPI_distribute_ranks(t_inc, MPIatom, MPIadapt, timings_1a, timings_1b, load_imbalance, nranks, myrank, initial=0)
    if (mpiadapt==1 .and. t_inc%i_iteration>1) then
      call check_communication_pattern(mpiatom, mpiadapt, timings_1a, timings_1b, load_imbalance, t_inc%nkmesh, t_inc%kmesh_ie)
    end if
    ! adapt MPI communicator grid to tackle load imbalance better
    if (mpiadapt>0) then
      ! create_subcomms_2d: first find maximal dimensions
      call find_dims_2d(nranks, t_inc%natyp, t_inc%ielast, dims, mpiatom)
      ! save in dims
      t_mpi_c_grid%dims = dims

      ! create communicator for atom/energy matrix (load_imbalance instead of t_inc%kmesh in callig list)
      call create_newcomms_group_ie(master, nranks, myrank, dims(1), dims(2), t_inc%nkmesh, load_imbalance, mympi_comm_ie, myrank_ie, nranks_ie, mympi_comm_at, myrank_at, &
        nranks_at, myrank_atcomm, nranks_atcomm)
      ! save grid info in type 't_mpi_c_grid'
      call save_t_mpi_c_grid(t_mpi_c_grid, dims, mympi_comm_ie, mympi_comm_at, myrank_ie, myrank_at, myrank_atcomm, nranks_ie, nranks_at, nranks_atcomm)
    end if

    if (myrank==master) call timing_stop('MPI 2')
#endif

    call timing_stop('Time in Iteration')

    if (myrank==master) call print_time_and_date('Iteration finished')

  end do                           ! scf-iteration
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< SCF-ITERATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! close allocated arrays and finalize MPI >>>>>>>>>>>>>>>
  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  i_all = -product(shape(t_inc%kmesh))*kind(t_inc%kmesh)
  deallocate (t_inc%kmesh, stat=i_stat)
  call memocc(i_stat, i_all, 't_inc%kmesh', 'main_all')

  ! deallocate arrays from t_params
  i_all = -product(shape(t_params%a))*kind(t_params%a)
  deallocate (t_params%a, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%A', 'main_all')
  i_all = -product(shape(t_params%b))*kind(t_params%b)
  deallocate (t_params%b, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%B', 'main_all')
  i_all = -product(shape(t_params%rmesh))*kind(t_params%rmesh)
  deallocate (t_params%rmesh, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%RMESH', 'main_all')
  i_all = -product(shape(t_params%rr))*kind(t_params%rr)
  deallocate (t_params%rr, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%RR', 'main_all')
  i_all = -product(shape(t_params%eu))*kind(t_params%eu)
  deallocate (t_params%eu, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%EU', 'main_all')
  i_all = -product(shape(t_params%ez))*kind(t_params%ez)
  deallocate (t_params%ez, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%EZ', 'main_all')
  i_all = -product(shape(t_params%rc))*kind(t_params%rc)
  deallocate (t_params%rc, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%RC', 'main_all')
  i_all = -product(shape(t_params%rmt))*kind(t_params%rmt)
  deallocate (t_params%rmt, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%RMT', 'main_all')
  i_all = -product(shape(t_params%cls))*kind(t_params%cls)
  deallocate (t_params%cls, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%CLS', 'main_all')
  i_all = -product(shape(t_params%ish))*kind(t_params%ish)
  deallocate (t_params%ish, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%ISH', 'main_all')
  i_all = -product(shape(t_params%jsh))*kind(t_params%jsh)
  deallocate (t_params%jsh, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%JSH', 'main_all')
  i_all = -product(shape(t_params%edc))*kind(t_params%edc)
  deallocate (t_params%edc, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%EDC', 'main_all')
  i_all = -product(shape(t_params%noq))*kind(t_params%noq)
  deallocate (t_params%noq, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%NOQ', 'main_all')
  i_all = -product(shape(t_params%rws))*kind(t_params%rws)
  deallocate (t_params%rws, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%RWS', 'main_all')
  i_all = -product(shape(t_params%gsh))*kind(t_params%gsh)
  deallocate (t_params%gsh, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%GSH', 'main_all')
  i_all = -product(shape(t_params%zat))*kind(t_params%zat)
  deallocate (t_params%zat, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%ZAT', 'main_all')
  i_all = -product(shape(t_params%wez))*kind(t_params%wez)
  deallocate (t_params%wez, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%WEZ', 'main_all')
  i_all = -product(shape(t_params%vbc))*kind(t_params%vbc)
  deallocate (t_params%vbc, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%VBC', 'main_all')
  i_all = -product(shape(t_params%imt))*kind(t_params%imt)
  deallocate (t_params%imt, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%IMT', 'main_all')
  i_all = -product(shape(t_params%irc))*kind(t_params%irc)
  deallocate (t_params%irc, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%IRC', 'main_all')
  i_all = -product(shape(t_params%nfu))*kind(t_params%nfu)
  deallocate (t_params%nfu, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%NFU', 'main_all')
  i_all = -product(shape(t_params%ilm_map))*kind(t_params%ilm_map)
  deallocate (t_params%ilm_map, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%ILM_MAP', 'main_all')
  i_all = -product(shape(t_params%phi))*kind(t_params%phi)
  deallocate (t_params%phi, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%PHI', 'main_all')
  i_all = -product(shape(t_params%txc))*kind(t_params%txc)
  deallocate (t_params%txc, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%TXC', 'main_all')
  i_all = -product(shape(t_params%ipan))*kind(t_params%ipan)
  deallocate (t_params%ipan, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%IPAN', 'main_all')
  i_all = -product(shape(t_params%jend))*kind(t_params%jend)
  deallocate (t_params%jend, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%JEND', 'main_all')
  i_all = -product(shape(t_params%rnew))*kind(t_params%rnew)
  deallocate (t_params%rnew, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%RNEW', 'main_all')
  i_all = -product(shape(t_params%iqat))*kind(t_params%iqat)
  deallocate (t_params%iqat, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%IQAT', 'main_all')
  i_all = -product(shape(t_params%icpa))*kind(t_params%icpa)
  deallocate (t_params%icpa, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%ICPA', 'main_all')
  i_all = -product(shape(t_params%ezoa))*kind(t_params%ezoa)
  deallocate (t_params%ezoa, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%EZOA', 'main_all')
  i_all = -product(shape(t_params%atom))*kind(t_params%atom)
  deallocate (t_params%atom, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%ATOM', 'main_all')
  i_all = -product(shape(t_params%lopt))*kind(t_params%lopt)
  deallocate (t_params%lopt, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%LOPT', 'main_all')
  i_all = -product(shape(t_params%zrel))*kind(t_params%zrel)
  deallocate (t_params%zrel, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%ZREL', 'main_all')
  i_all = -product(shape(t_params%ueff))*kind(t_params%ueff)
  deallocate (t_params%ueff, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%UEFF', 'main_all')
  i_all = -product(shape(t_params%jeff))*kind(t_params%jeff)
  deallocate (t_params%jeff, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%JEFF', 'main_all')
  i_all = -product(shape(t_params%espv))*kind(t_params%espv)
  deallocate (t_params%espv, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%ESPV', 'main_all')
  i_all = -product(shape(t_params%rhoc))*kind(t_params%rhoc)
  deallocate (t_params%rhoc, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%RHOC', 'main_all')
  i_all = -product(shape(t_params%conc))*kind(t_params%conc)
  deallocate (t_params%conc, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%CONC', 'main_all')
  i_all = -product(shape(t_params%rrot))*kind(t_params%rrot)
  deallocate (t_params%rrot, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%RROT', 'main_all')
  i_all = -product(shape(t_params%vref))*kind(t_params%vref)
  deallocate (t_params%vref, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%VREF', 'main_all')
  i_all = -product(shape(t_params%cleb))*kind(t_params%cleb)
  deallocate (t_params%cleb, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%CLEB', 'main_all')
  i_all = -product(shape(t_params%rcls))*kind(t_params%rcls)
  deallocate (t_params%rcls, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%RCLS', 'main_all')
  i_all = -product(shape(t_params%vins))*kind(t_params%vins)
  deallocate (t_params%vins, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%VINS', 'main_all')
  i_all = -product(shape(t_params%visp))*kind(t_params%visp)
  deallocate (t_params%visp, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%VISP', 'main_all')
  i_all = -product(shape(t_params%cscl))*kind(t_params%cscl)
  deallocate (t_params%cscl, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%CSCL', 'main_all')
  i_all = -product(shape(t_params%drdi))*kind(t_params%drdi)
  deallocate (t_params%drdi, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%DRDI', 'main_all')
  i_all = -product(shape(t_params%nsh1))*kind(t_params%nsh1)
  deallocate (t_params%nsh1, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%NSH1', 'main_all')
  i_all = -product(shape(t_params%nsh2))*kind(t_params%nsh2)
  deallocate (t_params%nsh2, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%NSH2', 'main_all')
  i_all = -product(shape(t_params%crel))*kind(t_params%crel)
  deallocate (t_params%crel, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%CREL', 'main_all')
  i_all = -product(shape(t_params%rrel))*kind(t_params%rrel)
  deallocate (t_params%rrel, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%RREL', 'main_all')
  i_all = -product(shape(t_params%irns))*kind(t_params%irns)
  deallocate (t_params%irns, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%IRNS', 'main_all')
  i_all = -product(shape(t_params%lmsp))*kind(t_params%lmsp)
  deallocate (t_params%lmsp, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%LMSP', 'main_all')
  i_all = -product(shape(t_params%optc))*kind(t_params%optc)
  deallocate (t_params%optc, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%OPTC', 'main_all')
  i_all = -product(shape(t_params%bzkp))*kind(t_params%bzkp)
  deallocate (t_params%bzkp, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%BZKP', 'main_all')
  i_all = -product(shape(t_params%irws))*kind(t_params%irws)
  deallocate (t_params%irws, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%IRWS', 'main_all')
  i_all = -product(shape(t_params%ecore))*kind(t_params%ecore)
  deallocate (t_params%ecore, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%ECORE', 'main_all')
  i_all = -product(shape(t_params%qmtet))*kind(t_params%qmtet)
  deallocate (t_params%qmtet, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%QMTET', 'main_all')
  i_all = -product(shape(t_params%qmphi))*kind(t_params%qmphi)
  deallocate (t_params%qmphi, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%QMPHI', 'main_all')
  i_all = -product(shape(t_params%rmrel))*kind(t_params%rmrel)
  deallocate (t_params%rmrel, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%RMREL', 'main_all')
  i_all = -product(shape(t_params%vtrel))*kind(t_params%vtrel)
  deallocate (t_params%vtrel, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%VTREL', 'main_all')
  i_all = -product(shape(t_params%btrel))*kind(t_params%btrel)
  deallocate (t_params%btrel, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%BTREL', 'main_all')
  i_all = -product(shape(t_params%srrel))*kind(t_params%srrel)
  deallocate (t_params%srrel, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%SRREL', 'main_all')
  i_all = -product(shape(t_params%drotq))*kind(t_params%drotq)
  deallocate (t_params%drotq, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%DROTQ', 'main_all')
  i_all = -product(shape(t_params%ratom))*kind(t_params%ratom)
  deallocate (t_params%ratom, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%RATOM', 'main_all')
  i_all = -product(shape(t_params%lcore))*kind(t_params%lcore)
  deallocate (t_params%lcore, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%LCORE', 'main_all')
  i_all = -product(shape(t_params%ncore))*kind(t_params%ncore)
  deallocate (t_params%ncore, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%NCORE', 'main_all')
  i_all = -product(shape(t_params%ircut))*kind(t_params%ircut)
  deallocate (t_params%ircut, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%IRCUT', 'main_all')
  i_all = -product(shape(t_params%icleb))*kind(t_params%icleb)
  deallocate (t_params%icleb, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%ICLEB', 'main_all')
  i_all = -product(shape(t_params%nacls))*kind(t_params%nacls)
  deallocate (t_params%nacls, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%NACLS', 'main_all')
  i_all = -product(shape(t_params%loflm))*kind(t_params%loflm)
  deallocate (t_params%loflm, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%LOFLM', 'main_all')
  i_all = -product(shape(t_params%kaoez))*kind(t_params%kaoez)
  deallocate (t_params%kaoez, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%KAOEZ', 'main_all')
  i_all = -product(shape(t_params%uldau))*kind(t_params%uldau)
  deallocate (t_params%uldau, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%ULDAU', 'main_all')
  i_all = -product(shape(t_params%wldau))*kind(t_params%wldau)
  deallocate (t_params%wldau, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%WLDAU', 'main_all')
  i_all = -product(shape(t_params%kmesh))*kind(t_params%kmesh)
  deallocate (t_params%kmesh, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%KMESH', 'main_all')
  i_all = -product(shape(t_params%r2nef))*kind(t_params%r2nef)
  deallocate (t_params%r2nef, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%R2NEF', 'main_all')
  i_all = -product(shape(t_params%mvevi))*kind(t_params%mvevi)
  deallocate (t_params%mvevi, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%MVEVI', 'main_all')
  i_all = -product(shape(t_params%irrel))*kind(t_params%irrel)
  deallocate (t_params%irrel, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%IRREL', 'main_all')
  i_all = -product(shape(t_params%nrrel))*kind(t_params%nrrel)
  deallocate (t_params%nrrel, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%NRREL', 'main_all')
  i_all = -product(shape(t_params%lmsp1))*kind(t_params%lmsp1)
  deallocate (t_params%lmsp1, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%LMSP1', 'main_all')
  i_all = -product(shape(t_params%ifunm))*kind(t_params%ifunm)
  deallocate (t_params%ifunm, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%IFUNM', 'main_all')
  i_all = -product(shape(t_params%llmsp))*kind(t_params%llmsp)
  deallocate (t_params%llmsp, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%LLSMP', 'main_all')
  i_all = -product(shape(t_params%irmin))*kind(t_params%irmin)
  deallocate (t_params%irmin, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%IRMIN', 'main_all')
  i_all = -product(shape(t_params%testc))*kind(t_params%testc)
  deallocate (t_params%testc, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%TESTC', 'main_all')
  i_all = -product(shape(t_params%volbz))*kind(t_params%volbz)
  deallocate (t_params%volbz, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%VOLBZ', 'main_all')
  i_all = -product(shape(t_params%nofks))*kind(t_params%nofks)
  deallocate (t_params%nofks, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%NOFKS', 'main_all')
  i_all = -product(shape(t_params%theta))*kind(t_params%theta)
  deallocate (t_params%theta, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%THETA', 'main_all')
  i_all = -product(shape(t_params%dsymll))*kind(t_params%dsymll)
  deallocate (t_params%dsymll, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%DSYMLL', 'main_all')
  i_all = -product(shape(t_params%itldau))*kind(t_params%itldau)
  deallocate (t_params%itldau, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%ITLDAU', 'main_all')
  i_all = -product(shape(t_params%jwsrel))*kind(t_params%jwsrel)
  deallocate (t_params%jwsrel, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%JWSREL', 'main_all')
  i_all = -product(shape(t_params%thetas))*kind(t_params%thetas)
  deallocate (t_params%thetas, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%THETAS', 'main_all')
  i_all = -product(shape(t_params%rmtnew))*kind(t_params%rmtnew)
  deallocate (t_params%rmtnew, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%RMTNEW', 'main_all')
  i_all = -product(shape(t_params%rhoorb))*kind(t_params%rhoorb)
  deallocate (t_params%rhoorb, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%RHOORB', 'main_all')
  i_all = -product(shape(t_params%nshell))*kind(t_params%nshell)
  deallocate (t_params%nshell, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%NSHELL', 'main_all')
  i_all = -product(shape(t_params%rho2ns))*kind(t_params%rho2ns)
  deallocate (t_params%rho2ns, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%RHO2NS', 'main_all')
  i_all = -product(shape(t_params%rmtref))*kind(t_params%rmtref)
  deallocate (t_params%rmtref, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%RMTREF', 'main_all')
  i_all = -product(shape(t_params%socscl))*kind(t_params%socscl)
  deallocate (t_params%socscl, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%SOCSCL', 'main_all')
  i_all = -product(shape(t_params%rbasis))*kind(t_params%rbasis)
  deallocate (t_params%rbasis, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%RBASIS', 'main_all')
  i_all = -product(shape(t_params%refpot))*kind(t_params%refpot)
  deallocate (t_params%refpot, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%REFPOT', 'main_all')
  i_all = -product(shape(t_params%ifunm1))*kind(t_params%ifunm1)
  deallocate (t_params%ifunm1, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%IFUNM1', 'main_all')
  i_all = -product(shape(t_params%ititle))*kind(t_params%ititle)
  deallocate (t_params%ititle, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%ITITLE', 'main_all')
  i_all = -product(shape(t_params%ntcell))*kind(t_params%ntcell)
  deallocate (t_params%ntcell, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%NTCELL', 'main_all')
  i_all = -product(shape(t_params%ixipol))*kind(t_params%ixipol)
  deallocate (t_params%ixipol, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%IXIPOL', 'main_all')
  i_all = -product(shape(t_params%imaxsh))*kind(t_params%imaxsh)
  deallocate (t_params%imaxsh, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%IMAXSH', 'main_all')
  i_all = -product(shape(t_params%nkcore))*kind(t_params%nkcore)
  deallocate (t_params%nkcore, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%NKCORE', 'main_all')
  i_all = -product(shape(t_params%volcub))*kind(t_params%volcub)
  deallocate (t_params%volcub, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%VOLCUB', 'main_all')
  i_all = -product(shape(t_params%iqcalc))*kind(t_params%iqcalc)
  deallocate (t_params%iqcalc, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%IQCALC', 'main_all')
  i_all = -product(shape(t_params%icheck))*kind(t_params%icheck)
  deallocate (t_params%icheck, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%ICHECK', 'main_all')
  i_all = -product(shape(t_params%phildau))*kind(t_params%phildau)
  deallocate (t_params%phildau, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%PHILDAU', 'main_all')
  i_all = -product(shape(t_params%drdirel))*kind(t_params%drdirel)
  deallocate (t_params%drdirel, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%DRDIREL', 'main_all')
  i_all = -product(shape(t_params%denefat))*kind(t_params%denefat)
  deallocate (t_params%denefat, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%DENEFAT', 'main_all')
  i_all = -product(shape(t_params%rclsimp))*kind(t_params%rclsimp)
  deallocate (t_params%rclsimp, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%RCLSIMP', 'main_all')
  i_all = -product(shape(t_params%mvevief))*kind(t_params%mvevief)
  deallocate (t_params%mvevief, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%MVEVIEF', 'main_all')
  i_all = -product(shape(t_params%irshift))*kind(t_params%irshift)
  deallocate (t_params%irshift, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%IRSHIFT', 'main_all')
  i_all = -product(shape(t_params%ijtabsh))*kind(t_params%ijtabsh)
  deallocate (t_params%ijtabsh, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%IJTABSH', 'main_all')
  i_all = -product(shape(t_params%atomimp))*kind(t_params%atomimp)
  deallocate (t_params%atomimp, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%ATOMIMP', 'main_all')
  i_all = -product(shape(t_params%vacflag))*kind(t_params%vacflag)
  deallocate (t_params%vacflag, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%VACFLAG', 'main_all')
  i_all = -product(shape(t_params%kapcore))*kind(t_params%kapcore)
  deallocate (t_params%kapcore, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%KAPCORE', 'main_all')
  i_all = -product(shape(t_params%hostimp))*kind(t_params%hostimp)
  deallocate (t_params%hostimp, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%HOSTIMP', 'main_all')
  i_all = -product(shape(t_params%ecorerel))*kind(t_params%ecorerel)
  deallocate (t_params%ecorerel, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%ECOREREL', 'main_all')
  i_all = -product(shape(t_params%ijtabsym))*kind(t_params%ijtabsym)
  deallocate (t_params%ijtabsym, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%IJTABSYM', 'main_all')
  i_all = -product(shape(t_params%erefldau))*kind(t_params%erefldau)
  deallocate (t_params%erefldau, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%EREFLDAU', 'main_all')
  i_all = -product(shape(t_params%qmphitab))*kind(t_params%qmphitab)
  deallocate (t_params%qmphitab, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%QMPHITAB', 'main_all')
  i_all = -product(shape(t_params%qmtettab))*kind(t_params%qmtettab)
  deallocate (t_params%qmtettab, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%QMTETTAB', 'main_all')
  i_all = -product(shape(t_params%qmgamtab))*kind(t_params%qmgamtab)
  deallocate (t_params%qmgamtab, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%QMGAMTAB', 'main_all')
  i_all = -product(shape(t_params%cmomhost))*kind(t_params%cmomhost)
  deallocate (t_params%cmomhost, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%CMOMHOST', 'main_all')
  i_all = -product(shape(t_params%socscale))*kind(t_params%socscale)
  deallocate (t_params%socscale, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%SOCSCALE', 'main_all')
  i_all = -product(shape(t_params%npan_tot))*kind(t_params%npan_tot)
  deallocate (t_params%npan_tot, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%NPAN_TOT', 'main_all')
  i_all = -product(shape(t_params%thetasnew))*kind(t_params%thetasnew)
  deallocate (t_params%thetasnew, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%THETASNEW', 'main_all')
  i_all = -product(shape(t_params%ijtabcalc))*kind(t_params%ijtabcalc)
  deallocate (t_params%ijtabcalc, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%IJTABCALC', 'main_all')
  i_all = -product(shape(t_params%r2drdirel))*kind(t_params%r2drdirel)
  deallocate (t_params%r2drdirel, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%R2DRDIREL', 'main_all')
  i_all = -product(shape(t_params%npan_eq_at))*kind(t_params%npan_eq_at)
  deallocate (t_params%npan_eq_at, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%NPAN_EQ_AT', 'main_all')
  i_all = -product(shape(t_params%symunitary))*kind(t_params%symunitary)
  deallocate (t_params%symunitary, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%SYMUNITARY', 'main_all')
  i_all = -product(shape(t_params%lefttinvll))*kind(t_params%lefttinvll)
  deallocate (t_params%lefttinvll, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%LEFTTINVLL', 'main_all')
  i_all = -product(shape(t_params%righttinvll))*kind(t_params%righttinvll)
  deallocate (t_params%righttinvll, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%RIGHTTTINVLL', 'main_all')
  i_all = -product(shape(t_params%npan_log_at))*kind(t_params%npan_log_at)
  deallocate (t_params%npan_log_at, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%NPAN_LOG_AT', 'main_all')
  i_all = -product(shape(t_params%ijtabcalc_i))*kind(t_params%ijtabcalc_i)
  deallocate (t_params%ijtabcalc_i, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%IJTABCALC_I', 'main_all')
  i_all = -product(shape(t_params%ipan_intervall))*kind(t_params%ipan_intervall)
  deallocate (t_params%ipan_intervall, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%IPAN_INTERVALL', 'main_all')
  i_all = -product(shape(t_params%rpan_intervall))*kind(t_params%rpan_intervall)
  deallocate (t_params%rpan_intervall, stat=i_stat)
  call memocc(i_stat, i_all, 't_params%RPAN_INTERVALL', 'main_all')

  ! delete temporary files
  if (myrank==master) then
    open (69, file='abvmad.unformatted')
    close (69, status='delete')
  end if

  ! deallocate arrays from t_wavefunctions
  if (t_wavefunctions%nwfsavemax>0) then

    i_all = -product(shape(t_wavefunctions%isave_wavefun))*kind(t_wavefunctions%isave_wavefun)
    deallocate (t_wavefunctions%isave_wavefun, stat=i_stat)
    call memocc(i_stat, i_all, 't_wavefunctions%isave_wavefun', 'main_all')

    if (t_wavefunctions%save_rll) then
      i_all = -product(shape(t_wavefunctions%rll))*kind(t_wavefunctions%rll)
      deallocate (t_wavefunctions%rll, stat=i_stat)
      call memocc(i_stat, i_all, 't_wavefunctions%rll', 'main_all')
    end if

    if (t_wavefunctions%save_sll) then
      i_all = -product(shape(t_wavefunctions%sll))*kind(t_wavefunctions%sll)
      deallocate (t_wavefunctions%sll, stat=i_stat)
      call memocc(i_stat, i_all, 't_wavefunctions%sll', 'main_all')
    end if

    if (t_wavefunctions%save_rllleft) then
      i_all = -product(shape(t_wavefunctions%rllleft))*kind(t_wavefunctions%rllleft)
      deallocate (t_wavefunctions%rllleft, stat=i_stat)
      call memocc(i_stat, i_all, 't_wavefunctions%rllleft', 'main_all')
    end if

    if (t_wavefunctions%save_sllleft) then
      i_all = -product(shape(t_wavefunctions%sllleft))*kind(t_wavefunctions%sllleft)
      deallocate (t_wavefunctions%sllleft, stat=i_stat)
      call memocc(i_stat, i_all, 't_wavefunctions%sllleft', 'main_all')
    end if

  end if                           ! (t_wavefunctions%Nwfsavemax>0)

#ifdef CPP_MPI
  ! deallocate arrays for MPIadapt
  if (mpiadapt>0) then
    i_all = -product(shape(timings_1a))*kind(timings_1a)
    deallocate (timings_1a, stat=i_stat)
    call memocc(i_stat, i_all, 'timings_1a', 'main_all')

    i_all = -product(shape(timings_1b))*kind(timings_1b)
    deallocate (timings_1b, stat=i_stat)
    call memocc(i_stat, i_all, 'timings_1b', 'main_all')

    i_all = -product(shape(load_imbalance))*kind(load_imbalance)
    deallocate (load_imbalance, stat=i_stat)
    call memocc(i_stat, i_all, 'load_imbalance', 'main_all')

  end if
#endif

  ! Deallocation of input arrays
  call allocate_cell(-1, naez, nemb, natyp, cls, imt, irws, irns, ntcell, refpot, kfg, kaoez, rmt, zat, rws, mtfac, rmtref, rmtrefat, rmtnew, rbasis, lmxc)
  call allocate_semi_inf_host(-1, nemb, tleft, tright)
  call allocate_potential(-1, irm, natyp, npotd, ipand, nfund, lmxspd, lmpot, irmind, nspotd, nfu, irc, ncore, irmin, lmsp, lmsp1, ircut, lcore, llmsp, ititle, fpradius, visp, &
    ecore, vins)
  call allocate_cpa(-1, naez, natyp, noq, icpa, iqat, hostimp, conc)
  call allocate_ldau(-1, natyp, lopt, ueff, jeff, erefldau)
  call allocate_ldau_potential(-1, irm, natyp, mmaxd, nspind, itldau, wldau, uldau, phildau)
  call allocate_magnetization(-1, naez, natyp, lmmaxd, inipol, ixipol, qmtet, qmphi, drotq)
  call allocate_soc(-1, krel, natyp, lmax, socscale, cscl, socscl)
  call allocate_energies(-1, iemxd, ez, dez, wez)
  call allocate_relativistic(-1, krel, irm, naez, natyp, zrel, jwsrel, irshift, vtrel, btrel, rmrel, drdirel, r2drdirel, qmgam, qmgamtab, qmphitab, qmtettab)
  call allocate_rel_transformations(-1, lmmaxd, nrrel, irrel, rc, crel, rrel, srrel)
  call allocate_clusters(-1, naez, lmax, ncleb, nclsd, nembd1, nsheld, naclsd, lmpot, natomimpd, nsh1, nsh2, nacls, nshell, atomimp, atom, ezoa, icleb, jend, ratom, rclsimp, &
    cmomhost, rcls)
  call allocate_expansion(-1, lm2d, irid, nfund, ntotd, ncleb, lassld, ncelld, ncheb, loflm, wg, cleb, yrg, thetas, thetasnew)
  call allocate_mesh(-1, irm, natyp, a, b, rmesh, drdi)
  call allocate_pannels(-1, natyp, ntotd, ipan, npan_tot, npan_eq_at, npan_log_at, ipan_intervall, rpan_intervall)
  call allocate_misc(-1, nrd, irm, irid, lmax, naez, natyp, nfund, nref, iemxd, ntotd, nsheld, lmmaxd, nembd1, ncheb, ncelld, lmxspd, nspindd, nsymaxd, nprincd, ifunm, ifunm1, &
    icheck, vref, s, rr, dror, rnew, rs, rrot, thesme, dsymll, dsymll1, lefttinvll, righttinvll)
  call allocate_green(-1, naez, iemxd, ngshd, nsheld, lmpot, nofgij, ish, jsh, kmesh, imaxsh, iqcalc, iofgij, jofgij, ijtabsh, ijtabsym, ijtabcalc, ijtabcalc_i, ilm_map, gsh)
  ! End of deallocation

  ! print memory report to stdout
  if (t_inc%i_write>0) call memocc(0, 0, 'count', 'stop')

#ifdef CPP_MPI
  ! finalize MPI
  call mpi_finalize(ierr)
#endif
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< close allocated arrays and finalize MPI !!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program kkrcode

