!-------------------------------------------------------------------------------
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
!< @note
!> - Jonathan Chico Jan. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
program kkrcode

   use Profiling
   use mod_main0
   use mod_main1a
   use mod_main1b
   use mod_main1c
   use mod_main2
   use mod_types
   use mod_timing
   use mod_md5sums
   use memoryhandling
   use mod_version_info
   use global_variables

#ifdef CPP_MPI
   use mod_mympi, only: mympi_init, myrank, nranks, master,find_dims_2d,      &
                        distribute_linear_on_tasks, create_newcomms_group_ie, &
                        MPIatom , MPIadapt, check_communication_pattern
   use mod_save_wavefun, only: t_wavefunctions, bcast_params_savewf
#else
   use mod_mympi, only: mympi_init, myrank, nranks, master, MPIatom, MPIadapt
   use mod_save_wavefun, only: t_wavefunctions
#endif

#ifdef CPP_MPI
   use mpi
#endif

   implicit none

   integer :: ierr,i_stat,i_all
#ifdef CPP_MPI
   integer :: myMPI_comm_at, myMPI_comm_ie, nranks_at, myrank_at
   integer :: myrank_ie, nranks_ie, nranks_atcomm, myrank_atcomm
   integer, dimension(2) :: dims
#endif
   character(len=3) :: ctemp !name for output file
   ! needed to use test('xxxxxxxx'):
   logical :: test
   external :: test

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! initialize MPI >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef CPP_MPI
   ! initialize MPI
   call MPI_Init ( ierr )
#endif
   ! set variables master, myrank and nranks for serial (myrank=master=0, nranks=1) as well as parallel execution
   call mympi_init()
   ! save myrank in ctemp, needed to open output unit 1337
   write(ctemp,'(I03.3)') myrank
   ! find serial number that is printed to files
   call construct_serialnr()
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< initialize MPI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! start KKR with main0 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   ! do all this only on master (read in parameters etc. in main0)
   if(myrank==master) then

      ! initialize timing
      call timing_init(myrank)

      ! start KKR program, first do main0, where everything is read in and initialized
      call timing_start('main0')

      ! open output files
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*) '!!! Most output written to output.myrank.txt files !!!'
      write(*,*) '!!! please check these files as well               !!!'
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      open(1337, file='output.'//trim(ctemp)//'.txt')
      call version_print_header(1337)

      ! default value on master (needed for writeout in main0)
      t_inc%i_write = 1

      ! run main0 only serially, then communicate everything that is in here
      call main0()
      call timing_stop('main0')

   end if ! myrank ==master
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< start KKR with main0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   ! without MPI (serial or openMP) something goes wrong if if files are not written out
   ! this seems to be only the case with the old solver
#ifndef CPP_MPI
   if(.not.t_inc%NEWSOSOL) then
      t_tgmat%tmat_to_file = .true.
      t_tgmat%gmat_to_file = .true.
      t_tgmat%gref_to_file = .true.
   end if
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! distribute stuff from main0 >>>>>>>>>>>>>>>>>>>>>>>>>>>
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef CPP_MPI
   ! now communicate type t_inc and t_tgmat switches (this has an implicit barrier, so that all other processes wait for master to finish with main0)
   if(myrank==master) call timing_start('MPI 1')
   call bcast_t_inc_tgmat(t_inc,t_tgmat,t_cpa)
   ! also communicate logicals from t_lloyd
   call bcast_t_lly_1(t_lloyd)

   ! communicate parameters that were written in wunfiles into t_params
   call bcast_t_params_scalars(t_params)
   if (myrank.ne.master) call init_t_params(t_params)
   call bcast_t_params_arrays(t_params)


   ! broadcast md5 sums (written to some output files (kkrflex_* etc.)
#ifdef CPP_MPI
   call myMPI_Bcast_md5sums(t_params%INS, myrank, master)
#endif

   ! in case of deci-out run everything only serially to write decifile correctly. This will be fixed in a later version when we get rid of the decifile
   if(t_inc%deci_out .and. t_mpi_c_grid%nranks_at>1) stop 'deci-out option chosen. Please run code serially in energy dimension!'

   ! call myMPI_distribute_ranks(MPIatom, MPIadapt, t_inc, timings_1a, timings_1b, load_imbalance, nranks, myrank, initial=1)
   ! create 2d matrix for processors so that more processors than energy points can be used.
   ! strategy here is to first parallelize the energy points ideally and then give all processors that are left to atom or k-point loop parallelization

   ! communicate logical MPIatom and MPIadapt which determine how the ranks are devided into groups
   call MPI_Bcast(MPIatom, 1, MPI_LOGICAL, master, MPI_COMM_WORLD, ierr)
   if(ierr/=MPI_SUCCESS) stop 'error broadcasting MPIatom in main_all'
   call MPI_Bcast(MPIadapt, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
   if(ierr/=MPI_SUCCESS) stop 'error broadcasting MPIadapt in main_all'
   ! allocate timing arrays
   if(MPIadapt>0) then
      allocate(timings_1a(t_inc%ielast,t_inc%natyp),stat=i_stat)
      call memocc(i_stat,product(shape(timings_1a))*kind(timings_1a),'timings_1a','main_all')
      allocate(timings_1b(t_inc%ielast),stat=i_stat)
      call memocc(i_stat,product(shape(timings_1b))*kind(timings_1b),'timings_1b','main_all')
      allocate(load_imbalance(t_inc%nkmesh), stat=i_stat)
      call memocc(i_stat,product(shape(load_imbalance))*kind(load_imbalance),'load_imbalance','main_all')
      timings_1a(:,:) = 0.0d0
      timings_1b(:) = 0.0d0
      load_imbalance(:) = 0
   end if
   ! create_subcomms_2d: first find maximal dimensions
   call find_dims_2d(nranks,t_inc%NATYP,t_inc%IELAST,dims,MPIatom)
   !save in dims
   t_mpi_c_grid%dims = dims

   ! create communicator for atom/energy matrix
   call create_newcomms_group_ie( nranks,myrank,dims(1),dims(2),t_inc%nkmesh, &
      t_inc%kmesh,mympi_comm_ie,myrank_ie,nranks_ie,mympi_comm_at,myrank_at,  &
      nranks_at, myrank_atcomm,nranks_atcomm)
   ! save grid info in type 't_mpi_c_grid'
   call save_t_mpi_c_grid(t_mpi_c_grid,dims, myMPI_comm_ie, myMPI_comm_at, &
      myrank_ie, myrank_at, myrank_atcomm, nranks_ie, nranks_at, nranks_atcomm)
   if(myrank==master) call timing_stop('MPI 1')
#endif

   !call set_writeout_timings(t_inc, myrank, master, ctemp)

   ! set verbosity levels for each rank set i_write and i_time to 0 or 1,2, depending on verbosity level specified in inputcard and rank
   ! take care of different ranks here for verbose0 and timings1,2
   ! convention: value of 0 mean do not write, value 1 write and reset file, value 2 write everthing
   if(t_inc%i_write==0 .and. myrank==master) t_inc%i_write = 1 ! only written to master, reset file after each iteration, output.init.txt for main0 output
   if(t_inc%i_time==0 .and. myrank==master) t_inc%i_time = 1  ! only master writes timings of current iteration
   if(t_inc%i_time==1) then
      t_inc%i_time = 0                          ! all ranks do not write timings but only the master
      if(myrank==master) t_inc%i_time = 2       ! master write all timings of all iterations
   endif !t_inc%i_time==1
   ! initialize timing files
   if(myrank.ne.master) call timing_init(myrank)

   ! set if(t_inc%i_write>0) in front of every write(1337) and in mod_timing for timing_stop writeout (if(t_inc%i_time>0))
   ! for i_write (or i_time) =2 do not reset files > here for output.*.txt, after main2, copy writeout after main0 to different file
   if (t_inc%i_write<2) then
      if(myrank==master) call SYSTEM('cp output.000.txt output.0.txt')
      if(myrank==master) close(1337, status='delete')
      if(t_inc%i_write>0) then
         open(1337, file='output.'//trim(ctemp)//'.txt')
         call version_print_header(1337)
      end if
   endif
   if(t_inc%i_write>0 .and. myrank.ne.master) then
      open(1337, file='output.'//trim(ctemp)//'.txt')
      call version_print_header(1337)
   end if

#ifdef CPP_MPI
   ! communicate parameters for save_wavefunctions
   call bcast_params_savewf(t_wavefunctions)
#endif
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< distribute stuff from main0 !!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SCF-ITERATION >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   ! Now start scf iterations and do all steps of the KKR formalism until convergence
   if(myrank==master) then
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,*) '+++            SCF ITERATIONS START                +++'
      write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      call print_time_and_date('started')
   end if
   do while ( (t_inc%i_iteration.lt.t_inc%N_iteration) .and. (t_inc%N_iteration.ne.-1) )

      ! reset files for t_inc%i_write<2
      ! first copy lat output to output.2.txt so that all information of the precious iteration can be accessed while the next iteration runs
      if (t_inc%i_write<2 .and. t_inc%i_write>0 .and. myrank==master .and. t_inc%i_iteration>1) call SYSTEM('cp output.000.txt output.2.txt')
      ! rewind output.xxx.txt
      if (t_inc%i_write<2 .and. t_inc%i_write>0) then
         rewind(1337)
         read(1337,*) ! skip first line to keep serial number
      end if
      ! rewind timing files if t_inc%i_time<2 (see mod_timing)
      if (t_inc%i_time<2 .and. t_inc%i_time>0) then
         rewind(43234059)
         read(43234059,*) ! skip first line to keep serial number
      end if

      call timing_start('Time in Iteration')

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate tmat and gref
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call timing_start('main1a')
      call main1a(INS,LLY,IRM,LM2D,ICST,IEND,NCLS,LMAX,NREF,NSRA,NEMB,NAEZ,NATYP,&
         NCLSD,NPOTD,ITSCF,NTOTD,MMAXD,LMPOT,NINEQ,NSPIN,NCHEB,LMMAXD,IELAST,    &
         NRMAXD,IRMIND,NATOMIMP,ALAT,R_LOG,TOLRDIF,DELTAE,CLS,IQAT,IRWS,NACLS,   &
         REFPOT,ATOM,ZAT,VREF,RMTREF,RCLS,SOLVER,SOCSCL,SOCSCALE,CSCL,NTLDAU,    &
         IDOLDAU,ITLDAU,UEFF,JEFF,IPAN,LOFLM,IRMIN,ATOMIMP,ICLEB,IRCUT,          &
         IPAN_INTERVALL,PHI,THETA,CLEB,VISP,DRDI,RNEW,RMESH,RPAN_INTERVALL,VINS, &
         ZREL,JWSREL,VTREL,BTREL,RMREL,DRDIREL,R2DRDIREL,ITRUNLDAU,LOPT,EREFLDAU,&
         WLDAU,ULDAU,PHILDAU)
      call timing_stop('main1a')

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate gmat
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call timing_start('main1b')
      call main1b(NR,LLY,ICC,IGF,INS,NPOL,LMAX,NREF,NSRA,NCLS,NCPA,NEMB,NAEZ,    &
         NATYP,NCLSD,NSPIN,KMROT,LMMAXD,IELAST,INVMOD,NSYMAT,NEMBD1,LMGF0D,      &
         NOFGIJ,NQCALC,NSPINDD,NLBASIS,NRBASIS,MAXMESH,NATOMIMP,ITCPAMAX,ALAT,   &
         CPATOL,NOQ,CLS,IQAT,NSH1,NSH2,ICPA,NACLS,NSHELL,REFPOT,ATOMIMP,EZOA,    &
         ATOM,KAOEZ,ICHECK,CONC,RMTREF,RR,RATOM,RBASIS,RROT,RCLS,SYMUNITARY,     &
         KMESH,IQCALC,IJTABSH,IJTABSYM,IJTABCALC,IJTABCALC_I,ISH,JSH,NRREL,IRREL,&
         VREF,RCLSIMP,RC,RREL,CREL,SRREL,DROTQ,DSYMLL,LEFTTINVLL,RIGHTTINVLL)
      call timing_stop('main1b')
      if(test('STOP1B  '))then
#ifdef CPP_MPI
         call MPI_Finalize(ierr)
#endif
         stop 'Stop after main1b'
      end if !test

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate density
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call timing_start('main1c')
      call main1c(INS,LLY,IRM,LM2D,ICST,NAEZ,NPOL,NSRA,LMAX,NTOTD,MMAXD,NATYP,   &
         NPOTD,KMROT,NSPIN,NCHEB,LMPOT,LMXSPD,IELAST,LMMAXD,NRMAXD,IRMIND,       &
         INTERVX,INTERVY,INTERVZ,IESEMICORE,TK,EMIN,EMAX,ALAT,EFERMI,SOLVER,IQAT,&
         ZREL,IPAN,IRWS,NCORE,JWSREL,NTCELL,ITITLE,CSCL,ZAT,CONC,SOCSCALE,NTLDAU,&
         IDOLDAU,ITRUNLDAU,ITLDAU,UEFF,JEFF,IEND,NFU,LOFLM,IRMIN,IRSHIFT,ICLEB,  &
         LCORE,IRCUT,IFUNM1,LMSP1,LLMSP,JEND,A,B,QMTET,QMPHI,CLEB,DRDI,ECORE,    &
         RMREL,SOCSCL,R2DRDIREL,VINS,VTREL,BTREL,DRDIREL,EZ,WEZ,LOPT,EREFLDAU,   &
         WLDAU,ULDAU,PHILDAU,R_LOG,NPAN_EQ,NPAN_LOG,NPAN_TOT,IPAN_INTERVALL,VISP,&
         RNEW,RPAN_INTERVALL,THETAS,THETASNEW)
      call timing_stop('main1c')

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate DFT stuff (potential from density, exc-potential, calculate total energy, ...)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call timing_start('main2')
      if (myrank==master) then
         call main2(LLY,ICC,INS,IPF,KTE,KXC,LPOT,IMIX,NAEZ,NSRA,LMAX,KPRE,NPOL,  &
            NPNT1,NPNT2,NPNT3,NATYP,NSPIN,ITSCF,KVMAD,LMPOT,NPOTD,NLEFT,NRIGHT,  &
            LMXSPD,IELAST,ISHIFT,ITDBRY,KSHAPE,KFORCE,IRMIND,NEMBD1,LMMAXD,      &
            IDOLDAU,NLBASIS,NRBASIS,SCFSTEPS,NATOMIMP,TK,FCM,ALAT,MIXING,QBOUND, &
            LAMBDA_XC,OPT,TEST,LRHOSYM,LINTERFACE,NOQ,IMT,IQAT,IPAN,ZREL,IRNS,   &
            IRWS,KAOEZ,JWSREL,NTCELL,ITITLE,NSHELL,ZAT,RMT,RWS,CONC,RMTNEW,TXC,  &
            EMIN,EMAX,TKSEMI,EMUSEMI,EBOTSEMI,FSEMICORE,IRC,NFU,LOPT,NCORE,IRMIN,&
            IXIPOL,IMAXSH,IRSHIFT,ATOMIMP,HOSTIMP,ILM,LMSP,LCORE,IRCUT,IFUNM,A,B,&
            VBC,GSH,FACT,QMGAM,QMTET,QMPHI,R,ECORE,DRDI,VISP,VTREL,BTREL,RMREL,  &
            DRDIREL,CMOMHOST,R2DRDIREL,VINS,THETAS,THESME,EZ,DEZ,WEZ)
      endif
      call timing_stop('main2')

      ! reset arrays for next iteration
      if(t_params%LLY/=0) then
         if(.not.t_lloyd%dtmat_to_file) t_lloyd%dtmat = (0.0d0, 0.0d0)
         if(.not.t_lloyd%tralpha_to_file) t_lloyd%tralpha = (0.0d0, 0.0d0)
         if(.not.t_lloyd%g0tr_to_file) t_lloyd%g0tr= (0.0d0, 0.0d0)
      end if

#ifdef CPP_MPI
      if(myrank==master) call timing_start('MPI 2')

      !   call myMPI_update_iteration()
      ! update i_iteration after check for convergence in main2:
      call MPI_Bcast(t_inc%i_iteration, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      if(ierr/=MPI_SUCCESS) stop 'error broadcasting i_iteration in main_all'
      ! broadcast parameter arrays from master (e.g. update nonco angles etc.)
      call bcast_t_params_scalars(t_params)
      call bcast_t_params_arrays(t_params)

      ! find out if MPI_communication pattern should be modified: (with test option 'MPIadapt' the program will be forced to change the communication grid after the first iteration and then compares the timings
      !   call myMPI_distribute_ranks(t_inc, MPIatom, MPIadapt, timings_1a, timings_1b, load_imbalance, nranks, myrank, initial=0)
      if(MPIadapt==1 .and. t_inc%i_iteration>1) then
         call check_communication_pattern(MPIatom, MPIadapt, timings_1a, timings_1b, load_imbalance, t_inc%nkmesh, t_inc%kmesh_ie)
      end if
      ! adapt MPI communicator grid to tackle load imbalance better
      if(MPIadapt>0) then
         ! create_subcomms_2d: first find maximal dimensions
         call find_dims_2d(nranks,t_inc%NATYP,t_inc%IELAST,dims,MPIatom)
         !save in dims
         t_mpi_c_grid%dims = dims

         ! create communicator for atom/energy matrix (load_imbalance instead of t_inc%kmesh in callig list)
         call create_newcomms_group_ie( nranks,myrank,dims(1),dims(2),t_inc%nkmesh,load_imbalance,mympi_comm_ie,  &
            myrank_ie,nranks_ie,mympi_comm_at,myrank_at,nranks_at, myrank_atcomm,nranks_atcomm)
         ! save grid info in type 't_mpi_c_grid'
         call save_t_mpi_c_grid(t_mpi_c_grid,dims, myMPI_comm_ie, myMPI_comm_at, myrank_ie, myrank_at,            &
            myrank_atcomm, nranks_ie, nranks_at, nranks_atcomm)
      end if

      if(myrank==master) call timing_stop('MPI 2')
#endif

      call timing_stop('Time in Iteration')

      if (myrank==master) call print_time_and_date('Iteration finished')

   end do ! scf-iteration
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< SCF-ITERATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! close allocated arrays and finalize MPI >>>>>>>>>>>>>>>
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   ! deallocate arrays from t_params
   i_all=-product(shape(t_params%A))*kind(t_params%A)
   deallocate(t_params%A,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%A','main_all')
   i_all=-product(shape(t_params%B))*kind(t_params%B)
   deallocate(t_params%B,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%B','main_all')
   i_all=-product(shape(t_params%R))*kind(t_params%R)
   deallocate(t_params%R,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%R','main_all')
   i_all=-product(shape(t_params%RR))*kind(t_params%RR)
   deallocate(t_params%RR,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%RR','main_all')
   i_all=-product(shape(t_params%EU))*kind(t_params%EU)
   deallocate(t_params%EU,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%EU','main_all')
   i_all=-product(shape(t_params%EZ))*kind(t_params%EZ)
   deallocate(t_params%EZ,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%EZ','main_all')
   i_all=-product(shape(t_params%RC))*kind(t_params%RC)
   deallocate(t_params%RC,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%RC','main_all')
   i_all=-product(shape(t_params%RMT))*kind(t_params%RMT)
   deallocate(t_params%RMT,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%RMT','main_all')
   i_all=-product(shape(t_params%CLS))*kind(t_params%CLS)
   deallocate(t_params%CLS,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%CLS','main_all')
   i_all=-product(shape(t_params%ISH))*kind(t_params%ISH)
   deallocate(t_params%ISH,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%ISH','main_all')
   i_all=-product(shape(t_params%JSH))*kind(t_params%JSH)
   deallocate(t_params%JSH,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%JSH','main_all')
   i_all=-product(shape(t_params%EDC))*kind(t_params%EDC)
   deallocate(t_params%EDC,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%EDC','main_all')
   i_all=-product(shape(t_params%NOQ))*kind(t_params%NOQ)
   deallocate(t_params%NOQ,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%NOQ','main_all')
   i_all=-product(shape(t_params%RWS))*kind(t_params%RWS)
   deallocate(t_params%RWS,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%RWS','main_all')
   i_all=-product(shape(t_params%GSH))*kind(t_params%GSH)
   deallocate(t_params%GSH,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%GSH','main_all')
   i_all=-product(shape(t_params%ZAT))*kind(t_params%ZAT)
   deallocate(t_params%ZAT,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%ZAT','main_all')
   i_all=-product(shape(t_params%WEZ))*kind(t_params%WEZ)
   deallocate(t_params%WEZ,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%WEZ','main_all')
   i_all=-product(shape(t_params%VBC))*kind(t_params%VBC)
   deallocate(t_params%VBC,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%VBC','main_all')
   i_all=-product(shape(t_params%IMT))*kind(t_params%IMT)
   deallocate(t_params%IMT,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%IMT','main_all')
   i_all=-product(shape(t_params%IRC))*kind(t_params%IRC)
   deallocate(t_params%IRC,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%IRC','main_all')
   i_all=-product(shape(t_params%NFU))*kind(t_params%NFU)
   deallocate(t_params%NFU,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%NFU','main_all')
   i_all=-product(shape(t_params%ILM))*kind(t_params%ILM)
   deallocate(t_params%ILM,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%ILM','main_all')
   i_all=-product(shape(t_params%PHI))*kind(t_params%PHI)
   deallocate(t_params%PHI,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%PHI','main_all')
   i_all=-product(shape(t_params%TXC))*kind(t_params%TXC)
   deallocate(t_params%TXC,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%TXC','main_all')
   i_all=-product(shape(t_params%IPAN))*kind(t_params%IPAN)
   deallocate(t_params%IPAN,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%IPAN','main_all')
   i_all=-product(shape(t_params%JEND))*kind(t_params%JEND)
   deallocate(t_params%JEND,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%JEND','main_all')
   i_all=-product(shape(t_params%RNEW))*kind(t_params%RNEW)
   deallocate(t_params%RNEW,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%RNEW','main_all')
   i_all=-product(shape(t_params%IQAT))*kind(t_params%IQAT)
   deallocate(t_params%IQAT,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%IQAT','main_all')
   i_all=-product(shape(t_params%ICPA))*kind(t_params%ICPA)
   deallocate(t_params%ICPA,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%ICPA','main_all')
   i_all=-product(shape(t_params%EZOA))*kind(t_params%EZOA)
   deallocate(t_params%EZOA,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%EZOA','main_all')
   i_all=-product(shape(t_params%ATOM))*kind(t_params%ATOM)
   deallocate(t_params%ATOM,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%ATOM','main_all')
   i_all=-product(shape(t_params%LOPT))*kind(t_params%LOPT)
   deallocate(t_params%LOPT,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%LOPT','main_all')
   i_all=-product(shape(t_params%ZREL))*kind(t_params%ZREL)
   deallocate(t_params%ZREL,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%ZREL','main_all')
   i_all=-product(shape(t_params%UEFF))*kind(t_params%UEFF)
   deallocate(t_params%UEFF,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%UEFF','main_all')
   i_all=-product(shape(t_params%JEFF))*kind(t_params%JEFF)
   deallocate(t_params%JEFF,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%JEFF','main_all')
   i_all=-product(shape(t_params%ESPV))*kind(t_params%ESPV)
   deallocate(t_params%ESPV,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%ESPV','main_all')
   i_all=-product(shape(t_params%RHOC))*kind(t_params%RHOC)
   deallocate(t_params%RHOC,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%RHOC','main_all')
   i_all=-product(shape(t_params%CONC))*kind(t_params%CONC)
   deallocate(t_params%CONC,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%CONC','main_all')
   i_all=-product(shape(t_params%RROT))*kind(t_params%RROT)
   deallocate(t_params%RROT,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%RROT','main_all')
   i_all=-product(shape(t_params%VREF))*kind(t_params%VREF)
   deallocate(t_params%VREF,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%VREF','main_all')
   i_all=-product(shape(t_params%CLEB))*kind(t_params%CLEB)
   deallocate(t_params%CLEB,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%CLEB','main_all')
   i_all=-product(shape(t_params%RCLS))*kind(t_params%RCLS)
   deallocate(t_params%RCLS,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%RCLS','main_all')
   i_all=-product(shape(t_params%VINS))*kind(t_params%VINS)
   deallocate(t_params%VINS,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%VINS','main_all')
   i_all=-product(shape(t_params%VISP))*kind(t_params%VISP)
   deallocate(t_params%VISP,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%VISP','main_all')
   i_all=-product(shape(t_params%CSCL))*kind(t_params%CSCL)
   deallocate(t_params%CSCL,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%CSCL','main_all')
   i_all=-product(shape(t_params%DRDI))*kind(t_params%DRDI)
   deallocate(t_params%DRDI,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%DRDI','main_all')
   i_all=-product(shape(t_params%NSH1))*kind(t_params%NSH1)
   deallocate(t_params%NSH1,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%NSH1','main_all')
   i_all=-product(shape(t_params%NSH2))*kind(t_params%NSH2)
   deallocate(t_params%NSH2,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%NSH2','main_all')
   i_all=-product(shape(t_params%CREL))*kind(t_params%CREL)
   deallocate(t_params%CREL,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%CREL','main_all')
   i_all=-product(shape(t_params%RREL))*kind(t_params%RREL)
   deallocate(t_params%RREL,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%RREL','main_all')
   i_all=-product(shape(t_params%IRNS))*kind(t_params%IRNS)
   deallocate(t_params%IRNS,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%IRNS','main_all')
   i_all=-product(shape(t_params%LMSP))*kind(t_params%LMSP)
   deallocate(t_params%LMSP,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%LMSP','main_all')
   i_all=-product(shape(t_params%OPTC))*kind(t_params%OPTC)
   deallocate(t_params%OPTC,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%OPTC','main_all')
   i_all=-product(shape(t_params%BZKP))*kind(t_params%BZKP)
   deallocate(t_params%BZKP,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%BZKP','main_all')
   i_all=-product(shape(t_params%IRWS))*kind(t_params%IRWS)
   deallocate(t_params%IRWS,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%IRWS','main_all')
   i_all=-product(shape(t_params%ECORE))*kind(t_params%ECORE)
   deallocate(t_params%ECORE,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%ECORE','main_all')
   i_all=-product(shape(t_params%QMTET))*kind(t_params%QMTET)
   deallocate(t_params%QMTET,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%QMTET','main_all')
   i_all=-product(shape(t_params%QMPHI))*kind(t_params%QMPHI)
   deallocate(t_params%QMPHI,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%QMPHI','main_all')
   i_all=-product(shape(t_params%RMREL))*kind(t_params%RMREL)
   deallocate(t_params%RMREL,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%RMREL','main_all')
   i_all=-product(shape(t_params%VTREL))*kind(t_params%VTREL)
   deallocate(t_params%VTREL,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%VTREL','main_all')
   i_all=-product(shape(t_params%BTREL))*kind(t_params%BTREL)
   deallocate(t_params%BTREL,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%BTREL','main_all')
   i_all=-product(shape(t_params%SRREL))*kind(t_params%SRREL)
   deallocate(t_params%SRREL,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%SRREL','main_all')
   i_all=-product(shape(t_params%DROTQ))*kind(t_params%DROTQ)
   deallocate(t_params%DROTQ,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%DROTQ','main_all')
   i_all=-product(shape(t_params%RATOM))*kind(t_params%RATOM)
   deallocate(t_params%RATOM,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%RATOM','main_all')
   i_all=-product(shape(t_params%LCORE))*kind(t_params%LCORE)
   deallocate(t_params%LCORE,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%LCORE','main_all')
   i_all=-product(shape(t_params%NCORE))*kind(t_params%NCORE)
   deallocate(t_params%NCORE,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%NCORE','main_all')
   i_all=-product(shape(t_params%IRCUT))*kind(t_params%IRCUT)
   deallocate(t_params%IRCUT,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%IRCUT','main_all')
   i_all=-product(shape(t_params%ICLEB))*kind(t_params%ICLEB)
   deallocate(t_params%ICLEB,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%ICLEB','main_all')
   i_all=-product(shape(t_params%NACLS))*kind(t_params%NACLS)
   deallocate(t_params%NACLS,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%NACLS','main_all')
   i_all=-product(shape(t_params%LOFLM))*kind(t_params%LOFLM)
   deallocate(t_params%LOFLM,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%LOFLM','main_all')
   i_all=-product(shape(t_params%KAOEZ))*kind(t_params%KAOEZ)
   deallocate(t_params%KAOEZ,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%KAOEZ','main_all')
   i_all=-product(shape(t_params%ULDAU))*kind(t_params%ULDAU)
   deallocate(t_params%ULDAU,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%ULDAU','main_all')
   i_all=-product(shape(t_params%WLDAU))*kind(t_params%WLDAU)
   deallocate(t_params%WLDAU,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%WLDAU','main_all')
   i_all=-product(shape(t_params%KMESH))*kind(t_params%KMESH)
   deallocate(t_params%KMESH,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%KMESH','main_all')
   i_all=-product(shape(t_params%R2NEF))*kind(t_params%R2NEF)
   deallocate(t_params%R2NEF,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%R2NEF','main_all')
   i_all=-product(shape(t_params%MVEVI))*kind(t_params%MVEVI)
   deallocate(t_params%MVEVI,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%MVEVI','main_all')
   i_all=-product(shape(t_params%IRREL))*kind(t_params%IRREL)
   deallocate(t_params%IRREL,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%IRREL','main_all')
   i_all=-product(shape(t_params%NRREL))*kind(t_params%NRREL)
   deallocate(t_params%NRREL,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%NRREL','main_all')
   i_all=-product(shape(t_params%LMSP1))*kind(t_params%LMSP1)
   deallocate(t_params%LMSP1,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%LMSP1','main_all')
   i_all=-product(shape(t_params%IFUNM))*kind(t_params%IFUNM)
   deallocate(t_params%IFUNM,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%IFUNM','main_all')
   i_all=-product(shape(t_params%LLMSP))*kind(t_params%LLMSP)
   deallocate(t_params%LLMSP,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%LLSMP','main_all')
   i_all=-product(shape(t_params%IRMIN))*kind(t_params%IRMIN)
   deallocate(t_params%IRMIN,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%IRMIN','main_all')
   i_all=-product(shape(t_params%TESTC))*kind(t_params%TESTC)
   deallocate(t_params%TESTC,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%TESTC','main_all')
   i_all=-product(shape(t_params%VOLBZ))*kind(t_params%VOLBZ)
   deallocate(t_params%VOLBZ,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%VOLBZ','main_all')
   i_all=-product(shape(t_params%NOFKS))*kind(t_params%NOFKS)
   deallocate(t_params%NOFKS,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%NOFKS','main_all')
   i_all=-product(shape(t_params%THETA))*kind(t_params%THETA)
   deallocate(t_params%THETA,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%THETA','main_all')
   i_all=-product(shape(t_params%DSYMLL))*kind(t_params%DSYMLL)
   deallocate(t_params%DSYMLL,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%DSYMLL','main_all')
   i_all=-product(shape(t_params%ITLDAU))*kind(t_params%ITLADU)
   deallocate(t_params%ITLDAU,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%ITLDAU','main_all')
   i_all=-product(shape(t_params%JWSREL))*kind(t_params%JWSREL)
   deallocate(t_params%JWSREL,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%JWSREL','main_all')
   i_all=-product(shape(t_params%THETAS))*kind(t_params%THETAS)
   deallocate(t_params%THETAS,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%THETAS','main_all')
   i_all=-product(shape(t_params%RMTNEW))*kind(t_params%RMTNEW)
   deallocate(t_params%RMTNEW,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%RMTNEW','main_all')
   i_all=-product(shape(t_params%RHOORB))*kind(t_params%RHOORB)
   deallocate(t_params%RHOORB,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%RHOORB','main_all')
   i_all=-product(shape(t_params%NSHELL))*kind(t_params%NSHELL)
   deallocate(t_params%NSHELL,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%NSHELL','main_all')
   i_all=-product(shape(t_params%RHO2NS))*kind(t_params%RHO2NS)
   deallocate(t_params%RHO2NS,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%RHO2NS','main_all')
   i_all=-product(shape(t_params%RMTREF))*kind(t_params%RMTREF)
   deallocate(t_params%RMTREF,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%RMTREF','main_all')
   i_all=-product(shape(t_params%SOCSCL))*kind(t_params%SOCSCL)
   deallocate(t_params%SOCSCL,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%SOCSCL','main_all')
   i_all=-product(shape(t_params%RBASIS))*kind(t_params%RBASIS)
   deallocate(t_params%RBASIS,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%RBASIS','main_all')
   i_all=-product(shape(t_params%REFPOT))*kind(t_params%REFPOT)
   deallocate(t_params%REFPOT,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%REFPOT','main_all')
   i_all=-product(shape(t_params%IFUNM1))*kind(t_params%IFUNM1)
   deallocate(t_params%IFUNM1,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%IFUNM1','main_all')
   i_all=-product(shape(t_params%ITITLE))*kind(t_params%ITITLE)
   deallocate(t_params%ITITLE,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%ITITLE','main_all')
   i_all=-product(shape(t_params%NTCELL))*kind(t_params%NTCELL)
   deallocate(t_params%NTCELL,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%NTCELL','main_all')
   i_all=-product(shape(t_params%IXIPOL))*kind(t_params%IXIPOL)
   deallocate(t_params%IXIPOL,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%IXIPOL','main_all')
   i_all=-product(shape(t_params%IMAXSH))*kind(t_params%IMAXSH)
   deallocate(t_params%IMAXSH,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%IMAXSH','main_all')
   i_all=-product(shape(t_params%NKCORE))*kind(t_params%NKCORE)
   deallocate(t_params%NKCORE,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%NKCORE','main_all')
   i_all=-product(shape(t_params%VOLCUB))*kind(t_params%VOLCUB)
   deallocate(t_params%VOLCUB,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%VOLCUB','main_all')
   i_all=-product(shape(t_params%IQCALC))*kind(t_params%IQCALC)
   deallocate(t_params%IQCALC,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%IQCALC','main_all')
   i_all=-product(shape(t_params%ICHECK))*kind(t_params%ICHECK)
   deallocate(t_params%ICHECK,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%ICHECK','main_all')
   i_all=-product(shape(t_params%PHILDAU))*kind(t_params%PHILDAU)
   deallocate(t_params%PHILDAU,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%PHILDAU','main_all')
   i_all=-product(shape(t_params%DRDIREL))*kind(t_params%DRDIREL)
   deallocate(t_params%DRDIREL,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%DRDIREL','main_all')
   i_all=-product(shape(t_params%DENEFAT))*kind(t_params%DENEFAT)
   deallocate(t_params%DENEFAT,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%DENEFAT','main_all')
   i_all=-product(shape(t_params%RCLSIMP))*kind(t_params%RCLSIMP)
   deallocate(t_params%RCLSIMP,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%RCLSIMP','main_all')
   i_all=-product(shape(t_params%MVEVIEF))*kind(t_params%MVEVIEF)
   deallocate(t_params%MVEVIEF,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%MVEVIEF','main_all')
   i_all=-product(shape(t_params%IRSHIFT))*kind(t_params%IRSHIFT)
   deallocate(t_params%IRSHIFT,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%IRSHIFT','main_all')
   i_all=-product(shape(t_params%IJTABSH))*kind(t_params%IJTABSH)
   deallocate(t_params%IJTABSH,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%IJTABSH','main_all')
   i_all=-product(shape(t_params%ATOMIMP))*kind(t_params%ATOMIMP)
   deallocate(t_params%ATOMIMP,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%ATOMIMP','main_all')
   i_all=-product(shape(t_params%VACFLAG))*kind(t_params%VACFLAG)
   deallocate(t_params%VACFLAG,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%VACFLAG','main_all')
   i_all=-product(shape(t_params%KAPCORE))*kind(t_params%KAPCORE)
   deallocate(t_params%KAPCORE,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%KAPCORE','main_all')
   i_all=-product(shape(t_params%HOSTIMP))*kind(t_params%HOSTIMP)
   deallocate(t_params%HOSTIMP,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%HOSTIMP','main_all')
   i_all=-product(shape(t_params%ECOREREL))*kind(t_params%ECOREREL)
   deallocate(t_params%ECOREREL,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%ECOREREL','main_all')
   i_all=-product(shape(t_params%IJTABSYM))*kind(t_params%IJTABSYM)
   deallocate(t_params%IJTABSYM,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%IJTABSYM','main_all')
   i_all=-product(shape(t_params%EREFLDAU))*kind(t_params%EREFLDAU)
   deallocate(t_params%EREFLDAU,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%EREFLDAU','main_all')
   i_all=-product(shape(t_params%QMPHITAB))*kind(t_params%QMPHITAB)
   deallocate(t_params%QMPHITAB,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%QMPHITAB','main_all')
   i_all=-product(shape(t_params%QMTETTAB))*kind(t_params%QMTETTAB)
   deallocate(t_params%QMTETTAB,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%QMTETTAB','main_all')
   i_all=-product(shape(t_params%QMGAMTAB))*kind(t_params%QMGAMTAB)
   deallocate(t_params%QMGAMTAB,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%QMGAMTAB','main_all')
   i_all=-product(shape(t_params%CMOMHOST))*kind(t_params%CMOMHOST)
   deallocate(t_params%CMOMHOST,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%CMOMHOST','main_all')
   i_all=-product(shape(t_params%SOCSCALE))*kind(t_params%SOCSCALE)
   deallocate(t_params%SOCSCALE,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%SOCSCALE','main_all')
   i_all=-product(shape(t_params%NPAN_TOT))*kind(t_params%NPAN_TOT)
   deallocate(t_params%NPAN_TOT,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%NPAN_TOT','main_all')
   i_all=-product(shape(t_params%THETASNEW))*kind(t_params%THETASNEW)
   deallocate(t_params%THETASNEW,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%THETASNEW','main_all')
   i_all=-product(shape(t_params%IJTABCALC))*kind(t_params%IJTABCALC)
   deallocate(t_params%IJTABCALC,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%IJTABCALC','main_all')
   i_all=-product(shape(t_params%R2DRDIREL))*kind(t_params%R2DRDIREL)
   deallocate(t_params%R2DRDIREL,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%R2DRDIREL','main_all')
   i_all=-product(shape(t_params%NPAN_EQNEW))*kind(t_params%NPAN_EQNEW)
   deallocate(t_params%NPAN_EQNEW,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%NPAN_EQNEW','main_all')
   i_all=-product(shape(t_params%SYMUNITARY))*kind(t_params%SYMUNITARY)
   deallocate(t_params%SYMUNITARY,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%SYMUNITARY','main_all')
   i_all=-product(shape(t_params%LEFTTINVLL))*kind(t_params%LEFTTINVLL)
   deallocate(t_params%LEFTTINVLL,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%LEFTTINVLL','main_all')
   i_all=-product(shape(t_params%RIGHTTINVLL))*kind(t_params%RIGHTTINVLL)
   deallocate(t_params%RIGHTTINVLL,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%RIGHTTTINVLL','main_all')
   i_all=-product(shape(t_params%NPAN_LOGNEW))*kind(t_params%NPAN_LOGNEW)
   deallocate(t_params%NPAN_LOGNEW,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%NPAN_LOGNEW','main_all')
   i_all=-product(shape(t_params%IJTABCALC_I))*kind(t_params%IJTABCALC_I)
   deallocate(t_params%IJTABCALC_I,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%IJTABCALC_I','main_all')
   i_all=-product(shape(t_params%IPAN_INTERVALL))*kind(t_params%IPAN_INTERVALL)
   deallocate(t_params%IPAN_INTERVALL,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%IPAN_INTERVALL','main_all')
   i_all=-product(shape(t_params%RPAN_INTERVALL))*kind(t_params%RPAN_INTERVALL)
   deallocate(t_params%RPAN_INTERVALL,stat=i_stat)
   call memocc(i_stat,i_all,'t_params%RPAN_INTERVALL','main_all')

   !delete temporary files
   if(myrank==master) then
      open(69, file='abvmad.unformatted')
      close(69, status='delete')
   end if

   ! deallocate arrays from t_wavefunctions
   if(t_wavefunctions%Nwfsavemax>0) then

      i_all=-product(shape(t_wavefunctions%isave_wavefun))*kind(t_wavefunctions%isave_wavefun)
      deallocate(t_wavefunctions%isave_wavefun, stat=i_stat)
      call memocc(i_stat,i_all,'t_wavefunctions%isave_wavefun','main_all')

      if(t_wavefunctions%save_rll) then
         i_all=-product(shape(t_wavefunctions%rll))*kind(t_wavefunctions%rll)
         deallocate(t_wavefunctions%rll, stat=i_stat)
         call memocc(i_stat,i_all,'t_wavefunctions%rll','main_all')
      endif

      if(t_wavefunctions%save_sll) then
         i_all=-product(shape(t_wavefunctions%sll))*kind(t_wavefunctions%sll)
         deallocate(t_wavefunctions%sll, stat=i_stat)
         call memocc(i_stat,i_all,'t_wavefunctions%sll','main_all')
      endif

      if(t_wavefunctions%save_rllleft) then
         i_all=-product(shape(t_wavefunctions%rllleft))*kind(t_wavefunctions%rllleft)
         deallocate(t_wavefunctions%rllleft, stat=i_stat)
         call memocc(i_stat,i_all,'t_wavefunctions%rllleft','main_all')
      endif

      if(t_wavefunctions%save_sllleft) then
         i_all=-product(shape(t_wavefunctions%sllleft))*kind(t_wavefunctions%sllleft)
         deallocate(t_wavefunctions%sllleft, stat=i_stat)
         call memocc(i_stat,i_all,'t_wavefunctions%sllleft','main_all')
      endif

   end if !(t_wavefunctions%Nwfsavemax>0)

#ifdef CPP_MPI
   ! deallocate arrays for MPIadapt
   if(MPIadapt>0) then
      i_all=-product(shape(timings_1a))*kind(timings_1a)
      deallocate(timings_1a,stat=i_stat)
      call memocc(i_stat,i_all,'timings_1a','main_all')

      i_all=-product(shape(timings_1b))*kind(timings_1b)
      deallocate(timings_1b,stat=i_stat)
      call memocc(i_stat,i_all,'timings_1b','main_all')

      i_all=-product(shape(load_imbalance))*kind(load_imbalance)
      deallocate(load_imbalance,stat=i_stat)
      call memocc(i_stat,i_all,'load_imbalance','main_all')

   end if
#endif

   ! Deallocation of input arrays
   call allocate_cell(-1,NAEZ,NEMB,NATYP,CLS,IMT,IRWS,IRNS,NTCELL,REFPOT,&
      KFG,KAOEZ,RMT,ZAT,RWS,MTFAC,RMTREF,RMTREFAT,RMTNEW,RBASIS)
   call allocate_semi_inf_host(-1,NEMB,TLEFT,TRIGHT)
   call allocate_potential(-1,NAEZ,NEMB,IRM,NATYP,NPOTD,IPAND,NFUND,LMXSPD,&
      LMPOT,IRMIND,NSPOTD,NFU,IRC,LMXC,NCORE,IRMIN,LMSP,LMSP1,IRCUT,LCORE,LLMSP,&
      ITITLE,FPRADIUS,VISP,ECORE,VINS)
   call allocate_cpa(-1,NAEZ,NEMB,NATYP,NOQ,ICPA,IQAT,HOSTIMP,CONC)
   call allocate_ldau(-1,NATYP,LOPT,UEFF,JEFF,EREFLDAU)
   call allocate_ldau_potential(-1,IRM,NATYP,MMAXD,NSPIND,ITLDAU,WLDAU,&
      ULDAU,PHILDAU)
   call allocate_magnetization(-1,NAEZ,NATYP,LMMAXD,INIPOL,IXIPOL,QMTET,&
      QMPHI,DROTQ)
   call allocate_SOC(-1,KREL,NATYP,LMAX,IMANSOC,SOCSCALE,CSCL,SOCSCL)
   call allocate_energies(-1,IEMXD,EZ,DEZ,WEZ)
   call allocate_relativistic(-1,KREL,IRM,NAEZ,NATYP,ZREL,JWSREL,IRSHIFT,&
      VTREL,BTREL,RMREL,DRDIREL,R2DRDIREL,QMGAM,QMGAMTAB,QMPHITAB,QMTETTAB)
   call allocate_rel_transformations(-1,LMMAXD,NRREL,IRREL,RC,CREL,RREL,SRREL)
   call allocate_clusters(-1,NAEZ,LMAX,NCLEB,NCLSD,NEMBD1,NSHELD,NACLSD,&
      LMPOT,NATOMIMPD,NSH1,NSH2,NACLS,NSHELL,ATOMIMP,ATOM,EZOA,ICLEB,JEND,RATOM,&
      RCLSIMP,CMOMHOST,RCLS)
   call allocate_expansion(-1,LM2D,IRID,NFUND,NTOTD,NCLEB,LASSLD,NCELLD,&
      NCHEBD,LOFLM,WG,CLEB,YRG,THETAS,THETASNEW)
   call allocate_mesh(-1,IRM,NATYP,A,B,R,DRDI)
   call allocate_pannels(-1,NATYP,NTOTD,IPAN,NPAN_TOT,NPAN_EQNEW,NPAN_LOGNEW,&
      IPAN_INTERVALL,RPAN_INTERVALL)
   call allocate_misc(-1,NR,IRM,IRID,LMAX,NAEZ,NATYP,NFUND,NREFD,IEMXD,&
      NTOTD,NSHELD,LMMAXD,NEMBD1,NCHEBD,NCELLD,LMXSPD,NSPINDD,NSYMAXD,NPRINCD,IFUNM,&
      IFUNM1,ICHECK,VREF,S,RR,DROR,RNEW,RS,RROT,THESME,DSYMLL,DSYMLL1,LEFTTINVLL,&
      RIGHTTINVLL)
   call allocate_green(-1,NAEZ,IEMXD,NGSHD,NSHELD,LMPOT,NOFGIJD,ISH,JSH,&
      KMESH,IMAXSH,IQCALC,IOFGIJ,JOFGIJ,IJTABSH,IJTABSYM,IJTABCALC,IJTABCALC_I,ILM,GSH)
   ! End of deallocation

#ifdef CPP_MPI
   ! finalize MPI
   call MPI_Finalize(ierr)
#endif
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< close allocated arrays and finalize MPI !!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



end program
