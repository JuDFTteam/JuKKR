program kkrcode

use mod_main0
use mod_main1a
use mod_main1b
use mod_main1c
use mod_main2
use mod_types
use mod_timing
use mod_version_info
use mod_md5sums

#ifdef CPP_MPI
 use mod_mympi, only: mympi_init, myrank, nranks, master,find_dims_2d, distribute_linear_on_tasks, &
                     &create_newcomms_group_ie, MPIatom , MPIadapt, check_communication_pattern
 use mod_save_wavefun, only: t_wavefunctions, bcast_params_savewf
#else
 use mod_mympi, only: mympi_init, myrank, nranks, master, MPIatom, MPIadapt
 use mod_save_wavefun, only: t_wavefunctions
#endif

#ifdef CPP_MPI
use mpi
#endif

implicit none

integer :: ierr
#ifdef CPP_MPI
integer :: myMPI_comm_at, myMPI_comm_ie, nranks_at, myrank_at, myrank_ie, nranks_ie, nranks_atcomm, myrank_atcomm
integer :: dims(2)
#endif
character(len=3) :: ctemp !name for output file
! needed to use test('xxxxxxxx'):
logical :: test, opt
external :: test, opt



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
   allocate(timings_1a(t_inc%ielast,t_inc%natyp), timings_1b(t_inc%ielast), load_imbalance(t_inc%nkmesh), stat=ierr)
   if(ierr/=0) stop '[main_all] Error allocating timing_1a, 1b arrays'
   timings_1a(:,:) = 0.0d0
   timings_1b(:) = 0.0d0
   load_imbalance(:) = 0
end if
! create_subcomms_2d: first find maximal dimensions
call find_dims_2d(nranks,t_inc%NATYP,t_inc%IELAST,dims,MPIatom)
!save in dims
t_mpi_c_grid%dims = dims

! create communicator for atom/energy matrix
call create_newcomms_group_ie( nranks,myrank,dims(1),dims(2),t_inc%nkmesh,t_inc%kmesh,mympi_comm_ie,  &
                               & myrank_ie,nranks_ie,mympi_comm_at,myrank_at,nranks_at, myrank_atcomm,nranks_atcomm)

if(t_inc%i_write>0 .and. myrank.ne.master) write(1337,*) 'after create_newcomms_ie', dims, myrank_ie, myrank_at, myrank_atcomm, nranks_ie, nranks_at, nranks_atcomm
! save grid info in type 't_mpi_c_grid'
call save_t_mpi_c_grid(t_mpi_c_grid,dims, myMPI_comm_ie, myMPI_comm_at, myrank_ie, myrank_at, myrank_atcomm, nranks_ie, nranks_at, nranks_atcomm)
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

  ! calculate tmat and gref
  call timing_start('main1a')
  call main1a()
  call timing_stop('main1a')

  ! calculate gmat
  call timing_start('main1b')
  call main1b()
  call timing_stop('main1b')
  if(test('STOP1B  '))then
#ifdef CPP_MPI
    call MPI_Finalize(ierr)
#endif
    if(.not. OPT('WRTGREEN')) write(*,*) 'done with WRTGREEN step'
    if(myrank==master) write(*,*) 'Stop after main1b'
    stop
  end if!test

  ! calculate density
  call timing_start('main1c')
  call main1c()
  call timing_stop('main1c')
  
  ! calculate DFT stuff (potential from density, exc-potential, calculate total energy, ...)
  call timing_start('main2')
  if (myrank==master) call main2()
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
                                   & myrank_ie,nranks_ie,mympi_comm_at,myrank_at,nranks_at, myrank_atcomm,nranks_atcomm)
    ! save grid info in type 't_mpi_c_grid'
    call save_t_mpi_c_grid(t_mpi_c_grid,dims, myMPI_comm_ie, myMPI_comm_at, myrank_ie, myrank_at, myrank_atcomm, nranks_ie, nranks_at, nranks_atcomm)
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
deallocate(t_params%EZ, t_params%WEZ, t_params%DROTQ, t_params%DSYMLL, t_params%LEFTTINVLL, t_params%RIGHTTINVLL, t_params%CREL, t_params%RC, t_params%RREL, t_params%SRREL, t_params%PHILDAU, t_params%VINS, t_params%VISP, t_params%VBC, t_params%VTREL, t_params%BTREL, t_params%SOCSCALE, t_params%DRDIREL, t_params%R2DRDIREL, t_params%RMREL, t_params%CMOMHOST       , t_params%ECORE, t_params%QMTET, t_params%QMPHI, t_params%QMPHITAB, t_params%QMTETTAB, t_params%QMGAMTAB, t_params%ZAT, t_params%R, t_params%DRDI, t_params%RMTREF, t_params%VREF, t_params%CLEB, t_params%RCLS, t_params%SOCSCL, t_params%CSCL, t_params%RBASIS, t_params%RR, t_params%CONC, t_params%RROT, t_params%RATOM, t_params%A, t_params%B, t_params%THETAS, t_params%RMT, t_params%RMTNEW, t_params%RWS, t_params%GSH, t_params%EREFLDAU, t_params%UEFF, t_params%JEFF, t_params%ULDAU, t_params%WLDAU, t_params%RPAN_INTERVALL, t_params%RNEW, t_params%MVEVI, t_params%MVEVIEF, t_params%THETASNEW, t_params%RHO2NS, t_params%R2NEF, t_params%RHOC, t_params%DENEFAT, t_params%ESPV, t_params%EDC, t_params%EU, t_params%RHOORB, t_params%ECOREREL, t_params%RCLSIMP, t_params%LOPT, t_params%ITLDAU , t_params%IRSHIFT, t_params%JWSREL , t_params%ZREL, t_params%LCORE, t_params%NCORE, t_params%IPAN , t_params%IRCUT, t_params%JEND , t_params%ICLEB, t_params%ATOM, t_params%CLS , t_params%NACLS, t_params%LOFLM, t_params%EZOA , t_params%KAOEZ, t_params%IQAT, t_params%ICPA, t_params%NOQ , t_params%KMESH , t_params%NSHELL, t_params%NSH1, t_params%NSH2, t_params%IJTABCALC, t_params%IJTABCALC_I, t_params%IJTABSYM, t_params%IJTABSH, t_params%ISH, t_params%JSH, t_params%IQCALC, t_params%ICHECK, t_params%ATOMIMP, t_params%REFPOT, t_params%IRREL, t_params%NRREL, t_params%IFUNM1, t_params%ITITLE, t_params%LMSP1, t_params%NTCELL, t_params%IXIPOL, t_params%IRNS  , t_params%IFUNM , t_params%LLMSP , t_params%LMSP, t_params%IMT , t_params%IRC , t_params%IRMIN, t_params%IRWS , t_params%NFU  , t_params%HOSTIMP, t_params%ILM    , t_params%IMAXSH , t_params%NPAN_LOGNEW, t_params%NPAN_EQNEW , t_params%NPAN_TOT, t_params%IPAN_INTERVALL, t_params%NKCORE , t_params%KAPCORE, t_params%SYMUNITARY, t_params%VACFLAG, t_params%TXC, t_params%TESTC, t_params%OPTC, t_params%BZKP, t_params%VOLCUB, t_params%VOLBZ, t_params%NOFKS, t_params%THETA, t_params%PHI, stat=ierr)
if(ierr/=0) stop '[main_all] Error deallocating arrays from t_params'

!delete temporary files
if(myrank==master) then
  open(69, file='abvmad.unformatted')
  close(69, status='delete')
end if

! deallocate arrays from t_wavefunctions
if(t_wavefunctions%Nwfsavemax>0) then

 deallocate(t_wavefunctions%isave_wavefun, stat=ierr)
 
 if(ierr/=0) stop '[main_all] Error deallocating arrays 1 from t_wavefunctions'
 if(t_wavefunctions%save_rll) then
   deallocate(t_wavefunctions%rll, stat=ierr)
   if(ierr/=0) stop '[main_all] Error deallocating arrays 2 from t_wavefunctions'
 endif
 
 if(t_wavefunctions%save_sll) then
   deallocate(t_wavefunctions%sll, stat=ierr)
   if(ierr/=0) stop '[main_all] Error deallocating arrays 3 from t_wavefunctions'
 endif
 
 if(t_wavefunctions%save_rllleft) then
   deallocate(t_wavefunctions%rllleft, stat=ierr)
   if(ierr/=0) stop '[main_all] Error deallocating arrays 4 from t_wavefunctions'
 endif
 
 if(t_wavefunctions%save_sllleft) then
   deallocate(t_wavefunctions%sllleft, stat=ierr)
   if(ierr/=0) stop '[main_all] Error deallocating arrays 5 from t_wavefunctions'
 endif
 
end if !(t_wavefunctions%Nwfsavemax>0)

#ifdef CPP_MPI
! deallocate arrays for MPIadapt
if(MPIadapt>0) then
   deallocate(timings_1a, timings_1b, load_imbalance, stat=ierr)
   if(ierr/=0) stop '[main_all] Error deallocating timing_1a, 1b arrays'
end if
#endif


#ifdef CPP_MPI
! finalize MPI  
call MPI_Finalize(ierr)
#endif
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< close allocated arrays and finalize MPI !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



end program
