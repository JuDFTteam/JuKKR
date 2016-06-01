program kkrcode

  use mod_main0
  use mod_main1a
  use mod_main1b
  use mod_main1c
  use mod_main2
  use mod_types
  use mod_timing

#ifdef CPP_MPI
  use mod_mympi, only: mympi_init, myrank, nranks, master,find_dims_2d, distribute_linear_on_tasks, create_newcomms_group_ie 
#else
  use mod_mympi, only: mympi_init, myrank, nranks, master
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

  integer :: ie,irec,i1,ispin
  character(len=3) :: ctemp

 
    
#ifdef CPP_MPI
  ! initialize MPI
  call MPI_Init ( ierr )
#endif
  ! set variables master, myrank and nranks for serial (myrank=master=0, nranks=1) as well as parallel execution
  call mympi_init()
  ! save myrank in ctemp, needed to open output unit 1337
  write(ctemp,'(I03.3)') myrank

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
     ! default value on master (needed for writeout in main0)
     t_inc%i_write = 1

     ! run main0 only serially, then communicate everything that is in here
     call main0()
     call timing_stop('main0')

  end if ! myrank ==master
  
  

! without MPI (serial or openMP) something goes wrong if if files are not written out
! this seems to be only the case with the old solver
#ifndef CPP_MPI
  if(.not.t_inc%NEWSOSOL) then
    t_tgmat%tmat_to_file = .true.
    t_tgmat%gmat_to_file = .true.
    t_tgmat%gref_to_file = .true.
  end if
#endif

#ifdef CPP_MPI
  ! now communicate type t_inc and t_tgmat switches (this has an implicit barrier, so that all other processes wait for master to finish with main0)
  if(myrank==master) call timing_start('MPI 1')
  call bcast_t_inc_tgmat(t_inc,t_tgmat,t_cpa)
  ! also communicate logicals from t_lloyd
  call bcast_t_lly_1(t_inc,t_lloyd)
  
  call bcast_t_params_scalars(t_params)
  if (myrank.ne.master) call init_t_params(t_params)
  call bcast_t_params_arrays(t_params)
  
  ! in case of deci-out run everything only serially to write decifile correctly. This will be fixed in a later version when we get rid of the decifile
  if(t_inc%deci_out .and. nranks>1) stop 'deci-out option chosen. Please run code serially!'
  

  ! create 2d matrix for processors so that more processors than energy points can be used. 
  ! strategy here is to first parallelize the energy points ideally and then give all processors that are left to atom or k-point loop parallelization
  
  ! create_subcomms_2d: first find maximal dimensions
  call find_dims_2d(nranks,t_inc%NATYP,t_inc%IELAST,dims)
  !save in dims
  t_mpi_c_grid%dims = dims

  ! create communicator for atom/energy matrix
  call create_newcomms_group_ie( nranks,myrank,dims(1),dims(2),t_inc%nkmesh,t_inc%kmesh,mympi_comm_ie,  &
                                 & myrank_ie,nranks_ie,mympi_comm_at,myrank_at,nranks_at, myrank_atcomm,nranks_atcomm)
  ! save grid info in type 't_mpi_c_grid'
  call save_t_mpi_c_grid(t_mpi_c_grid,dims, myMPI_comm_ie, myMPI_comm_at, myrank_ie, myrank_at, myrank_atcomm, nranks_ie, nranks_at, nranks_atcomm)
  if(myrank==master) call timing_stop('MPI 1')
#endif

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
       if(t_inc%i_write>0) open(1337, file='output.'//trim(ctemp)//'.txt')
    endif
    if(t_inc%i_write>0 .and. myrank.ne.master) open(1337, file='output.'//trim(ctemp)//'.txt')


  
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
    if (t_inc%i_write<2 .and. t_inc%i_write>0 .and. myrank==master) call SYSTEM('cp output.000.txt output.2.txt')
    ! rewind output.xxx.txt
    if (t_inc%i_write<2 .and. t_inc%i_write>0) rewind(1337)
    ! rewind timing files if t_inc%i_time<2
    if (t_inc%i_time<2 .and. t_inc%i_time>0) rewind(43234059) 
  
    call timing_start('Time in Iteration')
  
    ! calculate tmat and gref
    call timing_start('main1a')
    call main1a()
    call timing_stop('main1a')

    ! calculate gmat
    call timing_start('main1b')
    call main1b()
    call timing_stop('main1b')

    ! calculate density
    call timing_start('main1c')
    call main1c()
    call timing_stop('main1c')
    
    ! calculate DFT stuff (potential from density, exc-potential, calculate total energy, ...)
    call timing_start('main2')
    if (myrank==master) call main2()
    call timing_stop('main2')
    
    
#ifdef CPP_MPI
    ! update i_iteration after check for convergence in main2:
    call MPI_Bcast(t_inc%i_iteration, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error broadcasting i_iteration in main_all'
    ! broadcast parameter arrays from master (e.g. update nonco angles etc.)
    call bcast_t_params_scalars(t_params)
    call bcast_t_params_arrays(t_params)
#endif

    call timing_stop('Time in Iteration')
    
    if (myrank==master) call print_time_and_date('Iteration finished')
    
  end do ! scf-iteration
  
  
#ifdef CPP_MPI
  ! finalize MPI  
  call MPI_Finalize(ierr)
#endif


end program
