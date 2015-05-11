program kkrcode

  use mod_main0
  use mod_main1a
  use mod_main1b
  use mod_main1c
  use mod_main2
  use mod_types
  use mod_timing

#ifdef CPP_MPI
  use mod_mympi, only: mympi_init, myrank, nranks, master,find_dims_2d, distribute_linear_on_tasks, create_subarr_comm 
#else
  use mod_mympi, only: mympi_init, myrank, nranks, master,
#endif

  
#ifdef CPP_MPI
  use mpi
#endif

  implicit none
  
  integer :: ierr
  
#ifdef CPP_MPI
  integer :: myMPI_comm_grid, myMPI_comm_row, myMPI_comm_col, myrank_grid, myrank_row, myrank_col, nranks_row, nranks_col
  integer :: ntot1, dims(2)
#endif

  integer :: ie,irec,i1,ispin
  
  
#ifdef CPP_MPI
  ! initialize MPI
  call MPI_Init ( ierr )
#endif
  ! set variables master, myrank and nranks for serial (myrank=master=0, nranks=1) as well as parallel execution
  call mympi_init()

  ! initialize timing
  call timing_init(myrank)
    
  ! start KKR program, first do main0, where everything is read in and initialized
  call timing_start('main0')
  if(myrank==master) call main0()
  call timing_stop('main0')
  
  
#ifdef CPP_MPI
  !now communicate type t_inc and t_tgmat switches (this has an implicit barrier, so that all other processes wait for master to finish with main0)
  call timing_start('MPI 1')
  call bcast_t_inc_tgmat(t_inc,t_tgmat)
  !also communicate logicals from t_lloyd
  call bcast_t_lly_1(t_inc,t_lloyd)
  
  
  ! create 2d matrix for processors so that more processors than energy points can be used. 
  ! strategy here is to first parallelize the energy points ideally and then give all processors that are left to atom or k-point loop parallelization
  
  ! create_subcomms_2d:
  ntot1 = t_inc%NATYP
  call find_dims_2d(nranks,ntot1,t_inc%IELAST,dims)
  
  ! create cartesian communicator for atom/energy matrix
  call create_subarr_comm( dims, myMPI_comm_grid, myMPI_comm_row, myMPI_comm_col, myrank_grid, myrank_row, myrank_col, nranks_row, nranks_col)
  ! save grid info in type
  call save_t_mpi_c_grid(t_mpi_c_grid,dims, myMPI_comm_grid, myMPI_comm_row, myMPI_comm_col, myrank_grid, myrank_row, myrank_col, nranks_row, nranks_col)
  call timing_stop('MPI 1')
  
!     STOP
#endif
  
  
  ! Then start scf iterations and do all steps of the KKR formalism until convergence
  do while ( (t_inc%i_iteration.lt.t_inc%N_iteration) .and. (t_inc%N_iteration.ne.-1) )
  
  
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
    
    ! calculate do DFT stuff (potential from density, exc-potential, calculate total energy, ...)
    call timing_start('main2')
    if (myrank==master) call main2()
    call timing_stop('main2')
    
#ifdef CPP_MPI
    ! update i_iteration:
    call MPI_Bcast(t_inc%i_iteration, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting i_iteration in main_all'
#endif

    call timing_stop('Time in Iteration')
    
  end do ! scf-iteration
  
  
#ifdef CPP_MPI
  ! finalize MPI  
  call MPI_Finalize(ierr)
#endif


end program
