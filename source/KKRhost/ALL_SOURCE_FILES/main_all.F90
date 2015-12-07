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

  
  
!       integer,allocatable :: ntot_pT(:), ioff_pT(:)
!       integer :: ie_start, ie_end
    
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
  
  ! open output files
  if(myrank==master) then
     write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(*,*) '!!! Most output written to output.myrank.txt files !!!'
     write(*,*) '!!! please check these files as well               !!!'
     write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  end if
  write(ctemp,'(I03.3)') myrank
  open(1337, file='output.'//trim(ctemp)//'.txt')
  
!      write(*,*) 'nkmesh',myrank, t_inc%nkmesh, allocated(t_inc%kmesh)
  if(myrank==master) call main0()
!      write(*,*) 'nkmesh',myrank, t_inc%nkmesh, allocated(t_inc%kmesh)
  call timing_stop('main0')
  
#ifdef CPP_MPI
  !now communicate type t_inc and t_tgmat switches (this has an implicit barrier, so that all other processes wait for master to finish with main0)
  call timing_start('MPI 1')
  call bcast_t_inc_tgmat(t_inc,t_tgmat)
  !also communicate logicals from t_lloyd
  call bcast_t_lly_1(t_inc,t_lloyd)
  
  
  
  !in case of deci-out run everything only serially to write decifile correctly. This will be fixed in a later version when we get rid of the decifile
  if(t_inc%deci_out .and. nranks>1) stop 'deci-out option chosen. Please run code serially!'
  
  
  ! create 2d matrix for processors so that more processors than energy points can be used. 
  ! strategy here is to first parallelize the energy points ideally and then give all processors that are left to atom or k-point loop parallelization
  
  ! create_subcomms_2d:
  call find_dims_2d(nranks,t_inc%NATYP,t_inc%IELAST,dims)
!   write(*,*) 'find dims:', dims
  
  ! create communicator for atom/energy matrix
!   if(t_inc%nkmesh>1) then
    call create_newcomms_group_ie( nranks,myrank,dims(1),dims(2),t_inc%nkmesh,t_inc%kmesh,mympi_comm_ie,  &
                                 & myrank_ie,nranks_ie,mympi_comm_at,myrank_at,nranks_at, myrank_atcomm,nranks_atcomm)
  ! save grid info in type
    call save_t_mpi_c_grid(t_mpi_c_grid,dims, myMPI_comm_ie, myMPI_comm_at, myrank_ie, myrank_at, myrank_atcomm, nranks_ie, nranks_at, nranks_atcomm)
!   end if
  call timing_stop('MPI 1')
#endif

! bug found: for openmp run with SOC old files have to be written out, otherwise something goes wrong
! thus - for now - this is enforced here:
!#ifndef CPP_MPI
#ifdef CPP_OMP
  t_tgmat%tmat_to_file = .true.
  t_tgmat%gref_to_file = .true.
  t_tgmat%gmat_to_file = .true.
#endif
  
  
  ! Then start scf iterations and do all steps of the KKR formalism until convergence
  if(myrank==master) then
     write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     write(*,*) '+++            SCF ITERATIONS START                +++'
     write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     call print_time_and_date('started')
  end if
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

    if (myrank==master) call print_time_and_date('Iteration finished')
    
  end do ! scf-iteration
  
  
#ifdef CPP_MPI
  ! finalize MPI  
  call MPI_Finalize(ierr)
#endif


end program
