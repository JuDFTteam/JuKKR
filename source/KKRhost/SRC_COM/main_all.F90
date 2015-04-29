program kkrcode

  use mod_main0
  use mod_main1a
  use mod_main1b
  use mod_main1c
  use mod_main2
  use mod_types

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
  ! set variables master, myrank and nranks for serial as well as mpi-parallel execution
  call mympi_init()

    
  ! start KKR program, first do main0, where everything is read in and initialized
  if(myrank==master) call main0()
  
  
#ifdef CPP_MPI
  !now communicate type t_inc and t_tgmat switches (this has an implicit barrier, so that all other processes wait for master to finish with main0)
  call bcast_t_inc_tgmat(t_inc,t_tgmat)
  
  ! create 2d matrix for processors so that more processors than energy points can be used. 
  ! strategy here is to first parallelize the energy points ideally and then give all processors that are left to atom or k-point loop parallelization
  
  ! create_subcomms_2d:
  ntot1 = t_inc%NATYP
!   if(.not.t_inc%NEWSOSOL) ntot1 = t_inc%NSPIN*t_inc%NATYP
  call find_dims_2d(nranks,ntot1,t_inc%IELAST,dims)
  
!   write(*,*) 'finddims:',myrank,nranks,ntot1,t_inc%IELAST,dims

  ! create cartesian communicator for atom/energy matrix
  call create_subarr_comm( dims, myMPI_comm_grid, myMPI_comm_row, myMPI_comm_col, myrank_grid, myrank_row, myrank_col, nranks_row, nranks_col)
  ! save grid info in type
  call save_t_mpi_c_grid(t_mpi_c_grid,dims, myMPI_comm_grid, myMPI_comm_row, myMPI_comm_col, myrank_grid, myrank_row, myrank_col, nranks_row, nranks_col)
#endif
  
  ! Then start scf iterations and do all steps of the KKR formalism until convergence
  do while ( (t_inc%i_iteration.lt.t_inc%N_iteration) .and. (t_inc%N_iteration.ne.-1) )
  
      write(*,*) 'before main1a',myrank,myrank_col,myrank_row,master,nranks, nranks_col, nranks_row,t_inc%i_iteration,t_inc%N_iteration

    ! calculate tmat and gref
    !if (myrank==master) call main1a()
    call main1a()

    write(*,*) myrank,'finished main1a'
    !test      
    if(myrank==master) then
    write(*,*) t_inc%natyp, t_inc%NSPIN, t_inc%ielast
    do i1=1,t_inc%natyp
    do ispin=1,t_inc%NSPIN
    do ie=1,t_inc%ielast
!        write(*,*) myrank,i1,ie,irec,shape(t_tgmat%TMAT)
!        IREC = IE + t_inc%IELAST*(I1-1)
         IREC = ie + t_inc%IELAST*(ISPIN-1) + t_inc%IELAST*t_inc%NSPIN* (I1-1)
!          write(*,*) irec
       write(696969,*) t_tgmat%TMAT(:,:,irec)
    end do
    end do
    end do
    
    do ie=1,t_inc%ielast
!        write(*,*) myrank,ie,shape(t_tgmat%gref)
       write(686868,*) t_tgmat%gref(:,:,:,ie)
    end do
    endif

    
!     write(*,*) t_tgmat%gmat_to_file,t_tgmat%gref_to_file,t_tgmat%tmat_to_file
    
    
    ! calculate gmat
!     if (myrank==master) call main1b()
    call main1b()

   
    ! calculate density
!     if (myrank==master) call main1c()
    call main1c()
    
        

   
    ! calculate do DFT stuff (potential from density, exc-potential, calculate total energy, ...)
    if (myrank==master) call main2()
    
    
    call bcast_t_inc_tgmat(t_inc,t_tgmat)
!     call MPI_Barrier(MPI_COMM_WORLD,ie)
!     STOP
    
!     write(*,*) myrank,myrank_col,myrank_row,master,nranks, nranks_col, nranks_row,t_inc%i_iteration,t_inc%N_iteration
    
  end do
  
  
#ifdef CPP_MPI
  ! finalize MPI  
  call MPI_Finalize(ierr)
#endif


end program
