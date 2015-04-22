program kkrcode

  use mod_main0
  use mod_main1a
  use mod_main1b
  use mod_main1c
  use mod_main2
  use mod_types
  use mod_mympi, only: mympi_init, myrank, nranks, master
  
#ifdef CPP_MPI
  use mpi
#endif

  implicit none
  
  integer :: ierr
  
  
  !initialize MPI
#ifdef CPP_MPI
    call MPI_Init ( ierr )
#endif
    call mympi_init()

    
  ! KKR program, first do main0, where everything is read in and initialized.
  if(myrank==master) then
    call main0()
    call init_tgmat(t_inc,t_tgmat)
  end if
  
  !now communicate type t_inc and t_tgmat (this has an implicit barrier, so that all other processes wait for master to finish with main0)
#ifdef CPP_MPI
    call bcast_t_inc_tgmat(t_inc,t_tgmat)
#endif

  
  ! Then start scf iterations and do all steps of the KKR formalism until convergence
  if (myrank==master) then !test
  do while ( (type0%i_iteration.lt.type0%N_iteration) .and. (type0%N_iteration.ne.-1) )
  
    call main1a()
   
    call main1b()
   
    call main1c()
   
    call main2()
    
  end do
  end if!test
  
  
#ifdef CPP_MPI
    call MPI_Finalize(ierr)
#endif

end program
