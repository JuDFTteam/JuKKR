module mod_mympi
!ruess: taken from Pkkr_sidebranch2D_2014_12_16, created by Bernd Zimmermann

implicit none

  private
  public :: myrank, nranks, master, mympi_init

  integer, save :: myrank = -1
  integer, save :: nranks = -1
  integer, save :: master = -1

contains

  subroutine mympi_init()

#ifdef CPP_MPI
    use mpi
#endif

    integer :: ierr

    master = 0

#ifdef CPP_MPI
    call MPI_Comm_rank ( MPI_COMM_WORLD, myrank, ierr )
    call MPI_Comm_size ( MPI_COMM_WORLD, nranks, ierr )
#else
    myrank = master
    nranks = 1
#endif

  end subroutine mympi_init


end module mod_mympi
