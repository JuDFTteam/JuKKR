program Amatprecalc
  use type_inc,  only: inc_type
  use mod_mympi, only: mympi_init
  use mod_read,  only: read_inc
  use mod_scattering, only: read_scattmat, impcls_type
  use mpi
  implicit none

  type(inc_type) :: inc
  integer        :: ierr

  type(impcls_type), allocatable :: impcls(:)
  double complex,    allocatable :: Amat(:,:,:)

  !initialize MPI
  call MPI_Init ( ierr )
  call mympi_init()

  call read_inc(inc)

  call read_scattmat(inc, impcls, Amat, .true.)


  call MPI_Finalize(ierr)

end program Amatprecalc
