!-----------------------------------------------------------------------------------------!
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of kk-prime@juKKR and available as free software under the conditions !
! of the MIT license as expressed in the LICENSE file in more detail.                     !
!-----------------------------------------------------------------------------------------!


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
