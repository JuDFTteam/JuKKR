  subroutine init_arrays_tgmat(my_rank,mpi_size,mpi_iebounds)
! Allocate the storage according to the susc input
  use global

  implicit none
  
  integer(kind=i4b) :: is, il, im, ib, i, lm, ia, ilms, jlms, my_rank, mpi_size
  real(kind=r8b)    :: ram
! --> mpi_bounds
  integer(kind=i4b), intent(in) :: mpi_iebounds(2,0:mpi_size-1) 

  if (noparameters) stop 'init_arrays_tgmat: run init_param first!'
! -----------------------------------------------------------------------
!               Storage for t-matrices andstructural GF
! -----------------------------------------------------------------------
  if (my_rank == 0) then
    write(*,'(" init_arrays_tgmat: t-matrix  RAM=",f16.3," GB")') 16.d0*nasusc*nlms*nlms*nesusc/1024.d0**3
  end if ! my_rank
  allocate(tmatrix(nlms,nlms,nasusc,mpi_iebounds(1,my_rank):mpi_iebounds(2,my_rank)))
  tmatrix = 0.d0
  if (my_rank == 0) then
    write(*,'(" init_arrays_tgmat: struct GF RAM=",f16.3," GB")') 16.d0*nalms*nalms*nesusc/1024.d0**3
  end if ! my_rank
  allocate(gstruct(nalms,nalms,mpi_iebounds(1,my_rank):mpi_iebounds(2,my_rank)))
  gstruct = 0.d0
! All done!
  end subroutine init_arrays_tgmat
