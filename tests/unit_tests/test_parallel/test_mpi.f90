program test

  use mpi
  implicit none
  integer :: myrank, nranks, ierr

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myrank, ierr)
  call mpi_comm_size(mpi_comm_world, nranks, ierr)
  write(*,'(A,I5,A,I5,A)') 'myrank=', myrank, ' of ', nranks, ' ranks'
  call mpi_finalize(ierr)

end program test
