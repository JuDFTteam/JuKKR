program test

  use mpi
  use omp_lib
  implicit none
  integer :: myrank, nranks, ierr
  integer :: ith, nth

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myrank, ierr)
  call mpi_comm_size(mpi_comm_world, nranks, ierr)
  !$omp parallel
  nth = omp_get_num_threads()
  ith = omp_get_thread_num()
  write(*,'(A,2I5,A,2I5,A)') 'myrank, mythread=', myrank, ith, ' of ', nranks, nth, ' ranks,threads'
  !$omp end parallel
  call mpi_finalize(ierr)

end program test
