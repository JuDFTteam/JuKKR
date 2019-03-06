program test

  use omp_lib
  implicit none
  integer :: ith, nth

  !$omp parallel
  nth = omp_get_num_threads()
  ith = omp_get_thread_num()
  write(*,'(A,I5,A,I5,A)') 'mythread=', ith, ' of ', nth, ' threads'
  !$omp end parallel

end program test
