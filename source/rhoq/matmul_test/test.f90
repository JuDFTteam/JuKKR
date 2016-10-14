program test

implicit none

double complex, allocatable :: A(:,:), B(:,:), C(:,:)
integer nblock, i, j, N, lmax, nl1, nl2 
integer clock_rate, start_time, stop_time
double precision time

nl1 = 9000 
nl2 = 9000 
CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) ! Find the rate

write(*,*) '# nl1, nl2:', nl1, nl2
write(*,*) '# nblock, lmax, N, time'

do lmax=1,3
  do nblock=1,5
    N = 2*(lmax+1)**2 * nblock
    !write(*,*), nblock, lmax, N
   
    CALL SYSTEM_CLOCK(COUNT=start_time) ! Start timing
    allocate(A(N,N), B(N,N), C(N,N))
    do i=1,N
      do j=1,N
        A(i,j) = i*(1.0d0,0.0d0) + j*(0.0d0, 1.0d0)
        B(i,j) = j*(1.0d0,0.0d0) + i*(0.0d0, 1.0d0)
        C(i,j) = (i+j)*(1.0d0,0.0d0) + (i-j)*(0.0d0, 1.0d0)
      end do
    end do
    CALL SYSTEM_CLOCK(COUNT=stop_time) ! Stop timing
    time = (stop_time-start_time)/real(clock_rate)
    !write(*,*) 'time for initi:', time
   
    CALL SYSTEM_CLOCK(COUNT=start_time) ! Start timing
    do i=1,nl1/nblock/(2*(lmax+1)**2)
      do j=1,nl2/nblock/(2*(lmax+1)**2)
        CALL ZGEMM('N','N',N,N,N,(1.0d0, 0.0d0),A,N,B,N,(0.0d0, 0.0d0),C,N)
      end do
    end do
    CALL SYSTEM_CLOCK(COUNT=stop_time) ! Stop timing
    time = (stop_time-start_time)/real(clock_rate)
    !write(*,*) 'time for Zgemms:',time
    write(*,'(5I9,ES16.7)') nblock, lmax, nl1/nblock/(2*(lmax+1)**2), (nl1/nblock/(2*(lmax+1)**2)*N)**2, N, time
   
    deallocate(A,B,C)
  end do
end do


end program test
