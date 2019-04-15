program test

implicit none
integer, parameter :: m=10000
integer :: i, j, n, k, l
integer :: clock_rate
integer :: start_time
integer :: stop_time
real*8 :: timing, timing2, timing3
real*8, allocatable :: tmp(:,:), tmp2(:,:)
complex*16, allocatable :: a(:,:), b(:,:), c(:,:), b2(:,:)
complex*16 :: temp
complex*16, external :: zdotu

call system_clock(count_rate=clock_rate) ! Find the rate

do n=1,100
  allocate(a(n,n), b(n,n), c(n,n), tmp(n,n), tmp2(n,n))
  call random_number(tmp)
  call random_number(tmp2)
  a = cmplx(tmp, tmp2)
  call random_number(tmp)
  call random_number(tmp2)
  b = cmplx(tmp, tmp2)
 
  call system_clock(count=start_time)
  do j=1,m
    call zgemm('N','N', n,n,n, (1.0,0.0), a,n,b,n,(0.0,0.0),c,n)
  end do
  call system_clock(count=stop_time) ! Stop timing
  timing = (stop_time-start_time)/real(clock_rate)
 
  call system_clock(count=start_time)
  do j=1,m
    c = matmul(a, b)
  end do
  call system_clock(count=stop_time) ! Stop timing
  timing2 = (stop_time-start_time)/real(clock_rate)
 
  call system_clock(count=start_time)
  do j=1,m
    c = (0.0, 0.0)
    !$omp parallel do private(i, k) default(shared)
    do i=1,n
      do k=1,n
        temp = b(k,i)
        do l=1,n
          c(l,i)  = c(l,i) + temp*a(l,k)     
        end do
      end do
    end do
    !$omp end parallel do
  end do
  call system_clock(count=stop_time) ! Stop timing
  timing3 = (stop_time-start_time)/real(clock_rate)
 
  !write(*,*) 'c', c
  write(*,'(i5,3es25.16)') n, timing, timing2, timing3
 
  deallocate(a, b, c, tmp, tmp2)
end do

end program test
