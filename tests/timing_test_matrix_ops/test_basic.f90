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
  c = b
  temp = cmplx(tmp(1,n),tmp2(n,1))
 
  call system_clock(count=start_time)
  do j=1,m
    call zcopy(n*n, a,1,b,1)
  end do
  call system_clock(count=stop_time) ! Stop timing
  timing = (stop_time-start_time)/real(clock_rate)

  b = c
  call system_clock(count=start_time)
  do j=1,m
    b(:,:) = a(:,:)
  end do
  call system_clock(count=stop_time) ! Stop timing
  timing2 = (stop_time-start_time)/real(clock_rate)
 
  b = c
  call system_clock(count=start_time)
  do j=1,m
    !$omp parallel do private(i, k) default(shared)
    do i=1,n
      do k=1,n
        b(k,i) = a(k,i)
      end do
    end do
    !$omp end parallel do
  end do
  call system_clock(count=stop_time) ! Stop timing
  timing3 = (stop_time-start_time)/real(clock_rate)
 
  write(*,'(i5,3es25.16)') n, timing, timing2, timing3

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call system_clock(count=start_time)
  do j=1,m
    call zscal(n*n, temp, a, 1, b, 1)
  end do
  call system_clock(count=stop_time) ! Stop timing
  timing = (stop_time-start_time)/real(clock_rate)

  b = c
  call system_clock(count=start_time)
  do j=1,m
    b(:,:) = temp*a(:,:)
  end do
  call system_clock(count=stop_time) ! Stop timing
  timing2 = (stop_time-start_time)/real(clock_rate)

  b = c
  call system_clock(count=start_time)
  do j=1,m
    !$omp parallel do private(i, k) default(shared)
    do i=1,n
      do k=1,n
        b(k,i) = temp*a(k,i)
      end do
    end do
    !$omp end parallel do
  end do
  call system_clock(count=stop_time) ! Stop timing
  timing3 = (stop_time-start_time)/real(clock_rate)
 
  write(*,'(i5,3es25.16)') n, timing, timing2, timing3

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  call system_clock(count=start_time)
  do j=1,m
    call zaxpy(n*n, temp, a, 1, b, 1)
  end do
  call system_clock(count=stop_time) ! Stop timing
  timing = (stop_time-start_time)/real(clock_rate)

  call system_clock(count=start_time)
  do j=1,m
    b(:,:) = b(:,:) + temp*a(:,:)
  end do
  call system_clock(count=stop_time) ! Stop timing
  timing2 = (stop_time-start_time)/real(clock_rate)
 
  call system_clock(count=start_time)
  do j=1,m
    !$omp parallel do private(i, k) default(shared)
    do i=1,n
      do k=1,n
        b(k,i) = b(k,i)+ temp*a(k,i)
      end do
    end do
    !$omp end parallel do
  end do
  call system_clock(count=stop_time) ! Stop timing
  timing3 = (stop_time-start_time)/real(clock_rate)
 
  write(*,'(i5,3es25.16)') n, timing, timing2, timing3
 
  deallocate(a, b, c, tmp, tmp2)
end do

end program test
