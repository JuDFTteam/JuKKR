! this is a test program to compare the time a call of dot_product takes
! compared to a direct implementation of the loop or lapacks zdotu
! this kind of loop structure appears in rhoin
program test

implicit none
integer, parameter :: m=10000
integer :: i, j, n, k
integer :: clock_rate
integer :: start_time
integer :: stop_time
real*8 :: timing, timing2, timing3
real*8, allocatable :: tmp(:,:), tmp2(:,:)
complex*16, allocatable :: a(:,:), b(:,:), c(:), c2(:), a2(:), b2(:)
complex*16, external :: zdotu

call system_clock(count_rate=clock_rate) ! Find the rate

do n=1,100
  allocate(a(n,n), b(n,n), c(n), c2(n), a2(n*n), b2(n*n), tmp(n,n), tmp2(n,n))
  call random_number(tmp)
  call random_number(tmp2)
  a = cmplx(tmp, tmp2)
  call random_number(tmp)
  call random_number(tmp2)
  b = cmplx(tmp, tmp2)
 
  call system_clock(count=start_time)
  do j=1,m
    do i=1,n
      c(i) = zdotu(n,a(i,1),n,b(i,1),n)
    end do
  end do
  call system_clock(count=stop_time) ! Stop timing
  timing = (stop_time-start_time)/real(clock_rate)
 
  call system_clock(count=start_time)
  do j=1,m
    a2 = reshape(a, [n**2])
    b2 = reshape(b, [n**2])
    do i=1,n
      c2(i) = dot_product(conjg(a2(i::n)), b2(i::n))
    end do
  end do
  call system_clock(count=stop_time) ! Stop timing
  timing2 = (stop_time-start_time)/real(clock_rate)
 
  call system_clock(count=start_time)
  do j=1,m
    a2 = reshape(a, [n**2])
    b2 = reshape(b, [n**2])
    c2(:) = (0.0, 0.0)
    !$omp parallel do private(i, k) default(shared)
    do i=1,n
      do k=0,n-1
        c2(i) = c2(i) + a2(i+n*k)*b2(i+n*k)
      end do
    end do
    !$omp end parallel do
  end do
  call system_clock(count=stop_time) ! Stop timing
  timing3 = (stop_time-start_time)/real(clock_rate)
 
  !write(*,*) 'c', c
  !write(*,*) 'c2', c2
  !write(*,*) n, 'c-c2', sum(c-c2)
  write(*,'(i5,3es25.16)') n, timing, timing2, timing3
 
  deallocate(a, b, c, c2, a2, b2, tmp, tmp2)
end do

end program test
