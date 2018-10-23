program test
use omp_lib
implicit none

integer N, i,j, mythread, ierr
double complex, allocatable :: sum(:)
double complex, allocatable, save :: tmp(:,:)

!$omp threadprivate(tmp)

N = 2

!$omp parallel private(ierr, mythread)
mythread = omp_get_thread_num()
allocate(tmp(N,N), stat=ierr)
if(ierr/=0) then
   write(*,*) mythread,'error allocating tmp'
   stop
end if
!$omp end parallel

allocate(sum(N), stat=ierr)
if(ierr/=0) stop 'Error alloc sum'



!$omp parallel do private(i,j)
do i=1,N
  do j=1,N
    tmp(j,i) = i*N * (1.0d0, 0.0d0) + j*N**2*(-1)**i * (0.0d0, 1.0d0)
  end do
end do
!$omp end parallel do


!$omp parallel do private(i,j) reduction(+:sum)
do i=1,N
  sum(i) = (0.0d0, 0.0d0)
  do j=1,N
    sum(i) = sum(i) + tmp(j,i)
  end do
end do
!$omp end parallel do



!$omp parallel private(mythread, ierr)
mythread = omp_get_thread_num()
!$omp critical
write(*,'(I9,1000ES15.6)') mythread, tmp(:,:)
!$omp end critical

!$omp barrier
write(*,*)
!$omp barrier

!$omp critical
write(*,'(I9,1000ES15.6)') mythread, sum(:)
!$omp end critical

!$omp barrier
deallocate(tmp, stat=ierr)

! !$omp critical
if(ierr/=0) then
   write(*,*) mythread,'error deallocating tmp'
   stop
end if
! !$omp end critial
!$omp end parallel

deallocate(sum)


write(*,*)
write(*,*) '=================================='
write(*,*) 'TEST 2:'


allocate(sum(1))

!$omp parallel default(shared) private(mythread, i)

mythread = omp_get_thread_num()

sum = (0.0d0, 0.0d0)

!$omp do reduction(+:sum)
do i=1, 10
   sum(1) = sum(1) + (mythread+1)**2*(1.0d0, 0.0d0)
   write(*,*) 'in loop', mythread, dreal(sum(1))
end do
!$omp end do
!$omp barrier


write(*,*) 'after loop', mythread, dreal(sum(1))

!$omp end parallel

write(*,*) 'after omp', mythread, dreal(sum(1))


deallocate(sum)

end program test
