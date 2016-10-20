program test

use omp_lib
implicit none

double complex, allocatable :: A(:,:), B(:,:), C(:,:), D(:,:)
integer iaverage,Naverage,nblock, i, j, Nk, Nq 
integer clock_rate, start_time, stop_time
double precision time, t(3)

integer mythread, nthreads

nblock = 32 
Nk = 10000 !10000 !100x100
Nq = 10000 !10000
Naverage= 1  !50

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! init !!!!!!!!!!!!!!!!!!
CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) ! Find the clock rate
   
CALL SYSTEM_CLOCK(COUNT=start_time) ! Start timing
!allocate and initialize matrices
allocate(C(nblock,Nk*nblock), A(Nk*nblock,nblock), B(nblock,nblock), D(nblock, nblock), stat=i)
if(i/=0) stop 'ERROR allocating A,B,C,D'

do i=1,nblock*Nk
  do j=1,Nblock
    A(i,j) = i*(1.0d0,0.0d0)/Nblock/Nk + j*(0.0d0, 1.0d0)/Nk/Nblock
    C(j,i) = j*(1.0d0,0.0d0)/Nblock/Nk + i*(0.0d0, 1.0d0)/Nblock/Nk
    if(i<=nblock) B(i,j) = (i+j)*(1.0d0,0.0d0)/Nk/Nblock + (i-j)/Nk*(0.0d0, 1.0d0)
  end do
end do
D(:,:) = (0.0d0, 0.0d0)

CALL SYSTEM_CLOCK(COUNT=stop_time) ! Stop timing
time = (stop_time-start_time)/real(clock_rate)
!write(*,*) 'time for init:', time
t(1) = time 
!!!!!!!!!!!!!!!!!! init !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! big matrix !!!!!!!!!!!!!!
time = 0.0d0
!write(*,'("Loop over points:|",5(1X,I2,"%",5X,"|"),1X,I3,"%")') 0, 20, 40, 60, 80, 100
!write(*,FMT=190) !beginning of statusbar
do iaverage=1, Naverage
   CALL SYSTEM_CLOCK(COUNT=start_time) ! Start timing
   !A*b = A' (32*Nkx32.32x32=32*nkx32)
   CALL ZGEMM('N','N',Nk*nblock,nblock,nblock,(1.0d0, 0.0d0),A,Nk*nblock,B,nblock,(0.0d0, 0.0d0),A,Nk*nblock)
   !C*A'=D
   CALL ZGEMM('N','N',nblock,nblock,Nk*nblock,(1.0d0, 0.0d0),C,nblock,A,Nk*nblock,(0.0d0, 0.0d0),D,nblock)
   CALL SYSTEM_CLOCK(COUNT=stop_time) ! Stop timing
   time = time + (stop_time-start_time)/real(clock_rate)
   !update statusbar
   !if(mod(iaverage,1)==0) write(*,FMT=200)
end do
!write(*,*) ! status bar
!write(*,*) 'time for Zgemms:',time/Naverage, Naverage
!write(*,*) 'total time',time+t0 
t(2) = time
!!!!!!!!!!!!!!!! big matrix !!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! small matrices !!!!!!!!!!!!
time = 0.0d0
!write(*,'("Loop over points:|",5(1X,I2,"%",5X,"|"),1X,I3,"%")') 0, 20, 40, 60, 80, 100
!write(*,FMT=190) !beginning of statusbar
!$omp parallel default(shared) private(mythread, nthreads, i,B,D)
mythread = omp_get_thread_num()
nthreads= omp_get_num_threads()
if(.not.allocated(D)) allocate(D(nblock, nblock), stat=i)
if(.not.allocated(B)) allocate(B(nblock, nblock), stat=i)
!if(mythread==0) write(*,*) 'Using',nthreads,'OpenMP threads'

do iaverage=1, Naverage
   CALL SYSTEM_CLOCK(COUNT=start_time) ! Start timing
   !A*b = A' (32*Nkx32.32x32=32*nkx32)
   CALL ZGEMM('N','N',Nk*nblock,nblock,nblock,(1.0d0, 0.0d0),A,Nk*nblock,B,nblock,(0.0d0, 0.0d0),A,Nk*nblock)
   !C*A'=D
   !$omp do
   do i=1,Nk
      B = A(nblock*(i-1)+1:nblock*i,:)
      CALL ZGEMM('N','N',nblock,nblock,nblock,(1.0d0, 0.0d0),C,nblock,B,nblock,(0.0d0, 0.0d0),D,nblock)
   end do
   !$omp end do
   CALL SYSTEM_CLOCK(COUNT=stop_time) ! Stop timing
   time = time + (stop_time-start_time)/real(clock_rate)
   !update statusbar
   !if(mod(iaverage,1)==0 .and. mythread==0) write(*,FMT=200)
end do
!$omp end parallel
!write(*,*) ! status bar
!write(*,*) 'time for Zgemm2:',time/Naverage, Naverage
!write(*,*) 'total time',time+t0 
t(3) = time
!!!!!!!!!!!!!! small matrices !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!write(*,*)
!write(*,'(A)') '----------------------------------------'
 
!write(*,'(A3ESF16.7)') 'OMP,Nk, total/big/small per qpt ', t1, t1-time, time
!$omp parallel default(shared) private(mythread, nthreads)
mythread = omp_get_thread_num()
nthreads= omp_get_num_threads()
if(mythread==0.and.nthreads==1) write(*,*) '#  nthreads, Nk, total time, one big matrix/q, lots small/q ', omp_get_nested()
if(mythread==0) write(*,'(2I9,3ES16.7)') nthreads, Nk, sum(t), t(2)/Naverage, t(3)/Naverage
!$omp end parallel

!write(*,*) 'total run time: ', t1
!write(*,*) 'big matrix:     ', t1-time
!write(*,*) 'small:          ', time
!write(*,*) 'ratio:',(t1-time)/t1
!write(*,*)
!write(*,*) 'for',Nq,'qkpts:', (t1-time)*Nq/Naverage, t1*Nq/Naverage
!write(*,*) 'in mins:', (t1-time)*Nq/Naverage/60, t1*Nq/Naverage/60
!write(*,*) 'in hrs: ', (t1-time)*Nq/Naverage/3600, t1*Nq/Naverage/3600

deallocate(A,B,C,D)


190      FORMAT('                 |'$)   ! status bar
200      FORMAT('|'$)                    ! status bar

end program test
