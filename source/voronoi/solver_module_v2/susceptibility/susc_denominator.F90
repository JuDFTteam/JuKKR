  subroutine susc_denominator
! Analysis of 1 - KS susc * xc kernel
  use global

  implicit none

!  integer(kind=i4b), parameter :: itermax = 1
  complex(kind=c8b), parameter :: czero = (0.d0,0.d0), cminus = (-1.d0,0.d0)
  integer(kind=i4b) :: iq, jq, info, iter
  complex(kind=c8b) :: w(ndensum), wl, wr, work(2*ndensum)
  real(kind=r8b)    :: start, finish, rwork(2*ndensum), rw(ndensum), minlambda

! KS susc eigendecomposition
  call cpu_time(start)
!  kssusc0 = 0.5d0*(kssusc0 + transpose(kssusc0))
!  call zcopy(ndensum*ndensum,kssusc0,1,denominator,1)
!  call zgeev('N','N',ndensum,denominator,ndensum,w,wl,1,wr,1,work,2*ndensum,rwork,info)
!  if (info /= 0) stop 'susc_denominator: failure in zgeev: kssusc0'
!  call cpu_time(finish)
!  write(iodb,'(" denominator eigen time=",f10.3," s")') finish - start
!  rw = real(w)
! sort them into descending order (MKL utility)
!  call dlasrt('D',nlmsb*nlmsb,rw,info)
!  do iq=1,ndensum
!    if (abs(w(iq)) > atol) write(iodb,'("kssusc0 eval=",i4,2es16.8)') iq, w(iq)
!  end do
! xc kernel eigendecomposition
!  call cpu_time(start)
!  kxcsusc = real(kxcsusc)
!  call zcopy(ndensum*ndensum,kernel,1,denominator,1)
!  call zgeev('N','N',ndensum,denominator,ndensum,w,wl,1,wr,1,work,2*ndensum,rwork,info)
!  if (info /= 0) stop 'susc_denominator: failure in zgeev: kernel'
!  call cpu_time(finish)
!  write(*,'(" denominator eigen time=",f10.3," s")') finish - start
!  rw = real(w)
! sort them into descending order (MKL utility)
!  call dlasrt('D',nlmsb*nlmsb,rw,info)
!  do iq=1,ndensum
!    if (abs(w(iq)) > atol) write(iodb,'("kxcsusc eval=",i4,2es16.8)') iq, w(iq)
!  end do
! enhancement factor
  call zgemm('N','N',ndensum,ndensum,ndensum,cminus,kssusc0,ndensum,kernel,ndensum,czero,denominator,ndensum)
  do iq=1,ndensum
    denominator(iq,iq) = denominator(iq,iq) + 1.d0
!    do jq=0,ndensum-10,10
!      write(*,'(2i4,20es7.0)') iq, jq, denominator(jq+1:jq+10,iq)
!    end do
  end do
! Eigendecomposition
  call cpu_time(start)
  call zgeev('N','N',ndensum,denominator,ndensum,w,wl,1,wr,1,work,2*ndensum,rwork,info)
  if (info /= 0) stop 'susc_denominator: failure in zgeev: denominator'
  call cpu_time(finish)
  write(iodb,'(" denominator eigen time=",f10.3," s")') finish - start
!  rw = real(w)
! sort them into descending order (MKL utility)
!  call dlasrt('D',nlmsb*nlmsb,rw,info)
  minlambda = 1.d0
  do iq=1,ndensum
    if (abs(w(iq)) < abs(minlambda)) minlambda = w(iq)
    if (abs(w(iq)) < 0.999d0 .or. abs(w(iq)) > 1.001d0) write(iodb,'("denominator eval=",i4,2es16.8)') iq, w(iq)
  end do
  write(*,'("susc_denominator: minlambda=",es8.1)') minlambda
! *******************
  if (.not.lsoc) then
! *******************
    do iter=1,itermax
!     Correction
      if (minlambda < 0.d0) kernel = kernel*(1.d0 + (1.d0 + lambdamix)*minlambda)
      if (minlambda > 0.d0) kernel = kernel*(1.d0 + (1.d0 - lambdamix)*minlambda)
!     enhancement factor
      call zgemm('N','N',ndensum,ndensum,ndensum,cminus,kssusc0,ndensum,kernel,ndensum,czero,denominator,ndensum)
      do iq=1,ndensum
        denominator(iq,iq) = denominator(iq,iq) + 1.d0
      end do
!     Eigendecomposition
      call cpu_time(start)
      call zgeev('N','N',ndensum,denominator,ndensum,w,wl,1,wr,1,work,2*ndensum,rwork,info)
      if (info /= 0) stop 'susc_denominator: failure in zgeev'
      call cpu_time(finish)
      write(iodb,'(" denominator eigen time=",f10.3," s")') finish - start
!      rw = real(w)
!     sort them into descending order (MKL utility)
!      call dlasrt('D',nlmsb*nlmsb,rw,info)
      minlambda = 1.d0
      do iq=1,ndensum
        if (abs(w(iq)) < abs(minlambda)) minlambda = w(iq)
        if (abs(w(iq)) < 0.999d0 .or. abs(w(iq)) > 1.001d0) write(iodb,'("denominator eval=",i4,2es16.8)') iq, w(iq)
      end do
      write(*,'("susc_denominator: minlambda=",es8.1)') minlambda
    end do
! ******
  end if
! ******
  end subroutine susc_denominator
