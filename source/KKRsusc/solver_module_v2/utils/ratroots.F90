  subroutine rat_roots(zero,numd,dend,fit,numroots,denroots)
! Roots of a rational function's numerator and denominator
  use global, only: i4b, r8b, c8b

  implicit none

! Shift of origin
  complex(kind=c8b), intent(in)  :: zero
! Degree of numerator and denominator of rational function in fit
  integer(kind=i4b), intent(in)  :: numd, dend
! Fitting coefficients
  complex(kind=c8b), intent(in)  :: fit(numd+dend)
! Roots of numerator and denominator
  complex(kind=c8b), intent(out) :: numroots(numd), denroots(dend)
! -----------------------------------------
  real(kind=r8b), parameter :: tol = 1.d-12
! -----------------------------------------
  integer(kind=i4b) :: i, info
  complex(kind=c8b) :: w(numd+dend), vl, vr, work(2*(numd+dend)), a(numd+dend,numd+dend)
  real(kind=r8b)    :: rwork(2*(numd+dend))

! numerator loop
  numroots = 0.d0
  if (numd > 1) then
!   fill in upper Hessenber matrix
    a = 0.d0
    a(1,1) = -fit(numd-1)/fit(numd)
    do i=2,numd-1
      a(i,i-1) = 1.d0
      a(1,i) = -fit(numd-i)/fit(numd)
    end do
    call zgeev('N','N',numd-1,a,numd+dend,w,vl,1,vr,1,work,2*(numd+dend),rwork,info)
    numroots = w(1:numd-1)
!    write(iodb,'("num roots: ",2es16.8)') numroots
  end if
! denominator loop
  denroots = 0.d0
  if (dend > 1) then
! fill in upper Hessenber matrix
    a = 0.d0
    a(1,1) = -fit(numd+dend-1)/fit(numd+dend)
    do i=2,dend-1
      a(i,i-1) = 1.d0
      a(1,i) = -fit(numd+dend-i)/fit(numd+dend)
    end do
    call zgeev('N','N',dend-1,a,numd+dend,w,vl,1,vr,1,work,2*(numd+dend),rwork,info)
    denroots = w(1:dend-1)
!    write(iodb,'("den roots: ",2es16.8)') denroots
  end if
! All done!
  end subroutine rat_roots

