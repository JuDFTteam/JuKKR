  subroutine ratval_gf(ia,ja,e0,gf)
! Rational function evaluation
! TEST VERSION
  use global

  implicit none

! Which block of the GF
  integer(kind=i4b), intent(in)  :: ia, ja
! Desired point
  complex(kind=c8b), intent(in)  :: e0
! Result of fit
  complex(kind=c8b), intent(out) :: gf(nlmsb,nlmsb)
! -----------------------------------------
  real(kind=r8b), parameter :: tol = 1.d-16
! -----------------------------------------
  integer(kind=i4b) :: i, j

! Are we in the upper complex plane?
  if (aimag(e0) < 0.d0) stop 'ratval_gf: Im E < 0'
  gf = 0.d0
! Get coefficients
  do j=1,nlmsba(ja)
    do i=1,nlmsba(ia)
!      if (sum(abs(gffit(:,i,j,ia,ja))) < tol) cycle
      call ratval(numd,dend,gffit(:,i,j,ia,ja),eshift,e0,gf(i,j))
    end do
  end do
! All done!
  end subroutine ratval_gf

