  subroutine ratval2(numd,dend,fit,zero,x0,y,dy)
! Rational function evaluation
  use global, only: i4b, r8b, c8b

  implicit none

! Degree of numerator and denominator of rational function in fit
  integer(kind=i4b), intent(in)  :: numd, dend
! Fitting coefficients
  complex(kind=c8b), intent(in)  :: fit(numd+dend)
! Shift of origin
  complex(kind=c8b), intent(in)  :: zero
! Desired point
  complex(kind=c8b), intent(in)  :: x0
! Result of fit, gradient wrt to parameters
  complex(kind=c8b), intent(out) :: y, dy(numd+dend)
! -----------------------------------------
  real(kind=r8b), parameter :: tol = 1.d-12
! -----------------------------------------
  integer(kind=i4b) :: p
  complex(kind=c8b) :: xfac, num, den

! Rational function:
! numerator loop, degree numd-1
  xfac = 1.d0; num = 0.d0
  do p=1,numd
    num  = num + xfac*fit(p)
    xfac = xfac*(x0 - zero)
  end do
! denominator loop, degree dend-1
  xfac = 1.d0; den = 0.d0
  do p=1,dend
    den  = den + xfac*fit(numd+p)
    xfac = xfac*(x0 - zero)
  end do
  den = den + tol
  y = num/den
! ---------------------------------
! Gradient wrt to parameters:
! numerator loop, degree numd-1
  xfac = 1.d0
  do p=1,numd
    dy(p) = xfac/den
    xfac = xfac*(x0 - zero)
  end do
! denominator loop, degree dend-1
  xfac = 1.d0
  do p=1,dend
    dy(numd+p)  = -y*(xfac/den)
    xfac = xfac*(x0 - zero)
  end do
! All done!
  end subroutine ratval2

