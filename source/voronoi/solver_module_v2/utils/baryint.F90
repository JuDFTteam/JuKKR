  subroutine baryint(n,x,y,w,x0,y0)
! Rational barycentric interpolation
! Works also with equidistant points, as long as order is kept low

  implicit none

! Number of points
  integer(kind=i4b), intent(in)  :: n
! List of abcissas, function values and interpolation weights
  complex(kind=c8b), intent(in)  :: x(n), y(n), w(n)
! Desired point
  complex(kind=c8b), intent(in)  :: x0
! Result of interpolation
  complex(kind=c8b), intent(out) :: y0
! -----------------------------------------
  real(kind=r8b), parameter :: tol = 1.d-8
! -----------------------------------------
  integer(kind=i4b) :: k
  complex(kind=c8b) :: hstep, num, den

  num = 0.d0; den = 0.d0
  do k=1,n
    hstep = x0 - x(k)
    if (abs(hstep) < tol) then
      y0 = y(k)
      return
    end if
    hstep = w(k)/hstep
    num = num + hstep*y(k)
    den = den + hstep
  end do
  y0 = num/den
! All done!
  end subroutine baryint

