  subroutine baryint_gf(n,x,w,ia,ja,x0,gf0,onsite,struct)
! Rational barycentric interpolation
! Works also with equidistant points, as long as order is kept low

  implicit none

! Number of points
  integer(kind=i4b), intent(in)  :: n
! List of abcissas and interpolation weights
  complex(kind=c8b), intent(in)  :: x(n), w(n)
! Block of GF
  integer(kind=i4b), intent(in)  :: ia, ja
! Desired point
  complex(kind=c8b), intent(in)  :: x0
! Result of interpolation
  complex(kind=c8b), intent(out) :: gf0(nlmsb,nlmsb)
! To add onsite or structural parts of GF
  logical,           intent(in)  :: onsite, struct
! -----------------------------------------
  real(kind=r8b), parameter :: tol = 1.d-8
! -----------------------------------------
  integer(kind=i4b) :: k
  complex(kind=c8b) :: hstep, gf(nlmsb,nlmsb), den

  gf0 = 0.d0; den = 0.d0
  do k=1,n
    call projected_gf(k,ia,ja,gf,onsite,struct)
    hstep = x0 - x(k)
    if (abs(hstep) < tol) then
      gf0 = gf
      return
    end if
    hstep = w(k)/hstep
    gf0 = gf0 + hstep*gf
    den = den + hstep
  end do
  gf0 = gf0/den
! All done!
  end subroutine baryint_gf

