!***********************************************************************
subroutine gaunt2(w, yr, n)
! ************************************************************************
!     sets up values needed for gaunt
!        m. weinert  january 1982

!     changed for calculating with real spherical harmonics
!                                           b.drittler  july 1987

!     W(N)        integration weights on 4*LMAXD points in the intervall
!                 (-1,0) (from routine GRULE)

!     YR(N,L,M)   spherical harmonics on 4*LMAXD points to angular
!                 momentum indices (l,m) scaled with a factor
!                 of RF=(4*pi)**(1/3)

!-----------------------------------------------------------------------
  implicit none
!.. Arguments
  integer :: n
  double precision :: w(*), yr(n, 0:n, 0:n)
!..
!.. Local Scalars ..
  double precision :: a, cd, cth, fac, fpi, rf, sth, t
  integer :: k, l, lomax, m
!..
!.. Local Arrays ..
  double precision :: p(0:n+1, 0:n), x(n)
!..
!.. External Subroutines ..
  external :: grule
!..
!.. Intrinsic Functions ..
  intrinsic :: atan, sqrt
!     ..
  fpi = 16d0*atan(1d0)
  rf = fpi**(1d0/3d0)
  lomax = n

!--->    obtain gauss-legendre points and weights

  call grule(2*n, x, w)

!--->    generate associated legendre functions for m.ge.0

  do k = 1, n
    cth = x(k)
    sth = sqrt(1.d0-cth*cth)
    fac = 1.d0

!--->    loop over m values

    do m = 0, lomax
      fac = -dble(2*m-1)*fac
      p(m, m) = fac
      p(m+1, m) = dble(2*m+1)*cth*fac

!--->    recurse upward in l

      do l = m + 2, lomax
        p(l, m) = (dble(2*l-1)*cth*p(l-1,m)-dble(l+m-1)*p(l-2,m))/dble(l-m)
      end do

      fac = fac*sth
    end do

!--->    multiply in the normalization factors

    do l = 0, lomax
      a = rf*sqrt((2*l+1)/fpi)
      cd = 1.d0
      yr(k, l, 0) = a*p(l, 0)

      do m = 1, l
        t = dble((l+1-m)*(l+m))
        cd = cd/t
        yr(k, l, m) = a*sqrt(2.d0*cd)*p(l, m)
      end do
    end do
  end do
end subroutine
