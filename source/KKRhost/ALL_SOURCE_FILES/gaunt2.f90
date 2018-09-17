module mod_gaunt2
  use :: mod_datatypes, only: dp
  private :: dp

contains

  ! ***********************************************************************
  subroutine gaunt2(w, yr, n)
    use :: mod_grule
    ! ************************************************************************
    ! sets up values needed for gaunt
    ! m. weinert  january 1982

    ! changed for calculating with real spherical harmonics
    ! b.drittler  july 1987

    ! W(N)        integration weights on 4*LMAXD points in the intervall
    ! (-1,0) (from routine GRULE)

    ! YR(N,L,M)   spherical harmonics on 4*LMAXD points to angular
    ! momentum indices (l,m) scaled with a factor
    ! of RF=(4*pi)**(1/3)

    ! -----------------------------------------------------------------------
    implicit none
    ! .. Arguments
    integer :: n
    real (kind=dp) :: w(*), yr(n, 0:n, 0:n)
    ! ..
    ! .. Local Scalars ..
    real (kind=dp) :: a, cd, cth, fac, fpi, rf, sth, t
    integer :: k, l, lomax, m
    ! ..
    ! .. Local Arrays ..
    real (kind=dp) :: p(0:n+1, 0:n), x(n)
    ! ..
    fpi = 16e0_dp*atan(1e0_dp)
    rf = fpi**(1e0_dp/3e0_dp)
    lomax = n

    ! --->    obtain gauss-legendre points and weights

    call grule(2*n, x, w)

    ! --->    generate associated legendre functions for m.ge.0

    do k = 1, n
      cth = x(k)
      sth = sqrt(1.e0_dp-cth*cth)
      fac = 1.e0_dp

      ! --->    loop over m values

      do m = 0, lomax
        fac = -real(2*m-1, kind=dp)*fac
        p(m, m) = fac
        p(m+1, m) = real(2*m+1, kind=dp)*cth*fac

        ! --->    recurse upward in l

        do l = m + 2, lomax
          p(l, m) = (real(2*l-1,kind=dp)*cth*p(l-1,m)-real(l+m-1,kind=dp)*p(l-2,m))/real(l-m, kind=dp)
        end do

        fac = fac*sth
      end do

      ! --->    multiply in the normalization factors

      do l = 0, lomax
        a = rf*sqrt((2*l+1)/fpi)
        cd = 1.e0_dp
        yr(k, l, 0) = a*p(l, 0)

        do m = 1, l
          t = real((l+1-m)*(l+m), kind=dp)
          cd = cd/t
          yr(k, l, m) = a*sqrt(2.e0_dp*cd)*p(l, m)
        end do
      end do
    end do
  end subroutine gaunt2

end module mod_gaunt2
