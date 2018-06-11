!***********************************************************************
    Subroutine gaunt2(w, yr, n)
      Use mod_datatypes, Only: dp
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
      Implicit None
!.. Arguments
      Integer :: n
      Real (Kind=dp) :: w(*), yr(n, 0:n, 0:n)
!..
!.. Local Scalars ..
      Real (Kind=dp) :: a, cd, cth, fac, fpi, rf, sth, t
      Integer :: k, l, lomax, m
!..
!.. Local Arrays ..
      Real (Kind=dp) :: p(0:n+1, 0:n), x(n)
!..
!.. External Subroutines ..
      External :: grule
!..
!.. Intrinsic Functions ..
      Intrinsic :: atan, sqrt
!     ..
      fpi = 16E0_dp*atan(1E0_dp)
      rf = fpi**(1E0_dp/3E0_dp)
      lomax = n

!--->    obtain gauss-legendre points and weights

      Call grule(2*n, x, w)

!--->    generate associated legendre functions for m.ge.0

      Do k = 1, n
        cth = x(k)
        sth = sqrt(1.E0_dp-cth*cth)
        fac = 1.E0_dp

!--->    loop over m values

        Do m = 0, lomax
          fac = -real(2*m-1, kind=dp)*fac
          p(m, m) = fac
          p(m+1, m) = real(2*m+1, kind=dp)*cth*fac

!--->    recurse upward in l

          Do l = m + 2, lomax
            p(l, m) = (real(2*l-1,kind=dp)*cth*p(l-1,m)-real(l+m-1,kind=dp)*p( &
              l-2,m))/real(l-m, kind=dp)
          End Do

          fac = fac*sth
        End Do

!--->    multiply in the normalization factors

        Do l = 0, lomax
          a = rf*sqrt((2*l+1)/fpi)
          cd = 1.E0_dp
          yr(k, l, 0) = a*p(l, 0)

          Do m = 1, l
            t = real((l+1-m)*(l+m), kind=dp)
            cd = cd/t
            yr(k, l, m) = a*sqrt(2.E0_dp*cd)*p(l, m)
          End Do
        End Do
      End Do
    End Subroutine
