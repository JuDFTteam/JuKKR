! **********************************************************************
subroutine grule(n, x, w)

!***********************************************************************

!     determines the (n+1)/2 nonnegative points x(i) and
!     the corresponding weights w(i) of the n-point
!     gauss-legendre integration rule, normalized to the
!     interval [-1,1]. the x(i) appear in descending order.

!     this routine is from 'methods of numerical integration',
!     p.j. davis and p. rabinowitz, page 369.

!***********************************************************************
      Use mod_datatypes, Only: dp
implicit none

!.. Scalar Arguments ..
  integer :: n
!..
!.. Array Arguments ..
  real (kind=dp) :: w(*), x(*)
!..
!.. Local Scalars ..
  real (kind=dp) :: d1, d2pn, d3pn, d4pn, den, d_p, d_pn, e1, fx, h, p, pi, &
    pk, pkm1, pkp1, t, t1, u, v, x0
  integer :: i, it, k, m
!..
!.. Intrinsic Functions ..
  intrinsic :: cos, atan

  pi = 4.d0*atan(1.d0)
  m = (n+1)/2
  e1 = n*(n+1)
  do i = 1, m
    t = (4*i-1)*pi/(4*n+2)
    x0 = (1.0d0-(1.0d0-1.0d0/n)/(8.0d0*n*n))*cos(t)

!--->    iterate on the value  (m.w. jan. 1982)

    do it = 1, 2
      pkm1 = 1.
      pk = x0
      do k = 2, n
        t1 = x0*pk
        pkp1 = t1 - pkm1 - (t1-pkm1)/k + t1
        pkm1 = pk
        pk = pkp1
      end do
      den = 1. - x0*x0
      d1 = n*(pkm1-x0*pk)
      d_pn = d1/den
      d2pn = (2.*x0*d_pn-e1*pk)/den
      d3pn = (4.*x0*d2pn+(2.-e1)*d_pn)/den
      d4pn = (6.*x0*d3pn+(6.-e1)*d2pn)/den
      u = pk/d_pn
      v = d2pn/d_pn
      h = -u*(1.+.5*u*(v+u*(v*v-u*d3pn/(3.*d_pn))))
      p = pk + h*(d_pn+.5*h*(d2pn+h/3.*(d3pn+.25*h*d4pn)))
      d_p = d_pn + h*(d2pn+.5*h*(d3pn+h*d4pn/3.))
      h = h - p/d_p
      x0 = x0 + h
    end do
    x(i) = x0
    fx = d1 - h*e1*(pk+.5*h*(d_pn+h/3.*(d2pn+.25*h*(d3pn+.2*h*d4pn))))
    w(i) = 2.*(1.-x(i)*x(i))/(fx*fx)
  end do
  if (m+m>n) x(m) = 0.
end subroutine
