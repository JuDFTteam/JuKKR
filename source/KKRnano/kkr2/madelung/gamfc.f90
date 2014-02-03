subroutine gamfc(alpha,glh,lmax,r)
implicit none
!----------------------------------------------------------------------
!
!      calculation of convergence function
!
!       glh = i(alpha,l)/r**(l+1)*sqrt(pi)
!
!      with
!            alpha = r times the splitting paramter lamda
!      and
!            i(x,l) = erfc(x) + exp(-x*x)/sqrt(pi) *
!
!                                sum ( 2**i * x**(2i-1) / (2i-1)!! )
!                              1..i..l
!
!
! Note: gamfc( alpha -> 0, ... ) => glh = sqrt(pi) for all l  E.R.
!
!-----------------------------------------------------------------------
!     .. scalar arguments ..
double precision alpha,r
integer lmax
!     ..
!     .. array arguments ..
double precision glh(0:lmax)
!     ..
!     .. local scalars ..
double precision arg,facl,fex
integer l
!     ..
!     .. external functions ..
double precision erfcex
external erfcex
!     ..
!     .. intrinsic functions ..
intrinsic exp,real
!     ..
arg = alpha*alpha
glh(0) = erfcex(alpha)
facl = 2.0d0*alpha
!
!---> recursion
!
do 10 l = 1,lmax
   glh(l) = glh(l-1) + facl
   facl = facl*arg/ (real(l)+0.5d0)
10 continue

!     changed 21/10/99 
!     if arg is to big then cannot calculate 1/exp(arg) !!

fex = exp(-arg)

do 20 l = 0,lmax
   fex = fex/r
   glh(l) = glh(l)*fex
20 continue

end
