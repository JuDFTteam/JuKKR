SUBROUTINE gamfc(alpha,glh,lmax,r)
!----------------------------------------------------------------------

!      calculation of convergence function

!       glh = i(alpha,l)/r**(l+1)*sqrt(pi)

!      with
!            alpha = r times the splitting paramter lamda
!      and
!            i(x,l) = erfc(x) + exp(-x*x)/sqrt(pi) *

!                                sum ( 2**i * x**(2i-1) / (2i-1)!! )
!                              1..i..l

!-----------------------------------------------------------------------
implicit none
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
facl = 2.0D0*alpha

!---> recursion

DO  l = 1,lmax
  glh(l) = glh(l-1) + facl
  facl = facl*arg/ (REAL(l)+0.5D0)
END DO

!     changed 21/10/99
!     if arg is to big then cannot calculate 1/exp(arg) !!

fex = EXP(-arg)

DO  l = 0,lmax
  fex = fex/r
  glh(l) = glh(l)*fex
END DO

END SUBROUTINE gamfc
