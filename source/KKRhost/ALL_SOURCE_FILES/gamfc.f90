    Subroutine gamfc(alpha, glh, lmax, r)
      Use mod_datatypes, Only: dp
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
      Implicit None
!     .. scalar arguments ..
      Real (Kind=dp) :: alpha, r
      Integer :: lmax
!     ..
!     .. array arguments ..
      Real (Kind=dp) :: glh(0:lmax)
!     ..
!     .. local scalars ..
      Real (Kind=dp) :: arg, facl, fex
      Integer :: l
!     ..
!     .. external functions ..
      Real (Kind=dp) :: erfcex
      External :: erfcex
!     ..
!     .. intrinsic functions ..
      Intrinsic :: exp, real
!     ..
      arg = alpha*alpha
      glh(0) = erfcex(alpha)
      facl = 2.0E0_dp*alpha

!---> recursion

      Do l = 1, lmax
        glh(l) = glh(l-1) + facl
        facl = facl*arg/(real(l,kind=dp)+0.5E0_dp)
      End Do

!     changed 21/10/99
!     if arg is to big then cannot calculate 1/exp(arg) !!

      fex = exp(-arg)

      Do l = 0, lmax
        fex = fex/r
        glh(l) = glh(l)*fex
      End Do

    End Subroutine
