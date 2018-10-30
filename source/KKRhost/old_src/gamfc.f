      subroutine gamfc(alpha,glh,lmax,r)
      implicit none
c----------------------------------------------------------------------
c
c      calculation of convergence function
c
c       glh = i(alpha,l)/r**(l+1)*sqrt(pi)
c
c      with
c            alpha = r times the splitting paramter lamda
c      and
c            i(x,l) = erfc(x) + exp(-x*x)/sqrt(pi) *
c
c                                sum ( 2**i * x**(2i-1) / (2i-1)!! )
c                              1..i..l
c
c-----------------------------------------------------------------------
c     .. scalar arguments ..
      double precision alpha,r
      integer lmax
c     ..
c     .. array arguments ..
      double precision glh(0:lmax)
c     ..
c     .. local scalars ..
      double precision arg,facl,fex
      integer l
c     ..
c     .. external functions ..
      double precision erfcex
      external erfcex
c     ..
c     .. intrinsic functions ..
      intrinsic exp,real
c     ..
      arg = alpha*alpha
      glh(0) = erfcex(alpha)
      facl = 2.0d0*alpha
c
c---> recursion
c
      do 10 l = 1,lmax
         glh(l) = glh(l-1) + facl
         facl = facl*arg/ (real(l)+0.5d0)
   10 continue

c     changed 21/10/99 
c     if arg is to big then cannot calculate 1/exp(arg) !!

      fex = exp(-arg)

      do 20 l = 0,lmax
         fex = fex/r
         glh(l) = glh(l)*fex
 20   continue
 
      end
