c***********************************************************************
      SUBROUTINE splint_real(xa,ya,y2a,n,x,y,yderiv)
      implicit none
      INTEGER n
      REAL*8         x,y,yderiv,xa(*),ya(*),y2a(*)
c Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a
c function (with the xai's in order), and given the array y2a(1:n), which
c is the output from spline above, and given a value of x, this routine
c returns a cubic-spline interpolated value y and the derivative yderiv.
c Taken from "Numerical Recipes in Fortran 77", W.H.Press et al.
      INTEGER k,khi,klo
      REAL*8         a,b,h
c We will  nd the right place in the table by means of bisection.
c This is optimal if sequential calls to this routine are at random
c values of x. If sequential calls are in order, and closely
c spaced, one would do better to store previous values of
c klo and khi and test if they remain appropriate on the
c next call.
      klo=1
      khi=n
 1    if (khi-klo.gt.1) then
         k=(khi+klo)/2
         if(xa(k).gt.x)then
            khi=k
         else
            klo=k
         endif
         goto 1
      endif
c klo and khi now bracket the input value of x.
      h=xa(khi)-xa(klo)
c The xa's must be distinct.
      if (h.eq.0.d0) pause 'bad xa input in splint_real'
c Cubic spline polynomial is now evaluated.
      a = (xa(khi)-x)/h
      b = (x-xa(klo))/h
      y = a*ya(klo) + b*ya(khi) +
     &     ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi)) * (h**2)/6.d0
      yderiv = (ya(khi)-ya(klo))/h - 
     &     ((3.d0*a*a-1.d0)*y2a(klo) - (3.d0*b*b-1.d0)*y2a(khi))*h/6.d0

      return
      END
