c***********************************************************************
      SUBROUTINE spline(NMAX,x,y,n,yp1,ypn,y2) 
      implicit none
      INTEGER n,NMAX 
C      REAL*8          yp1,ypn,x(NMAX),y(NMAX),y2(NMAX) 
      REAL*8           x(NMAX)
      COMPLEX          yp1,ypn,y(NMAX),y2(NMAX) 
c Given arrays x(1:n) and  y(1:n) containing a tabulated function, 
c i.e., y i = f(xi), with x1<x2<...<xN , and given values yp1 and ypn 
c for the 1rst derivative of the interpolating function at points 
c 1 and n, respectively, this routine returns an array y2(1:n) of 
c length n which contains the second derivatives of the interpolating 
c function at the tabulated points xi. 
c If yp1 and/or ypn are equal to 1.e30 or larger, the routine is
c signaled to set the corresponding boundary condition for a natural
c spline, with zero second derivative on that boundary. 
c Parameter: NMAX is the largest anticipated value of n. 
c Taken from "Numerical Recipes in Fortran 77", W.H.Press et al.
      INTEGER i,k 
c      REAL*8          p,qn,sig,un,u(NMAX) 
      COMPLEX          p,qn,sig,un,u(NMAX) 

      if (n.gt.nmax) stop 'SPLINE: n > NMAX.'
      if (ABS(yp1).gt.0.99d30) then
c The lower boundary condition is set either to be "natural" 
         y2(1) = 0.d0
         u(1) = 0.d0
      else
c or else to have a specified first derivative. 
         y2(1) = -0.5d0
         u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1) 
      endif 

      do i = 2,n-1  
c This is the decomposition loop of the tridiagonal algorithm. y2 and u 
c are used for temporary storage of the decomposed factors. 
         sig = (x(i)-x(i-1)) / (x(i+1)-x(i-1))
         p = sig * y2(i-1) + 2.d0 
         y2(i) = (sig-1.d0)/p
         u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) 
     &        /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1)) / p 
      enddo

      if (ABS(ypn) .gt.0.99d30) then   
c The upper boundary condition is set either to be "natural"
         qn = 0.d0
         un = 0.d0
      else
c or else to have a specified 1rst derivative. 
         qn = 0.5d0
         un = (3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n) = (un-qn*u(n-1)) / (qn*y2(n-1)+1.d0) 
      do k = n-1,1,-1 
c This is the backsubstitution loop of the tridiagonal algorithm. 
         y2(k)=y2(k)*y2(k+1)+u(k) 
      enddo

      return 
      END
