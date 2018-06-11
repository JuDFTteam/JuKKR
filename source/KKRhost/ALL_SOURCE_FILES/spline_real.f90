!***********************************************************************
    Subroutine spline_real(nmax, x, y, n, yp1, ypn, y2)
      Use mod_datatypes, Only: dp
! Given arrays x(1:n) and  y(1:n) containing a tabulated function,
! i.e., y i = f(xi), with x1<x2<...<xN , and given values yp1 and ypn
! for the 1rst derivative of the interpolating function at points
! 1 and n, respectively, this routine returns an array y2(1:n) of
! length n which contains the second derivatives of the interpolating
! function at the tabulated points xi.
! If yp1 and/or ypn are equal to 1.e30 or larger, the routine is
! signaled to set the corresponding boundary condition for a natural
! spline, with zero second derivative on that boundary.
! Parameter: NMAX is the largest anticipated value of n.
! Taken from "Numerical Recipes in Fortran 77", W.H.Press et al.
      Implicit None
      Integer :: n, nmax
      Real (Kind=dp) :: yp1, ypn, x(nmax), y(nmax), y2(nmax)
      Integer :: i, k
      Real (Kind=dp) :: p, qn, sig, un, u(nmax)

      If (n>nmax) Stop 'SPLINE: n > NMAX.'
      If (abs(yp1)>0.99E30_dp) Then
! The lower boundary condition is set either to be "natural"
        y2(1) = 0.E0_dp
        u(1) = 0.E0_dp
      Else
! or else to have a specified first derivative.
        y2(1) = -0.5E0_dp
        u(1) = (3.E0_dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      End If

      Do i = 2, n - 1
! This is the decomposition loop of the tridiagonal algorithm. y2 and u
! are used for temporary storage of the decomposed factors.
        sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p = sig*y2(i-1) + 2.E0_dp
        y2(i) = (sig-1.E0_dp)/p
        u(i) = (6.E0_dp*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)- &
          x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      End Do

      If (abs(ypn)>0.99E30_dp) Then
! The upper boundary condition is set either to be "natural"
        qn = 0.E0_dp
        un = 0.E0_dp
      Else
! or else to have a specified 1rst derivative.
        qn = 0.5E0_dp
        un = (3.E0_dp/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      End If
      y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.E0_dp)
      Do k = n - 1, 1, -1
! This is the backsubstitution loop of the tridiagonal algorithm.
        y2(k) = y2(k)*y2(k+1) + u(k)
      End Do

      Return
    End Subroutine
