!***********************************************************************
    Subroutine splint_real(xa, ya, y2a, n, x, y, yderiv)
      Use mod_datatypes, Only: dp
! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a
! function (with the xai's in order), and given the array y2a(1:n), which
! is the output from spline above, and given a value of x, this routine
! returns a cubic-spline interpolated value y and the derivative yderiv.
! Taken from "Numerical Recipes in Fortran 77", W.H.Press et al.
! We will  nd the right place in the table by means of bisection.
! This is optimal if sequential calls to this routine are at random
! values of x. If sequential calls are in order, and closely
! spaced, one would do better to store previous values of
! klo and khi and test if they remain appropriate on the
! next call.
      Implicit None
      Integer :: n
      Real (Kind=dp) :: x, y, yderiv, xa(*), ya(*), y2a(*)
      Integer :: k, khi, klo
      Real (Kind=dp) :: a, b, h

      klo = 1
      khi = n
100   If (khi-klo>1) Then
        k = (khi+klo)/2
        If (xa(k)>x) Then
          khi = k
        Else
          klo = k
        End If
        Go To 100
      End If
! klo and khi now bracket the input value of x.
      h = xa(khi) - xa(klo)
! The xa's must be distinct.
      If (h==0.E0_dp) Stop 'bad xa input in splint'
! Cubic spline polynomial is now evaluated.
      a = (xa(khi)-x)/h
      b = (x-xa(klo))/h
      y = a*ya(klo) + b*ya(khi) + ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2) &
        /6.E0_dp
      yderiv = (ya(khi)-ya(klo))/h - ((3.E0_dp*a*a-1.E0_dp)*y2a(klo)-(3.E0_dp* &
        b*b-1.E0_dp)*y2a(khi))*h/6.E0_dp

      Return
    End Subroutine
