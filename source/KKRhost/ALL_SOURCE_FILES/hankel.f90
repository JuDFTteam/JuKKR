subroutine hankel(h, l, arg)
!  this subroutine uses the explicit formulas for the hankel
!  functions. for higher l-values these formulas may lead to
!  loss of significant figures. This subroutine should be used
!  only for core states.
  implicit none
!.. Scalar Arguments ..
  double complex :: arg
  integer :: l
!..
!.. Array Arguments ..
  double complex :: h(*)
!..
!.. Local Scalars ..
  double complex :: a1, a2, a3, a4
!..
!.. Intrinsic Functions ..
  intrinsic :: exp
!..
!.. Parameters ..
  double complex :: ci
  parameter (ci=(0.0d0,1.0d0))
!     ..
  h(1) = -exp(arg*ci)/arg
  if (l/=1) then
    a1 = (1.d0, 0.d0) - arg*ci
    h(2) = h(1)*a1/arg
    if (l/=2) then
      a1 = 3.d0*a1
      a2 = arg*arg
      h(3) = h(1)*(a1-a2)/a2
      if (l/=3) then
        a1 = 5.d0*a1
        a3 = a2*arg*ci
        a4 = a2*arg
        a2 = 6.d0*a2
        h(4) = h(1)*(a1-a2+a3)/a4
        if (l/=4) then
          a1 = 7.d0*a1
          a2 = 7.5d0*a2
          a3 = 10.d0*a3
          a4 = a4*arg
          h(5) = h(1)*(a1-a2+a3+a4)/a4
          if (l/=5) then
            h(6) = (9.0d0, 0.0d0)*h(5)/arg - h(4)
            if (l/=6) then
              write (6, fmt=100) l
              stop 'HANKEL'

            end if

          end if

        end if

      end if

    end if

  end if

  return


100 format (2x, ' hankel :  l=', i2, ' is too large')
end subroutine
