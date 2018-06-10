SUBROUTINE fplaneg(lamda,g,pref,lmax,ga,vol)
! **************************************************************************
! This sub calculates the derivatives of the inverce
! space contribution to the ewald sum


!     l     (                                                      )
!    d   pi*( exp(gz)*erfc(g/lam/lam + 2*z)*lam/2                  )  /
!           (             +   exp(-gz)*erfc(g/lam/lam - 2*z)*lam/2 ) / (g Vol)
!    ---    (                                                      )
!      l
!    dz

!   And the limit z -> 0 is taken (lam is the lamda parameter

! *********************************************************************
      implicit none
      integer lmax
      double precision alpha,g(0:4),pref(0:lmax)
      integer l
      double precision derfc,lamda,er,ex,pi,ga,vol
      double precision sqpi

DO l=0,4
  g(l) = 0.d0
END DO
DO l=0,lmax
  pref(l) = 0.d0
END DO
pi = 4.d0*ATAN(1.d0)
sqpi=SQRT(pi)
alpha = ga/2.d0/lamda
er = derfc(alpha)
ex = EXP(-alpha*alpha)

IF (ABS(ga) > 1.d-6) THEN
  g(0) = 2.d0*pi/vol*er/ga
ELSE
  g(0) = 0.d0
END IF

pref(2) = SQRT(5.d0/pi)/2.d0/2.d0
g(2) = pref(2)/vol*( 2.d0*pi*ga*er - ex*4.d0*sqpi*lamda)

pref(4) = 3.d0*SQRT(9.d0/pi)/16.d0/9.d0
g(4) = pref(4)/vol*(2.d0*pi*ga*ga*ga*er - ex*4.d0*sqpi*  &
    (lamda*ga*ga - 2.d0*lamda**3))

!      pref(6) = sqrt(13.d0/pi)/32.d0/13.d0
!      g(6) = pref(6)/vol*(2.d0*pi*ga**5*er - ex*4.d0*sqpi*
!     &                               (lamda*ga**4
!     &                               -2.d0*ga*ga*lamda*lamda
!     &                               +12.d0*lamda**5))
END SUBROUTINE fplaneg
