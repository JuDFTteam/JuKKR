SUBROUTINE fplaner(alpha,g,r)
! ************************************************
! This sub calculates the derivatives of the real
! space contribution to the ewald sum .

!              l
!             d     erfc(lamda*sqrt(d*d+z*z))
!      lim    --   ------------------------
!      z->0     l        sqrt(d*d+z*z)
!             dz

!  Up to l = 4 (l=1,3,5,7 etc vanish)



! ************************************************

implicit none
double precision alpha,g(0:4),r
integer l
double precision derfc,lamda,er,ex,pi,pref,sqpi

DO l=0,4
  g(l) = 0.d0
END DO
pi = 4.d0*ATAN(1.d0)
sqpi=SQRT(pi)
er = derfc(alpha)
ex = EXP(-alpha*alpha)
lamda = alpha/r

g(0) = er/r

pref = SQRT(5.d0/pi)/4.d0
g(2) = -pref*( er/r/r/r  +  ex*2.d0*lamda/r/r/sqpi)

pref = 3.d0*SQRT(9.d0/pi)/16.d0/9.d0
g(4) = pref*(9.d0*er + ex*(12.d0*alpha**3+18.d0*alpha)/sqpi) /r/r/r/r/r

END SUBROUTINE fplaner
