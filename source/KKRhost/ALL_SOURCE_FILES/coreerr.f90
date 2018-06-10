SUBROUTINE coreerr(ERR,var,s,nsol,pow,qow,piw,qiw)
!   ********************************************************************
!   *                                                                  *
!   *   CALCULATE THE MISMATCH OF THE RADIAL WAVE FUNCTIONS AT THE     *
!   *   POINT  NMATCH  FOR OUT- AND INWARD INTEGRATION                 *
!   *                                                                  *
!   ********************************************************************
      IMPLICIT NONE
! Dummy arguments
      INTEGER NSOL,S
      REAL*8 ERR(4),PIW(2,2),POW(2,2),QIW(2,2),QOW(2,2),VAR(4)

! Local variables
      INTEGER T

ERR(1) = pow(s,s) - piw(s,s)*var(2)
ERR(2) = qow(s,s) - qiw(s,s)*var(2)

IF ( nsol == 1 ) RETURN

t = 3 - s

ERR(1) = ERR(1) + pow(s,t)*var(3) - piw(s,t)*var(2)*var(4)
ERR(2) = ERR(2) + qow(s,t)*var(3) - qiw(s,t)*var(2)*var(4)
ERR(3) = pow(t,s) - piw(t,s)*var(2) + pow(t,t)*var(3) - piw(t,t)  &
    *var(2)*var(4)
ERR(4) = qow(t,s) - qiw(t,s)*var(2) + qow(t,t)*var(3) - qiw(t,t)  &
    *var(2)*var(4)

END SUBROUTINE coreerr
