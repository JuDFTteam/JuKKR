SUBROUTINE csinwd(f,fint,lmmsqd,irmind,irmd,irmin,ipan,ircut)
!-----------------------------------------------------------------------
!     this subroutine does an inwards integration of llmax
!     functions f with an extended 3-point-simpson :


!                               irmax
!                   fint(ll,i) = { f(ll,i') di'
!                                ir

!  the starting value for this integration at ist - 1 is determined by
!    a 4 point lagrangian integration  , coefficients given by
!    m. abramowitz and i.a. stegun, handbook of mathematical functions,
!    nbs applied mathematics series 55 (1968)

!  attention in case of radial integration :
!       the weights drdi have to be multiplied before calling this
!       subroutine .

!                                     b. drittler mar. 1989

!    modified for functions with kinks - at each kink the integration
!      is restarted

!    attention : it is supposed that irmin + 3 is less than imt !


!                                     b. drittler july 1989
!    modified by m. ogura, june 2015
!-----------------------------------------------------------------------
!     ..
!.. Parameters ..
      DOUBLE PRECISION A1,A2,A3
      PARAMETER (A1=5.D0/12.D0,A2=8.D0/12.D0,A3=-1.D0/12.D0)
!..
!.. Scalar Arguments ..
      INTEGER IPAN,IRMD,IRMIND,LMMSQD,IRMIN
!..
!.. Array Arguments ..
      DOUBLE COMPLEX F(LMMSQD,IRMIND:IRMD),FINT(LMMSQD,IRMIND:IRMD)
      INTEGER IRCUT(0:IPAN)
!..
!.. Local Scalars ..
      INTEGER I,IEN,IP,IST,LL

!---> loop over kinks

DO  ip = ipan,1,-1
  ist = ircut(ip)
  ien = ircut(ip-1) + 1
  IF (ip == 1) ien = irmin
  
  IF (ip == ipan) THEN
    DO  ll = 1,lmmsqd
      fint(ll,ist) = 0.0D0
    END DO
    
  ELSE
    DO  ll = 1,lmmsqd
      fint(ll,ist) = fint(ll,ist+1)
    END DO
  END IF
  
!---> calculate fint with an extended 3-point-simpson
  
  DO  i = ist,ien+2,-2
    DO  ll = 1,lmmsqd
      fint(ll,i-1)=fint(ll,i) +f(ll,i)*a1+f(ll,i-1)*a2+f(ll,i-2)*a3
      fint(ll,i-2)=fint(ll,i-1) +f(ll,i)*a3+f(ll,i-1)*a2+f(ll,i-2)*a1
    END DO
  END DO
  IF(MOD(ist-ien,2) == 1)THEN
    DO  ll=1,lmmsqd
      fint(ll,ien)=fint(ll,ien+1) +f(ll,ien)*a1+f(ll,ien+1)*a2+f(ll,ien+2)*a3
    END DO
  END IF
END DO

END SUBROUTINE csinwd
