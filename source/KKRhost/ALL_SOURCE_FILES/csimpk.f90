SUBROUTINE csimpk(cf,cfint,ipan,ircut,drdi)
!-----------------------------------------------------------------------
!     this subroutine does an integration up to rcut of an
!     complex function cf with an extended 3-point-simpson :

!                             rcut
!                      cfint = { cf(r') dr'
!                              0

!     modified for functions with kinks - at each kink the
!     integration is restarted .

!     attention : input cf is destroyed !

!-----------------------------------------------------------------------
!.. Scalar Arguments ..
      DOUBLE COMPLEX CFINT
      INTEGER IPAN
!..
!.. Array Arguments ..
      DOUBLE COMPLEX CF(*)
      DOUBLE PRECISION DRDI(*)
      INTEGER IRCUT(0:IPAN)
!..
!.. Local Scalars ..
      DOUBLE PRECISION A1,A2
      INTEGER I,IEN,IP,IST,N
!..
!.. External Functions ..
      DOUBLE COMPLEX CSUM
      EXTERNAL CSUM
!..
!.. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
a1 = 4.0D0/3.0D0
a2 = 2.0D0/3.0D0
cfint = 0.0D0

DO  ip = 1,ipan
  
!---> loop over kinks
  
  ist = ircut(ip-1) + 1
  ien = ircut(ip)
  
  DO  i = ist,ien
    cf(i) = cf(i)*drdi(i)
  END DO
  
  IF (MOD(ien-ist,2) == 0) THEN
    cfint = cfint + (cf(ist)-cf(ien))/3.0D0
    ist = ist + 1
    n = (ien-ist+1)/2
    
  ELSE
!---> four point lagrange integration for the first step
    cfint = cfint + (9.0D0*cf(ist)+19.0D0*cf(ist+1)-  &
        5.0D0*cf(ist+2)+cf(ist+3))/24.0D0 + (cf(ist+1)-cf(ien))/3.0D0
    ist = ist + 2
    n = (ien-ist+1)/2
  END IF
  
!---> calculate with an extended 3-point-simpson
  
  cfint = cfint + a1*csum(n,cf(ist),2) + a2*csum(n,cf(ist+1),2)
END DO

END SUBROUTINE csimpk
