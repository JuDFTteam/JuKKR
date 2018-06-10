FUNCTION cnlz(l,z)
!   ********************************************************************
!   *                                                                  *
!   *     von NEUMANN  - FUNCTION  N(L,Z)  FOR COMPLEX ARGUMENT  Z     *
!   *                  see:  e.g. MERZBACHER EQ. (10.34)               *
!   *                                                                  *
!   ********************************************************************
IMPLICIT NONE

! PARAMETER definitions
COMPLEX*16 C1
PARAMETER (C1=(1.0D0,0.0D0))
INTEGER LP2MAX
PARAMETER (LP2MAX=25)

! Dummy arguments
INTEGER L
COMPLEX*16 Z
COMPLEX*16 CNLZ

! Local variables
COMPLEX*16 CJLZ

REAL*8 DFAC
COMPLEX*16 DT,S(LP2MAX),T,ZSQ
INTEGER I,K,LLP1

IF ( l < 0 ) THEN
  cnlz = cjlz(l+1,z)
  RETURN
END IF

zsq = z*z
llp1 = l + l + 1

IF ( ABS(zsq/DBLE(llp1)) <= 10.d0 ) THEN
  
  dfac = 1.0D0
  DO k = 3,llp1,2
    dfac = dfac*DBLE(k)
  END DO
  
  dt = c1
  t = c1
  
  DO i = 2,400,2
    dt = -dt*zsq/DBLE(i*(i-llp1))
    t = t + dt
    IF ( ABS(dt) < 1.0D-10 ) GO TO 50
  END DO
  
  50      CONTINUE
  cnlz = -t*z**(-l-1)*dfac/DBLE(llp1)
  
ELSE
  IF ( l > 23 ) STOP '<cnlz>: l too large'
  s(2) = COS(z)
  IF ( l <= 0 ) THEN
    cnlz = -s(2)*z**(-l-1)
    RETURN
  END IF
  
  s(1) = -SIN(z)/z
  DO i = 3,l + 2
    s(i) = s(i-1)*(2*i-5) - zsq*s(i-2)
  END DO
  cnlz = -s(l+2)*z**(-l-1)
  
END IF

END FUNCTION cnlz
