SUBROUTINE sphere_gga(lmax,yr,wtyr,rij,ijd,lmmaxd,thet,ylm,  &
        dylmt1,dylmt2,dylmf1,dylmf2,dylmtf)
!-----------------------------------------------------------------------
!     generate an angular mesh and spherical harmonics at those
!     mesh points. For an angular integration the weights are ge-
!     rated .

!     R. Zeller      Feb. 1996
!     Small change for GGA implementation
!     Nikos          Dec. 1996
!-----------------------------------------------------------------------
IMPLICIT NONE
!.. Scalar Arguments ..
      INTEGER IJD,LMAX,LMMAXD
!..
!.. Local Scalars ..
      DOUBLE PRECISION DX1,DX2,DX3,F0,PI,R,R1,R2,R3
      INTEGER IJ,LM1
!..
!.. External Subroutines ..
      EXTERNAL CYLM02,YMY
!..
!.. Array Arguments ..
DOUBLE PRECISION DYLMF1(IJD,LMMAXD),DYLMF2(IJD,LMMAXD), &
                 DYLMT1(IJD,LMMAXD),DYLMT2(IJD,LMMAXD), &
                 DYLMTF(IJD,LMMAXD),RIJ(IJD,3),THET(IJD), &
                 WTYR(IJD,*),YLM(IJD,LMMAXD),YR(IJD,*)
!..
!.. Local Arrays ..

DOUBLE PRECISION COSFI(IJD),COSX(IJD),FAI(IJD),ND(3,3),SINFI(IJD), &
                 WGHT,Y(1000)
!..
!.. Intrinsic Functions ..
      INTRINSIC ABS,ACOS,ATAN,COS,SIN,SQRT
!     ..
pi = 4.d0*ATAN(1.d0)
WRITE (6,*) 'SPHERE for GGA: read LEBEDEV mesh'
IF (ijd > 1000) STOP 'SPHERE'


DO  ij = 1,ijd
  CALL lebedev (ij,r1,r2,r3,wght)
  
!      make a small rotation
  
  f0 = 0.08D0
  nd(1,1) = COS(f0)
  nd(1,2) = 0D0
  nd(1,3) = SIN(f0)
  nd(2,1) = 0D0
  nd(2,2) = 1D0
  nd(2,3) = 0D0
  nd(3,1) = -SIN(f0)
  nd(3,2) = 0D0
  nd(3,3) = COS(f0)
  
  dx1 = nd(1,1)*r1 + nd(2,1)*r2 + nd(3,1)*r3
  dx2 = nd(1,2)*r1 + nd(2,2)*r2 + nd(3,2)*r3
  dx3 = nd(1,3)*r1 + nd(2,3)*r2 + nd(3,3)*r3
  
  r1 = dx1
  r2 = dx2
  r3 = dx3
  
  rij(ij,1) = r1
  rij(ij,2) = r2
  rij(ij,3) = r3
  
  CALL ymy(r1,r2,r3,r,y,lmax)
  DO  lm1 = 1, (lmax+1)**2
    yr(ij,lm1) = y(lm1)
  END DO
  
!---> multiply the spherical harmonics with the weights
  
  DO  lm1 = 1, (lmax+1)**2
    wtyr(ij,lm1) = yr(ij,lm1)*wght*pi*4.d0
  END DO
  
!---> produce what is needed for GGA
  
  cosx(ij) = r3
  IF (ABS(r3) /= 1.d0) THEN
    cosfi(ij) = r1/SQRT(1.d0-r3*r3)
    sinfi(ij) = r2/SQRT(1.d0-r3*r3)
    IF (ABS(cosfi(ij)) > 1.d0) cosfi(ij) = cosfi(ij)/ ABS(cosfi(ij))
    IF (ABS(sinfi(ij)) > 1.d0) sinfi(ij) = sinfi(ij)/ ABS(sinfi(ij))
    fai(ij) = ACOS(cosfi(ij))
  ELSE IF (sinfi(ij) == 0.d0) THEN
    fai(ij) = pi/2.d0
  ELSE
    cosfi(ij) = r1
    sinfi(ij) = r2
    IF (ABS(cosfi(ij)) > 1.d0) cosfi(ij) = cosfi(ij)/ ABS(cosfi(ij))
  END IF
  fai(ij) = ACOS(cosfi(ij))
END DO

CALL cylm02(lmax,cosx,fai,2*lmax+1,lmmaxd,thet,ylm,dylmt1,dylmt2,  &
    dylmf1,dylmf2,dylmtf)

END SUBROUTINE sphere_gga
