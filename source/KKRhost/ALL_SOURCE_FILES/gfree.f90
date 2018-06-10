! **********************************************************************
SUBROUTINE gfree13(rdiff,e0,gmll,dgmll,cleb,icleb,loflm,iend)
! **********************************************************************
!     .. Parameters ..
      IMPLICIT NONE
      include 'inc.p'
      INTEGER LMAX
      PARAMETER (LMAX=LMAXD)
      INTEGER LMGF0D
      PARAMETER (LMGF0D= (LMAX+1)**2)
      INTEGER LMAX2P,LMX2SQ
      PARAMETER (LMAX2P=LMAX*2+1,LMX2SQ=LMAX2P**2)
      DOUBLE COMPLEX CZERO,CI
      PARAMETER (CZERO= (0.0D0,0.0D0),CI= (0.0D0,1.0D0))
!     ..
!     .. Scalar Arguments ..
      DOUBLE COMPLEX E0
      INTEGER IEND
!     ..
!     .. Array Arguments ..
      DOUBLE COMPLEX GMLL(LMGF0D,LMGF0D)
      DOUBLE COMPLEX DGMLL(LMGF0D,LMGF0D) ! LLY
      DOUBLE PRECISION CLEB(NCLEB),RDIFF(*)
      INTEGER ICLEB(NCLEB,4),LOFLM(*)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION FPI,PI,RABS,RFPI,X,Y,Z
      INTEGER IFAC,J,LM1,LM2,LM3
!     ..
!     .. Local Arrays ..
      DOUBLE COMPLEX BL(LMAX2P),HL(LMAX2P),HYL(LMX2SQ),NL(LMAX2P)
      DOUBLE COMPLEX DHL(LMAX2P),DHYL(LMX2SQ) ! LLY
      DOUBLE PRECISION YL(LMX2SQ)
      INTEGER LF(LMX2SQ)
!     ..
!     .. External Subroutines ..
      EXTERNAL BESHAN,YMY
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ATAN,SQRT
!     ..
pi = 4.d0*ATAN(1.d0)
fpi = 4.d0*pi
rfpi = SQRT(fpi)
!-----------------------------------------------------------------------
!---- CALCULATION OF FREE ELECTRON GREEN'S FUNCTION :  G(M)LL'(E0)
!-----------------------------------------------------------------------
!     Also:
!---- derivative of free electron green function matrix elements
!     returned in array DGMLL=dGMLL / dE
!
!     the analytical formula for the derivative of spherical Hankel
!     functions is used:
!
!     d                     l+1
!     --  h (x) = h   (x) - --- h (x)
!     dx   l       l-1       x   l
!
!     which for x = sqrt(E0)*r leads to
!
!      d                       r           rl
!     --- ( sqrt(E0) h (x) ) = - h   (x) - -- h (x) )
!     dE0             l        2  l-1      2x  l
!
!     Ported from KKRnano by Phivos Mavropoulos 10.10.2013
!-----------------------------------------------------------------------
DO  lm1 = 1,lmx2sq
  lf(lm1) = loflm(lm1) + 1
END DO
x = rdiff(1)
y = rdiff(2)
z = rdiff(3)
CALL ymy(x,y,z,rabs,yl,lmax*2)
CALL beshan(hl,bl,nl,SQRT(e0)*rabs,lmax*2)

! Derivatives of Hankel functions ! LLY
dhl(1) = 0.5D0*ci*rabs*hl(1)
DO lm1 = 2,lmax2p
  dhl(lm1) = 0.5D0*(rabs*hl(lm1-1)-(lm1-1)*hl(lm1)/SQRT(e0))
END DO

DO  lm1 = 1,lmx2sq
  hyl(lm1) = -fpi*ci*SQRT(e0)*yl(lm1)*hl(lf(lm1))
  dhyl(lm1) = -fpi*ci*yl(lm1)*dhl(lf(lm1))   ! LLY
END DO


DO  lm1 = 1,lmgf0d
  gmll(lm1,lm1) = hyl(1)/rfpi
  dgmll(lm1,lm1) = dhyl(1)/rfpi  ! LLY
  DO  lm2 = 1,lm1 - 1
    gmll(lm1,lm2) = czero
    dgmll(lm1,lm2) = czero    ! LLY
  END DO
END DO

DO  j = 1,iend
  lm1 = icleb(j,1)
  lm2 = icleb(j,2)
  lm3 = icleb(j,3)
  gmll(lm1,lm2) = gmll(lm1,lm2) + cleb(j)*hyl(lm3)
  dgmll(lm1,lm2) = dgmll(lm1,lm2) + cleb(j)*dhyl(lm3)  ! LLY
END DO
DO  lm1 = 1,lmgf0d
  DO  lm2 = 1,lm1 - 1
    ifac = (-1)** (loflm(lm1)+loflm(lm2))
    gmll(lm2,lm1) = ifac*gmll(lm1,lm2)
    dgmll(lm2,lm1) = ifac*dgmll(lm1,lm2)  ! LLY
  END DO
END DO
RETURN

END SUBROUTINE gfree13
