SUBROUTINE dirbsrad(xbs,y,dydx,drdi,b,v,r,nmesh)
!   ********************************************************************
!   *                                                                  *
!   *   supply the derivatives for the coupled set of                  *
!   *   radial dirac equation in case of a spin-dependent potential    *
!   *                                                                  *
!   ********************************************************************

IMPLICIT NONE

INCLUDE 'sprkkr_rmesh.dim'


! PARAMETER definitions
      INTEGER NLAG
      PARAMETER (NLAG=3)

! COMMON variables
      REAL*8 CGD(2),CGMD(2),CGO,KAP(2)
      REAL*8 CSQR
      COMPLEX*16 EBS
      INTEGER NRADBS,NSOLBS
      COMMON /COMMBS/ EBS,CSQR,CGD,CGMD,CGO,KAP,NSOLBS,NRADBS

! Dummy arguments
      INTEGER NMESH
      REAL*8 XBS
      REAL*8 B(NRMAX),DRDI(NRMAX),R(NRMAX),V(NRMAX)
      COMPLEX*16 DYDX(NCFMAX),Y(NCFMAX)

! Local variables
      REAL*8 BBS,BPP,BQQ,DRDIBS,RBS,VBS
      COMPLEX*16 EMVPP,EMVQQ
      INTEGER I,J,K,M

CALL dirbslag(xbs,vbs,bbs,rbs,drdibs,v,b,r,drdi,nradbs,nlag,nmesh)

emvqq = (ebs-vbs+csqr)*drdibs/csqr
emvpp = -(ebs-vbs)*drdibs
bqq = bbs*drdibs/csqr
bpp = bbs*drdibs


m = 0
DO j = 1,nsolbs
  DO i = 1,nsolbs
    k = 3 - 4*(i-1)
    dydx(m+1) = -kap(i)*y(m+1)/rbs*drdibs + (emvqq+bqq*cgmd(i)) *y(m+2)
    dydx(m+2) = kap(i)*y(m+2)/rbs*drdibs + (emvpp+bpp*cgd(i))  &
        *y(m+1) + bpp*cgo*y(m+k)
    m = m + 2
  END DO
END DO
END SUBROUTINE dirbsrad
