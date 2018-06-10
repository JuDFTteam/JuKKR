SUBROUTINE dirbsmid(y,dydx,nv,xs,htot,nstep,yout,b,v,r,drdi,nmesh)
!   ********************************************************************
!   *                                                                  *
!   *   modified midpoint step to support the  Burlisch-Stoer method   *
!   *   on exit:  the incremented variable is in   YOUT                *
!   *                                                                  *
!   *   see: numerical recipes chapter 15.3                            *
!   *                                                                  *
!   ********************************************************************

IMPLICIT NONE

INCLUDE 'sprkkr_rmesh.dim'

! Dummy arguments
REAL*8 HTOT,XS
INTEGER NMESH,NSTEP,NV
REAL*8 B(NRMAX),DRDI(NRMAX),R(NRMAX),V(NRMAX)
COMPLEX*16 DYDX(NV),Y(NV),YOUT(NV)

! Local variables
REAL*8 H,H2,X
INTEGER I,N
COMPLEX*16 SWAP,YM(NCFMAX),YN(NCFMAX)

h = htot/nstep
DO i = 1,nv
  ym(i) = y(i)
  yn(i) = y(i) + h*dydx(i)
END DO
x = xs + h

CALL dirbsrad(x,yn,yout,drdi,b,v,r,nmesh)

h2 = 2.d0*h
DO n = 2,nstep
  DO i = 1,nv
    swap = ym(i) + h2*yout(i)
    ym(i) = yn(i)
    yn(i) = swap
  END DO
  x = x + h
  CALL dirbsrad(x,yn,yout,drdi,b,v,r,nmesh)
END DO
DO i = 1,nv
  yout(i) = 0.5D0*(ym(i)+yn(i)+h*yout(i))
  
END DO

END SUBROUTINE dirbsmid
