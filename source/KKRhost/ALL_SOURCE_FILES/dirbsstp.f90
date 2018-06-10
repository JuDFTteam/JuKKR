SUBROUTINE dirbsstp(y,dydx,nv,x,htry,eps,yscal,b,v,r,drdi,nmesh)
!   ********************************************************************
!   *                                                                  *
!   *   Burlisch-Stoer step with monitoring of local truncation error  *
!   *   on entry: X,Y,DXDY  for last mesh-point                        *
!   *   on exit:  X,Y,DXDY  updated for X = X(last) + HTRY             *
!   *                                                                  *
!   *   see: numerical recipes chapter 15.4                            *
!   *                                                                  *
!   *   note: don't set NUSE    > NUSEMAX in <DIRBSRZE>                *
!   *         don't set ISEQMAX > ISEQMAX in <DIRBSRZE>                *
!   *         no step size adjusted in case of no convergency > STOP   *
!   *                                                                  *
!   ********************************************************************

use mod_types, only: t_inc
IMPLICIT NONE

INCLUDE 'sprkkr_rmesh.dim'

! PARAMETER definitions
INTEGER ISEQMAX,NUSE
PARAMETER (ISEQMAX=30,NUSE=7)
COMPLEX*16 TINY
! Bereshad:
!      parameter (  tiny  = (1.0d-20,1.0d-20)  )
! dadurch wird der Imaginaerteil fuer die Skalierungsfunktion nicht 00 klein.
! Ich hatte da spruenge in dem errmax; das verlaeuft jetzt stetig, hat aber mit
! den kleinen Werten doch so seine Konvergenzprobleme.

PARAMETER (TINY=(1.0D-20,1.0D-20))

! Dummy arguments
REAL*8 EPS,HTRY,X
INTEGER NMESH,NV
REAL*8 B(NRMAX),DRDI(NRMAX),R(NRMAX),V(NRMAX)
COMPLEX*16 DYDX(NV),Y(NV),YSCAL(NV)

! Local variables
COMPLEX*16 DYSAV(NCFMAX),YERR(NCFMAX),YSAV(NCFMAX),YSEQ(NCFMAX)
REAL*8 ERRMAX,H,XEST,XSAV
INTEGER I,J,NSEQ(ISEQMAX)

DATA nseq/2,4,6,8,12,16,24,32,48,64,96,128,192,256,384,512,768,  &
    1024,1536,2048,3072,4096, 6144, 8192,12288,16384,24576,32768,49152,65536/

h = htry
xsav = x

DO i = 1,nv
  ysav(i) = y(i)
  dysav(i) = dydx(i)
END DO
DO i = 1,iseqmax
  
  CALL dirbsmid(ysav,dysav,nv,xsav,h,nseq(i),yseq, b,v,r,drdi,nmesh)
  xest = (h/nseq(i))**2
  
  CALL dirbsrze(i,xest,yseq,y,yerr,nv,nuse)
  
  errmax = 0.0D0
  DO j = 1,nv
    errmax = DBLE(MAX(errmax,ABS(yerr(j)/(yscal(j)+tiny))))
  END DO
  errmax = errmax/eps
  IF ( errmax < 1.0D0 ) THEN
    x = x + h
    
    CALL dirbsrad(x,y,dydx,drdi,b,v,r,nmesh)
    
    
    RETURN
  END IF
END DO

IF ( errmax < 1000D0 ) THEN
  x = x + h
  CALL dirbsrad(x,y,dydx,drdi,b,v,r,nmesh)
  IF(t_inc%i_write>0) THEN
    WRITE (1337,*) '<DIRBSSTP>  not converged after ',iseqmax, ' refinements'
    WRITE (1337,*) 'step size will not be adjusted !!!!!!'
    WRITE (1337,*) 'max. relative error : ',errmax*eps
    WRITE (1337,*) 'tolerance             ',eps
    WRITE (1337,*) 'grid position  X      ',x
  END IF
ELSE
  IF(t_inc%i_write>0) THEN
    WRITE (6,*) '<DIRBSSTP>  not converged after ',iseqmax, ' refinements'
    WRITE (6,*) 'step size will not be adjusted !!!!!!'
    WRITE (6,*) 'max. relative error : ',errmax*eps
    WRITE (6,*) 'tolerance             ',eps
    WRITE (6,*) 'grid position  X      ',x
  END IF
  STOP
END IF


RETURN
END SUBROUTINE dirbsstp
