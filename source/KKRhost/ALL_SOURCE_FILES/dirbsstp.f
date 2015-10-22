      SUBROUTINE DIRBSSTP(Y,DYDX,NV,X,HTRY,EPS,YSCAL,B,V,R,DRDI,NMESH)
C
C   ********************************************************************
C   *                                                                  *
C   *   Burlisch-Stoer step with monitoring of local truncation error  *
C   *   on entry: X,Y,DXDY  for last mesh-point                        *
C   *   on exit:  X,Y,DXDY  updated for X = X(last) + HTRY             *
C   *                                                                  *
C   *   see: numerical recipes chapter 15.4                            *
C   *                                                                  *
C   *   note: don't set NUSE    > NUSEMAX in <DIRBSRZE>                *
C   *         don't set ISEQMAX > ISEQMAX in <DIRBSRZE>                *
C   *         no step size adjusted in case of no convergency > STOP   *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE

      INCLUDE 'sprkkr_rmesh.dim'
C
C PARAMETER definitions
C
      INTEGER ISEQMAX,NUSE
      PARAMETER (ISEQMAX=30,NUSE=7)
      COMPLEX*16 TINY
C Bereshad:
C      parameter (  tiny  = (1.0d-20,1.0d-20)  )
C dadurch wird der Imaginaerteil fuer die Skalierungsfunktion nicht 00 klein.
C Ich hatte da spruenge in dem errmax; das verlaeuft jetzt stetig, hat aber mit
C den kleinen Werten doch so seine Konvergenzprobleme.
C
      PARAMETER (TINY=(1.0D-20,1.0D-20))
C
C Dummy arguments
C
      REAL*8 EPS,HTRY,X
      INTEGER NMESH,NV
      REAL*8 B(NRMAX),DRDI(NRMAX),R(NRMAX),V(NRMAX)
      COMPLEX*16 DYDX(NV),Y(NV),YSCAL(NV)
C
C Local variables
C
      
      COMPLEX*16 DYSAV(NCFMAX),YERR(NCFMAX),YSAV(NCFMAX),YSEQ(NCFMAX)
      REAL*8 ERRMAX,H,XEST,XSAV
      INTEGER I,J,NSEQ(ISEQMAX)
C
      DATA NSEQ/2,4,6,8,12,16,24,32,48,64,96,128,192,256,384,512,768,
     &     1024,1536,2048,3072,4096,
     &     6144, 8192,12288,16384,24576,32768,49152,65536/
C
      H = HTRY
      XSAV = X
C
      DO I = 1,NV
         YSAV(I) = Y(I)
         DYSAV(I) = DYDX(I)
      END DO
      DO I = 1,ISEQMAX
C
         CALL DIRBSMID(YSAV,DYSAV,NV,XSAV,H,NSEQ(I),YSEQ,
     &                 B,V,R,DRDI,NMESH)
         XEST = (H/NSEQ(I))**2
C
         CALL DIRBSRZE(I,XEST,YSEQ,Y,YERR,NV,NUSE)
C
         ERRMAX = 0.0D0
         DO J = 1,NV
            ERRMAX = DBLE(MAX(ERRMAX,ABS(YERR(J)/(YSCAL(J)+TINY))))
         END DO
         ERRMAX = ERRMAX/EPS
         IF ( ERRMAX.LT.1.0D0 ) THEN
            X = X + H
C
            CALL DIRBSRAD(X,Y,DYDX,DRDI,B,V,R,NMESH)
C
C
            RETURN
         END IF
      END DO

      IF ( ERRMAX.LT.1000D0 ) THEN
          X = X + H
          CALL DIRBSRAD(X,Y,DYDX,DRDI,B,V,R,NMESH)
          WRITE (1337,*) '<DIRBSSTP>  not converged after ',ISEQMAX,
     &         ' refinements'
          WRITE (1337,*) 'step size will not be adjusted !!!!!!'
          WRITE (1337,*) 'max. relative error : ',ERRMAX*EPS
          WRITE (1337,*) 'tolerance             ',EPS
          WRITE (1337,*) 'grid position  X      ',X
      ELSE 
          WRITE (6,*) '<DIRBSSTP>  not converged after ',ISEQMAX,
     &         ' refinements'
          WRITE (6,*) 'step size will not be adjusted !!!!!!'
          WRITE (6,*) 'max. relative error : ',ERRMAX*EPS
          WRITE (6,*) 'tolerance             ',EPS
          WRITE (6,*) 'grid position  X      ',X
          STOP 
      END IF
C
      
      RETURN
c     STOP
      END 
