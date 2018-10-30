      SUBROUTINE DIRBSMID(Y,DYDX,NV,XS,HTOT,NSTEP,YOUT,B,V,R,DRDI,NMESH)
C
C   ********************************************************************
C   *                                                                  *
C   *   modified midpoint step to support the  Burlisch-Stoer method   *
C   *   on exit:  the incremented variable is in   YOUT                *
C   *                                                                  *
C   *   see: numerical recipes chapter 15.3                            *
C   *                                                                  *
C   ********************************************************************
C

      IMPLICIT NONE

      INCLUDE 'sprkkr_rmesh.dim'

C
C Dummy arguments
C
      REAL*8 HTOT,XS
      INTEGER NMESH,NSTEP,NV
      REAL*8 B(NRMAX),DRDI(NRMAX),R(NRMAX),V(NRMAX)
      COMPLEX*16 DYDX(NV),Y(NV),YOUT(NV)
C
C Local variables
C
      REAL*8 H,H2,X
      INTEGER I,N
      COMPLEX*16 SWAP,YM(NCFMAX),YN(NCFMAX)
C
      SAVE
C
      H = HTOT/NSTEP
      DO I = 1,NV
         YM(I) = Y(I)
         YN(I) = Y(I) + H*DYDX(I)
      END DO
      X = XS + H
C
      CALL DIRBSRAD(X,YN,YOUT,DRDI,B,V,R,NMESH)
C
      H2 = 2.D0*H
      DO N = 2,NSTEP
         DO I = 1,NV
            SWAP = YM(I) + H2*YOUT(I)
            YM(I) = YN(I)
            YN(I) = SWAP
         END DO
         X = X + H
         CALL DIRBSRAD(X,YN,YOUT,DRDI,B,V,R,NMESH)
      END DO
      DO I = 1,NV
         YOUT(I) = 0.5D0*(YM(I)+YN(I)+H*YOUT(I))
C
      END DO
C
      END 
