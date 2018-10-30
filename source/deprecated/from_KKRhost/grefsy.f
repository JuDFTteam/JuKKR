      SUBROUTINE GREFSY13(GTMAT,GMAT,DGTDE,LLY_G0TR,IPVT,
     &                    NDIM,LLY,LMGF0D,NGD)
! **********************************************************************
!   Solve the Dyson equation to get reference Green function            
! Calculate also (1-gt)^-1 * d(1-gt)/dE and the trace LLY_G0TR for      
! Lloyds formula (ported from KKRnano by Phivos Mavropoulos 11.10.2013) 
! **********************************************************************

      IMPLICIT NONE
!     .. PARAMETERS ..
      DOUBLE COMPLEX CZERO,CONE
      PARAMETER (CZERO=(0.D0,0.D0),CONE= (1.D0,0.D0))
!     ..
!     .. SCALAR ARGUMENTS ..
      INTEGER NDIM,NGD,LMGF0D
      INTEGER LLY ! LLY .ne. 0: calculate GF derivative; input
      DOUBLE COMPLEX LLY_G0TR ! LLY Trace of (1-gt)^-1 * d(1-gt)/dE ; output
!     ..
!     .. ARRAY ARGUMENTS ..
      DOUBLE COMPLEX GMAT(NGD,LMGF0D),GTMAT(NGD,NGD)
!     On input, DGTDE contains the source term -dg/dE * t - g * dt/dE
!     (where g is the free GF, t the ref.-sys. t-matrix); 
!     on output it contains (1-gt)^-1 * d(1-gt)/dE
!     (PhD thesis Alex Thiess, eq. 5.28)
      DOUBLE COMPLEX DGTDE(NGD,LMGF0D) ! LLY input/output
!     ..
!     .. LOCAL SCALARS ..
      INTEGER II,INFO
!     ..
!     .. LOCAL ARRAYS ..
      INTEGER IPVT(NGD)
!     ..
!     .. EXTERNAL SUBROUTINES ..
      EXTERNAL ZGETRF,ZGETRS
!     ..

! GTMAT =  -g*t
! GMAT = g
! DGTDE = -dg/dE * t - g * dt/dE (on input)


      DO 10 II = 1,NDIM
        GTMAT(II,II) = CONE + GTMAT(II,II) ! GTMAT= 1 - g * t
   10 CONTINUE
!
!---> SOLVE THE SYSTEM OF LINEAR EQUATIONS
!
      CALL ZGETRF(NDIM,NDIM,GTMAT,NGD,IPVT,INFO)
      CALL ZGETRS('N',NDIM,LMGF0D,GTMAT,NGD,IPVT,GMAT,NGD,INFO)

      LLY_G0TR = CZERO

      IF (LLY.EQ.0) RETURN

      ! LLY Calculate GF derivative and trace for Lloyd
      CALL ZGETRS('N',NDIM,LMGF0D,GTMAT,NGD,IPVT,DGTDE,NGD,INFO) ! LLY
      ! DGTDE now contains (1-gt)^-1 * d(1-gt)/dE
      ! (Thiess PhD eq. 5.28)

      DO II = 1,LMGF0D
         LLY_G0TR = LLY_G0TR - DGTDE(II,II)                      ! LLY
      ENDDO

      ! LLY_G0TR contains  -Trace[ (1-gt)^-1 * d(1-gt)/dE ]


      END
! **********************************************************************

      ! Obsolete, replaced by grefsy13 returning also the derivative on demand.
      SUBROUTINE GREFSY(GTMAT,GMAT,NDIM,LMGF0D,NGD)
C **********************************************************************
C * Solve the Dyson equation to get reference Green function           *
C **********************************************************************
C
      IMPLICIT NONE
!     .. PARAMETERS ..
      DOUBLE COMPLEX CONE
      PARAMETER (CONE= (1.D0,0.D0))
!     ..
!     .. SCALAR ARGUMENTS ..
      INTEGER NDIM,NGD,LMGF0D
!     ..
!     .. ARRAY ARGUMENTS ..
      DOUBLE COMPLEX GMAT(NGD,LMGF0D),GTMAT(NGD,NGD)
!     ..
!     .. LOCAL SCALARS ..
      INTEGER I,INFO
!     ..
!     .. LOCAL ARRAYS ..
      INTEGER IPVT(NGD)
!     ..
!     .. EXTERNAL SUBROUTINES ..
      EXTERNAL ZGETRF,ZGETRS
!     ..
 
      DO 10 I = 1,NDIM
        GTMAT(I,I) = CONE + GTMAT(I,I) ! GTMAT= 1 - G * T
   10 CONTINUE
 
!---> SOLVE THE SYSTEM OF LINEAR EQUATIONS

      CALL ZGETRF(NDIM,NDIM,GTMAT,NGD,IPVT,INFO)
      CALL ZGETRS('N',NDIM,LMGF0D,GTMAT,NGD,IPVT,GMAT,NGD,INFO)
      RETURN

      END
