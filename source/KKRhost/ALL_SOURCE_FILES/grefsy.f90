SUBROUTINE grefsy13(gtmat,gmat,dgtde,lly_g0tr,ipvt,  &
        ndim,lly,lmgf0d,ngd)
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


DO  ii = 1,ndim
  gtmat(ii,ii) = cone + gtmat(ii,ii) ! GTMAT= 1 - g * t
END DO
!
!---> SOLVE THE SYSTEM OF LINEAR EQUATIONS
!
CALL zgetrf(ndim,ndim,gtmat,ngd,ipvt,info)
CALL zgetrs('N',ndim,lmgf0d,gtmat,ngd,ipvt,gmat,ngd,info)

lly_g0tr = czero

IF (lly == 0) RETURN

! LLY Calculate GF derivative and trace for Lloyd
CALL zgetrs('N',ndim,lmgf0d,gtmat,ngd,ipvt,dgtde,ngd,info) ! LLY
! DGTDE now contains (1-gt)^-1 * d(1-gt)/dE
! (Thiess PhD eq. 5.28)

DO ii = 1,lmgf0d
  lly_g0tr = lly_g0tr - dgtde(ii,ii)                      ! LLY
END DO

! LLY_G0TR contains  -Trace[ (1-gt)^-1 * d(1-gt)/dE ]


END SUBROUTINE grefsy13
! **********************************************************************

! Obsolete, replaced by grefsy13 returning also the derivative on demand.

SUBROUTINE grefsy(gtmat,gmat,ndim,lmgf0d,ngd)
! **********************************************************************
! * Solve the Dyson equation to get reference Green function           *
! **********************************************************************
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

DO  i = 1,ndim
  gtmat(i,i) = cone + gtmat(i,i) ! GTMAT= 1 - G * T
END DO

!---> SOLVE THE SYSTEM OF LINEAR EQUATIONS

CALL zgetrf(ndim,ndim,gtmat,ngd,ipvt,info)
CALL zgetrs('N',ndim,lmgf0d,gtmat,ngd,ipvt,gmat,ngd,info)
RETURN

END SUBROUTINE grefsy
