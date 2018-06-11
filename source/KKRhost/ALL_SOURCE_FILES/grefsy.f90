    Subroutine grefsy13(gtmat, gmat, dgtde, lly_g0tr, ipvt, ndim, lly, lmgf0d, &
      ngd)
      Use mod_datatypes, Only: dp
! **********************************************************************
!   Solve the Dyson equation to get reference Green function
! Calculate also (1-gt)^-1 * d(1-gt)/dE and the trace LLY_G0TR for
! Lloyds formula (ported from KKRnano by Phivos Mavropoulos 11.10.2013)
! **********************************************************************
      Implicit None
!     .. PARAMETERS ..
      Complex (Kind=dp) :: czero, cone
      Parameter (czero=(0.E0_dp,0.E0_dp), cone=(1.E0_dp,0.E0_dp))
!     ..
!     .. SCALAR ARGUMENTS ..
      Integer :: ndim, ngd, lmgf0d
      Integer :: lly ! LLY .ne. 0: calculate GF derivative; input
      Complex (Kind=dp) :: lly_g0tr ! LLY Trace of (1-gt)^-1 * d(1-gt)/dE ; output
!     ..
!     .. ARRAY ARGUMENTS ..
      Complex (Kind=dp) :: gmat(ngd, lmgf0d), gtmat(ngd, ngd)
!     On input, DGTDE contains the source term -dg/dE * t - g * dt/dE
!     (where g is the free GF, t the ref.-sys. t-matrix); 
!     on output it contains (1-gt)^-1 * d(1-gt)/dE
!     (PhD thesis Alex Thiess, eq. 5.28)
      Complex (Kind=dp) :: dgtde(ngd, lmgf0d) ! LLY input/output
!     ..
!     .. LOCAL SCALARS ..
      Integer :: ii, info
!     ..
!     .. LOCAL ARRAYS ..
      Integer :: ipvt(ngd)
!     ..
!     .. EXTERNAL SUBROUTINES ..
      External :: zgetrf, zgetrs
!     ..

! GTMAT =  -g*t
! GMAT = g
! DGTDE = -dg/dE * t - g * dt/dE (on input)


      Do ii = 1, ndim
        gtmat(ii, ii) = cone + gtmat(ii, ii) ! GTMAT= 1 - g * t
      End Do
!
!---> SOLVE THE SYSTEM OF LINEAR EQUATIONS
!
      Call zgetrf(ndim, ndim, gtmat, ngd, ipvt, info)
      Call zgetrs('N', ndim, lmgf0d, gtmat, ngd, ipvt, gmat, ngd, info)

      lly_g0tr = czero

      If (lly==0) Return

! LLY Calculate GF derivative and trace for Lloyd
      Call zgetrs('N', ndim, lmgf0d, gtmat, ngd, ipvt, dgtde, ngd, info) ! LLY
! DGTDE now contains (1-gt)^-1 * d(1-gt)/dE
! (Thiess PhD eq. 5.28)

      Do ii = 1, lmgf0d
        lly_g0tr = lly_g0tr - dgtde(ii, ii) ! LLY
      End Do

! LLY_G0TR contains  -Trace[ (1-gt)^-1 * d(1-gt)/dE ]


    End Subroutine
! **********************************************************************

! Obsolete, replaced by grefsy13 returning also the derivative on demand.

    Subroutine grefsy(gtmat, gmat, ndim, lmgf0d, ngd)
      Use mod_datatypes, Only: dp
! **********************************************************************
! * Solve the Dyson equation to get reference Green function           *
! **********************************************************************
      Implicit None
!     .. PARAMETERS ..
      Complex (Kind=dp) :: cone
      Parameter (cone=(1.E0_dp,0.E0_dp))
!     ..
!     .. SCALAR ARGUMENTS ..
      Integer :: ndim, ngd, lmgf0d
!     ..
!     .. ARRAY ARGUMENTS ..
      Complex (Kind=dp) :: gmat(ngd, lmgf0d), gtmat(ngd, ngd)
!     ..
!     .. LOCAL SCALARS ..
      Integer :: i, info
!     ..
!     .. LOCAL ARRAYS ..
      Integer :: ipvt(ngd)
!     ..
!     .. EXTERNAL SUBROUTINES ..
      External :: zgetrf, zgetrs
!     ..

      Do i = 1, ndim
        gtmat(i, i) = cone + gtmat(i, i) ! GTMAT= 1 - G * T
      End Do

!---> SOLVE THE SYSTEM OF LINEAR EQUATIONS

      Call zgetrf(ndim, ndim, gtmat, ngd, ipvt, info)
      Call zgetrs('N', ndim, lmgf0d, gtmat, ngd, ipvt, gmat, ngd, info)
      Return

    End Subroutine
