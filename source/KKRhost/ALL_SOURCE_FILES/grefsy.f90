subroutine grefsy13(gtmat, gmat, dgtde, lly_g0tr, ipvt, ndim, lly, lmgf0d, &
  ngd)
! **********************************************************************
!   Solve the Dyson equation to get reference Green function
! Calculate also (1-gt)^-1 * d(1-gt)/dE and the trace LLY_G0TR for
! Lloyds formula (ported from KKRnano by Phivos Mavropoulos 11.10.2013)
! **********************************************************************
  implicit none
!     .. PARAMETERS ..
  double complex :: czero, cone
  parameter (czero=(0.d0,0.d0), cone=(1.d0,0.d0))
!     ..
!     .. SCALAR ARGUMENTS ..
  integer :: ndim, ngd, lmgf0d
  integer :: lly ! LLY .ne. 0: calculate GF derivative; input
  double complex :: lly_g0tr ! LLY Trace of (1-gt)^-1 * d(1-gt)/dE ; output
!     ..
!     .. ARRAY ARGUMENTS ..
  double complex :: gmat(ngd, lmgf0d), gtmat(ngd, ngd)
!     On input, DGTDE contains the source term -dg/dE * t - g * dt/dE
!     (where g is the free GF, t the ref.-sys. t-matrix); 
!     on output it contains (1-gt)^-1 * d(1-gt)/dE
!     (PhD thesis Alex Thiess, eq. 5.28)
  double complex :: dgtde(ngd, lmgf0d) ! LLY input/output
!     ..
!     .. LOCAL SCALARS ..
  integer :: ii, info
!     ..
!     .. LOCAL ARRAYS ..
  integer :: ipvt(ngd)
!     ..
!     .. EXTERNAL SUBROUTINES ..
  external :: zgetrf, zgetrs
!     ..

! GTMAT =  -g*t
! GMAT = g
! DGTDE = -dg/dE * t - g * dt/dE (on input)


  do ii = 1, ndim
    gtmat(ii, ii) = cone + gtmat(ii, ii) ! GTMAT= 1 - g * t
  end do
!
!---> SOLVE THE SYSTEM OF LINEAR EQUATIONS
!
  call zgetrf(ndim, ndim, gtmat, ngd, ipvt, info)
  call zgetrs('N', ndim, lmgf0d, gtmat, ngd, ipvt, gmat, ngd, info)

  lly_g0tr = czero

  if (lly==0) return

! LLY Calculate GF derivative and trace for Lloyd
  call zgetrs('N', ndim, lmgf0d, gtmat, ngd, ipvt, dgtde, ngd, info) ! LLY
! DGTDE now contains (1-gt)^-1 * d(1-gt)/dE
! (Thiess PhD eq. 5.28)

  do ii = 1, lmgf0d
    lly_g0tr = lly_g0tr - dgtde(ii, ii) ! LLY
  end do

! LLY_G0TR contains  -Trace[ (1-gt)^-1 * d(1-gt)/dE ]


end subroutine
! **********************************************************************

! Obsolete, replaced by grefsy13 returning also the derivative on demand.

subroutine grefsy(gtmat, gmat, ndim, lmgf0d, ngd)
! **********************************************************************
! * Solve the Dyson equation to get reference Green function           *
! **********************************************************************
  implicit none
!     .. PARAMETERS ..
  double complex :: cone
  parameter (cone=(1.d0,0.d0))
!     ..
!     .. SCALAR ARGUMENTS ..
  integer :: ndim, ngd, lmgf0d
!     ..
!     .. ARRAY ARGUMENTS ..
  double complex :: gmat(ngd, lmgf0d), gtmat(ngd, ngd)
!     ..
!     .. LOCAL SCALARS ..
  integer :: i, info
!     ..
!     .. LOCAL ARRAYS ..
  integer :: ipvt(ngd)
!     ..
!     .. EXTERNAL SUBROUTINES ..
  external :: zgetrf, zgetrs
!     ..

  do i = 1, ndim
    gtmat(i, i) = cone + gtmat(i, i) ! GTMAT= 1 - G * T
  end do

!---> SOLVE THE SYSTEM OF LINEAR EQUATIONS

  call zgetrf(ndim, ndim, gtmat, ngd, ipvt, info)
  call zgetrs('N', ndim, lmgf0d, gtmat, ngd, ipvt, gmat, ngd, info)
  return

end subroutine
