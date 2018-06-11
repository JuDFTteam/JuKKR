! **********************************************************************
    Subroutine gfree13(rdiff, e0, gmll, dgmll, cleb, icleb, loflm, iend)
      Use mod_datatypes, Only: dp
! **********************************************************************
!     .. Parameters ..
      Implicit None
      Include 'inc.p'
      Integer :: lmax
      Parameter (lmax=lmaxd)
      Integer :: lmgf0d
      Parameter (lmgf0d=(lmax+1)**2)
      Integer :: lmax2p, lmx2sq
      Parameter (lmax2p=lmax*2+1, lmx2sq=lmax2p**2)
      Complex (Kind=dp) :: czero, ci
      Parameter (czero=(0.0E0_dp,0.0E0_dp), ci=(0.0E0_dp,1.0E0_dp))
! LLY
!     ..
      Complex (Kind=dp) :: e0
      Integer :: iend
!     .. Local Scalars ..
!     ..
      Complex (Kind=dp) :: gmll(lmgf0d, lmgf0d)
      Complex (Kind=dp) :: dgmll(lmgf0d, lmgf0d) !     .. Local Arrays ..
      Real (Kind=dp) :: cleb(ncleb), rdiff(*)
      Integer :: icleb(ncleb, 4), loflm(*)
! LLY
!     ..
      Real (Kind=dp) :: fpi, pi, rabs, rfpi, x, y, z
      Integer :: ifac, j, lm1, lm2, lm3
!     .. External Subroutines ..
!     ..
      Complex (Kind=dp) :: bl(lmax2p), hl(lmax2p), hyl(lmx2sq), nl(lmax2p)
      Complex (Kind=dp) :: dhl(lmax2p), dhyl(lmx2sq) !     .. Intrinsic Functions ..
      Real (Kind=dp) :: yl(lmx2sq)
      Integer :: lf(lmx2sq)
!     ..
!-----------------------------------------------------------------------
      External :: beshan, ymy
!---- CALCULATION OF FREE ELECTRON GREEN'S FUNCTION :  G(M)LL'(E0)
!-----------------------------------------------------------------------
      Intrinsic :: atan, sqrt
!     Also:
      pi = 4.E0_dp*atan(1.E0_dp)
      fpi = 4.E0_dp*pi
      rfpi = sqrt(fpi)
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

! Derivatives of Hankel functions ! LLY

! LLY
      Do lm1 = 1, lmx2sq
        lf(lm1) = loflm(lm1) + 1
      End Do
      x = rdiff(1)
      y = rdiff(2)
      z = rdiff(3)
      Call ymy(x, y, z, rabs, yl, lmax*2)
      Call beshan(hl, bl, nl, sqrt(e0)*rabs, lmax*2)


      dhl(1) = 0.5E0_dp*ci*rabs*hl(1)
      Do lm1 = 2, lmax2p
        dhl(lm1) = 0.5E0_dp*(rabs*hl(lm1-1)-(lm1-1)*hl(lm1)/sqrt(e0))
      End Do
! LLY
      Do lm1 = 1, lmx2sq
        hyl(lm1) = -fpi*ci*sqrt(e0)*yl(lm1)*hl(lf(lm1))
        dhyl(lm1) = -fpi*ci*yl(lm1)*dhl(lf(lm1)) ! LLY
      End Do

! LLY
      Do lm1 = 1, lmgf0d
        gmll(lm1, lm1) = hyl(1)/rfpi
        dgmll(lm1, lm1) = dhyl(1)/rfpi ! LLY
        Do lm2 = 1, lm1 - 1
          gmll(lm1, lm2) = czero
          dgmll(lm1, lm2) = czero
        End Do
      End Do
! **********************************************************************
      Do j = 1, iend
        lm1 = icleb(j, 1)
        lm2 = icleb(j, 2)
        lm3 = icleb(j, 3)
        gmll(lm1, lm2) = gmll(lm1, lm2) + cleb(j)*hyl(lm3)
        dgmll(lm1, lm2) = dgmll(lm1, lm2) + cleb(j)*dhyl(lm3) ! **********************************************************************
      End Do
      Do lm1 = 1, lmgf0d
        Do lm2 = 1, lm1 - 1
          ifac = (-1)**(loflm(lm1)+loflm(lm2))
          gmll(lm2, lm1) = ifac*gmll(lm1, lm2)
          dgmll(lm2, lm1) = ifac*dgmll(lm1, lm2) !     .. Parameters ..
        End Do
      End Do
      Return
! set to 1 if NEWSOSOL under RUNOPT, otherwise 0
    End Subroutine
