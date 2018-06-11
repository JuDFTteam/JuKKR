    Subroutine dirbsrad(xbs, y, dydx, drdi, b, v, r, nmesh)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   *   supply the derivatives for the coupled set of                  *
!   *   radial dirac equation in case of a spin-dependent potential    *
!   *                                                                  *
!   ********************************************************************

      Implicit None

      Include 'sprkkr_rmesh.dim'


! PARAMETER definitions
      Integer :: nlag
      Parameter (nlag=3)

! COMMON variables
      Real (Kind=dp) :: cgd(2), cgmd(2), cgo, kap(2)
      Real (Kind=dp) :: csqr
      Complex (Kind=dp) :: ebs
      Integer :: nradbs, nsolbs
      Common /commbs/ebs, csqr, cgd, cgmd, cgo, kap, nsolbs, nradbs

! Dummy arguments
      Integer :: nmesh
      Real (Kind=dp) :: xbs
      Real (Kind=dp) :: b(nrmax), drdi(nrmax), r(nrmax), v(nrmax)
      Complex (Kind=dp) :: dydx(ncfmax), y(ncfmax)

! Local variables
      Real (Kind=dp) :: bbs, bpp, bqq, drdibs, rbs, vbs
      Complex (Kind=dp) :: emvpp, emvqq
      Integer :: i, j, k, m

      Call dirbslag(xbs, vbs, bbs, rbs, drdibs, v, b, r, drdi, nradbs, nlag, &
        nmesh)

      emvqq = (ebs-vbs+csqr)*drdibs/csqr
      emvpp = -(ebs-vbs)*drdibs
      bqq = bbs*drdibs/csqr
      bpp = bbs*drdibs


      m = 0
      Do j = 1, nsolbs
        Do i = 1, nsolbs
          k = 3 - 4*(i-1)
          dydx(m+1) = -kap(i)*y(m+1)/rbs*drdibs + (emvqq+bqq*cgmd(i))*y(m+2)
          dydx(m+2) = kap(i)*y(m+2)/rbs*drdibs + (emvpp+bpp*cgd(i))*y(m+1) + &
            bpp*cgo*y(m+k)
          m = m + 2
        End Do
      End Do
    End Subroutine
