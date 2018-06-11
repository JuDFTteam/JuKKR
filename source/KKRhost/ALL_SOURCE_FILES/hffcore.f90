    Subroutine hffcore(rnuc, jtop, kap1, kap2, nsol, mj, gc, fc, nrc, shf, s, &
      nmemax, nkmmax, r, drdi, sdia, smdia, soff, smoff, qdia, qoff, qmdia, &
      qmoff, nucleus, jlim)
      Use mod_datatypes, Only: dp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Calculates matrix elements of several hyperfine interaction
!     connected quantities in the core.
!     All the related arrays have a counting index as
!     the last index of the array indicates the corresponding physical
!     property.
!     Index-list
!     1      electron-Spin-electron-Spin Hyperfine field
!     2      nuclear-spin-electron-orbit hyperfine field
!     3      electron-spin-nulceus-spin-contact hyperfine field
!     4      expectation value of (1/r)^3
!     5      Total Hyperfine Field (see Rose (1961))
!     called by core
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Implicit None


! PARAMETER definitions
      Real (Kind=dp) :: mb, a0, f1, f2

!BOHR-MAGNETON       IN ERG/GAUSS
!CONVERSION FACTOR FOR HYPERFINE FIELDS FROM A.U. TO GAUSS
!ELECTRON CHARGE     IN ESU

      Parameter (mb=9.274078E-21_dp, a0=0.52917706E-08_dp, f1=1.0E0_dp, &
        f2=2.0E0_dp*mb/(a0*a0*a0))

! Dummy arguments
      Integer :: jlim, jtop, kap1, kap2, nkmmax, nmemax, nrc, nsol, nucleus, s
      Real (Kind=dp) :: mj, rnuc
      Real (Kind=dp) :: drdi(nrc), fc(2, 2, nrc), gc(2, 2, nrc), qdia(nkmmax), &
        qmdia(nkmmax), qmoff(nkmmax), qoff(nkmmax), r(nrc), sdia(nkmmax), &
        shf(2, 2, nmemax), smdia(nkmmax), smoff(nkmmax), soff(nkmmax)

! Local variables
      Real (Kind=dp) :: ame(2, 2), cff(2, 2), cfg(2, 2), cgf(2, 2), cgg(2, 2), &
        cqf(2, 2), cqg(2, 2), csf(2, 2), csg(2, 2), dovr(nrc), drovrn(nrc), &
        drovrn1(nrc), f(nrc, 2), ff(2, 2), ff1(2, 2), ff2(2, 2), fg(2, 2), &
        fg1(2, 2), fg2(2, 2), g(nrc, 2), gf(2, 2), gf1(2, 2), gf2(2, 2), &
        gg(2, 2), gg1(2, 2), gg2(2, 2)
      Real (Kind=dp) :: dble, dsqrt
      Integer :: i, ikm1, ikm2, j, k, k1, k2, n
      Integer :: ikapmue
      Integer :: nint

      If (kap2==0) kap2 = kap1

      Call rinit(4, gg)
      Call rinit(4, ff)
      Call rinit(4, gg1)
      Call rinit(4, ff1)
      Call rinit(4, gg2)
      Call rinit(4, ff2)
      Call rinit(4, gf)
      Call rinit(4, fg)
      Call rinit(4, gf1)
      Call rinit(4, fg1)
      Call rinit(4, gf2)
      Call rinit(4, fg2)
      Call rinit(2*nrc, g)
      Call rinit(2*nrc, f)

      Do k = 1, 2
        Do n = 1, jtop
          g(n, k) = gc(k, s, n)
          f(n, k) = fc(k, s, n)
        End Do
      End Do
!     prepare meshes for finite nucleus calculation
      Do i = 1, jtop
        dovr(i) = drdi(i)/r(i)
        If (nucleus/=0) Then
          drovrn1(i) = (r(i)/rnuc)**3*drdi(i)
          drovrn(i) = drovrn1(i)/r(i)
        End If
      End Do
      ikm1 = ikapmue(kap1, nint(mj-0.5E0_dp))
      ikm2 = ikapmue(kap2, nint(mj-0.5E0_dp))
!     angular hyperfine matrix elements   see e.g.  E.M.Rose
!     the factor  i  has been omitted
      Call rinit(4, ame)
      ame(1, 1) = 4.0E0_dp*kap1*mj/(4.0E0_dp*kap1*kap1-1.0E0_dp)
      If (nsol==2) Then
        ame(2, 2) = 4.0E0_dp*kap2*mj/(4.0E0_dp*kap2*kap2-1.0E0_dp)
        ame(1, 2) = dsqrt(0.25E0_dp-(mj/dble(kap1-kap2))**2)
        ame(2, 1) = ame(1, 2)
      End If
!     coefficients for the spin-dipolar matrix elements
      Call rinit(4, csf)
      Call rinit(4, csg)
      csg(1, 1) = sdia(ikm1)
      csf(1, 1) = smdia(ikm1)
      If (nsol==2) Then
        csg(2, 2) = sdia(ikm2)
        csg(1, 2) = soff(ikm1)
        csg(2, 1) = csg(1, 2)
        csf(2, 2) = smdia(ikm2)
        csf(1, 2) = smoff(ikm1)
        csf(2, 1) = smoff(ikm1)
      End If
!     COEFFICIENTS FOR THE QUADRUPOLAR MATRIX ELEMENTS
      cqg(1, 1) = qdia(ikm1)
      cqg(2, 2) = qdia(ikm2)
      cqg(1, 2) = qoff(ikm1)
      cqg(2, 1) = cqg(1, 2)
      Call rinit(4, cqf)
      cqf(1, 1) = qmdia(ikm1)
      cqf(2, 2) = qmdia(ikm2)
      cqf(1, 2) = qmoff(ikm1)
      cqf(2, 1) = cqf(1, 2)
!     coefficients to calculate the spin-spin field
      Call rinit(4, cgg)
      Call rinit(4, cgf)
      cgg(1, 1) = -mj/(kap1+0.5E0_dp)
      cgf(1, 1) = -mj/(-kap1+0.5E0_dp)
      If (nsol==2) Then
        cgg(1, 2) = -dsqrt(1.0E0_dp-(mj/(kap1+0.5E0_dp))**2)
        cgg(2, 1) = cgg(1, 2)
        cgg(2, 2) = -mj/(kap2+0.5E0_dp)
        cgf(2, 2) = -mj/(-kap2+0.5E0_dp)
!     CGF(1,2) = -DSQRT( 1.0D0 - (MJ/(- KAP1+0.5D0))**2 )
!     CGF(2,1) = CGF(1,2)
      End If
!     coefficients to calculate the orbital field
      Call rinit(4, cfg)
      Call rinit(4, cff)
      cfg(1, 1) = mj*(kap1+1.0E0_dp)/(kap1+0.5E0_dp)
      cff(1, 1) = mj*(-kap1+1.0E0_dp)/(-kap1+0.5E0_dp)
      If (nsol==2) Then
        cfg(2, 2) = mj*(kap2+1.0E0_dp)/(kap2+0.5E0_dp)
        cfg(1, 2) = 0.5E0_dp*dsqrt(1.0E0_dp-(mj/(kap1+0.5E0_dp))**2)
        cfg(2, 1) = cfg(1, 2)
        cff(2, 2) = mj*(-kap2+1.0E0_dp)/(-kap2+0.5E0_dp)
!     CFF(1,2) = 0.5D0 * DSQRT( 1.0D0 - (MJ/(- KAP1 + 0.5D0))**2 )
!     CFF(2,1) = CFF(1,2)
      End If
! Calculates integrals from 0.0 to jtop
      Call hffint(gg, g, g, dovr, r, 0.0E0_dp, nsol, jtop, nrc)
      Call hffint(ff, f, f, dovr, r, 0.0E0_dp, nsol, jtop, nrc)
      Call hffint(gf, g, f, drdi, r, 0.0E0_dp, nsol, jtop, nrc)
      Call hffint(fg, f, g, drdi, r, 0.0E0_dp, nsol, jtop, nrc)
      Call rsumupint(shf(1,1,5), f1, gg, cqg, f1, ff, cqf, nsol)
      If (nucleus/=0) Then
!     calculates integrals inside nucleus at RNUC in order to get
!     contribution outside the nucleus
        Call hffint(gg1, g, g, dovr, r, rnuc, nsol, jlim, nrc)
        Call hffint(ff1, f, f, dovr, r, rnuc, nsol, jlim, nrc)
        Call hffint(gf1, g, f, drdi, r, rnuc, nsol, jlim, nrc)
        Call hffint(fg1, f, g, drdi, r, rnuc, nsol, jlim, nrc)
!     calculates contribution from RNUC to jtop
        Do i = 1, nsol
          Do j = 1, nsol
            gg(i, j) = gg(i, j) - gg1(i, j)
            ff(i, j) = ff(i, j) - ff1(i, j)
            gf(i, j) = gf(i, j) - gf1(i, j)
            fg(i, j) = fg(i, j) - fg1(i, j)
          End Do
        End Do
      End If !end of nucleus.eq.0
!     calculates B_sp which is zero inside the nucleus
      Call rsumupint(shf(1,1,1), f1, gg, csg, -f1, ff, csf, nsol)
!     calculates hyperfine integrals from 0.0 to RNUC which are added
!     external integrals
      If (nucleus/=0) Then
        Call hffint(gg2, g, g, drovrn, r, rnuc, nsol, jlim, nrc)
        Call hffint(ff2, f, f, drovrn, r, rnuc, nsol, jlim, nrc)
        Call hffint(gf2, g, f, drovrn1, r, rnuc, nsol, jlim, nrc)
        Call hffint(fg2, f, g, drovrn1, r, rnuc, nsol, jlim, nrc)
        Do i = 1, nsol
          Do j = 1, nsol
            gg(i, j) = gg(i, j) + gg2(i, j)
            ff(i, j) = ff(i, j) + ff2(i, j)
            gf(i, j) = gf(i, j) + gf2(i, j)
            fg(i, j) = fg(i, j) + fg2(i, j)
          End Do
        End Do
      End If
!     calculates B_nseo and B_tot
      Call rsumupint(shf(1,1,2), f2, gg, cfg, -f2, ff, cff, nsol)
!      CALL RSUMUPINT(SHF(1,1,5),CAUTOG,GF,AME,CAUTOG,FG,AME,NSOL)
!     modifications for B_ssc which is zero outside the nucleus
      If (nucleus/=0) Then
        Do i = 1, nsol
          Do j = 1, nsol
            gg(i, j) = gg2(i, j)
            ff(i, j) = ff2(i, j)
          End Do
        End Do
      End If
!     calculates B_ssc
      Call rsumupint(shf(1,1,3), f2, gg, cgg, -f2, ff, cgf, nsol)
!     for testing purposes write in 4 the sum of 1,2,3
      Do k1 = 1, nsol
        Do k2 = 1, nsol
          shf(k1, k2, 4) = shf(k1, k2, 1) + shf(k1, k2, 2) + shf(k1, k2, 3)
        End Do
      End Do

    End Subroutine

    Subroutine rsumupint(sum, vg, g, wg, vf, f, wf, n)
      Use mod_datatypes, Only: dp
      Implicit None

! Dummy arguments
      Integer :: n
      Real (Kind=dp) :: vf, vg
      Real (Kind=dp) :: f(n, n), g(n, n), sum(n, n), wf(n, n), wg(n, n)

! Local variables
      Integer :: i, j

      Do j = 1, n
        Do i = 1, n
          sum(i, j) = vg*g(i, j)*wg(i, j) + vf*f(i, j)*wf(i, j)
        End Do
      End Do
    End Subroutine

    Subroutine hffint(gg, ga, gb, dr, r, rnuc, nsol, jtop, nrc)
      Use mod_datatypes, Only: dp
!     Calculates Hyperfine integrals, extrapolates to zero and
!     intrapolates to exact nuclear radius RNUC
      Implicit None

! Dummy arguments
      Integer :: jtop, nrc, nsol
      Real (Kind=dp) :: rnuc
      Real (Kind=dp) :: dr(nrc), ga(nrc, 2), gb(nrc, 2), gg(2, 2), r(nrc)

! Local variables
      Integer :: i, k1, k2
      Real (Kind=dp) :: x(5), y(5), yi(nrc), zi(nrc)
      Real (Kind=dp) :: ylag

      Do k1 = 1, nsol
        Do k2 = 1, nsol
          Do i = 1, jtop
            yi(i) = ga(i, k1)*gb(i, k2)*dr(i)
          End Do
          Call rint4pts(yi, jtop, zi)
!     Intrapolation
          If (rnuc/=0.0E0_dp) Then
            Do i = 1, 5
              x(i) = r(jtop-5+i)
              y(i) = zi(jtop-5+i)
            End Do
            zi(jtop) = ylag(rnuc, x, y, 0, 4, 5)
          End If
!     Extrapolation
          x(1) = 1.0E0_dp
          x(2) = 6.0E0_dp
          x(3) = 11.0E0_dp
          y(1) = zi(jtop) - zi(1)
          y(2) = zi(jtop) - zi(5)
          y(3) = zi(jtop) - zi(9)
          gg(k1, k2) = ylag(0.0E0_dp, x, y, 0, 2, 3)
        End Do
      End Do
    End Subroutine
