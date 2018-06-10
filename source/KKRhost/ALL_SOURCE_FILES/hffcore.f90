subroutine hffcore(rnuc, jtop, kap1, kap2, nsol, mj, gc, fc, nrc, shf, s, &
  nmemax, nkmmax, r, drdi, sdia, smdia, soff, smoff, qdia, qoff, qmdia, qmoff, &
  nucleus, jlim)
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
  implicit none


! PARAMETER definitions
  double precision :: mb, a0, f1, f2

!BOHR-MAGNETON       IN ERG/GAUSS
!CONVERSION FACTOR FOR HYPERFINE FIELDS FROM A.U. TO GAUSS
!ELECTRON CHARGE     IN ESU

  parameter (mb=9.274078d-21, a0=0.52917706d-08, f1=1.0d0, &
    f2=2.0d0*mb/(a0*a0*a0))

! Dummy arguments
  integer :: jlim, jtop, kap1, kap2, nkmmax, nmemax, nrc, nsol, nucleus, s
  double precision :: mj, rnuc
  double precision :: drdi(nrc), fc(2, 2, nrc), gc(2, 2, nrc), qdia(nkmmax), &
    qmdia(nkmmax), qmoff(nkmmax), qoff(nkmmax), r(nrc), sdia(nkmmax), &
    shf(2, 2, nmemax), smdia(nkmmax), smoff(nkmmax), soff(nkmmax)

! Local variables
  double precision :: ame(2, 2), cff(2, 2), cfg(2, 2), cgf(2, 2), cgg(2, 2), &
    cqf(2, 2), cqg(2, 2), csf(2, 2), csg(2, 2), dovr(nrc), drovrn(nrc), &
    drovrn1(nrc), f(nrc, 2), ff(2, 2), ff1(2, 2), ff2(2, 2), fg(2, 2), &
    fg1(2, 2), fg2(2, 2), g(nrc, 2), gf(2, 2), gf1(2, 2), gf2(2, 2), gg(2, 2), &
    gg1(2, 2), gg2(2, 2)
  double precision :: dble, dsqrt
  integer :: i, ikm1, ikm2, j, k, k1, k2, n
  integer :: ikapmue
  integer :: nint

  if (kap2==0) kap2 = kap1

  call rinit(4, gg)
  call rinit(4, ff)
  call rinit(4, gg1)
  call rinit(4, ff1)
  call rinit(4, gg2)
  call rinit(4, ff2)
  call rinit(4, gf)
  call rinit(4, fg)
  call rinit(4, gf1)
  call rinit(4, fg1)
  call rinit(4, gf2)
  call rinit(4, fg2)
  call rinit(2*nrc, g)
  call rinit(2*nrc, f)

  do k = 1, 2
    do n = 1, jtop
      g(n, k) = gc(k, s, n)
      f(n, k) = fc(k, s, n)
    end do
  end do
!     prepare meshes for finite nucleus calculation
  do i = 1, jtop
    dovr(i) = drdi(i)/r(i)
    if (nucleus/=0) then
      drovrn1(i) = (r(i)/rnuc)**3*drdi(i)
      drovrn(i) = drovrn1(i)/r(i)
    end if
  end do
  ikm1 = ikapmue(kap1, nint(mj-0.5d0))
  ikm2 = ikapmue(kap2, nint(mj-0.5d0))
!     angular hyperfine matrix elements   see e.g.  E.M.Rose
!     the factor  i  has been omitted
  call rinit(4, ame)
  ame(1, 1) = 4.0d0*kap1*mj/(4.0d0*kap1*kap1-1.0d0)
  if (nsol==2) then
    ame(2, 2) = 4.0d0*kap2*mj/(4.0d0*kap2*kap2-1.0d0)
    ame(1, 2) = dsqrt(0.25d0-(mj/dble(kap1-kap2))**2)
    ame(2, 1) = ame(1, 2)
  end if
!     coefficients for the spin-dipolar matrix elements
  call rinit(4, csf)
  call rinit(4, csg)
  csg(1, 1) = sdia(ikm1)
  csf(1, 1) = smdia(ikm1)
  if (nsol==2) then
    csg(2, 2) = sdia(ikm2)
    csg(1, 2) = soff(ikm1)
    csg(2, 1) = csg(1, 2)
    csf(2, 2) = smdia(ikm2)
    csf(1, 2) = smoff(ikm1)
    csf(2, 1) = smoff(ikm1)
  end if
!     COEFFICIENTS FOR THE QUADRUPOLAR MATRIX ELEMENTS
  cqg(1, 1) = qdia(ikm1)
  cqg(2, 2) = qdia(ikm2)
  cqg(1, 2) = qoff(ikm1)
  cqg(2, 1) = cqg(1, 2)
  call rinit(4, cqf)
  cqf(1, 1) = qmdia(ikm1)
  cqf(2, 2) = qmdia(ikm2)
  cqf(1, 2) = qmoff(ikm1)
  cqf(2, 1) = cqf(1, 2)
!     coefficients to calculate the spin-spin field
  call rinit(4, cgg)
  call rinit(4, cgf)
  cgg(1, 1) = -mj/(kap1+0.5d0)
  cgf(1, 1) = -mj/(-kap1+0.5d0)
  if (nsol==2) then
    cgg(1, 2) = -dsqrt(1.0d0-(mj/(kap1+0.5d0))**2)
    cgg(2, 1) = cgg(1, 2)
    cgg(2, 2) = -mj/(kap2+0.5d0)
    cgf(2, 2) = -mj/(-kap2+0.5d0)
!     CGF(1,2) = -DSQRT( 1.0D0 - (MJ/(- KAP1+0.5D0))**2 )
!     CGF(2,1) = CGF(1,2)
  end if
!     coefficients to calculate the orbital field
  call rinit(4, cfg)
  call rinit(4, cff)
  cfg(1, 1) = mj*(kap1+1.0d0)/(kap1+0.5d0)
  cff(1, 1) = mj*(-kap1+1.0d0)/(-kap1+0.5d0)
  if (nsol==2) then
    cfg(2, 2) = mj*(kap2+1.0d0)/(kap2+0.5d0)
    cfg(1, 2) = 0.5d0*dsqrt(1.0d0-(mj/(kap1+0.5d0))**2)
    cfg(2, 1) = cfg(1, 2)
    cff(2, 2) = mj*(-kap2+1.0d0)/(-kap2+0.5d0)
!     CFF(1,2) = 0.5D0 * DSQRT( 1.0D0 - (MJ/(- KAP1 + 0.5D0))**2 )
!     CFF(2,1) = CFF(1,2)
  end if
! Calculates integrals from 0.0 to jtop
  call hffint(gg, g, g, dovr, r, 0.0d0, nsol, jtop, nrc)
  call hffint(ff, f, f, dovr, r, 0.0d0, nsol, jtop, nrc)
  call hffint(gf, g, f, drdi, r, 0.0d0, nsol, jtop, nrc)
  call hffint(fg, f, g, drdi, r, 0.0d0, nsol, jtop, nrc)
  call rsumupint(shf(1,1,5), f1, gg, cqg, f1, ff, cqf, nsol)
  if (nucleus/=0) then
!     calculates integrals inside nucleus at RNUC in order to get
!     contribution outside the nucleus
    call hffint(gg1, g, g, dovr, r, rnuc, nsol, jlim, nrc)
    call hffint(ff1, f, f, dovr, r, rnuc, nsol, jlim, nrc)
    call hffint(gf1, g, f, drdi, r, rnuc, nsol, jlim, nrc)
    call hffint(fg1, f, g, drdi, r, rnuc, nsol, jlim, nrc)
!     calculates contribution from RNUC to jtop
    do i = 1, nsol
      do j = 1, nsol
        gg(i, j) = gg(i, j) - gg1(i, j)
        ff(i, j) = ff(i, j) - ff1(i, j)
        gf(i, j) = gf(i, j) - gf1(i, j)
        fg(i, j) = fg(i, j) - fg1(i, j)
      end do
    end do
  end if !end of nucleus.eq.0
!     calculates B_sp which is zero inside the nucleus
  call rsumupint(shf(1,1,1), f1, gg, csg, -f1, ff, csf, nsol)
!     calculates hyperfine integrals from 0.0 to RNUC which are added
!     external integrals
  if (nucleus/=0) then
    call hffint(gg2, g, g, drovrn, r, rnuc, nsol, jlim, nrc)
    call hffint(ff2, f, f, drovrn, r, rnuc, nsol, jlim, nrc)
    call hffint(gf2, g, f, drovrn1, r, rnuc, nsol, jlim, nrc)
    call hffint(fg2, f, g, drovrn1, r, rnuc, nsol, jlim, nrc)
    do i = 1, nsol
      do j = 1, nsol
        gg(i, j) = gg(i, j) + gg2(i, j)
        ff(i, j) = ff(i, j) + ff2(i, j)
        gf(i, j) = gf(i, j) + gf2(i, j)
        fg(i, j) = fg(i, j) + fg2(i, j)
      end do
    end do
  end if
!     calculates B_nseo and B_tot
  call rsumupint(shf(1,1,2), f2, gg, cfg, -f2, ff, cff, nsol)
!      CALL RSUMUPINT(SHF(1,1,5),CAUTOG,GF,AME,CAUTOG,FG,AME,NSOL)
!     modifications for B_ssc which is zero outside the nucleus
  if (nucleus/=0) then
    do i = 1, nsol
      do j = 1, nsol
        gg(i, j) = gg2(i, j)
        ff(i, j) = ff2(i, j)
      end do
    end do
  end if
!     calculates B_ssc
  call rsumupint(shf(1,1,3), f2, gg, cgg, -f2, ff, cgf, nsol)
!     for testing purposes write in 4 the sum of 1,2,3
  do k1 = 1, nsol
    do k2 = 1, nsol
      shf(k1, k2, 4) = shf(k1, k2, 1) + shf(k1, k2, 2) + shf(k1, k2, 3)
    end do
  end do

end subroutine

subroutine rsumupint(sum, vg, g, wg, vf, f, wf, n)
  implicit none

! Dummy arguments
  integer :: n
  double precision :: vf, vg
  double precision :: f(n, n), g(n, n), sum(n, n), wf(n, n), wg(n, n)

! Local variables
  integer :: i, j

  do j = 1, n
    do i = 1, n
      sum(i, j) = vg*g(i, j)*wg(i, j) + vf*f(i, j)*wf(i, j)
    end do
  end do
end subroutine

subroutine hffint(gg, ga, gb, dr, r, rnuc, nsol, jtop, nrc)
!     Calculates Hyperfine integrals, extrapolates to zero and
!     intrapolates to exact nuclear radius RNUC
  implicit none

! Dummy arguments
  integer :: jtop, nrc, nsol
  double precision :: rnuc
  double precision :: dr(nrc), ga(nrc, 2), gb(nrc, 2), gg(2, 2), r(nrc)

! Local variables
  integer :: i, k1, k2
  double precision :: x(5), y(5), yi(nrc), zi(nrc)
  double precision :: ylag

  do k1 = 1, nsol
    do k2 = 1, nsol
      do i = 1, jtop
        yi(i) = ga(i, k1)*gb(i, k2)*dr(i)
      end do
      call rint4pts(yi, jtop, zi)
!     Intrapolation
      if (rnuc/=0.0d0) then
        do i = 1, 5
          x(i) = r(jtop-5+i)
          y(i) = zi(jtop-5+i)
        end do
        zi(jtop) = ylag(rnuc, x, y, 0, 4, 5)
      end if
!     Extrapolation
      x(1) = 1.0d0
      x(2) = 6.0d0
      x(3) = 11.0d0
      y(1) = zi(jtop) - zi(1)
      y(2) = zi(jtop) - zi(5)
      y(3) = zi(jtop) - zi(9)
      gg(k1, k2) = ylag(0.0d0, x, y, 0, 2, 3)
    end do
  end do
end subroutine

