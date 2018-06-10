subroutine vxcspo(exc, fpirho, vxc, kxc, ijend, ijd)
!-----------------------------------------------------------------------
!     calculate the spin-polarized exchange-correlation potential
!     and the spin-polarized exchange-correlation energy .
!     kxc=0 means : spin-polarized exchange-correlation potential
!                   u. von barth and l.hedin, j.phys.c5,1629 (1972)
!                   with parametrization of moruzzi,janak,williams
!     kxc=1 means : spin-polarized exchange-correlation potential
!                   u. von barth and l.hedin, j.phys.c5,1629 (1972)
!                   with parametrization of von barth,hedin

!     use as input the density generated on an angular mesh (see
!     subroutine vxclm) . fpirho(.,1) contains the charge density
!     times 4 pi and fpirho(.,2) the spin density times 4 pi .
!     then the ex.-cor. potential and the ex.-cor. energy on those
!     mesh points is calculated .
!     the spin-down potential is stored in vxc(.,1) .

!                                  b.drittler    june 1987
!-----------------------------------------------------------------------
!..
!.. Scalar Arguments ..
  integer :: ijend, kxc, ijd
!..
!.. Array Arguments ..
  double precision :: exc(*), fpirho(ijd, 2), vxc(ijd, 2)
!..
!.. Local Scalars ..
  double precision :: cex, cf, cfln, cfmjw, cfvbh, cp, cpln, cpmjw, cpvbh, d1, &
    d2, dcfx, excfrs, excprs, exfrs, exprs, fac, ff, onthrd, rf, rfmjw, rfvbh, &
    rp, rpmjw, rpvbh, rs, smag, te1b3, vxcc, x, xfac
  integer :: ij
!..
!.. Intrinsic Functions ..
  intrinsic :: abs, dsign, log, max, min
!..
!.. Statement Functions ..
  double precision :: f
!..
!.. Save statement ..
  save :: cpmjw, cfmjw, rpmjw, rfmjw, cpvbh, cfvbh, rpvbh, rfvbh, ff, cex, &
    onthrd, te1b3
!..
!.. Data statements ..
!
!---> ff=1/(2**(1/3)-1) , cex=2*(3/(2*pi))**(2/3) , te1b3=2**(1/3)
!
  data cpmjw, cfmjw, rpmjw, rfmjw/0.045d0, 0.0225d0, 21.d0, 52.916684096d0/
  data cpvbh, cfvbh, rpvbh, rfvbh/0.0504d0, 0.0254d0, 30.d0, 75.d0/
  data ff, cex/3.847322101863d0, 1.221774115422d0/
  data onthrd, te1b3/0.333333333333d0, 1.259921049899d0/
!     ..
!     .. Statement Function definitions ..

  f(x) = (1.d0+x*x*x)*log(1.d0+1.d0/x) + 0.5d0*x - x*x - 1.0d0/3.0d0
!     ..

!---> get key dependent the right parameters

  if (kxc==1) then
    cp = cpvbh
    cf = cfvbh
    rp = rpvbh
    rf = rfvbh

  else
    cp = cpmjw
    cf = cfmjw
    rp = rpmjw
    rf = rfmjw
  end if

!---> loop over the angular mesh points

  do ij = 1, ijend
    fpirho(ij, 1) = max(1.0d-10, fpirho(ij,1))
    smag = dsign(1.0d0, fpirho(ij,2))
    fpirho(ij, 2) = smag*min(fpirho(ij,1)-1.0d-10, abs(fpirho(ij,2)))
    rs = (3.d0/fpirho(ij,1))**onthrd
    cpln = cp*log(1.d0+rp/rs)
    cfln = cf*log(1.d0+rf/rs)
    dcfx = (cf*f(rs/rf)-cp*f(rs/rp))*4.d0*onthrd
    d1 = (1.d0+fpirho(ij,2)/fpirho(ij,1))**onthrd
    d2 = (1.d0-fpirho(ij,2)/fpirho(ij,1))**onthrd
    fac = (d1**4+d2**4-2.d0)*0.5d0

!---> calculate ex.-cor. energy

    exprs = -0.75d0*cex/rs
    exfrs = exprs*te1b3
    excprs = exprs - cp*f(rs/rp)
    excfrs = exfrs - cf*f(rs/rf)
    exc(ij) = excprs + (excfrs-excprs)*fac*ff

!---> calculate ex.-cor. potential

    vxcc = -cpln + (fac*(cpln-cfln+dcfx)+dcfx)*ff
    xfac = -cex/rs - dcfx*ff
    vxc(ij, 2) = vxcc + d1*xfac
    vxc(ij, 1) = vxcc + d2*xfac
  end do
end subroutine
