subroutine vosko(exc, fpirho, vxc, ijend, ijd)
!-----------------------------------------------------------------------
!     calculate the spin-polarized exchange-correlation potential
!     and the spin-polarized exchange-correlation energy from
!     ceperley-alder ( parametrization of vosko, wilk and nusair )
!                                            ( m. manninen )
!     use as input the density generated on an angular mesh (see
!     subroutine vxclm) . fpirho(.,1) contains the charge density
!     times 4 pi and fpirho(.,2) the spin density times 4 pi .
!     then the ex.-cor. potential and the ex.-cor. energy on those
!     mesh points is calculated .
!     the spin-down potential is stored in vxc(.,1) .

!                                  b.drittler    june 1987
!-----------------------------------------------------------------------
  implicit none
!..
!.. Scalar Arguments ..
  integer :: ijend, ijd
!..
!.. Array Arguments ..
  double precision :: exc(*), fpirho(ijd, 2), vxc(ijd, 2)
!..
!.. Local Scalars ..
  double precision :: af, ap, atnf, atnp, beta, bf, bp, cbrt1, cbrt2, cf, cf1, &
    cf2, cf3, cp, cp1, cp2, cp3, dbeta, dfs, duc, duc1, duc2, ec, ecf, ecp, &
    fs, onthrd, qf, qp, rs, s, s4, smag, tf1, tp1, uc0, uc1, uc10, uc2, uc20, &
    ucf, ucp, x, xf0, xfx, xp0, xpx
  integer :: ij
!..
!.. Intrinsic Functions ..
  intrinsic :: abs, atan, log, max, min, sign, sqrt
!..
!.. Save statement ..
  save :: ap, xp0, bp, cp, qp, cp1, cp2, cp3, af, xf0, bf, cf, qf, cf1, cf2, &
    cf3
!..
!.. Data statements ..
  data ap, xp0, bp, cp, qp, cp1, cp2, cp3/0.0621814d0, -0.10498d0, 3.72744d0, &
    12.9352d0, 6.1519908d0, 1.2117833d0, 1.1435257d0, -0.031167608d0/
  data af, xf0, bf, cf, qf, cf1, cf2, cf3/0.0310907d0, -0.32500d0, 7.06042d0, &
    18.0578d0, 4.7309269d0, 2.9847935d0, 2.7100059d0, -0.1446006d0/
!     ..

  onthrd = 1.0d0/3.0d0

!---> loop over the angular mesh points

  do ij = 1, ijend
    fpirho(ij, 1) = max(1.0d-10, fpirho(ij,1))
    smag = sign(1.0d0, fpirho(ij,2))
    fpirho(ij, 2) = smag*min(fpirho(ij,1)-1.0d-10, abs(fpirho(ij,2)))
    rs = (3.d0/fpirho(ij,1))**onthrd
    s = fpirho(ij, 2)/fpirho(ij, 1)
    x = sqrt(rs)
    xpx = x*x + bp*x + cp
    xfx = x*x + bf*x + cf
    s4 = s**4 - 1.d0
    cbrt1 = (1.d0+s)**(1.d0/3.d0)
    cbrt2 = (1.d0-s)**(1.d0/3.d0)
    fs = ((1.d0+s)**(4.d0/3.d0)+(1.d0-s)**(4.d0/3.d0)-2.d0)/ &
      (2.d0**(4.d0/3.d0)-2.d0)
    beta = 1.d0/(2.74208d0+3.182d0*x+0.09873d0*x*x+0.18268d0*x**3)
    dfs = 4.d0/3.d0*(cbrt1-cbrt2)/(2.d0**(4.d0/3.d0)-2.d0)
    dbeta = -(0.27402d0*x+0.09873d0+1.591d0/x)*beta**2
    atnp = atan(qp/(2.d0*x+bp))
    atnf = atan(qf/(2.d0*x+bf))
    ecp = ap*(log(x*x/xpx)+cp1*atnp-cp3*(log((x-xp0)**2/xpx)+cp2*atnp))
    ecf = af*(log(x*x/xfx)+cf1*atnf-cf3*(log((x-xf0)**2/xfx)+cf2*atnf))
    ec = ecp + fs*(ecf-ecp)*(1.d0+s4*beta)

!---> calculate ex.-cor. energy

    exc(ij) = ec - 0.9163306d0/rs - 0.2381735d0/rs*fs
    tp1 = (x*x+bp*x)/xpx
    tf1 = (x*x+bf*x)/xfx
    ucp = ecp - ap/3.d0*(1.d0-tp1-cp3*(x/(x-xp0)-tp1-xp0*x/xpx))
    ucf = ecf - af/3.d0*(1.d0-tf1-cf3*(x/(x-xf0)-tf1-xf0*x/xfx))
    uc0 = ucp + (ucf-ucp)*fs
    uc10 = uc0 - (ecf-ecp)*(s-1.d0)*dfs
    uc20 = uc0 - (ecf-ecp)*(s+1.d0)*dfs
    duc = (ucf-ucp)*beta*s4*fs + (ecf-ecp)*(-rs/3.d0)*dbeta*s4*fs
    duc1 = duc - (ecf-ecp)*beta*(s-1.d0)*(4.d0*s**3*fs+s4*dfs)
    duc2 = duc - (ecf-ecp)*beta*(s+1.d0)*(4.d0*s**3*fs+s4*dfs)
    uc1 = uc10 + duc1
    uc2 = uc20 + duc2

!---> calculate exc.-cor. potential

    vxc(ij, 2) = uc1 - 1.221774d0/rs*cbrt1
    vxc(ij, 1) = uc2 - 1.221774d0/rs*cbrt2
    if (abs(fpirho(ij,1))<=1.0d-10) then
      vxc(ij, 1) = 0.0d0
      vxc(ij, 2) = 0.0d0
    end if

  end do
end subroutine
