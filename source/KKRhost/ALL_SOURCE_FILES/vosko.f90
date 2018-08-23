module mod_vosko

contains

subroutine vosko(exc, fpirho, vxc, ijend, ijd)
  use :: mod_datatypes, only: dp
  ! -----------------------------------------------------------------------
  ! calculate the spin-polarized exchange-correlation potential
  ! and the spin-polarized exchange-correlation energy from
  ! ceperley-alder ( parametrization of vosko, wilk and nusair )
  ! ( m. manninen )
  ! use as input the density generated on an angular mesh (see
  ! subroutine vxclm) . fpirho(.,1) contains the charge density
  ! times 4 pi and fpirho(.,2) the spin density times 4 pi .
  ! then the ex.-cor. potential and the ex.-cor. energy on those
  ! mesh points is calculated .
  ! the spin-down potential is stored in vxc(.,1) .

  ! b.drittler    june 1987
  ! -----------------------------------------------------------------------
  implicit none
  ! ..
  ! .. Scalar Arguments ..
  integer :: ijend, ijd
  ! ..
  ! .. Array Arguments ..
  real (kind=dp) :: exc(*), fpirho(ijd, 2), vxc(ijd, 2)
  ! ..
  ! .. Local Scalars ..
  real (kind=dp) :: af, ap, atnf, atnp, beta, bf, bp, cbrt1, cbrt2, cf, cf1, &
    cf2, cf3, cp, cp1, cp2, cp3, dbeta, dfs, duc, duc1, duc2, ec, ecf, ecp, &
    fs, onthrd, qf, qp, rs, s, s4, smag, tf1, tp1, uc0, uc1, uc10, uc2, uc20, &
    ucf, ucp, x, xf0, xfx, xp0, xpx
  integer :: ij
  ! ..
  ! .. Intrinsic Functions ..
  intrinsic :: abs, atan, log, max, min, sign, sqrt
  ! ..
  ! .. Save statement ..
  save :: ap, xp0, bp, cp, qp, cp1, cp2, cp3, af, xf0, bf, cf, qf, cf1, cf2, &
    cf3
  ! ..
  ! .. Data statements ..
  data ap, xp0, bp, cp, qp, cp1, cp2, cp3/0.0621814e0_dp, -0.10498e0_dp, &
    3.72744e0_dp, 12.9352e0_dp, 6.1519908e0_dp, 1.2117833e0_dp, &
    1.1435257e0_dp, -0.031167608e0_dp/
  data af, xf0, bf, cf, qf, cf1, cf2, cf3/0.0310907e0_dp, -0.32500e0_dp, &
    7.06042e0_dp, 18.0578e0_dp, 4.7309269e0_dp, 2.9847935e0_dp, &
    2.7100059e0_dp, -0.1446006e0_dp/
  ! ..

  onthrd = 1.0e0_dp/3.0e0_dp

  ! ---> loop over the angular mesh points

  do ij = 1, ijend
    fpirho(ij, 1) = max(1.0e-10_dp, fpirho(ij,1))
    smag = sign(1.0e0_dp, fpirho(ij,2))
    fpirho(ij, 2) = smag*min(fpirho(ij,1)-1.0e-10_dp, abs(fpirho(ij,2)))
    rs = (3.e0_dp/fpirho(ij,1))**onthrd
    s = fpirho(ij, 2)/fpirho(ij, 1)
    x = sqrt(rs)
    xpx = x*x + bp*x + cp
    xfx = x*x + bf*x + cf
    s4 = s**4 - 1.e0_dp
    cbrt1 = (1.e0_dp+s)**(1.e0_dp/3.e0_dp)
    cbrt2 = (1.e0_dp-s)**(1.e0_dp/3.e0_dp)
    fs = ((1.e0_dp+s)**(4.e0_dp/3.e0_dp)+(1.e0_dp-s)**(4.e0_dp/3.e0_dp)- &
      2.e0_dp)/(2.e0_dp**(4.e0_dp/3.e0_dp)-2.e0_dp)
    beta = 1.e0_dp/(2.74208e0_dp+3.182e0_dp*x+0.09873e0_dp*x*x+0.18268e0_dp*x &
      **3)
    dfs = 4.e0_dp/3.e0_dp*(cbrt1-cbrt2)/(2.e0_dp**(4.e0_dp/3.e0_dp)-2.e0_dp)
    dbeta = -(0.27402e0_dp*x+0.09873e0_dp+1.591e0_dp/x)*beta**2
    atnp = atan(qp/(2.e0_dp*x+bp))
    atnf = atan(qf/(2.e0_dp*x+bf))
    ecp = ap*(log(x*x/xpx)+cp1*atnp-cp3*(log((x-xp0)**2/xpx)+cp2*atnp))
    ecf = af*(log(x*x/xfx)+cf1*atnf-cf3*(log((x-xf0)**2/xfx)+cf2*atnf))
    ec = ecp + fs*(ecf-ecp)*(1.e0_dp+s4*beta)

    ! ---> calculate ex.-cor. energy

    exc(ij) = ec - 0.9163306e0_dp/rs - 0.2381735e0_dp/rs*fs
    tp1 = (x*x+bp*x)/xpx
    tf1 = (x*x+bf*x)/xfx
    ucp = ecp - ap/3.e0_dp*(1.e0_dp-tp1-cp3*(x/(x-xp0)-tp1-xp0*x/xpx))
    ucf = ecf - af/3.e0_dp*(1.e0_dp-tf1-cf3*(x/(x-xf0)-tf1-xf0*x/xfx))
    uc0 = ucp + (ucf-ucp)*fs
    uc10 = uc0 - (ecf-ecp)*(s-1.e0_dp)*dfs
    uc20 = uc0 - (ecf-ecp)*(s+1.e0_dp)*dfs
    duc = (ucf-ucp)*beta*s4*fs + (ecf-ecp)*(-rs/3.e0_dp)*dbeta*s4*fs
    duc1 = duc - (ecf-ecp)*beta*(s-1.e0_dp)*(4.e0_dp*s**3*fs+s4*dfs)
    duc2 = duc - (ecf-ecp)*beta*(s+1.e0_dp)*(4.e0_dp*s**3*fs+s4*dfs)
    uc1 = uc10 + duc1
    uc2 = uc20 + duc2

    ! ---> calculate exc.-cor. potential

    vxc(ij, 2) = uc1 - 1.221774e0_dp/rs*cbrt1
    vxc(ij, 1) = uc2 - 1.221774e0_dp/rs*cbrt2
    if (abs(fpirho(ij,1))<=1.0e-10_dp) then
      vxc(ij, 1) = 0.0e0_dp
      vxc(ij, 2) = 0.0e0_dp
    end if

  end do
end subroutine vosko

end module mod_vosko
