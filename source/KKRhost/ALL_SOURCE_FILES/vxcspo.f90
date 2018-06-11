    Subroutine vxcspo(exc, fpirho, vxc, kxc, ijend, ijd)
      Use mod_datatypes, Only: dp
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
      Integer :: ijend, kxc, ijd
!..
!.. Array Arguments ..
      Real (Kind=dp) :: exc(*), fpirho(ijd, 2), vxc(ijd, 2)
!..
!.. Local Scalars ..
      Real (Kind=dp) :: cex, cf, cfln, cfmjw, cfvbh, cp, cpln, cpmjw, cpvbh, &
        d1, d2, dcfx, excfrs, excprs, exfrs, exprs, fac, ff, onthrd, rf, &
        rfmjw, rfvbh, rp, rpmjw, rpvbh, rs, smag, te1b3, vxcc, x, xfac
      Integer :: ij
!..
!.. Intrinsic Functions ..
      Intrinsic :: abs, sign, log, max, min
!..
!.. Statement Functions ..
      Real (Kind=dp) :: f
!..
!.. Save statement ..
      Save :: cpmjw, cfmjw, rpmjw, rfmjw, cpvbh, cfvbh, rpvbh, rfvbh, ff, cex, &
        onthrd, te1b3
!..
!.. Data statements ..
!
!---> ff=1/(2**(1/3)-1) , cex=2*(3/(2*pi))**(2/3) , te1b3=2**(1/3)
!
      Data cpmjw, cfmjw, rpmjw, rfmjw/0.045E0_dp, 0.0225E0_dp, 21.E0_dp, &
        52.916684096E0_dp/
      Data cpvbh, cfvbh, rpvbh, rfvbh/0.0504E0_dp, 0.0254E0_dp, 30.E0_dp, &
        75.E0_dp/
      Data ff, cex/3.847322101863E0_dp, 1.221774115422E0_dp/
      Data onthrd, te1b3/0.333333333333E0_dp, 1.259921049899E0_dp/
!     ..
!     .. Statement Function definitions ..

      f(x) = (1.E0_dp+x*x*x)*log(1.E0_dp+1.E0_dp/x) + 0.5E0_dp*x - x*x - &
        1.0E0_dp/3.0E0_dp
!     ..

!---> get key dependent the right parameters

      If (kxc==1) Then
        cp = cpvbh
        cf = cfvbh
        rp = rpvbh
        rf = rfvbh

      Else
        cp = cpmjw
        cf = cfmjw
        rp = rpmjw
        rf = rfmjw
      End If

!---> loop over the angular mesh points

      Do ij = 1, ijend
        fpirho(ij, 1) = max(1.0E-10_dp, fpirho(ij,1))
        smag = sign(1.0E0_dp, fpirho(ij,2))
        fpirho(ij, 2) = smag*min(fpirho(ij,1)-1.0E-10_dp, abs(fpirho(ij,2)))
        rs = (3.E0_dp/fpirho(ij,1))**onthrd
        cpln = cp*log(1.E0_dp+rp/rs)
        cfln = cf*log(1.E0_dp+rf/rs)
        dcfx = (cf*f(rs/rf)-cp*f(rs/rp))*4.E0_dp*onthrd
        d1 = (1.E0_dp+fpirho(ij,2)/fpirho(ij,1))**onthrd
        d2 = (1.E0_dp-fpirho(ij,2)/fpirho(ij,1))**onthrd
        fac = (d1**4+d2**4-2.E0_dp)*0.5E0_dp

!---> calculate ex.-cor. energy

        exprs = -0.75E0_dp*cex/rs
        exfrs = exprs*te1b3
        excprs = exprs - cp*f(rs/rp)
        excfrs = exfrs - cf*f(rs/rf)
        exc(ij) = excprs + (excfrs-excprs)*fac*ff

!---> calculate ex.-cor. potential

        vxcc = -cpln + (fac*(cpln-cfln+dcfx)+dcfx)*ff
        xfac = -cex/rs - dcfx*ff
        vxc(ij, 2) = vxcc + d1*xfac
        vxc(ij, 1) = vxcc + d2*xfac
      End Do
    End Subroutine
