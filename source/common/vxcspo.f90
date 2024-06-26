!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Calculate the spin-polarized exchange-correlation potential
!> Author: B. Drittler
!> Calculate the spin-polarized exchange-correlation potential and the 
!> spin-polarized exchange-correlation energy kxc=0 means : spin-polarized 
!> exchange-correlation potential U. Von Barth and l.hedin, J. Phys. C5,1629 (1972)
!> with parametrization of Moruzzi, Janak, Williams kxc=1 means : spin-polarized 
!> exchange-correlation potential U. Von Barth and L. Hedin, J. Phys. C5,1629 (1972)
!> with parametrization of Von Barth, Hedin
!> use as input the density generated on an angular mesh (see subroutine `vxclm`).
!> `fpirho(.,1)` contains the charge density times \(4\pi\) and `fpirho(.,2)` 
!> the spin density times \(4 \pi\).
!> Then the ex.-cor. potential and the ex.-cor. energy on those mesh points is calculated .
!> the spin-down potential is stored in `vxc(.,1)`.
!------------------------------------------------------------------------------------
module mod_vxcspo
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculate the spin-polarized exchange-correlation potential and the spin-polarized exchange-correlation energy .
  !> Author: B. Drittler
  !> Category: xc-potential, KKRhost
  !> Deprecated: False 
  !> Calculate the spin-polarized exchange-correlation potential and the 
  !> spin-polarized exchange-correlation energy kxc=0 means : spin-polarized 
  !> exchange-correlation potential U. Von Barth and l.hedin, J. Phys. C5,1629 (1972)
  !> with parametrization of Moruzzi, Janak, Williams kxc=1 means : spin-polarized 
  !> exchange-correlation potential U. Von Barth and L. Hedin, J. Phys. C5,1629 (1972)
  !> with parametrization of Von Barth, Hedin
  !> use as input the density generated on an angular mesh (see subroutine `vxclm`).
  !> `fpirho(.,1)` contains the charge density times \(4\pi\) and `fpirho(.,2)` 
  !> the spin density times \(4 \pi\).
  !> Then the ex.-cor. potential and the ex.-cor. energy on those mesh points is calculated .
  !> the spin-down potential is stored in `vxc(.,1)`.
  !-------------------------------------------------------------------------------
  subroutine vxcspo(exc, fpirho, vxc, kxc, ijend, ijd)

    ! .. Input variables
    integer, intent(in) :: kxc !! Type of xc-potential 0=vBH 1=MJW 2=VWN 3=PW91
    integer, intent(in) :: ijd
    integer, intent(in) :: ijend
    ! .. Output variables
    real (kind=dp), dimension(*), intent(out) :: exc !! xc-energy
    real (kind=dp), dimension(ijd, 2), intent(out) :: vxc !! spin polarized xc-potential
    ! .. In/Out variables
    real (kind=dp), dimension(ijd, 2), intent(inout) :: fpirho !! spin resolved density
    ! ..
    ! .. Local Scalars ..
    integer :: ij
    real (kind=dp) :: cex, cf, cfln, cfmjw, cfvbh, cp, cpln, cpmjw, cpvbh, d1, d2
    real (kind=dp) :: dcfx, excfrs, excprs, exfrs, exprs, fac, ff, onthrd, rf, rfmjw
    real (kind=dp) :: rfvbh, rp, rpmjw, rpvbh, rs, smag, te1b3, vxcc, x, xfac
    ! ..
    ! .. Intrinsic Functions ..
    intrinsic :: abs, sign, log, max, min
    ! ..
    ! .. Statement Functions ..
    real (kind=dp) :: f
    ! ..
    ! .. Save statement ..
    save :: cpmjw, cfmjw, rpmjw, rfmjw, cpvbh, cfvbh, rpvbh, rfvbh, ff, cex, onthrd, te1b3
    ! ..
    ! .. Data statements ..

    ! ---> ff=1/(2**(1/3)-1) , cex=2*(3/(2*pi))**(2/3) , te1b3=2**(1/3)

    data cpmjw, cfmjw, rpmjw, rfmjw/0.045e0_dp, 0.0225e0_dp, 21.e0_dp, 52.916684096e0_dp/
    data cpvbh, cfvbh, rpvbh, rfvbh/0.0504e0_dp, 0.0254e0_dp, 30.e0_dp, 75.e0_dp/
    data ff, cex/3.847322101863e0_dp, 1.221774115422e0_dp/
    data onthrd, te1b3/0.333333333333e0_dp, 1.259921049899e0_dp/
    ! ..
    ! .. Statement Function definitions ..

    f(x) = (1.e0_dp+x*x*x)*log(1.e0_dp+1.e0_dp/x) + 0.5e0_dp*x - x*x - 1.0e0_dp/3.0e0_dp
    ! ..

    ! ---> get key dependent the right parameters

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

    ! ---> loop over the angular mesh points

    do ij = 1, ijend
      fpirho(ij, 1) = max(1.0e-10_dp, fpirho(ij,1))
      smag = sign(1.0e0_dp, fpirho(ij,2))
      fpirho(ij, 2) = smag*min(fpirho(ij,1)-1.0e-10_dp, abs(fpirho(ij,2)))
      rs = (3.e0_dp/fpirho(ij,1))**onthrd
      cpln = cp*log(1.e0_dp+rp/rs)
      cfln = cf*log(1.e0_dp+rf/rs)
      dcfx = (cf*f(rs/rf)-cp*f(rs/rp))*4.e0_dp*onthrd
      d1 = (1.e0_dp+fpirho(ij,2)/fpirho(ij,1))**onthrd
      d2 = (1.e0_dp-fpirho(ij,2)/fpirho(ij,1))**onthrd
      fac = (d1**4+d2**4-2.e0_dp)*0.5e0_dp

      ! ---> calculate ex.-cor. energy

      exprs = -0.75e0_dp*cex/rs
      exfrs = exprs*te1b3
      excprs = exprs - cp*f(rs/rp)
      excfrs = exfrs - cf*f(rs/rf)
      exc(ij) = excprs + (excfrs-excprs)*fac*ff

      ! ---> calculate ex.-cor. potential

      vxcc = -cpln + (fac*(cpln-cfln+dcfx)+dcfx)*ff
      xfac = -cex/rs - dcfx*ff
      vxc(ij, 2) = vxcc + d1*xfac
      vxc(ij, 1) = vxcc + d2*xfac
    end do
  end subroutine vxcspo

end module mod_vxcspo
