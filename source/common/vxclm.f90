!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Add the exchange-correlation-potential to the given potential and if total energies should be calculated (`kte=1`) the exchange-correlation-energies are calculated.
!> Author: B. Drittler
!> Add the exchange-correlation-potential to the given potential and if total 
!> energies should be calculated (`kte=1`) the exchange-correlation-energies are calculated.
!> Use as input the charge density times \(r^2\) (`rho2ns(...,1)`) and in the 
!> spin-polarized case (`nspin=2`) the spin density times \(r^2\) (`rho2ns(...,2)`) .
!> The density times \(4 \pi\) is generated at an angular mesh. tThe exchange-correlation 
!> potential and the exchange-correlation energy are calculated at those mesh 
!> points with a subroutine. In the paramagnetic case the _spin-density_ is set 
!> equal zero. After that the exchange-correlation potential and in the case of
!> total energies (kte=1) the exchange-correlation energy are expanded into spherical 
!> harmonics. The ex.-cor. potential is added to the given potential.
!> The expansion into spherical harmonics uses the orthogonality of these harmonics.
!> Therefore a gauss-legendre integration for \(\theta\) and a gauss-tschebyscheff 
!> integration for \(\phi\) is used .
!> All needed values for the angular mesh and angular integration are generate in 
!> the subroutine sphere. The ex.-cor. potential is extrapolated to the origin only
!> for the lm=1 value .
!------------------------------------------------------------------------------------
!> @note 
!> - Modified for shape functions B. Drittler oct. 1989
!> - Simplified and modified for Paragon X/PS R. Zeller Nov. 1993
!> - Cor error 23/6/1996
!> @endnote
!------------------------------------------------------------------------------------
module mod_vxclm

contains

  !-------------------------------------------------------------------------------
  !> Summary: Add the exchange-correlation-potential to the given potential and if total energies should be calculated (`kte=1`) the exchange-correlation-energies are calculated.
  !> Author: B. Drittler
  !> Category: xc-potential, KKRhost
  !> Deprecated: False 
  !> Add the exchange-correlation-potential to the given potential and if total 
  !> energies should be calculated (`kte=1`) the exchange-correlation-energies are calculated.
  !> Use as input the charge density times \(r^2\) (`rho2ns(...,1)`) and in the 
  !> spin-polarized case (`nspin=2`) the spin density times \(r^2\) (`rho2ns(...,2)`) .
  !> The density times \(4 \pi\) is generated at an angular mesh. tThe exchange-correlation 
  !> potential and the exchange-correlation energy are calculated at those mesh 
  !> points with a subroutine. In the paramagnetic case the _spin-density_ is set 
  !> equal zero. After that the exchange-correlation potential and in the case of
  !> total energies (kte=1) the exchange-correlation energy are expanded into spherical 
  !> harmonics. The ex.-cor. potential is added to the given potential.
  !> The expansion into spherical harmonics uses the orthogonality of these harmonics.
  !> Therefore a gauss-legendre integration for \(\theta\) and a gauss-tschebyscheff 
  !> integration for \(\phi\) is used .
  !> All needed values for the angular mesh and angular integration are generate in 
  !> the subroutine sphere. The ex.-cor. potential is extrapolated to the origin only
  !> for the lm=1 value .
  !-------------------------------------------------------------------------------
  !> @note 
  !> - Modified for shape functions B. Drittler oct. 1989
  !> - Simplified and modified for Paragon X/PS R. Zeller Nov. 1993
  !> - Cor error 23/6/1996
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine vxclm(exc,kte,kxc,lmax,nspin,iatyp,rho2ns,v,r,drdi,irws,ircut,ipan,    &
    kshape,gsh,ilm,imaxsh,ifunm,thetas,yr,wtyr,ijend,lmsp)

    use :: mod_datatypes, only: dp
    use :: global_variables, only: lmxspd, ipand, lmpotd, ngshd, irmd, irid, nfund, lpotd
    use :: mod_vosko, only: vosko
    use :: mod_vxcspo, only: vxcspo
    use :: mod_simpk, only: simpk
    use :: mod_simp3, only: simp3
    use :: mod_constants, only: pi
    implicit none
    ! ..
    ! .. Scalar Arguments ..
    integer, intent(in) :: kte    !! Calculation of the total energy On/Off (1/0)
    integer, intent(in) :: kxc    !! Type of xc-potential 0=vBH 1=MJW 2=VWN 3=PW91
    integer, intent(in) :: ipan   !! Number of panels in non-MT-region
    integer, intent(in) :: irws   !! R point at WS radius
    integer, intent(in) :: lmax   !! Maximum l component in wave function expansion
    integer, intent(in) :: nspin  !! Counter for spin directions
    integer, intent(in) :: iatyp
    integer, intent(in) :: ijend
    integer, intent(in) :: kshape !! Exact treatment of WS cell
    integer, dimension(lmxspd), intent(in) :: lmsp !! 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
    integer, dimension(lmxspd), intent(in) :: ifunm
    integer, dimension(0:ipand), intent(in) :: ircut !! R points of panel borders
    integer, dimension(0:lmpotd), intent(in) :: imaxsh
    integer, dimension(ngshd, 3), intent(in) :: ilm
    real (kind=dp), dimension(irmd), intent(in) :: r    !! Set of real space vectors (in a.u.)
    real (kind=dp), dimension(*), intent(in)    :: gsh
    real (kind=dp), dimension(irmd), intent(in) :: drdi !! Derivative dr/di
    real (kind=dp), dimension(ijend, *), intent(in) :: yr
    real (kind=dp), dimension(ijend, *), intent(in) :: wtyr
    real (kind=dp), dimension(irid, nfund), intent(in) :: thetas !! shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
    real (kind=dp), dimension(irmd, lmpotd, 2), intent(in) :: rho2ns  !! radial density
    ! .. Int/Out variables
    real (kind=dp), dimension(0:lpotd, *), intent(inout) :: exc !! xc energy
    real (kind=dp), dimension(irmd, lmpotd, 2), intent(inout)  :: v
    ! ..
    ! .. Local Scalars ..
    integer :: ifun, ij, ipot, ir, irc1, irh, irs1, is, ispin, j, l, lm, lm2, lmmax0d, m
    real (kind=dp) :: elmxc, fpi, fpipr2, vlmxc, vxc1, vxc2, vxc3, factor
    ! ..
    ! .. Local Arrays ..
    real (kind=dp), dimension(ijend)          :: excij
    real (kind=dp), dimension(irmd, 0:lpotd)  :: er
    real (kind=dp), dimension(ijend, 2)       :: vxc
    real (kind=dp), dimension(2:3, 2)         :: vxcr
    real (kind=dp), dimension(ijend, 2)       :: fprho
    real (kind=dp), dimension(irmd, lmpotd)   :: estor
    ! ..
    ! .. External Functions ..
    real (kind=dp) :: ddot
    external :: ddot

    write (1337, *) 'Including cutoff of vxc for small density'
    fpi =4.0_dp*pi 
    lmmax0d = (lmax+1)*(lmax+1)

    ! loop over given representive atoms

    if (kshape/=0) then
      irc1 = ircut(ipan)
      irs1 = ircut(1)
    else
      irc1 = irws
      irs1 = irc1
    end if

    do ispin = 1, nspin
      vxcr(2, ispin) = 0.0_dp
      vxcr(3, ispin) = 0.0_dp
    end do

    ! initialize for ex.-cor. energy
    if (kte==1) then
      do l = 0, lmax
        exc(l, iatyp) = 0.0_dp
        do ir = 1, irc1
          er(ir, l) = 0.0_dp
        end do
      end do

      do lm = 1, lmmax0d
        do ir = 1, irc1
          estor(ir, lm) = 0.0_dp
        end do
      end do
    end if

    ! loop over radial mesh
    do ir = 2, irc1
      ! generate the densities on an angular mesh
      do is = 1, 2
        do ij = 1, ijend
          fprho(ij, is) = 0._dp
        end do
      end do

      fpipr2 = fpi/r(ir)**2
      do ispin = 1, nspin
        do lm = 1, lmmax0d
          call daxpy(ijend, rho2ns(ir,lm,ispin)*fpipr2, yr(1,lm), 1, fprho(1,ispin), 1)
        end do
      end do

      ! calculate the ex.-cor. potential

      if (kxc<=1) then
        call vxcspo(excij, fprho, vxc, kxc, ijend, ijend)
      else
        call vosko(excij, fprho, vxc, ijend, ijend)
      end if

      do ij = 1, ijend
        factor = (1.0_dp-exp(-abs(fprho(ij,1))*1000.0_dp))
        do ispin = 1, nspin
          vxc(ij, ispin) = vxc(ij, ispin)*factor ! cutoff
        end do
      end do
      ! expand the ex.-cor. potential into spherical harmonics ,
      ! using the orthogonality
      do ispin = 1, nspin
        ! determine the corresponding potential number
        ipot = ispin
        do lm = 1, lmmax0d
          vlmxc = ddot(ijend, vxc(1,ispin), 1, wtyr(1,lm), 1)
          v(ir, lm, ipot) = v(ir, lm, ipot) + vlmxc
          ! store the ex.-c. potential of ir=2 and =3 for the extrapolation
          if (lm==1 .and. (ir==2 .or. ir==3)) vxcr(ir, ispin) = vlmxc
        end do
      end do
      ! file er in case of total energies
      if (kte==1) then
        ! expand ex.-cor. energy into spherical harmonics
        ! using the orthogonality
        do l = 0, lmax
          do m = -l, l
            lm = l*l + l + m + 1
            elmxc = ddot(ijend, excij, 1, wtyr(1,lm), 1)
            ! multiply the lm-component of the ex.-cor. energy with the same
            ! lm-component of the charge density times r**2 and sum over lm
            ! this corresponds to a integration over the angular .
            if ((kshape/=0) .and. (ir>irs1)) then
              estor(ir, lm) = elmxc
            else
              er(ir, l) = er(ir, l) + rho2ns(ir, lm, 1)*elmxc
            end if
          end do
        end do
      end if
    end do

    ! integrate er in case of total energies to get exc
    if (kte==1) then
      if (kshape==0) then
        do l = 0, lmax
          call simp3(er(1,l), exc(l,iatyp), 1, irs1, drdi)
        end do
      else
        do l = 0, lmax
          do m = -l, l
            lm = l*l + l + m + 1
            ! convolute with shape function
            do j = imaxsh(lm-1) + 1, imaxsh(lm)
              lm2 = ilm(j, 2)
              if (lmsp(ilm(j,3))>0) then
                ifun = ifunm(ilm(j,3))
                do ir = irs1 + 1, irc1
                  irh = ir - irs1
                  er(ir, l) = er(ir, l) + rho2ns(ir, lm, 1)*gsh(j)*thetas(irh, ifun)*estor(ir, lm2)
                end do
              end if
            end do
          end do
          call simpk(er(1,l), exc(l,iatyp), ipan, ircut, drdi)
        end do
      end if
    end if

    ! extrapolate ex.-cor potential to the origin only for lm=1
    do ispin = 1, nspin
      ipot = ispin
      vxc2 = vxcr(2, ispin)
      vxc3 = vxcr(3, ispin)
      vxc1 = vxc2 - r(2)*(vxc3-vxc2)/(r(3)-r(2))
      v(1, 1, ipot) = v(1, 1, ipot) + vxc1
    end do

  end subroutine vxclm

end module mod_vxclm
