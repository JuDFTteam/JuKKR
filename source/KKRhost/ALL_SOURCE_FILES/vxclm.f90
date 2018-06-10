!-------------------------------------------------------------------------------
! SUBROUTINE: VXCLM
!> @brief Add the exchange-correlation-potential to the given potential
!> and if total energies should be calculated (KTE=1) the
!> exchange-correlation-energies are calculated.
!> @details Use as input the charge density times \f$r^2\f$ (rho2ns(...,1)) and
!> in the spin-polarized case (NSPIN=2) the spin density times \f$r^2\f$
!> (rho2ns(...,2)) .
!> The density times \f$4\pi\f$ is generated at an angular mesh.
!> The exchange-correlation potential and the exchange-correlation
!> energy are calculated at those mesh points with a subroutine.
!> In the paramagnetic case the "spin-density" is set equal zero.
!> After that the exchange-correlation potential and in the case of
!> total energies (KTE=1) the exchange-correlation energy are
!> expanded into spherical harmonics.
!> The ex.-cor. potential is added to the given potential.
!> The expansion into spherical harmonics uses the orthogonality
!> of these harmonics.
!> - Therefore a gauss-legendre integration for \f$\theta\f$ and a
!> gauss-tschebyscheff integration for \f$\phi\f$ is used.
!>
!> All needed values for the angular mesh and angular integration
!> are generate in the subroutine sphere.
!> The ex.-cor. potential is extrapolated to the origin only
!> for the lm=1 value .
!> @author B. Drittler, R. Zeller
!> @note
!> - B. Drittler Oct. 1989: Modified for shape functions
!> - R. Zeller Nov. 1993: simplified and modified for Paragon X/PS
!> - R. Zeller 23/6/1996: cor error
!> - Jonathan Chico: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine vxclm(exc, kte, kxc, lmax, nspin, iatyp, rho2ns, v, r, drdi, irws, &
  ircut, ipan, kshape, gsh, ilm_map, imaxsh, ifunm, thetas, yr, wtyr, ijend, &
  lmsp, lmpot, lmxspd, lmmax, irm, lpot, natyp)

  use :: constants
  use :: global_variables

  implicit none

! .. Scalar Arguments
  integer, intent (in) :: irm !< Maximum number of radial points
  integer, intent (in) :: kte !< Calculation of the total energy On/Off (1/0)
  integer, intent (in) :: kxc !< Type of xc-potential 0=vBH 1=MJW 2=VWN 3=PW91
  integer, intent (in) :: lmax !< Maximum l component in wave function expansion
  integer, intent (in) :: irws !< IATYP Entry in the IRWS array with the R point at WS radius
  integer, intent (in) :: ipan !< IATYP Entry in the IPAN array with the number of panels in non-MT-region
  integer, intent (in) :: lpot !< Maximum l component in potential expansion
  integer, intent (in) :: natyp !< Number of kinds of atoms in unit cell
  integer, intent (in) :: iatyp
  integer, intent (in) :: ijend
  integer, intent (in) :: lmpot !< (LPOT+1)**2
  integer, intent (in) :: nspin !< Counter for spin directions
  integer, intent (in) :: lmmax !< (LMAX+1)^2
  integer, intent (in) :: kshape !< Exact treatment of WS cell
  integer, intent (in) :: lmxspd !< (2*LPOT+1)**2
! .. Array Arguments
  integer, dimension (lmxspd), intent (in) :: lmsp !< 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
  integer, dimension (0:ipand), intent (in) :: ircut !< R points of panel borders
  integer, dimension (lmxspd), intent (in) :: ifunm
  integer, dimension (0:lmpot), intent (in) :: imaxsh
  integer, dimension (ngshd, 3), intent (in) :: ilm_map
  double precision, dimension (irm), intent (in) :: r !< Radial mesh ( in units a Bohr)
  double precision, dimension (ijend, lmpot), intent (in) :: yr
  double precision, dimension (ngshd), intent (in) :: gsh
  double precision, dimension (irm), intent (in) :: drdi !< Derivative dr/di
  double precision, dimension (ijend, lmpot), intent (in) :: wtyr
  double precision, dimension (irid, nfund), intent (in) :: thetas !< shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
  double precision, dimension (irm, lmpot, 2), intent (in) :: rho2ns !< radial density
! .. Input/Output variables
  double precision, dimension (0:lpot, natyp), intent (inout) :: exc !< exchange correlation energy
  double precision, dimension (irm, lmpot, 2), intent (inout) :: v
! .. Local Scalars
  integer :: ifun, ij, ipot, ir, irc1, irh, irs1, is, ispin, j, l, lm, lm2, m
  double precision :: elmxc, fpi, fpipr2, vlmxc, vxc1, vxc2, vxc3, factor
! .. Local Arrays
  double precision, dimension (ijend) :: excij
  double precision, dimension (irm, 0:lpot) :: er
  double precision, dimension (ijend, 2) :: vxc
  double precision, dimension (2:3, 2) :: vxcr
  double precision, dimension (ijend, 2) :: fprho
  double precision, dimension (irm, lmpot) :: estor
! .. External Functions
  double precision :: ddot
  external :: ddot
! .. External Subroutines ..
  external :: daxpy, simp3, simpk, vosko, vxcspo
! ..
  write (1337, *) 'Including cutoff of vxc for small density'
  fpi = 4.0d0*pi

!----------------------------------------------------------------------------
! Loop over given representive atoms
!----------------------------------------------------------------------------
  if (kshape/=0) then
    irc1 = ircut(ipan)
    irs1 = ircut(1)
  else
    irc1 = irws
    irs1 = irc1
  end if

  do ispin = 1, nspin
    vxcr(2, ispin) = 0.0d0
    vxcr(3, ispin) = 0.0d0
  end do ! ISPIN
!----------------------------------------------------------------------------
! Initialize for ex.-cor. energy
!----------------------------------------------------------------------------
  if (kte==1) then
    do l = 0, lmax
      exc(l, iatyp) = 0.0d0
      do ir = 1, irc1
        er(ir, l) = 0.0d0
      end do ! IR
    end do ! L
!
    do lm = 1, lmmax
      do ir = 1, irc1
        estor(ir, lm) = 0.0d0
      end do ! IR
    end do ! LM
  end if
!----------------------------------------------------------------------------
! Loop over radial mesh
!----------------------------------------------------------------------------
  do ir = 2, irc1
!-------------------------------------------------------------------------
! Generate the densities on an angular mesh
!-------------------------------------------------------------------------
    do is = 1, 2
      do ij = 1, ijend
        fprho(ij, is) = 0.d0
      end do ! IJ
    end do ! IS

    fpipr2 = fpi/r(ir)**2
    do ispin = 1, nspin
      do lm = 1, lmmax
        call daxpy(ijend, rho2ns(ir,lm,ispin)*fpipr2, yr(1,lm), 1, &
          fprho(1,ispin), 1)
      end do ! LM
    end do
!-------------------------------------------------------------------------
! Calculate the ex.-cor. potential
!-------------------------------------------------------------------------
    if (kxc<=1) then
      call vxcspo(excij, fprho, vxc, kxc, ijend, ijend)
    else
      call vosko(excij, fprho, vxc, ijend, ijend)
    end if

    do ij = 1, ijend
      factor = (1.d0-dexp(-dabs(fprho(ij,1))*1000.d0))
      do ispin = 1, nspin
        vxc(ij, ispin) = vxc(ij, ispin)*factor !cutoff
      end do ! ISPIN
    end do !IJ
!-------------------------------------------------------------------------
! Expand the ex.-cor. potential into spherical harmonics ,
!   using the orthogonality
!-------------------------------------------------------------------------
    do ispin = 1, nspin
!----------------------------------------------------------------------
! Determine the corresponding potential number
!----------------------------------------------------------------------
      ipot = ispin
      do lm = 1, lmmax
        vlmxc = ddot(ijend, vxc(1,ispin), 1, wtyr(1,lm), 1)
        v(ir, lm, ipot) = v(ir, lm, ipot) + vlmxc
!-------------------------------------------------------------------
! store the ex.-c. potential of ir=2 and =3 for the extrapolation
!-------------------------------------------------------------------
        if (lm==1 .and. (ir==2 .or. ir==3)) vxcr(ir, ispin) = vlmxc
      end do ! LM
    end do ! ISPIN
!-------------------------------------------------------------------------
! File er in case of total energies
!-------------------------------------------------------------------------
    if (kte==1) then
!----------------------------------------------------------------------
! expand ex.-cor. energy into spherical harmonics
!   using the orthogonality
!----------------------------------------------------------------------
      do l = 0, lmax
        do m = -l, l
          lm = l*l + l + m + 1
          elmxc = ddot(ijend, excij, 1, wtyr(1,lm), 1)
!----------------------------------------------------------------
! multiply the lm-component of the ex.-cor. energy with the same
! lm-component of the charge density times r**2 and sum over lm
! this corresponds to a integration over the angular .
!----------------------------------------------------------------
          if ((kshape/=0) .and. (ir>irs1)) then
            estor(ir, lm) = elmxc
          else
            er(ir, l) = er(ir, l) + rho2ns(ir, lm, 1)*elmxc
          end if
        end do !M
      end do !L
    end if
  end do !IR

!----------------------------------------------------------------------------
! Integrate er in case of total energies to get exc
!----------------------------------------------------------------------------
  if (kte==1) then
    if (kshape==0) then
      do l = 0, lmax
        call simp3(er(1,l), exc(l,iatyp), 1, irs1, drdi)
      end do !L
    else
      do l = 0, lmax
        do m = -l, l
          lm = l*l + l + m + 1
!----------------------------------------------------------------
! Convolute with shape function
!----------------------------------------------------------------
          do j = imaxsh(lm-1) + 1, imaxsh(lm)
            lm2 = ilm_map(j, 2)
            if (lmsp(ilm_map(j,3))>0) then
              ifun = ifunm(ilm_map(j,3))
              do ir = irs1 + 1, irc1
                irh = ir - irs1
                er(ir, l) = er(ir, l) + rho2ns(ir, lm, 1)*gsh(j)*thetas(irh, &
                  ifun)*estor(ir, lm2)
              end do ! IR
            end if
          end do ! J
        end do ! M
        call simpk(er(1,l), exc(l,iatyp), ipan, ircut, drdi)
      end do ! L
    end if
  end if
!----------------------------------------------------------------------------
! Extrapolate ex.-cor potential to the origin only for lm=1
!----------------------------------------------------------------------------
  do ispin = 1, nspin
    ipot = ispin
    vxc2 = vxcr(2, ispin)
    vxc3 = vxcr(3, ispin)
    vxc1 = vxc2 - r(2)*(vxc3-vxc2)/(r(3)-r(2))
    v(1, 1, ipot) = v(1, 1, ipot) + vxc1
  end do ! ISPIN
!
end subroutine
