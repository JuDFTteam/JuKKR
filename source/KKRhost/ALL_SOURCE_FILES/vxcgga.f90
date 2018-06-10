!-------------------------------------------------------------------------------
! SUBROUTINE: VXCGGA
!> @brief Add the exchange-correlation-potential given by GGA to the potential
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
subroutine vxcgga(exc, kte, kxc, lmax, nspin, iatyp, rho2ns, v, r, drdi, a, &
  irws, ircut, ipan, kshape, gsh, ilm_map, imaxsh, ifunm, thetas, wtyr, ijend, &
  lmsp, thet, ylm, dylmt1, dylmt2, dylmf1, dylmf2, dylmtf, lmpot, lmxspd, &
  lmmax, irm, lpot, natyp)

  use :: constants
  use :: global_variables

  implicit none

! .. Input variables
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
  double precision, intent (in) :: a !< IATYP entry for the array A with the constants for exponential R mesh
! .. Array Arguments
  integer, dimension (lmxspd), intent (in) :: lmsp !< 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
  integer, dimension (0:ipand), intent (in) :: ircut !< R points of panel borders
  integer, dimension (lmxspd), intent (in) :: ifunm
  integer, dimension (0:lmpot), intent (in) :: imaxsh
  integer, dimension (ngshd, 3), intent (in) :: ilm_map
  double precision, dimension (irm), intent (in) :: r !< IATYP entry of the radial mesh ( in units a Bohr)
  double precision, dimension (ngshd), intent (in) :: gsh
  double precision, dimension (irm), intent (in) :: drdi !< IATYP entry of the derivative dr/di
  double precision, dimension (ijend), intent (in) :: thet
  double precision, dimension (ijend, lmpot), intent (in) :: ylm
  double precision, dimension (ijend, lmpot), intent (in) :: wtyr
  double precision, dimension (irid, nfund), intent (in) :: thetas !< IATYP entry of the shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
  double precision, dimension (ijend, lmpot), intent (in) :: dylmf1
  double precision, dimension (ijend, lmpot), intent (in) :: dylmf2
  double precision, dimension (ijend, lmpot), intent (in) :: dylmt1
  double precision, dimension (ijend, lmpot), intent (in) :: dylmt2
  double precision, dimension (ijend, lmpot), intent (in) :: dylmtf
  double precision, dimension (irm, lmpot, 2), intent (in) :: rho2ns !< radial density
! .. Input/Output variables
  double precision, dimension (0:lpot, natyp), intent (inout) :: exc !< exchange correlation energy
  double precision, dimension (irm, lmpot, 2), intent (inout) :: v
! .. Local Scalars
  double precision :: vxc1, vxc2, vxc3, zero, zero1
  double precision :: chgden, dx, elmxc, fpi, r1, r2, rpoint, spiden, vlmxc
  integer :: lm2, m, mesh, nspin2
  integer :: ifun, ipan1, ipot, ir, irc0, irc1, irh, irs1, ispin, j, l, l1max, &
    lm
! .. Local Arrays
  double precision, dimension (ijend) :: excij
  double precision, dimension (irm, 0:lpot) :: er
  double precision, dimension (ijend, 2) :: vxc
  double precision, dimension (irm, lmpot) :: drrl
  double precision, dimension (2:3, 2) :: vxcr
  double precision, dimension (irm, lmpot) :: estor
  double precision, dimension (lmpot, 2) :: rholm
  double precision, dimension (irm, lmpot) :: drrul
  double precision, dimension (irm, lmpot) :: ddrrl
  double precision, dimension (irm, lmpot) :: ddrrul
  double precision, dimension (irm, 2, lmpot) :: rhol
! .. External Functions
  double precision :: ddot
  external :: ddot
! .. External Subroutines ..
  external :: gradrl, mkxcpe, simp3, simpk, mkxcpe2
! .. Intrinsic Functions ..
  intrinsic :: abs, atan, mod
! .. Data statements ..
  data zero, zero1/0.d0, 1.d-12/
!     ..
  write (1337, fmt=*) ' GGA CALCULATION '
  fpi = 4.0d0*pi
!----------------------------------------------------------------------------
! Loop over given representive atoms
!----------------------------------------------------------------------------
  if (kshape/=0) then
    ipan1 = ipan
    irc1 = ircut(ipan)
    irs1 = ircut(1)
    irc0 = 2
    if (krel==1) stop ' REL + FULL POTENTIAL N/A '
  else
    irc1 = irws
    irs1 = irc1
    ipan1 = 1
    irc0 = 2
    if (krel==1) irc0 = 2 + mod(ircut(1), 2)
  end if

  do ispin = 1, nspin
    vxcr(2, ispin) = 0.0d0
    vxcr(3, ispin) = 0.0d0
  end do
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
!
  l1max = lmax + 1
  mesh = irws
  dx = a
!
  if (nspin==2) then
    do lm = 1, lmmax
      do ir = 2, mesh
        r1 = r(ir)
        r2 = r1*r1
        chgden = rho2ns(ir, lm, 1)/r2
        spiden = rho2ns(ir, lm, 2)/r2
        if (abs(chgden)<=zero1) chgden = zero
        if (abs(spiden)<=zero1) spiden = zero
        rhol(ir, 2, lm) = (chgden+spiden)/2.d0
        rhol(ir, 1, lm) = (chgden-spiden)/2.d0
      end do ! IR
! extrapolate
      rhol(1, 1, lm) = rhol(2, 1, lm)
      rhol(1, 2, lm) = rhol(2, 2, lm)
    end do ! LM
  else
!
    do lm = 1, lmmax
      do ir = 2, mesh
        r1 = r(ir)
        r2 = r1*r1
!
        chgden = rho2ns(ir, lm, 1)/r2
        if (abs(chgden)<=zero1) chgden = zero
        rhol(ir, 1, lm) = chgden/2.d0
        rhol(ir, 2, lm) = chgden/2.d0
      end do ! IR
! extrapolate
      rhol(1, 1, lm) = rhol(2, 1, lm)
      rhol(1, 2, lm) = rhol(2, 2, lm)
    end do ! LM
  end if

  call gradrl(nspin, mesh, l1max, dx, rhol, r, drdi, ipan1, ipand, ircut, &
    drrl, ddrrl, drrul, ddrrul, irm, lmpot)

!----------------------------------------------------------------------------
! Loop over radial mesh
!----------------------------------------------------------------------------

  do ir = irc0, irc1
    rpoint = r(ir)
!-------------------------------------------------------------------------
! Calculate the ex.-cor. potential
!-------------------------------------------------------------------------
    nspin2 = 2

    do ispin = 1, nspin2
      do lm = 1, lmmax
        rholm(lm, ispin) = rhol(ir, ispin, lm)
      end do
    end do
!    only for spin-polarized
!
! PW91 functional
    if (kxc==3) then
      call mkxcpe(nspin2, ir, ijend, l1max, rpoint, rholm, vxc, excij, thet, &
        ylm, dylmt1, dylmt2, dylmf1, dylmf2, dylmtf, drrl, ddrrl, drrul, &
        ddrrul, irm, lmpot)
! PBE functional
    else if (kxc==4) then
      call mkxcpe2(ir, ijend, rpoint, rholm, vxc, excij, ylm, dylmt1, dylmf1, &
        dylmf2, dylmtf, drrl, ddrrl, drrul, ddrrul, irm, lmpot, lmmax, &
        .false.)
! PBEsol functional
    else if (kxc==5) then
      call mkxcpe2(ir, ijend, rpoint, rholm, vxc, excij, ylm, dylmt1, dylmf1, &
        dylmf2, dylmtf, drrl, ddrrl, drrul, ddrrul, irm, lmpot, lmmax, .true.)
    else
      write (1337, *) ' KXC ???'
      stop
    end if
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
! Store the ex.-c. potential of ir=2 and =3 for the extrapolation
!-------------------------------------------------------------------
        if (lm==1 .and. (ir==2 .or. ir==3)) vxcr(ir, ispin) = vlmxc
      end do ! LM
    end do ! ISPIN
!-------------------------------------------------------------------------
! File er in case of total energies
!-------------------------------------------------------------------------
    if (kte==1) then
!----------------------------------------------------------------------
! Expand ex.-cor. energy into spherical harmonics
!   using the orthogonality
!----------------------------------------------------------------------
      do l = 0, lmax
        do m = -l, l
          lm = l*l + l + m + 1
          elmxc = ddot(ijend, excij, 1, wtyr(1,lm), 1)
!----------------------------------------------------------------
! Multiply the lm-component of the ex.-cor. energy with the same
! lm-component of the charge density times r**2 and sum over lm
! this corresponds to a integration over the angular .
!----------------------------------------------------------------
          if ((kshape/=0) .and. (ir>irs1)) then
            estor(ir, lm) = elmxc
          else
            er(ir, l) = er(ir, l) + rho2ns(ir, lm, 1)*elmxc
          end if
        end do ! M
      end do ! L
    end if
  end do !IR
!----------------------------------------------------------------------------
! Integrate er in case of total energies to get exc
!----------------------------------------------------------------------------
  if (kte==1) then
    if (kshape==0) then
      do l = 0, lmax
        call simp3(er(1,l), exc(l,iatyp), 1, irs1, drdi)
      end do
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
        call simpk(er(1,l), exc(l,iatyp), ipan1, ircut, drdi)
      end do ! L
    end if
  end if
!----------------------------------------------------------------------------
! Extrapolate ex.-cor potential to the origin only for lm=1
!----------------------------------------------------------------------------
  do ispin = 1, nspin
    ipot = ispin
!
    vxc2 = vxcr(2, ispin)
    vxc3 = vxcr(3, ispin)
    vxc1 = vxc2 - r(2)*(vxc3-vxc2)/(r(3)-r(2))
!
    v(1, 1, ipot) = v(1, 1, ipot) + vxc1
  end do
!
end subroutine
