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
    Subroutine vxcgga(exc, kte, kxc, lmax, nspin, iatyp, rho2ns, v, r, drdi, &
      a, irws, ircut, ipan, kshape, gsh, ilm_map, imaxsh, ifunm, thetas, wtyr, &
      ijend, lmsp, thet, ylm, dylmt1, dylmt2, dylmf1, dylmf2, dylmtf, lmpot, &
      lmmax, lpot, natyp)

      Use constants
      Use global_variables
      Use mod_datatypes, Only: dp

      Implicit None

! .. Input variables
      Integer, Intent (In) :: kte !< Calculation of the total energy On/Off (1/0)
      Integer, Intent (In) :: kxc !< Type of xc-potential 0=vBH 1=MJW 2=VWN 3=PW91
      Integer, Intent (In) :: lmax !< Maximum l component in wave function expansion
      Integer, Intent (In) :: irws !< IATYP Entry in the IRWS array with the R point at WS radius
      Integer, Intent (In) :: ipan !< IATYP Entry in the IPAN array with the number of panels in non-MT-region
      Integer, Intent (In) :: lpot !< Maximum l component in potential expansion
      Integer, Intent (In) :: natyp !< Number of kinds of atoms in unit cell
      Integer, Intent (In) :: iatyp
      Integer, Intent (In) :: ijend
      Integer, Intent (In) :: lmpot !< (LPOT+1)**2
      Integer, Intent (In) :: nspin !< Counter for spin directions
      Integer, Intent (In) :: lmmax !< (LMAX+1)^2
      Integer, Intent (In) :: kshape !< Exact treatment of WS cell
      Real (Kind=dp), Intent (In) :: a !< IATYP entry for the array A with the constants for exponential R mesh
! .. Array Arguments
      Integer, Dimension (lmxspd), Intent (In) :: lmsp !< 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
      Integer, Dimension (0:ipand), Intent (In) :: ircut !< R points of panel borders
      Integer, Dimension (lmxspd), Intent (In) :: ifunm
      Integer, Dimension (0:lmpot), Intent (In) :: imaxsh
      Integer, Dimension (ngshd, 3), Intent (In) :: ilm_map
      Real (Kind=dp), Dimension (irmd), Intent (In) :: r !< IATYP entry of the radial mesh ( in units a Bohr)
      Real (Kind=dp), Dimension (ngshd), Intent (In) :: gsh
      Real (Kind=dp), Dimension (irmd), Intent (In) :: drdi !< IATYP entry of the derivative dr/di
      Real (Kind=dp), Dimension (ijend), Intent (In) :: thet
      Real (Kind=dp), Dimension (ijend, lmpot), Intent (In) :: ylm
      Real (Kind=dp), Dimension (ijend, lmpot), Intent (In) :: wtyr
      Real (Kind=dp), Dimension (irid, nfund), Intent (In) :: thetas !< IATYP entry of the shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
      Real (Kind=dp), Dimension (ijend, lmpot), Intent (In) :: dylmf1
      Real (Kind=dp), Dimension (ijend, lmpot), Intent (In) :: dylmf2
      Real (Kind=dp), Dimension (ijend, lmpot), Intent (In) :: dylmt1
      Real (Kind=dp), Dimension (ijend, lmpot), Intent (In) :: dylmt2
      Real (Kind=dp), Dimension (ijend, lmpot), Intent (In) :: dylmtf
      Real (Kind=dp), Dimension (irmd, lmpot, 2), Intent (In) :: rho2ns !< radial density
! .. Input/Output variables
      Real (Kind=dp), Dimension (0:lpot, natyp), Intent (Inout) :: exc !< exchange correlation energy
      Real (Kind=dp), Dimension (irmd, lmpot, 2), Intent (Inout) :: v
! .. Local Scalars
      Real (Kind=dp) :: vxc1, vxc2, vxc3, zero, zero1
      Real (Kind=dp) :: chgden, dx, elmxc, fpi, r1, r2, rpoint, spiden, vlmxc
      Integer :: lm2, m, mesh, nspin2
      Integer :: ifun, ipan1, ipot, ir, irc0, irc1, irh, irs1, ispin, j, l, &
        l1max, lm
! .. Local Arrays
      Real (Kind=dp), Dimension (ijend) :: excij
      Real (Kind=dp), Dimension (irmd, 0:lpot) :: er
      Real (Kind=dp), Dimension (ijend, 2) :: vxc
      Real (Kind=dp), Dimension (irmd, lmpot) :: drrl
      Real (Kind=dp), Dimension (2:3, 2) :: vxcr
      Real (Kind=dp), Dimension (irmd, lmpot) :: estor
      Real (Kind=dp), Dimension (lmpot, 2) :: rholm
      Real (Kind=dp), Dimension (irmd, lmpot) :: drrul
      Real (Kind=dp), Dimension (irmd, lmpot) :: ddrrl
      Real (Kind=dp), Dimension (irmd, lmpot) :: ddrrul
      Real (Kind=dp), Dimension (irmd, 2, lmpot) :: rhol
! .. External Functions
      Real (Kind=dp) :: ddot
      External :: ddot
! .. External Subroutines ..
      External :: gradrl, mkxcpe, simp3, simpk, mkxcpe2
! .. Intrinsic Functions ..
      Intrinsic :: abs, atan, mod
! .. Data statements ..
      Data zero, zero1/0.E0_dp, 1.E-12_dp/
!     ..
      Write (1337, Fmt=*) ' GGA CALCULATION '
      fpi = 4.0E0_dp*pi
!----------------------------------------------------------------------------
! Loop over given representive atoms
!----------------------------------------------------------------------------
      If (kshape/=0) Then
        ipan1 = ipan
        irc1 = ircut(ipan)
        irs1 = ircut(1)
        irc0 = 2
        If (krel==1) Stop ' REL + FULL POTENTIAL N/A '
      Else
        irc1 = irws
        irs1 = irc1
        ipan1 = 1
        irc0 = 2
        If (krel==1) irc0 = 2 + mod(ircut(1), 2)
      End If

      Do ispin = 1, nspin
        vxcr(2, ispin) = 0.0E0_dp
        vxcr(3, ispin) = 0.0E0_dp
      End Do
!----------------------------------------------------------------------------
! Initialize for ex.-cor. energy
!----------------------------------------------------------------------------
      If (kte==1) Then
        Do l = 0, lmax
          exc(l, iatyp) = 0.0E0_dp
          Do ir = 1, irc1
            er(ir, l) = 0.0E0_dp
          End Do ! IR
        End Do ! L
!
        Do lm = 1, lmmax
          Do ir = 1, irc1
            estor(ir, lm) = 0.0E0_dp
          End Do ! IR
        End Do ! LM
      End If
!
      l1max = lmax + 1
      mesh = irws
      dx = a
!
      If (nspin==2) Then
        Do lm = 1, lmmax
          Do ir = 2, mesh
            r1 = r(ir)
            r2 = r1*r1
            chgden = rho2ns(ir, lm, 1)/r2
            spiden = rho2ns(ir, lm, 2)/r2
            If (abs(chgden)<=zero1) chgden = zero
            If (abs(spiden)<=zero1) spiden = zero
            rhol(ir, 2, lm) = (chgden+spiden)/2.E0_dp
            rhol(ir, 1, lm) = (chgden-spiden)/2.E0_dp
          End Do ! IR
! extrapolate
          rhol(1, 1, lm) = rhol(2, 1, lm)
          rhol(1, 2, lm) = rhol(2, 2, lm)
        End Do ! LM
      Else
!
        Do lm = 1, lmmax
          Do ir = 2, mesh
            r1 = r(ir)
            r2 = r1*r1
!
            chgden = rho2ns(ir, lm, 1)/r2
            If (abs(chgden)<=zero1) chgden = zero
            rhol(ir, 1, lm) = chgden/2.E0_dp
            rhol(ir, 2, lm) = chgden/2.E0_dp
          End Do ! IR
! extrapolate
          rhol(1, 1, lm) = rhol(2, 1, lm)
          rhol(1, 2, lm) = rhol(2, 2, lm)
        End Do ! LM
      End If

      Call gradrl(nspin, mesh, l1max, dx, rhol, r, drdi, ipan1, ipand, ircut, &
        drrl, ddrrl, drrul, ddrrul, irmd, lmpot)

!----------------------------------------------------------------------------
! Loop over radial mesh
!----------------------------------------------------------------------------

      Do ir = irc0, irc1
        rpoint = r(ir)
!-------------------------------------------------------------------------
! Calculate the ex.-cor. potential
!-------------------------------------------------------------------------
        nspin2 = 2

        Do ispin = 1, nspin2
          Do lm = 1, lmmax
            rholm(lm, ispin) = rhol(ir, ispin, lm)
          End Do
        End Do
!    only for spin-polarized
!
! PW91 functional
        If (kxc==3) Then
          Call mkxcpe(nspin2, ir, ijend, l1max, rpoint, rholm, vxc, excij, &
            thet, ylm, dylmt1, dylmt2, dylmf1, dylmf2, dylmtf, drrl, ddrrl, &
            drrul, ddrrul, irmd, lmpot)
! PBE functional
        Else If (kxc==4) Then
          Call mkxcpe2(ir, ijend, rpoint, rholm, vxc, excij, ylm, dylmt1, &
            dylmf1, dylmf2, dylmtf, drrl, ddrrl, drrul, ddrrul, irmd, lmpot, &
            lmmax, .False.)
! PBEsol functional
        Else If (kxc==5) Then
          Call mkxcpe2(ir, ijend, rpoint, rholm, vxc, excij, ylm, dylmt1, &
            dylmf1, dylmf2, dylmtf, drrl, ddrrl, drrul, ddrrul, irmd, lmpot, &
            lmmax, .True.)
        Else
          Write (1337, *) ' KXC ???'
          Stop
        End If
!-------------------------------------------------------------------------
! Expand the ex.-cor. potential into spherical harmonics ,
!   using the orthogonality
!-------------------------------------------------------------------------
        Do ispin = 1, nspin
!----------------------------------------------------------------------
! Determine the corresponding potential number
!----------------------------------------------------------------------
          ipot = ispin
          Do lm = 1, lmmax
            vlmxc = ddot(ijend, vxc(1,ispin), 1, wtyr(1,lm), 1)
            v(ir, lm, ipot) = v(ir, lm, ipot) + vlmxc
!-------------------------------------------------------------------
! Store the ex.-c. potential of ir=2 and =3 for the extrapolation
!-------------------------------------------------------------------
            If (lm==1 .And. (ir==2 .Or. ir==3)) vxcr(ir, ispin) = vlmxc
          End Do ! LM
        End Do ! ISPIN
!-------------------------------------------------------------------------
! File er in case of total energies
!-------------------------------------------------------------------------
        If (kte==1) Then
!----------------------------------------------------------------------
! Expand ex.-cor. energy into spherical harmonics
!   using the orthogonality
!----------------------------------------------------------------------
          Do l = 0, lmax
            Do m = -l, l
              lm = l*l + l + m + 1
              elmxc = ddot(ijend, excij, 1, wtyr(1,lm), 1)
!----------------------------------------------------------------
! Multiply the lm-component of the ex.-cor. energy with the same
! lm-component of the charge density times r**2 and sum over lm
! this corresponds to a integration over the angular .
!----------------------------------------------------------------
              If ((kshape/=0) .And. (ir>irs1)) Then
                estor(ir, lm) = elmxc
              Else
                er(ir, l) = er(ir, l) + rho2ns(ir, lm, 1)*elmxc
              End If
            End Do ! M
          End Do ! L
        End If
      End Do !IR
!----------------------------------------------------------------------------
! Integrate er in case of total energies to get exc
!----------------------------------------------------------------------------
      If (kte==1) Then
        If (kshape==0) Then
          Do l = 0, lmax
            Call simp3(er(1,l), exc(l,iatyp), 1, irs1, drdi)
          End Do
        Else
          Do l = 0, lmax
            Do m = -l, l
              lm = l*l + l + m + 1
!----------------------------------------------------------------
! Convolute with shape function
!----------------------------------------------------------------
              Do j = imaxsh(lm-1) + 1, imaxsh(lm)
                lm2 = ilm_map(j, 2)
                If (lmsp(ilm_map(j,3))>0) Then
                  ifun = ifunm(ilm_map(j,3))
                  Do ir = irs1 + 1, irc1
                    irh = ir - irs1
                    er(ir, l) = er(ir, l) + rho2ns(ir, lm, 1)*gsh(j)*thetas( &
                      irh, ifun)*estor(ir, lm2)
                  End Do ! IR
                End If
              End Do ! J
            End Do ! M
            Call simpk(er(1,l), exc(l,iatyp), ipan1, ircut, drdi)
          End Do ! L
        End If
      End If
!----------------------------------------------------------------------------
! Extrapolate ex.-cor potential to the origin only for lm=1
!----------------------------------------------------------------------------
      Do ispin = 1, nspin
        ipot = ispin
!
        vxc2 = vxcr(2, ispin)
        vxc3 = vxcr(3, ispin)
        vxc1 = vxc2 - r(2)*(vxc3-vxc2)/(r(3)-r(2))
!
        v(1, 1, ipot) = v(1, 1, ipot) + vxc1
      End Do
!
    End Subroutine
