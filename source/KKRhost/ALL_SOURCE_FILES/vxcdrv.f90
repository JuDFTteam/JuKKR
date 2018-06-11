!-------------------------------------------------------------------------------
! SUBROUTINE: VXCDRV
!> @brief Wrapper for the calculation of the exchange-correlation energy, for
!> the different treatments of the xc-potential.
!> @note
!> - Jonathan Chico Jan. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
    Subroutine vxcdrv(exc, kte, kxc, lpot, nspin, nstart, nend, rho2ns, vons, &
      r, drdi, a, irws, ircut, ipan, ntcell, kshape, gsh, ilm_map, imaxsh, &
      ifunm, thetas, lmsp, npotd, lmpot, lmxspd, irm, natyp, lmmaxd)

      Use global_variables
      Use mod_datatypes, Only: dp

      Implicit None

! .. Input variables
      Integer, Intent (In) :: irm !< Maximum number of radial points
      Integer, Intent (In) :: kte !< Calculation of the total energy On/Off (1/0)
      Integer, Intent (In) :: kxc !< Type of xc-potential 0=vBH 1=MJW 2=VWN 3=PW91
      Integer, Intent (In) :: nend
      Integer, Intent (In) :: lpot !< Maximum l component in potential expansion
      Integer, Intent (In) :: natyp !< Number of kinds of atoms in unit cell
      Integer, Intent (In) :: npotd !< (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP)
      Integer, Intent (In) :: lmpot !< (LPOT+1)**2
      Integer, Intent (In) :: nspin !< Counter for spin directions
      Integer, Intent (In) :: lmmaxd !< (KREL+KORBIT+1)(LMAX+1)^2
      Integer, Intent (In) :: nstart
      Integer, Intent (In) :: kshape !< Exact treatment of WS cell
      Integer, Intent (In) :: lmxspd !< (2*LPOT+1)**2
! .. Array Arguments ..
      Real (Kind=dp), Dimension (natyp), Intent (In) :: a !< Constants for exponential R mesh
      Real (Kind=dp), Dimension (ngshd), Intent (In) :: gsh
      Real (Kind=dp), Dimension (irm, natyp), Intent (In) :: r !< Radial mesh ( in units a Bohr)
      Real (Kind=dp), Dimension (irm, natyp), Intent (In) :: drdi !< Derivative dr/di
      Real (Kind=dp), Dimension (irid, nfund, ncelld), Intent (In) :: thetas !< shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
      Real (Kind=dp), Dimension (irm, lmpot, natyp, 2), Intent (In) :: rho2ns !< radial density
      Integer, Dimension (natyp), Intent (In) :: irws !< R point at WS radius
      Integer, Dimension (natyp), Intent (In) :: ipan !< Number of panels in non-MT-region
      Integer, Dimension (natyp), Intent (In) :: ntcell !< index for WS cell
      Integer, Dimension (0:lmpot), Intent (In) :: imaxsh
      Integer, Dimension (ngshd, 3), Intent (In) :: ilm_map
      Integer, Dimension (natyp, lmxspd), Intent (In) :: lmsp !< 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
      Integer, Dimension (natyp, lmxspd), Intent (In) :: ifunm
      Integer, Dimension (0:ipand, natyp), Intent (In) :: ircut !< r points of panel borders
! .. In/Out variables
      Real (Kind=dp), Dimension (0:lpot, natyp), Intent (Inout) :: exc !< exchange correlation energy
! .. Output variables
      Real (Kind=dp), Dimension (irm, lmpot, npotd), Intent (Out) :: vons !< output potential (nonspherical VONS)
! ..
! .. External Subroutines
      External :: dcopy, sphere_gga, sphere_nogga, vxcgga, vxclm
! .. Local Scalars
      Integer :: iatyp, icell, ipot, lmx1
      Integer :: ijd
! .. Parameters
      Parameter (ijd=434)
! .. Local Arrays ..
      Integer, Dimension (lmxspd) :: lmspiat
      Integer, Dimension (lmxspd) :: ifunmiat
      Real (Kind=dp), Dimension (ijd) :: thet
      Real (Kind=dp), Dimension (ijd, lmpot) :: yr
      Real (Kind=dp), Dimension (ijd, 3) :: rij
      Real (Kind=dp), Dimension (ijd, lmpot) :: ylm
      Real (Kind=dp), Dimension (ijd, lmpot) :: wtyr
      Real (Kind=dp), Dimension (ijd, lmpot) :: dylmf1
      Real (Kind=dp), Dimension (ijd, lmpot) :: dylmf2
      Real (Kind=dp), Dimension (ijd, lmpot) :: dylmt1
      Real (Kind=dp), Dimension (ijd, lmpot) :: dylmt2
      Real (Kind=dp), Dimension (ijd, lmpot) :: dylmtf
      Real (Kind=dp), Dimension (irm, lmpot, 2) :: rho2iat
!     ..
      If (kxc<3) Then
        Call sphere_nogga(lpot, yr, wtyr, rij, ijd)
      Else
        Call sphere_gga(lpot, yr, wtyr, rij, ijd, lmpot, thet, ylm, dylmt1, &
          dylmt2, dylmf1, dylmf2, dylmtf)
      End If
      Do iatyp = nstart, nend
        icell = ntcell(iatyp)
        ipot = nspin*(iatyp-1) + 1
        Do lmx1 = 1, lmxspd
          ifunmiat(lmx1) = ifunm(icell, lmx1)
          lmspiat(lmx1) = lmsp(icell, lmx1)
        End Do
        Call dcopy(irm*lmpot, rho2ns(1,1,iatyp,1), 1, rho2iat(1,1,1), 1)
        If (nspin==2 .Or. krel==1) Then
          Call dcopy(irm*lmpot, rho2ns(1,1,iatyp,2), 1, rho2iat(1,1,2), 1)
        End If
        If (kxc<3) Then
          Call vxclm(exc, kte, kxc, lpot, nspin, iatyp, rho2iat, &
            vons(1,1,ipot), r(1,iatyp), drdi(1,iatyp), irws(iatyp), &
            ircut(0,iatyp), ipan(iatyp), kshape, gsh, ilm_map, imaxsh, &
            ifunmiat, thetas(1,1,icell), yr, wtyr, ijd, lmspiat, lmpot, &
            lmxspd, lmmaxd, irm, lpot, natyp)
        Else
!----------------------------------------------------------------------
! GGA EX-COR POTENTIAL
!----------------------------------------------------------------------
          Call vxcgga(exc, kte, kxc, lpot, nspin, iatyp, rho2iat, &
            vons(1,1,ipot), r(1,iatyp), drdi(1,iatyp), a(iatyp), irws(iatyp), &
            ircut(0,iatyp), ipan(iatyp), kshape, gsh, ilm_map, imaxsh, &
            ifunmiat, thetas(1,1,icell), wtyr, ijd, lmspiat, thet, ylm, &
            dylmt1, dylmt2, dylmf1, dylmf2, dylmtf, lmpot, lmxspd, lmmaxd, &
            irm, lpot, natyp)
        End If
      End Do
    End Subroutine
