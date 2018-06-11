!-------------------------------------------------------------------------------
!> @brief add core and valence density expanded in spherical harmonics
!>         ( convention see subroutine rholm )
!> @details In the paramagnetic case (nspin=1) the core valence charge times
!> r**2 is added to the valence charge density times r**2
!> then only rho2ns(irmd,lmxtsq,natypd,1) is used .
!> In the spin-polarized case (nspin=2) the spin-splitted core
!> charge density times r**2 is converted into core charge
!> density times r**2 and core spin density times r**2 .
!> then these parts are added to corresponding parts of
!> the valence densities times r**2 , that are rho2ns(...,1)
!> which contains the charge density  and rho2ns(...,2) which
!> contains in that case the spin density .
!> (see notes by b.drittler)
!>
!> Attention : the core density is spherically averaged and multiplied by 4 pi.
!> therefore the core density is only added to l=0 part .
!
!> @author B. Drittler
!> @date   Nov. 1989
!
!> @note -V. Popescu March 2002: Total orbital moment within the WS sphere is also calculated
!> in the relativistic case; orbital density is normalised in the
!> same way as the charge density.
!> @note -Jonathan Chico Apr. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
    Subroutine rhototb(ipf, natyp, naez, nspin, rho2ns, rhoc, rhoorb, z, drdi, &
      irws, ircut, lpot, nfu, llmsp, thetas, ntcell, kshape, ipan, chrgnt, &
      itc, nshell, noq, conc, kaoez, catom, irm, nemb, lmpot)

      Use global_variables
      Use mod_datatypes, Only: dp

      Implicit None

!     .. Parameters ..
! .. Input variables
      Integer, Intent (In) :: itc
      Integer, Intent (In) :: ipf
      Integer, Intent (In) :: irm !< Maximum number of radial points
      Integer, Intent (In) :: nemb !< Number of 'embedding' positions
      Integer, Intent (In) :: lpot !< Maximum l component in potential expansion
      Integer, Intent (In) :: naez !< Number of atoms in unit cell
      Integer, Intent (In) :: natyp !< Number of kinds of atoms in unit cell
      Integer, Intent (In) :: nspin !< Counter for spin directions
      Integer, Intent (In) :: lmpot !< (LPOT+1)**2
      Integer, Intent (In) :: kshape !< Exact treatment of WS cell
!     .. Array Arguments ..
      Integer, Dimension (naez), Intent (In) :: noq !< Number of diff. atom types located
      Integer, Dimension (*), Intent (In) :: nfu !< number of shape function components in cell 'icell'
      Integer, Dimension (*), Intent (In) :: ipan !< Number of panels in non-MT-region
      Integer, Dimension (*), Intent (In) :: irws !< R point at WS radius
      Integer, Dimension (0:nsheld), Intent (In) :: nshell !< Index of atoms/pairs per shell (ij-pairs); nshell(0) = number of shells
      Integer, Dimension (*), Intent (In) :: ntcell !< Index for WS cell

      Integer, Dimension (0:ipand, *), Intent (In) :: ircut !< R points of panel borders
      Integer, Dimension (natyp, *), Intent (In) :: llmsp !< lm=(l,m) of 'nfund'th nonvanishing component of non-spherical pot.
      Integer, Dimension (natyp, naez+nemb), Intent (In) :: kaoez !< Kind of atom at site in elem. cell
      Real (Kind=dp), Dimension (*), Intent (In) :: z
      Real (Kind=dp), Dimension (natyp), Intent (In) :: conc !< Concentration of a given atom
      Real (Kind=dp), Dimension (irm, *), Intent (In) :: drdi !< Derivative dr/di
      Real (Kind=dp), Dimension (irm, *), Intent (In) :: rhoc !< core charge density
      Real (Kind=dp), Dimension (irm*krel+(1-krel), natyp), &
        Intent (In) :: rhoorb !< Orbital density
      Real (Kind=dp), Dimension (irid, nfund, *), Intent (In) :: thetas !< shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion

! .. In/Out variables
      Real (Kind=dp), Dimension (irm, lmpot, natyp, *), &
        Intent (Inout) :: rho2ns
! .. Output variables
      Real (Kind=dp), Intent (Out) :: chrgnt
      Real (Kind=dp), Dimension (natyp, 2*krel+(1-krel)*nspin), &
        Intent (Out) :: catom
! .. Local variables
      Integer :: i, i1, iatyp, icell, ifun, ipan1, ipotd, ipotu, irc1, irs1, &
        ispin, lm, iqez, ioez
      Real (Kind=dp) :: diff, factor, rfpi, sum, totsmom, totomom, sumo
      Real (Kind=dp), Dimension (natyp) :: omom !< Orbital moment
      Real (Kind=dp), Dimension (irm) :: rho
      Real (Kind=dp), Dimension (naez, 2*krel+(1-krel)*nspin) :: csite
      Real (Kind=dp), Dimension (krel*naez+(1-krel)) :: muosite
!
      Logical :: opt
! .. External Subroutines
      External :: simp3, simpk, opt
! .. Intrinsic Functions
      Intrinsic :: atan, sqrt
! .. Save statement
      Save
!     ..
      rfpi = sqrt(16.0E0_dp*atan(1.0E0_dp))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Loop over atomic sites
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Do iqez = 1, naez
!-------------------------------------------------------------------------
! Loop over atoms located on IQEZ
!-------------------------------------------------------------------------
        Do ispin = 1, nspin
          csite(iqez, ispin) = 0.0E0_dp
        End Do

        If (krel==1) muosite(iqez) = 0.0E0_dp

        Do ioez = 1, noq(iqez)

          iatyp = kaoez(ioez, iqez)
!---------------------------------------------------------------------
! Determine the right potential numbers for rhoc
!---------------------------------------------------------------------
          If (nspin==2) Then
            ipotd = 2*iatyp - 1
            ipotu = 2*iatyp
            factor = 1.0E0_dp
          Else
            ipotd = iatyp
            ipotu = iatyp
            factor = 0.5E0_dp
          End If

          If (kshape/=0) Then
            ipan1 = ipan(iatyp)
            irs1 = ircut(1, iatyp)
            irc1 = ircut(ipan1, iatyp)
          Else
            irs1 = irws(iatyp)
            irc1 = irs1
          End If
!-----------------------------------------------------------------------
          Do i = 2, irs1
!-------------------------------------------------------------------
! Convert core density
!-------------------------------------------------------------------
            sum = (rhoc(i,ipotd)+rhoc(i,ipotu))*factor/rfpi
            diff = (rhoc(i,ipotu)-rhoc(i,ipotd))/rfpi
!-------------------------------------------------------------------
! Add this to the lm=1 component of rho2ns
!-------------------------------------------------------------------
            rho2ns(i, 1, iatyp, 1) = rho2ns(i, 1, iatyp, 1) + sum
            rho2ns(i, 1, iatyp, nspin) = rho2ns(i, 1, iatyp, nspin) + diff
          End Do
!----------------------------------------------------------------------
! Calculate  charge and moment of the atom
!----------------------------------------------------------------------
          Do ispin = 1, nspin
!
            If (kshape==0) Then
!----------------------------------------------------------------
! Integrate over wigner seitz sphere - no shape correction
!----------------------------------------------------------------
              Call simp3(rho2ns(1,1,iatyp,ispin), sum, 1, irs1, drdi(1,iatyp))
!----------------------------------------------------------------
! The result has to be multiplied by sqrt(4 pi)
! (4 pi for integration over angle and 1/sqrt(4 pi) for
! the spherical harmonic y(l=0))
!----------------------------------------------------------------
              sum = sum*rfpi
            Else ! (KSHAPE.EQ.0)
!----------------------------------------------------------------
! convolute charge density with shape function to get the
! charge in the exact cell - if kshape .gt. 0
!----------------------------------------------------------------
              icell = ntcell(iatyp)

              Do i = 1, irs1
                rho(i) = rho2ns(i, 1, iatyp, ispin)*rfpi
              End Do
!
              Do i = irs1 + 1, irc1
                rho(i) = 0.0E0_dp
              End Do
!
              Do ifun = 1, nfu(icell)
                lm = llmsp(icell, ifun)
                If (lm<=lmpot) Then
                  Do i = irs1 + 1, irc1
                    rho(i) = rho(i) + rho2ns(i, lm, iatyp, ispin)*thetas(i- &
                      irs1, ifun, icell)
                  End Do
                End If
              End Do
!----------------------------------------------------------------
! Integrate over circumscribed sphere
!----------------------------------------------------------------
              Call simpk(rho, sum, ipan1, ircut(0,iatyp), drdi(1,iatyp))
            End If ! (KSHAPE.EQ.0)

            catom(iatyp, ispin) = sum
            csite(iqez, ispin) = csite(iqez, ispin) + &
              catom(iatyp, ispin)*conc(iatyp)

            If (ispin/=1) Then
!----------------------------------------------------------------
! Calculate orbital moment (ASA) and add it to the total
!----------------------------------------------------------------
              If ((krel==1) .And. (kshape==0)) Then
                Call simp3(rhoorb(1,iatyp), sumo, 1, irs1, drdi(1,iatyp))
                sumo = sumo*rfpi
                omom(iatyp) = sumo
                muosite(iqez) = muosite(iqez) + omom(iatyp)*conc(iatyp)
              End If

              If (kshape/=0) Then
                Write (ipf, Fmt=110) sum
              Else
                Write (ipf, Fmt=130) sum
                If (krel==1) Then
                  Write (ipf, Fmt=140) omom(iatyp)
                  Write (ipf, Fmt=150) sum + omom(iatyp)
                End If
              End If
            Else ! (ISPIN.NE.1)
              If (kshape/=0) Then
                Write (ipf, Fmt=100) iatyp, sum
              Else
                Write (ipf, Fmt=120) iatyp, sum
              End If
            End If ! (ISPIN.NE.1)
          End Do ! ISPIN = 1,NSPIN
!----------------------------------------------------------------------
          If (ioez/=noq(iqez)) Write (ipf, '(2X,77("-"))')
        End Do
!-------------------------------------------------------------------------
! IOEZ = 1, NOQ(IQEZ)
!-------------------------------------------------------------------------
        If (noq(iqez)>1) Then
          Write (ipf, '(2X,77("="))')
          Write (ipf, Fmt=200) iqez, csite(iqez, 1)
          If (nspin==2) Then
            Write (ipf, Fmt=210) csite(iqez, nspin)
            If (krel==1) Then
              Write (ipf, Fmt=220) muosite(iqez)
              Write (ipf, Fmt=230) csite(iqez, nspin) + muosite(iqez)
            End If
          End If
          If (iqez/=naez) Write (ipf, '(2X,77("="))')
        Else
          If (iqez/=naez) Write (ipf, '(2X,77("="))')
        End If
      End Do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! IQEZ = 1, NAEZ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Write (ipf, *)

      chrgnt = 0.0E0_dp
      Do i1 = 1, natyp
        chrgnt = chrgnt + real(nshell(i1), kind=dp)*(catom(i1,1)-z(i1))*conc( &
          i1)
      End Do

      Write (ipf, '(79("+"))')
      Write (ipf, Fmt=160) itc, chrgnt
      Write (6, Fmt=160) itc, chrgnt

      If (nspin==2) Then
        totsmom = 0.0E0_dp
        If (krel==1) totomom = 0.0E0_dp
        Do i1 = 1, natyp
          totsmom = totsmom + real(nshell(i1), kind=dp)*catom(i1, nspin)*conc( &
            i1)
          If (krel==1) totomom = totomom + real(nshell(i1), kind=dp)*omom(i1)* &
            conc(i1)
        End Do

        If (krel==0) Then
          Write (ipf, Fmt=170) totsmom
          Write (6, Fmt=170) totsmom
        Else
          Write (ipf, Fmt=170) totsmom + totomom
          Write (ipf, Fmt=180) totsmom
          Write (ipf, Fmt=190) totomom
          Write (6, Fmt=170) totsmom + totomom
          Write (6, Fmt=180) totsmom
          Write (6, Fmt=190) totomom
        End If
      End If
      Write (ipf, *)

      Return

100   Format ('  Atom ', I4, ' charge in wigner seitz cell =', F10.6)
110   Format (7X, 'spin moment in wigner seitz cell =', F10.6)
120   Format ('  Atom ', I4, ' charge in wigner seitz sphere =', F10.6)
130   Format (7X, 'spin moment in wigner seitz sphere =', F10.6)
140   Format (7X, 'orb. moment in wigner seitz sphere =', F10.6)
150   Format (7X, 'total magnetic moment in WS sphere =', F10.6)
160   Format ('      ITERATION', I4, ' charge neutrality in unit cell = ', &
        F12.6)
170   Format ('                   ', ' TOTAL mag. moment in unit cell = ', &
        F12.6)
180   Format ('                   ', '           spin magnetic moment = ', &
        F12.6)
190   Format ('                   ', '        orbital magnetic moment = ', &
        F12.6)
200   Format ('      Site ', I3, ' total charge =', F10.6)
210   Format ('         ', ' total spin moment =', F10.6)
220   Format ('         ', ' total orb. moment =', F10.6)
230   Format ('      total magnetic moment =', F10.6)

    End Subroutine
