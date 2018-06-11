!-------------------------------------------------------------------------------
! SUBROUTINE: VINTRAS
!> @brief Calculate the electron-intracell-potentials and the charge-moments
!> of given charge densities. ( For each spin-direction the potential is the
!> same in the polarized case.)
!> @details Initialize the potential \f$ V\f$ with the electron-intracell-potentials
!> the intracell-potential is expanded into spherical harmonics .
!> the lm-term of the intracell-potential of the representive atom i is given by
!>
!> \f$V\left(r,lm,i\right)=\frac{8\pi}{2l+1} \left(\int_{0}^{r}dr'\frac{r'^{l}}{r^{l+1}} rho2ns(r',lm,i,1) + \int_{r}^{r_{cut}}dr' \frac{r^l}{r'^{l+1}}rho2ns(r',lm,i,1) ) \right)\f$
!>
!> the lm contribution of the charge moment of the representive atom i is given by
!>
!> \f$ cmom\left(lm,i\right)=\int_{0}^{r_{cut}} dr' r'^{l}rho2ns(r',lm,i,1)\f$
!>
!>  (see notes by b.drittler and u.klemradt) \f$ r_{cut}\f$ is muffin tin or
!> Wigner-Seitz sphere radius, depending on kshape turned on or off
!
!> @note Attention : \f$ rho2ns(...,1)\f$ is the real charge density times \f$r^2\f$
!> developed into spherical harmonics . (see deck rholm)
!>
!> - Jonathan Chico Jan. 2018: Removed inc.p dependencies and rewrote to Fortran90
!> @author B. Drittler
!> @date May 1987
!-----------------------------------------------------------------------

    Subroutine vintras(cmom, cminst, lmax, nspin, nstart, nend, rho2ns, v, r, &
      drdi, irws, ircut, ipan, kshape, ntcell, ilm_map, ifunm, imaxsh, gsh, &
      thetas, lmsp, lmpot, natyp)

      Use constants
      Use global_variables
      Use mod_datatypes, Only: dp

      Implicit None

! .. Input Variables
      Integer, Intent (In) :: lmax !< Maximum l component in wave function expansion
      Integer, Intent (In) :: nend
      Integer, Intent (In) :: nspin !< Counter for spin directions
      Integer, Intent (In) :: lmpot !< (LPOT+1)**2
      Integer, Intent (In) :: natyp !< Number of kinds of atoms in unit cell
      Integer, Intent (In) :: nstart
      Integer, Intent (In) :: kshape !< Exact treatment of WS cell
      Integer, Dimension (natyp), Intent (In) :: irws !< R point at WS radius
      Integer, Dimension (natyp), Intent (In) :: ipan !< Number of panels in non-MT-region
      Integer, Dimension (natyp), Intent (In) :: ntcell !< Index for WS cell
      Integer, Dimension (0:lmpot), Intent (In) :: imaxsh
      Integer, Dimension (ngshd, 3), Intent (In) :: ilm_map
      Integer, Dimension (natyp, lmxspd), Intent (In) :: lmsp !< 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
      Integer, Dimension (natyp, lmxspd), Intent (In) :: ifunm
      Integer, Dimension (0:ipand, natyp), Intent (In) :: ircut !< R points of panel borders
      Real (Kind=dp), Dimension (ngshd), Intent (In) :: gsh
      Real (Kind=dp), Dimension (irmd, natyp), Intent (In) :: r !< Radial mesh ( in units a Bohr)
      Real (Kind=dp), Dimension (irmd, natyp), Intent (In) :: drdi !< Derivative dr/di
      Real (Kind=dp), Dimension (irid, nfund, ncelld), Intent (In) :: thetas !< shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
      Real (Kind=dp), Dimension (irmd, lmpot, natyp, 2), Intent (In) :: rho2ns !< radial density
! .. Output variables
      Real (Kind=dp), Dimension (lmpot, natyp), Intent (Out) :: cmom !< LM moment of total charge
      Real (Kind=dp), Dimension (lmpot, natyp), Intent (Out) :: cminst
      Real (Kind=dp), Dimension (irmd, lmpot, npotd), Intent (Out) :: v
! .. Local Variables
      Real (Kind=dp) :: fac, rl
      Integer :: i, iatyp, icell, iend, ifun, ipot, irc1, irs1, istart, j, l, &
        lm, lm2, lm3, m
! .. Local Arrays
      Integer, Dimension (0:ipand) :: ircutm
      Real (Kind=dp), Dimension (irmd) :: v1
      Real (Kind=dp), Dimension (irmd) :: v2
      Real (Kind=dp), Dimension (irmd) :: vint1
      Real (Kind=dp), Dimension (irmd) :: vint2
! .. External Subroutines ..
      External :: sinwk, soutk
! .. Intrinsic Functions ..
      Intrinsic :: atan, real
!     ..
      Do iatyp = nstart, nend
        If (kshape/=0) Then
          irs1 = ircut(1, iatyp)
          irc1 = ircut(ipan(iatyp), iatyp)
          icell = ntcell(iatyp)
          Do i = 0, ipan(iatyp)
            ircutm(i) = ircut(i, iatyp)
          End Do ! I
        Else
          irs1 = irws(iatyp)
          irc1 = irs1
          ircutm(0) = ircut(0, iatyp)
          ircutm(1) = irc1
        End If
!-------------------------------------------------------------------------
! Determine the right potential numbers
!-------------------------------------------------------------------------
        ipot = nspin*iatyp

        Do l = 0, lmax
          fac = 8.0E0_dp*pi/real(2*l+1, kind=dp)
          Do m = -l, l
            lm = l*l + l + m + 1
!-------------------------------------------------------------------
! Set up of the integrands v1 and v2
!-------------------------------------------------------------------
            v1(1) = 0.0E0_dp
            v2(1) = 0.0E0_dp
            Do i = 2, irs1
              rl = r(i, iatyp)**l
              v1(i) = rho2ns(i, lm, iatyp, 1)*rl*drdi(i, iatyp)
              v2(i) = rho2ns(i, lm, iatyp, 1)/r(i, iatyp)/rl*drdi(i, iatyp)
            End Do ! I
!-------------------------------------------------------------------
! Convolute charge density of interstial with shape function if kshape.gt.0
!-------------------------------------------------------------------
            If (kshape/=0) Then
              Do i = irs1 + 1, irc1
                v1(i) = 0.0E0_dp
              End Do ! I
              istart = imaxsh(lm-1) + 1
              iend = imaxsh(lm)
              Do j = istart, iend
                lm2 = ilm_map(j, 2)
                lm3 = ilm_map(j, 3)
                If (lmsp(icell,lm3)>0) Then
                  ifun = ifunm(icell, lm3)
                  Do i = irs1 + 1, irc1
                    v1(i) = v1(i) + gsh(j)*rho2ns(i, lm2, iatyp, 1)*thetas(i- &
                      irs1, ifun, icell)
                  End Do ! I
                End If
              End Do ! J

              Do i = irs1 + 1, irc1
                rl = r(i, iatyp)**l
                v2(i) = v1(i)/r(i, iatyp)/rl*drdi(i, iatyp)
                v1(i) = v1(i)*rl*drdi(i, iatyp)
              End Do ! I
            End If
!-------------------------------------------------------------------
! Now integrate v1 and v2
!-------------------------------------------------------------------
            Call soutk(v1, vint1, ipan(iatyp), ircutm)
            Call sinwk(v2, vint2, ipan(iatyp), ircutm)
!-------------------------------------------------------------------
! Gather all parts
!-------------------------------------------------------------------
            If (lm==1) Then
              v(1, lm, ipot) = fac*vint2(1)
            Else
              v(1, lm, ipot) = 0.0E0_dp
            End If

            Do i = 2, irc1
              rl = r(i, iatyp)**l
              v(i, lm, ipot) = fac*(vint1(i)/r(i,iatyp)/rl+vint2(i)*rl)
            End Do ! I
!-------------------------------------------------------------------
! Store charge moment - in case of kshape.gt.0 this is the moment
!      of the charge in the muffin tin sphere
!-------------------------------------------------------------------
            cmom(lm, iatyp) = vint1(irs1)
!-------------------------------------------------------------------
! Store charge moment of interstial in case of kshape.gt.0
!-------------------------------------------------------------------
            If (kshape/=0) cminst(lm, iatyp) = vint1(irc1) - vint1(irs1)
!
            If (nspin==2) Then
              Do i = 1, irc1
                v(i, lm, ipot-1) = v(i, lm, ipot)
              End Do ! I
            End If
          End Do ! M
        End Do ! L
      End Do ! IATYP
      Return

    End Subroutine
