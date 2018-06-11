!-------------------------------------------------------------------------------
! SUBROUTINE: RHOSYMM
!> @brief Symmetrize the charge densities and magnetic moments of
!> atoms which are magnetic 'antisymmetric'
!> (dependencies in IXIPOL(*))
!
!> @author P. Zahn
!> @date Aug. 1996
!> @note -Jonathan Chico Apr. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
    Subroutine rhosymm(lmpot, nspin, nstart, nend, rho2ns, ixipol, irws, &
      ircut, ipan, kshape, natyp, irm)

      Use global_variables
      Use mod_datatypes, Only: dp

      Implicit None
! .. Input variables

      Integer, Intent (In) :: irm !< Maximum number of radial points
      Integer, Intent (In) :: nend
      Integer, Intent (In) :: natyp !< Number of kinds of atoms in unit cell
      Integer, Intent (In) :: lmpot !< (LPOT+1)**2
      Integer, Intent (In) :: nspin !< Counter for spin directions
      Integer, Intent (In) :: kshape !< Exact treatment of WS cell
      Integer, Intent (In) :: nstart
      Integer, Dimension (*), Intent (In) :: ipan !< Number of panels in non-MT-region
      Integer, Dimension (*), Intent (In) :: irws !< R point at WS radius
      Integer, Dimension (*), Intent (In) :: ixipol !< Constraint of spin pol.
      Integer, Dimension (0:ipand, *), Intent (In) :: ircut !< R points of panel borders
! .. In/Out variables
      Real (Kind=dp), Dimension (irm, lmpot, natyp, *), &
        Intent (Inout) :: rho2ns !< radial density
! .. Local variables
      Integer :: i, iatyp, iatyp1, irc, irc1, lm
      Real (Kind=dp) :: fac
! .. Intrinsic Functions
      Intrinsic :: abs
!----------------------------------------------------------------------------
!
      Do iatyp = nstart, nend
!
        iatyp1 = abs(ixipol(iatyp))
!
        fac = 1.E0_dp
        If (ixipol(iatyp)<0) fac = -1.E0_dp
!
        If (iatyp1>=iatyp) Then
          Write (1337, *) 'Symmetrize atom ', iatyp, ' with ', iatyp1, '.'
          If (kshape/=0) Then
            irc = ircut(ipan(iatyp), iatyp)
            irc1 = ircut(ipan(iatyp1), iatyp1)
          Else
            irc = irws(iatyp)
            irc1 = irws(iatyp1)
          End If
!
          If (irc/=irc1) Then
            Write (6, *) 'Error in RHOSYMM : ***********************'
            Write (6, *) 'Radial mesh of atoms ', iatyp, ' and ', iatyp1, &
              ' are not equal.'
          End If
!
          Do lm = 1, lmpot
            Do i = 1, irc1
              rho2ns(i, lm, iatyp, 1) = (rho2ns(i,lm,iatyp,1)+rho2ns(i,lm, &
                iatyp1,1))/2.E0_dp
              rho2ns(i, lm, iatyp1, 1) = rho2ns(i, lm, iatyp, 1)
              If (nspin>1) Then
                rho2ns(i, lm, iatyp, 2) = (rho2ns(i,lm,iatyp,2)+fac*rho2ns(i, &
                  lm,iatyp1,2))/2.E0_dp
                rho2ns(i, lm, iatyp1, 2) = fac*rho2ns(i, lm, iatyp, 2)
              End If
            End Do ! I =1,IRC1
          End Do ! LM =1,LMPOT
        End If ! (IATYP1.GT.IATYP)
      End Do ! IATYP=NSTART,NEND

      Return

    End Subroutine
