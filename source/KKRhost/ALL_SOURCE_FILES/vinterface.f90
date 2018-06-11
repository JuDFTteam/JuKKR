!-------------------------------------------------------------------------------
! SUBROUTINE: VINTERFACE
!> @brief This is calculating the intra-atomic contibution of the potential in
!>  the case of an interface taking into account the bulk potential on
!>  the two sides.
!
!> @details It uses the structure dependent matrices AVMAD which are calculated
!>  once in the subroutine MADELUNG2D() and saved in the DA-file
!>  avmad.unformatted ( May 2004)
!>
!>  For each site in a layer the summation in all other layers is split
!>  into three parts: within the slab, over the NLEFT*NLBASIS left host
!>  sites and over the NRIGHT*NRBASIS right host sites, the last two
!>  steps only in case of decimation run
!
!-------------------------------------------------------------------------------
!> @note
!> - Adapted for the case of more atoms on the same site, summation is
!>  done over the occupants of that site, the charge is weighted with
!>  the appropriate concentration of the occupant  V. Popescu feb. 2002
!-------------------------------------------------------------------------------
!>
!> - Impurity-program adopted feb. 2004 (according to N. Papanikalou)
!>
!> - Jonathan Chico Feb. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
    Subroutine vinterface(cmom, cminst, lpot, nspin, nlayers, natyp, v, zat, &
      r, irws, ircut, ipan, kshape, noq, kaoez, iqat, conc, catom, icc, &
      hostimp, nlbasis, nleft, nrbasis, nright, cmomhost, chrgnt, vinters, &
      naez, lmpot)

      Use constants
      Use global_variables
      Use mod_datatypes, Only: dp

      Implicit None

! .. Input variables ..
      Integer, Intent (In) :: icc !< Enables the calculation of off-diagonal elements of the GF.(0=SCF/DOS; 1=cluster; -1=custom)
      Integer, Intent (In) :: lpot !< Maximum l component in potential expansion
      Integer, Intent (In) :: naez !< Number of atoms in unit cell
      Integer, Intent (In) :: lmpot !< (LPOT+1)**2
      Integer, Intent (In) :: nspin !< Counter for spin directions
      Integer, Intent (In) :: natyp !< Number of kinds of atoms in unit cell
      Integer, Intent (In) :: nleft !< Number of repeated basis for left host to get converged electrostatic potentials
      Integer, Intent (In) :: nright !< Number of repeated basis for right host to get converged electrostatic potentials
      Integer, Intent (In) :: kshape !< Exact treatment of WS cell
      Integer, Intent (In) :: nlayers
      Integer, Intent (In) :: nlbasis !< Number of basis layers of left host (repeated units)
      Integer, Intent (In) :: nrbasis !< Number of basis layers of right host (repeated units)
      Real (Kind=dp), Intent (In) :: chrgnt
      Integer, Dimension (naez), Intent (In) :: noq !< Number of diff. atom types located
      Integer, Dimension (natyp), Intent (In) :: irws !< R point at WS radius
      Integer, Dimension (natyp), Intent (In) :: ipan !< Number of panels in non-MT-region
      Integer, Dimension (natyp), Intent (In) :: iqat !< The site on which an atom is located on a given site
      Integer, Dimension (0:natyp), Intent (In) :: hostimp
      Integer, Dimension (0:ipand, natyp), Intent (In) :: ircut !< R points of panel borders
      Integer, Dimension (natyp, naez+nembd1-1), Intent (In) :: kaoez !< Kind of atom at site in elem. cell
      Real (Kind=dp), Dimension (natyp), Intent (In) :: zat !< Nuclear charge
      Real (Kind=dp), Dimension (natyp), Intent (In) :: conc !< Concentration of a given atom
      Real (Kind=dp), Dimension (natyp), Intent (In) :: catom
      Real (Kind=dp), Dimension (irmd, natyp), Intent (In) :: r !< Radial mesh ( in units a Bohr)
      Real (Kind=dp), Dimension (lmpot, natyp), Intent (In) :: cmom !< LM moment of total charge
      Real (Kind=dp), Dimension (lmpot, natyp), Intent (In) :: cminst !< charge moment of interstitial
      Real (Kind=dp), Dimension (lmpot, nembd1), Intent (In) :: cmomhost !< Charge moments of each atom of the (left/right) host
! .. In/out variables
      Real (Kind=dp), Dimension (irmd, lmpot, npotd), Intent (Inout) :: v

! .. Local variables
      Integer :: ileft, iright
      Integer :: i, iatom, ib, ih1, ilay1, ilay2, io2
      Integer :: ipot, irs1, ispin, it1, it2, l, lm, lm2, m
      Integer :: lrecamad, irec, nleftoff, nrightoff, nleftall, nrightall
      Real (Kind=dp) :: cm1, fpi
      Logical :: opt, test, lread
      Real (Kind=dp), Dimension (lmpot) :: ac
      Real (Kind=dp), Dimension (lmpot) :: cm
      Real (Kind=dp), Dimension (2) :: charge
      Real (Kind=dp), Dimension (naez) :: monopol
      Real (Kind=dp), Dimension (lmpot, lmpot) :: avmad
      Real (Kind=dp), Dimension (lmpot, naez) :: vinters
! .. Intrinsic Functions ..
      Intrinsic :: atan, sqrt
! .. External Functions/Subroutines
      External :: opt, test

      If (test('flow    ')) Write (1337, *) '>>>>>> Vinterface'

      Inquire (File='avmad.unformatted', Exist=lread) ! ewald2d

      If (lread) Then
        lrecamad = wlength*2*lmpot*lmpot
        Open (69, Access='direct', Recl=lrecamad, File='avmad.unformatted', &
          Form='unformatted')
      Else
        lrecamad = wlength*2*lmpot*lmpot + wlength*2*lmpot
        Open (69, Access='direct', Recl=lrecamad, File='abvmad.unformatted', &
          Form='unformatted')
      End If

      Write (1337, Fmt=100)
      Write (1337, Fmt=110)

      fpi = 4.E0_dp*pi

      If (opt('DECIMATE')) Then
!-------------------------------------------------------------------------
! Setup the charges to put in the ghost layers in the case of
! decimation technique to achieve charge neutrality
!-------------------------------------------------------------------------
        charge(1) = -chrgnt/(2.E0_dp*sqrt(fpi))
        charge(2) = -chrgnt/(2.E0_dp*sqrt(fpi))
!
        nleftoff = nlayers*nlayers ! record offsets
        nrightoff = nleftoff + nlayers*nleft*nlbasis ! left and right
        nleftall = nleft*nlbasis
        nrightall = nright*nrbasis
      End If
!----------------------------------------------------------------------------
!                   START CALCULATION IN THE LAYERS
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! Loop over atoms in slab
!----------------------------------------------------------------------------
!
      Do it1 = 1, natyp
!-------------------------------------------------------------------------
! Take a site occupied by IT1
!-------------------------------------------------------------------------
        ilay1 = iqat(it1)
!
        If (kshape/=0) Then
          irs1 = ircut(ipan(it1), it1)
        Else
          irs1 = irws(it1)
        End If
!
        Do lm = 1, lmpot
          ac(lm) = 0.E0_dp
        End Do
!-------------------------------------------------------------------------
! 1.  Summation in all layers in the slab
!-------------------------------------------------------------------------
        Do ilay2 = 1, nlayers
          irec = ilay2 + nlayers*(ilay1-1)
          Read (69, Rec=irec) avmad

!----------------------------------------------------------------------
! Keep the monopole term -- Hoshino is doing (SMONOPOL(I) -SMONOPOL(0))
!----------------------------------------------------------------------
          If (ilay1==ilay2) monopol(ilay1) = avmad(1, 1)
!----------------------------------------------------------------------
! Loop over all occupants of site ILAY2
!----------------------------------------------------------------------
          Do io2 = 1, noq(ilay2)
            it2 = kaoez(io2, ilay2)
!
            Do lm = 1, lmpot
              cm(lm) = cmom(lm, it2)
!----------------------------------------------------------------
! Add contribution of interstial in case of shapes
!----------------------------------------------------------------
              If (kshape/=0) cm(lm) = cm(lm) + cminst(lm, it2)
            End Do
            cm(1) = cm(1) - zat(it2)/sqrt(fpi)
!
            Do lm = 1, lmpot
              Do lm2 = 1, lmpot
                ac(lm) = ac(lm) + avmad(lm, lm2)*cm(lm2)*conc(it2)
              End Do
            End Do
          End Do
!----------------------------------------------------------------------
! Loop over all occupants of site ILAY2
!----------------------------------------------------------------------
        End Do ! ILAY2 loop in all interface planes
!-------------------------------------------------------------------------
        Do ilay2 = 1, nlayers
!----------------------------------------------------------------------
! Loop over all occupants of site ILAY2
!----------------------------------------------------------------------
          Do io2 = 1, noq(ilay2)
            it2 = kaoez(io2, ilay2)
!
            cm1 = cmom(1, it2)
            If (kshape/=0) cm1 = cm1 + cminst(1, it2)
!
            cm1 = cm1 - zat(it2)/sqrt(fpi)
            ac(1) = ac(1) - monopol(ilay1)*cm1*conc(it2)
!
          End Do
!----------------------------------------------------------------------
! Loop over all occupants of site ILAY2
!----------------------------------------------------------------------
        End Do
!-------------------------------------------------------------------------
! Correction: charge neutrality is imposed (see P. Lang)
!-------------------------------------------------------------------------
        If (opt('DECIMATE')) Then
!----------------------------------------------------------------------
! 2.  Summation in the LEFT bulk side
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Loop over all occupants of LEFT host
!----------------------------------------------------------------------
          ileft = 0
          Do ih1 = 1, nleft
            Do ib = 1, nlbasis
              ileft = ileft + 1
              irec = ileft + nleftall*(ilay1-1) + nleftoff
              Read (69, Rec=irec) avmad
!
              iatom = ib
              Do lm = 1, lmpot
                Do lm2 = 1, lmpot
                  ac(lm) = ac(lm) + avmad(lm, lm2)*cmomhost(lm2, iatom)
                End Do
              End Do
!
              If ((ih1==1) .And. (ib==1)) Then
                ac(1) = ac(1) + (avmad(1,1)-monopol(ilay1))*charge(1)
              End If
            End Do
          End Do
!----------------------------------------------------------------------
          If (ileft/=nleftall) Then
            Write (6, *) ' < VINTERFACE > : index error ', &
              'ILEFT <> NLEFT*NLBASIS'
            Stop
          End If
!----------------------------------------------------------------------
! 3.  Summation in the RIGHT bulk side
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Loop over all occupants of RIGHT host
!----------------------------------------------------------------------
          iright = 0
          Do ih1 = 1, nright
            Do ib = 1, nrbasis
              iright = iright + 1
              irec = iright + nrightall*(ilay1-1) + nrightoff
              Read (69, Rec=irec) avmad
!
              iatom = nlbasis + ib
              Do lm = 1, lmpot
                Do lm2 = 1, lmpot
                  ac(lm) = ac(lm) + avmad(lm, lm2)*cmomhost(lm2, iatom)
                End Do
              End Do
!
              If ((ih1==1) .And. (ib==1)) Then
                ac(1) = ac(1) + (avmad(1,1)-monopol(ilay1))*charge(2)
              End If
            End Do
          End Do
!----------------------------------------------------------------------
          If (iright/=nrightall) Then
            Write (6, *) ' < VINTERFACE > : index error ', &
              'IRIGHT <> NRIGHT*NRBASIS'
            Stop
          End If
        End If ! (OPT(DECIMATE)
!-------------------------------------------------------------------------
        Write (1337, Fmt=120) it1, (catom(it1)-zat(it1)), &
          (ac(1)/sqrt(4.E0_dp*pi)), (ac(3)/sqrt(4.E0_dp*pi))
!-------------------------------------------------------------------------
! Loop over spins of atom IT1
!-------------------------------------------------------------------------
        Do ispin = 1, nspin
!----------------------------------------------------------------------
! Determine the right potential number
!----------------------------------------------------------------------
          ipot = nspin*(it1-1) + ispin
!----------------------------------------------------------------------
! In the case of l=0 : r(1)**l is not defined
!----------------------------------------------------------------------
          v(1, 1, ipot) = v(1, 1, ipot) + ac(1)
!
          Do l = 0, lpot
            Do m = -l, l
              lm = l*l + l + m + 1
              Do i = 2, irs1
                v(i, lm, ipot) = v(i, lm, ipot) + (-r(i,it1))**l*ac(lm)
              End Do
            End Do
          End Do
        End Do
!-------------------------------------------------------------------------
! This part (ICC.GT.0) should be presumably reconsidered for impurity
! calculation in host-CPA case
!-------------------------------------------------------------------------
        If (icc>0 .Or. opt('KKRFLEX ')) Then
          Do l = 0, lpot
            Do m = -l, l
              lm = l*l + l + m + 1
              vinters(lm, ilay1) = ac(lm)
            End Do
          End Do
        End If
!-------------------------------------------------------------------------
      End Do
!----------------------------------------------------------------------------
      Close (69)
      Write (1337, '(15X,45("-"),/)')
      Write (1337, '(79("="))')
      If ((icc==0) .And. (.Not. opt('KKRFLEX '))) Return
!----------------------------------------------------------------------------
! Now Prepare output for Impurity calculation
!----------------------------------------------------------------------------
      Open (91, File='intercell_ref', Status='unknown', Form='formatted')
      Write (1337, *)
      Write (1337, *) '                     ', &
        'Writing intercell potential for impurity'
      Write (1337, '(/,20X,55("-"))')
      Write (1337, 130) hostimp(0), lmpot
      Write (1337, '(20X,55("-"),/,35X,"  i host lm  Vint")')
      Do i = 1, hostimp(0)
        Write (1337, *)
        lm = 1
        Write (1337, '(35X,I4,I4,I3,1X,F10.6)') i, hostimp(i), lm, &
          vinters(lm, hostimp(i))
        Do lm = 2, 9
          Write (1337, '(43X,I3,1X,F10.6)') lm, vinters(lm, hostimp(i))
        End Do
        Write (1337, '(20X,55("-"))')
      End Do
      Write (1337, '(79("="),/)')

      Write (91, 140) hostimp(0), lmpot
      Do i = 1, hostimp(0)
        Write (91, 150)(vinters(lm,hostimp(i)), lm=1, lmpot)
      End Do
      Close (91)
!
      Return
!
100   Format (79('='), /, 25X, ' INTERFACE MADELUNG POTENTIALS ')
110   Format (/, 15X, ' ATOM ', '  Delta_Q  ', '   MONOPOLE       DIPOLE', /, &
        15X, 45('-'))
120   Format (15X, I4, 2X, F10.6, 1X, 1P, D13.6, 1X, 1P, D13.6)
130   Format (22X, I4, ' host atoms, LMPOT = ', I2, ' output up to LM = 9')
140   Format (3I6)
150   Format (4D20.10)
    End Subroutine
