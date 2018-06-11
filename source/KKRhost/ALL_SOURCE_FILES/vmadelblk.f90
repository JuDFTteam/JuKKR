!-------------------------------------------------------------------------------
! SUBROUTINE: VMADELBLK
!> @brief Calculate the madelung potentials and add these to the potential \f$V\f$
!> (in he spin-polarized case for each spin-direction this is the same)
!> @details It uses the structure dependent matrices AVMAD and BVMAD which
!> are calculated once in the subroutine MADELUNG3D() and saved in
!> the DA-file abvmad.unformatted (May 2004)
!> The charge-moments are calculated in the subroutine vintras,
!> therefore vintras has to be called first.
!> The madelung-potential is expanded into spherical harmonics.
!> The lm-term of the potential \f$V\f$ of the atom \f$i\f$ is given by
!> \f$ V(r,lm,i) = \sum_{i2}^{N} \sum_{l'm'} (-r)^l * \left\{avmad(i,i2,lm,l'm')*cmom(i2,l'm') +bvmad(i,i2,lm)*z(i2)\right\}\f$
!> where \f$ N\f$ is the number of atoms
!> @author B. Drittler
!> @date Nov. 1989
!> @note
!> - V. Popescu Feb. 2002: Adopted for the case of more atoms on the same site, summation is done over the occupants of that site, the charge is weighted with the appropriate concentration of the occupant
!> - Impurity-program adopted feb. 2004 (according to N. Papanikalou)
!> - Jonathan Chico Jan. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
    Subroutine vmadelblk(cmom, cminst, lmax, nspin, naez, v, zat, r, irws, &
      ircut, ipan, kshape, noq, kaoez, conc, catom, icc, hostimp, vinters, &
      irm, nemb, lmpot, npotd, lmmaxd, natyp)

      Use constants
      Use global_variables
      Use mod_datatypes, Only: dp

      Implicit None

! .. Input variables
      Integer, Intent (In) :: icc !< Enables the calculation of off-diagonal elements of the GF.(0=SCF/DOS; 1=cluster; -1=custom)
      Integer, Intent (In) :: irm !< Maximum number of radial points
      Integer, Intent (In) :: naez !< Number of atoms in unit cell
      Integer, Intent (In) :: lmax !< Maximum l component in wave function expansion
      Integer, Intent (In) :: nemb !< Number of 'embedding' positions
      Integer, Intent (In) :: natyp !< Number of kinds of atoms in unit cell
      Integer, Intent (In) :: nspin !< Counter for spin directions
      Integer, Intent (In) :: lmpot !< (LPOT+1)**2
      Integer, Intent (In) :: npotd !< (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP)
      Integer, Intent (In) :: kshape !< Exact treatment of WS cell
      Integer, Intent (In) :: lmmaxd !< (KREL+KORBIT+1)(LMAX+1)^2
! .. Array Arguments
      Integer, Dimension (naez), Intent (In) :: noq !< Number of diff. atom types located
      Integer, Dimension (natyp), Intent (In) :: irws !< Position of atoms in the unit cell in units of bravais vectors
      Integer, Dimension (natyp), Intent (In) :: ipan !< Number of panels in non-MT-region
      Integer, Dimension (0:natyp), Intent (In) :: hostimp
      Integer, Dimension (0:ipand, natyp), Intent (In) :: ircut !< R points of panel borders
      Integer, Dimension (natyp, naez+nemb), Intent (In) :: kaoez !< Kind of atom at site in elem. cell
      Real (Kind=dp), Dimension (natyp), Intent (In) :: zat !< Nuclear charge
      Real (Kind=dp), Dimension (natyp), Intent (In) :: conc !< Concentration of a given atom
      Real (Kind=dp), Dimension (natyp), Intent (In) :: catom
      Real (Kind=dp), Dimension (irm, natyp), Intent (In) :: r !< Radial mesh ( in units a Bohr)
      Real (Kind=dp), Dimension (lmpot, natyp), Intent (In) :: cmom !< LM moment of total charge
      Real (Kind=dp), Dimension (lmpot, natyp), Intent (In) :: cminst !< charge moment of interstitial
! .. Input/Ouput variables
      Real (Kind=dp), Dimension (irm, lmpot, npotd), Intent (Inout) :: v
! .. Output variables
      Real (Kind=dp), Dimension (lmpot, naez), Intent (Out) :: vinters
! .. Local Scalars
      Integer :: lrecabmad, irec
      Integer :: i, l, lm, lm2, lmmax, m, io1, io2, ipot, iq1, iq2
      Integer :: irs1, ispin, it1, it2, noqval
      Real (Kind=dp) :: ac
! .. Local Arrays
      Real (Kind=dp), Dimension (lmpot) :: bvmad !< Structure dependent matrix
      Real (Kind=dp), Dimension (lmpot, lmpot) :: avmad !< Structure dependent matrix
      Logical :: opt
! .. Intrinsic Functions ..
      Intrinsic :: sqrt
!----------------------------------------------------------------------------
      Write (1337, Fmt=100)
      Write (1337, Fmt=110)
!
      lrecabmad = wlength*2*lmpot*lmpot + wlength*2*lmpot
      Open (69, Access='direct', Recl=lrecabmad, File='abvmad.unformatted', &
        Form='unformatted')
!
      lmmax = (lmax+1)*(lmax+1)
!
      If (icc/=0) Then
        Do iq1 = 1, naez
          Do lm = 1, lmpot
            vinters(lm, iq1) = 0E0_dp
          End Do
        End Do
      End If
!----------------------------------------------------------------------------
! Loop over all types in unit cell
!----------------------------------------------------------------------------
      Do iq1 = 1, naez ! added bauer 2/7/2012
        noqval = noq(iq1) ! added bauer 2/7/2012
        If (noqval<1) noqval = 1 ! added bauer 2/7/2012
        Do io1 = 1, noqval ! added bauer 2/7/2012
          it1 = kaoez(io1, iq1) ! added bauer 2/7/2012

!----------------------------------------------------------------------
! Take a site occupied by atom IT1
!----------------------------------------------------------------------
          If (it1/=-1) Then ! added bauer 2/7/2012
            If (kshape/=0) Then
              irs1 = ircut(ipan(it1), it1)
            Else
              irs1 = irws(it1)
            End If
          End If ! added bauer 2/7/2012
!----------------------------------------------------------------------
          Do l = 0, lmax
!-------------------------------------------------------------------
            Do m = -l, l
              lm = l*l + l + m + 1
              ac = 0.0E0_dp
!----------------------------------------------------------------
              If (naez==1) Then
                irec = iq1 + naez*(iq1-1)
                Read (69, Rec=irec) avmad, bvmad
!-------------------------------------------------------------
! Loop over all occupants of site IQ2=IQ1
!-------------------------------------------------------------
                Do io2 = 1, noq(iq1)
                  it2 = kaoez(io2, iq1)
!----------------------------------------------------------
! lm = 1 component disappears if there is only one host atom
! take moments of sphere
!----------------------------------------------------------
                  Do lm2 = 2, lmmax
                    ac = ac + avmad(lm, lm2)*cmom(lm2, it2)*conc(it2)
                  End Do
!----------------------------------------------------------
! Add contribution of interstial in case of shapes
!----------------------------------------------------------
                  If (kshape/=0) Then
                    Do lm2 = 2, lmmax
                      ac = ac + avmad(lm, lm2)*cminst(lm2, it2)*conc(it2)
                    End Do
                  End If
                End Do
!-------------------------------------------------------------
              Else
!-------------------------------------------------------------
! Loop over all sites
!-------------------------------------------------------------
                Do iq2 = 1, naez
                  irec = iq2 + naez*(iq1-1)
                  Read (69, Rec=irec) avmad, bvmad
!----------------------------------------------------------
! Loop over all occupants of site IQ2
!----------------------------------------------------------
                  Do io2 = 1, noq(iq2)
!
                    it2 = kaoez(io2, iq2)
                    ac = ac + bvmad(lm)*zat(it2)*conc(it2)
!-------------------------------------------------------
! Take moments of sphere
!-------------------------------------------------------
                    Do lm2 = 1, lmmax
                      ac = ac + avmad(lm, lm2)*cmom(lm2, it2)*conc(it2)
                    End Do
!-------------------------------------------------------
! Add contribution of interstial in case of shapes
!-------------------------------------------------------
                    If (kshape/=0) Then
                      Do lm2 = 1, lmmax
                        ac = ac + avmad(lm, lm2)*cminst(lm2, it2)*conc(it2)
                      End Do
                    End If
                  End Do ! IO2 = 1, NOQ(IQ2)
!----------------------------------------------------------
                End Do ! IQ2 = 1, NAEZ
!-------------------------------------------------------------
              End If ! NAEZ.GT.1
!----------------------------------------------------------------
              If (lm==1) Then
                Write (1337, Fmt=120) it1, (catom(it1)-zat(it1)), &
                  (ac/sqrt(4.E0_dp*pi))
              End If
!----------------------------------------------------------------
! Add to v the intercell-potential
!----------------------------------------------------------------
!----------------------------------------------------------------
! SPIN
!----------------------------------------------------------------
              Do ispin = 1, nspin
!-------------------------------------------------------------
! Determine the right potential number
!-------------------------------------------------------------
                ipot = nspin*(it1-1) + ispin
!-------------------------------------------------------------
! In the case of l=0 : r(1)**l is not defined
!-------------------------------------------------------------
                If (it1/=-1) Then ! added bauer 2/7/2012
                  If (l==0) v(1, 1, ipot) = v(1, 1, ipot) + ac
                  Do i = 2, irs1
                    v(i, lm, ipot) = v(i, lm, ipot) + (-r(i,it1))**l*ac
                  End Do
                End If
              End Do ! added bauer 2/7/2012
!----------------------------------------------------------------
! SPIN
!----------------------------------------------------------------
              If (icc/=0 .Or. opt('KKRFLEX ')) Then
                lm = l*l + l + m + 1
                Write (1337, *) 'ac', iq1, lm, ac
                vinters(lm, iq1) = ac
              End If
!
            End Do
!-------------------------------------------------------------------
          End Do
!----------------------------------------------------------------------
        End Do
      End Do
!----------------------------------------------------------------------------
      Close (69)
!----------------------------------------------------------------------------
      Write (1337, *) 'ICC in VMADELBLK', icc
      Write (1337, '(25X,30("-"),/)')
      Write (1337, '(79("="))')
!
      If ((icc==0) .And. (.Not. opt('KKRFLEX '))) Return
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now Prepare output for Impurity calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Open (91, File='intercell_ref', Status='unknown', Form='formatted')
      Write (1337, *)
      Write (1337, *) '                     ', &
        'Writing intercell potential for impurity'
      Write (1337, '(/,20X,55("-"))')
      Write (1337, 130) hostimp(0), lmmax
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
!
      Write (91, 140) hostimp(0), lmmax
      Do i = 1, hostimp(0)
        Write (91, 150)(vinters(lm,hostimp(i)), lm=1, lmmax)
      End Do
      Close (91)

      Return
!
100   Format (79('='), /, 18X, ' MADELUNG POTENTIALS ', &
        '(spherically averaged) ')
110   Format (/, 25X, ' ATOM ', '  Delta_Q  ', '     VMAD', /, 25X, 30('-'))
120   Format (25X, I4, 2X, F10.6, 1X, F12.6)
130   Format (22X, I4, ' host atoms, LMPOT = ', I2, ' output up to LM = 9')
140   Format (3I6)
150   Format (4D20.10)
    End Subroutine
