!-------------------------------------------------------------------------------
! SUBROUTINE: STARTB1
!> @brief Reads the input potentials
!> @details   units :
!> - Rydbergs - units for energy
!> - The lattice constant and all other lengths given in bohr units
!> - The planck constant \f$\frac{h}{2\pi}=1\f$
!> - The electron charge \f$e=\sqrt{2}\f$
!> - The electron mass \f$m=\frac{1}{2}\f$
!> - The speed of light \f$c = \frac{2}{\alpha} = 274.0720442\f$ with the fine structure constant \f$\alpha\f$
!>
!> In case of shape corrections this routine reads from unit 19 a suitable radial mesh 'xrn',its derivate 'drn' and the shape
!> functions 'thetas'. Thus, the region from the muffin-tin to the circumscribed sphere radii is divided  into 'npan'
!> pannels, each one containing 'nm(ipan)' points in order to take care of the  discontinuities of the shape-function  derivative.
!
!> @note remember that the input potentials do not include the electro-static contribution of the nucleus of the cell itself
!> this has to be added explicitly!
!> @note Modified for bandstructure code
!> - Jonathan Chico: Removed inc.p dependencies and rewrote to Fortran90
!> @author B. Drittler
!> @date Nov. 1989
!-------------------------------------------------------------------------------
    Subroutine startb1(ifile, ipf, ipfe, ipe, krel, kws, lmax, nbeg, nend, &
      alat, rmtnew, rmt, ititle, imt, irc, vconst, ins, irns, fpradius, lpot, &
      nspin, vins, irmin, kshape, ntcell, ircut, ipan, thetas, ifunm, nfu, &
      llmsp, lmsp, efermi, vbc, dror, rs, s, vm2z, rws, ecore, lcore, ncore, &
      drdi, r, zat, a, b, irws, iinfo, lmpot, irmind, irm, lmxspd, ipand, &
      irid, irnsd, natyp, ncelld, nfund, nspotd, ivshift, npotd)

      Use constants
      Use mod_datatypes, Only: dp

      Implicit None
! ..
! .. Input variables
      Integer, Intent (In) :: irm !< Maximum number of radial points
      Integer, Intent (In) :: kws !< 0 (MT), 1(ASA)
      Integer, Intent (In) :: ins !< 0 (MT), 1(ASA), 2(Full Potential)
      Integer, Intent (In) :: ipe !< Not real used, IPFE should be 0
      Integer, Intent (In) :: ipf !< Not real used, IPFE should be 0
      Integer, Intent (In) :: ipfe !< Not real used, IPFE should be 0
      Integer, Intent (In) :: irid !< Shape functions parameters in non-spherical part
      Integer, Intent (In) :: lmax !< Maximum l component in wave function expansion
      Integer, Intent (In) :: lpot !< Maximum l component in potential expansion
      Integer, Intent (In) :: nbeg !< Starting number for reading the potential
      Integer, Intent (In) :: nend !< Final number for reading the potential
      Integer, Intent (In) :: krel !< Switch for non-relativistic/relativistic (0/1) program. Attention: several other parameters depend explicitly on KREL, they are set automatically Used for Dirac solver in ASA
      Integer, Intent (In) :: npotd !< (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP)
      Integer, Intent (In) :: nspin !< Counter for spin directions
      Integer, Intent (In) :: ipand !< Number of panels in non-spherical part
      Integer, Intent (In) :: irnsd
      Integer, Intent (In) :: nfund !< Shape functions parameters in non-spherical part
      Integer, Intent (In) :: ifile !< Unit specifier for potential card
      Integer, Intent (In) :: iinfo
      Integer, Intent (In) :: natyp !< Number of kinds of atoms in unit cell
      Integer, Intent (In) :: lmpot !< (LPOT+1)**2
      Integer, Intent (In) :: irmind !< IRM-IRNSD
      Integer, Intent (In) :: lmxspd !< (2*LPOT+1)**2
      Integer, Intent (In) :: ncelld !< Number of cells (shapes) in non-spherical part
      Integer, Intent (In) :: nspotd !< Number of potentials for storing non-sph. potentials
      Integer, Intent (In) :: kshape !< Exact treatment of WS cell
      Integer, Intent (In) :: ivshift
      Real (Kind=dp), Intent (In) :: vconst !< Potential shift
      Integer, Dimension (natyp), Intent (In) :: ntcell !< Index for WS cell
      Real (Kind=dp), Dimension (natyp), Intent (In) :: fpradius !< R point at which full-potential treatment starts
! .. In/Out variables
      Real (Kind=dp), Intent (Inout) :: alat !< Lattice constant in a.u.
      Real (Kind=dp), Intent (Inout) :: efermi !< Fermi energy
      Integer, Dimension (natyp), Intent (Inout) :: nfu !< number of shape function components in cell 'icell'
      Integer, Dimension (natyp), Intent (Inout) :: imt !< R point at MT radius
      Integer, Dimension (natyp), Intent (Inout) :: irc !< R point for potential cutting
      Integer, Dimension (natyp), Intent (Inout) :: ipan !< Number of panels in non-MT-region
      Integer, Dimension (natyp), Intent (Inout) :: irns !< Position of atoms in the unit cell in units of bravais vectors
      Integer, Dimension (natyp), Intent (Inout) :: irws !< R point at WS radius
      Integer, Dimension (natyp), Intent (Inout) :: irmin !< Max R for spherical treatment
      Integer, Dimension (npotd), Intent (Inout) :: ncore !< Number of core states
      Integer, Dimension (natyp, lmxspd), Intent (Inout) :: lmsp !< 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
      Integer, Dimension (20, npotd), Intent (Inout) :: lcore !< Angular momentum of core states
      Integer, Dimension (natyp, nfund), Intent (Inout) :: llmsp !< lm=(l,m) of 'nfund'th nonvanishing component of non-spherical pot.
      Integer, Dimension (0:ipand, natyp), Intent (Inout) :: ircut !< R points of panel borders
      Integer, Dimension (natyp, lmxspd), Intent (Inout) :: ifunm
      Integer, Dimension (20, npotd), Intent (Inout) :: ititle !< Titles of the potential card
      Real (Kind=dp), Dimension (natyp), Intent (Inout) :: a !< Constants for exponential R mesh
      Real (Kind=dp), Dimension (natyp), Intent (Inout) :: b !< Constants for exponential R mesh
      Real (Kind=dp), Dimension (natyp), Intent (Inout) :: zat !< Nuclear charge
      Real (Kind=dp), Dimension (2), Intent (Inout) :: vbc !< Potential constants
      Real (Kind=dp), Dimension (natyp), Intent (Inout) :: rmt !< Muffin-tin radius of true system
      Real (Kind=dp), Dimension (natyp), Intent (Inout) :: rws !< Wigner Seitz radius
      Real (Kind=dp), Dimension (natyp), Intent (Inout) :: rmtnew !< Adapted muffin-tin radius
      Real (Kind=dp), Dimension (0:lmax, natyp), Intent (Inout) :: s
      Real (Kind=dp), Dimension (irm, natyp), Intent (Inout) :: r !< Radial mesh ( in units a Bohr)
      Real (Kind=dp), Dimension (irm, natyp), Intent (Inout) :: drdi !< Derivative dr/di
      Real (Kind=dp), Dimension (irm, natyp), Intent (Inout) :: dror
      Real (Kind=dp), Dimension (irm, npotd), Intent (Inout) :: vm2z
      Real (Kind=dp), Dimension (20, npotd), Intent (Inout) :: ecore !< Core energies
      Real (Kind=dp), Dimension (irm, 0:lmax, natyp), Intent (Inout) :: rs
      Real (Kind=dp), Dimension (irmind:irm, lmpot, nspotd), &
        Intent (Inout) :: vins !< Non-spherical part of the potential
      Real (Kind=dp), Dimension (irid, nfund, ncelld), &
        Intent (Inout) :: thetas !< shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
! .. Local Scalars
      Integer :: inslpd, lmshapemax
      Integer :: j, l, lm, lm1, lmpotp, n, ncell, nfun, nr
      Integer :: irminm, irminp, irns1p, irt1p, irws1, isave, ispin, isum
      Integer :: i, ia, icell, icore, ifun, ih, imt1, inew, io, ipan1, ir, &
        irc1, iri
      Real (Kind=dp) :: a1, b1, ea, efnew, s1, z1
      Logical :: test
! ..
! .. Local Arrays
      Integer, Dimension (ncelld) :: npan
      Integer, Dimension (ncelld) :: meshn
      Integer, Dimension (ipand, ncelld) :: nm
      Real (Kind=dp), Dimension (irm) :: u
      Real (Kind=dp), Dimension (ncelld) :: scale
      Real (Kind=dp), Dimension (irid) :: rdummy
      Real (Kind=dp), Dimension (irm) :: vspsme ! dummy for potcut IMPURITY-compatible
      Real (Kind=dp), Dimension (irid, ncelld) :: drn
      Real (Kind=dp), Dimension (irid, ncelld) :: xrn
! ..
! .. External Subroutines ..
      External :: calrmt, potcut, rinit, test
! ..
! .. Intrinsic Functions ..
      Intrinsic :: nint, exp, log, max, mod, real, sqrt
! ..
! .. Data statement ..
      Integer :: ishape
      Data ishape/0/
!----------------------------------------------------------------------------
! Output of radial mesh information
!----------------------------------------------------------------------------
      io = 0
      If (iinfo/=0 .And. test('RMESH   ')) io = 1
!----------------------------------------------------------------------------
! Set speed of light
!----------------------------------------------------------------------------
      inslpd = (irnsd+1)*lmpot*nspotd
      lmshapemax = (4*lmax+1)**2
      Call rinit(inslpd, vins(irmind,1,1))
!----------------------------------------------------------------------------
! Read radial mesh information of the shape functions and
! shape functions THETAS in the first iteration - if needed
!----------------------------------------------------------------------------
      If ((kshape/=0) .And. (ishape==0)) Then
        ishape = 1
        Read (19, Fmt=100) ncell
        Write (1337, Fmt=*) '  ncell : ', ncell, ncelld
!       check consistency with shape numbers from inputcard
        If (maxval(ntcell(1:natyp))>ncell) Then
          Write (*, *) 'Found ', ncell, 'shapes in shapefun file but need', &
            maxval(ntcell(1:natyp)), 'according to inputcard/default values'
          Write (*, *) 'Did you set <SHAPE> correctly in inputcard?'
          Stop 'Error consistency shapes from input/shapefun file'
        End If
!
        If (ncell>ncelld) Then
          Write (6, *) 'Please, change the parameter ncelld (', ncelld, &
            ') in the inputcard to', ncell
          Stop 'STARTB - NCELLD'
        End If
!
        Read (19, Fmt=110)(scale(icell), icell=1, ncell)
        Do icell = 1, ncell
          Read (19, Fmt=100) npan(icell), meshn(icell)
!
          If (npan(icell)+1>ipand) Then
            Write (6, *) 'Please, change the parameter ipand (', ipand, &
              ') in the inputcard to', npan(icell) + 1
            Stop 'STARTB - IPAND'
          End If
!
          If (meshn(icell)>irid) Then
            Write (6, *) 'Please, change the parameter irid (', irid, &
              ') in the inputcard to', meshn(icell)
            Stop 'STARTB - IRID'
          End If
!
          Read (19, Fmt=100)(nm(ipan1,icell), ipan1=2, npan(icell)+1)
          Read (19, Fmt=110)(xrn(ir,icell), drn(ir,icell), ir=1, meshn(icell))

          Read (19, Fmt=100) nfu(icell)
          nfun = nfu(icell)
          Write (1337, Fmt=*) '  nfun  : ', nfun, nfund
!
          If (nfun>nfund) Then
            Write (6, *) 'Please, change the parameter nfund (', nfund, &
              ') in the inputcard to', nfun
            Stop 'STARTB - NFUND'
          End If
!
          Do lm = 1, lmxspd
            lmsp(icell, lm) = 0
          End Do

          Do ifun = 1, nfun
            Read (19, Fmt=100) lm
            If (lm<=lmshapemax) Then
              llmsp(icell, ifun) = lm
              lmsp(icell, lm) = 1
              ifunm(icell, lm) = ifun
              Read (19, Fmt=110)(thetas(n,ifun,icell), n=1, meshn(icell))
            Else
              Read (19, Fmt=110)(rdummy(n), n=1, meshn(icell))
            End If
          End Do

        End Do
      End If ! ((KSHAPE.NE.0) .AND. (IFILE.NE.0))
!----------------------------------------------------------------------------
!LMPOT = (LPOT+1)* (LPOT+1)
      Do ih = nbeg, nend
        Do ispin = 1, nspin
          i = nspin*(ih-1) + ispin

          If (ifile/=0) Then
            ircut(0, ih) = 0
            If (ins/=0) Then
! p.z.            IF (KSHAPE.NE.0) THEN
              icell = ntcell(ih)
              ipan(ih) = 1 + npan(icell)
            Else
              ipan(ih) = 1
            End If
!-------------------------------------------------------------------
! Read title of potential card
!-------------------------------------------------------------------
            Read (ifile, Fmt=120)(ititle(ia,i), ia=1, 20)
            If (iinfo/=0) Then
              If (ins==0) Then
                Write (1337, Fmt=180)(ititle(ia,i), ia=1, 20)
              Else
                Write (1337, Fmt=190)(ititle(ia,i), ia=1, 20)
              End If
            End If
!
!-------------------------------------------------------------------
! Read muffin-tin radius , lattice constant and new muffin radius
! (new mt radius is adapted to the given radial mesh)
!-------------------------------------------------------------------
            Read (ifile, Fmt=*) rmt(ih), alat, rmtnew(ih)
!READ (IFILE,FMT=9030) RMT(IH),ALAT,RMTNEW(IH)
!
!-------------------------------------------------------------------
! Read nuclear charge , lmax of the core states ,
! wigner seitz radius , fermi energy and energy difference
! between electrostatic zero and muffin tin zero
!-------------------------------------------------------------------
!            READ (IFILE,FMT=9040) ZAT(IH),RWS(IH),EFNEW,VBC(ISPIN)
            Read (ifile, *) z1
            Read (ifile, *) rws(ih), efnew, vbc(ispin)

!             READ (IFILE,*) Z1,RWS(IH),EFNEW,VBC(ISPIN)
            If (zat(ih)<0.E0_dp) zat(ih) = z1
            If (z1/=zat(ih) .And. zat(ih)>=0.E0_dp) Then
              Write (*, *) 'Warning: For atom ', ih, &
                ': ZATOM different in inputcard and in potential.', zat(ih), &
                z1
            End If
!-------------------------------------------------------------------
! If efermi .eq. 0 use value from in5
!-------------------------------------------------------------------
            If (efnew/=0.0E0_dp .And. i==1) efermi = efnew
!-------------------------------------------------------------------
! Read : number of radial mesh points
! (in case of ws input-potential: last mesh point corresponds
! to ws-radius, in case of shape-corrected input-potential
! last mesh point of the exponential mesh corresponds to
! mt-radius/nevertheless this point is always in the array
! irws(ih)),number of points for the radial non-muffin-tin
! mesh  needed for shape functions, the constants a and b
! for the radial exponential mesh : r(i) = b*(exp(a*(i-1))-1)
! the no. of different core states and some other stuff
!-------------------------------------------------------------------
!READ (IFILE,FMT=9050) IRWS(IH),A(IH),B(IH),NCORE(I),INEW
            Read (ifile, Fmt=*) irws(ih)
            Read (ifile, Fmt=*) a(ih), b(ih)
            Read (ifile, Fmt=*) ncore(i), inew
!
            nr = irws(ih)

            If (nr>irm) Then
              Write (6, *) 'Increase parameter IRM in the inputcard ', &
                ' to a value .ge. ', nr, ' (= IRWS(', ih, ')).'
              Stop 'STARTB1 - IRWS'
            End If
!-------------------------------------------------------------------
! Read the different core states : l and energy
!-------------------------------------------------------------------
            If (ncore(i)>=1) Then
              Do icore = 1, ncore(i)
                Read (ifile, Fmt=170) lcore(icore, i), ecore(icore, i)
              End Do
            End If
!
            If (ins<1) Then
!----------------------------------------------------------------
! Read radial mesh points, its derivative, the spherically averaged
! charge density and the input potential without the nuclear pot.
!----------------------------------------------------------------
              If (inew==0) Then
                Read (ifile, Fmt=160)(r(ir,ih), drdi(ir,ih), vm2z(ir,i), ir=1, &
                  nr)
              Else
                Read (ifile, Fmt=*)(vm2z(ir,i), ir=1, nr)
              End If
            Else ! (INS.LT.1)
!-------------------------------------------------------------------
! Read full potential - the non spherical contribution from irmin
! to irt - remember that the lm = 1 contribution is multiplied by
! 1/sqrt(4 pi)
!-------------------------------------------------------------------
              Read (ifile, Fmt=200) irt1p, irns1p, lmpotp, isave
              irminp = irt1p - irns1p
              irminm = max(irminp, irmind)
              Read (ifile, Fmt=210)(vm2z(ir,i), ir=1, nr)
              If (lmpotp>1) Then
                lm1 = 2
                Do lm = 2, lmpotp
                  If (lm1/=1) Then
                    If (isave==1) Then
                      Read (ifile, Fmt=200) lm1
                    Else
                      lm1 = lm
                    End If

                    If (lm1>1) Then
                      Read (ifile, Fmt=210)(u(ir), ir=irminp, nr)
                      If (lm1<=lmpot) Then
                        Do ir = irminm, nr
                          vins(ir, lm1, i) = u(ir)
                        End Do
                      End If
                    End If
                  End If
                End Do
              End If
            End If ! (INS.LT.1)
!
            irws1 = irws(ih)
!----------------------------------------------------------------------
! Redefine new mt-radius in case of shape corrections
!----------------------------------------------------------------------
            If (ins/=0) Then
! p.z.      IF (KSHAPE.NE.0) THEN
              rmtnew(ih) = scale(icell)*alat*xrn(1, icell)
              imt1 = nint(log(rmtnew(ih)/b(ih)+1.0E0_dp)/a(ih)) + 1
!-------------------------------------------------------------------
! For proper core treatment imt must be odd
! shift potential by one mesh point if imt is even
!-------------------------------------------------------------------
              If (mod(imt1,2)==0) Then
                imt1 = imt1 + 1
                Do ir = imt1, 2, -1
                  vm2z(ir, i) = vm2z(ir-1, i)
                End Do
              End If
!
              imt(ih) = imt1
              b(ih) = rmtnew(ih)/(exp(a(ih)*real(imt1-1,kind=dp))-1.0E0_dp)
            End If ! (KSHAPE.NE.0)
!----------------------------------------------------------------------
! Generate radial mesh - potential only is stored in potential card
! INEW = 1
! p. zahn, jan. 99
!----------------------------------------------------------------------
            a1 = a(ih)
            b1 = b(ih)
            r(1, ih) = 0.0E0_dp
            drdi(1, ih) = a1*b1
            Do ir = 2, irws1
              ea = exp(a1*real(ir-1,kind=dp))
              r(ir, ih) = b1*(ea-1.0E0_dp)
              drdi(ir, ih) = a1*b1*ea
              dror(ir, ih) = a1/(1.0E0_dp-1.0E0_dp/ea)
            End Do
!----------------------------------------------------------------------
! Fill cell-type depending mesh points in the non-muffin-tin-region
!----------------------------------------------------------------------
            If (ins/=0) Then
! p.z.      IF (KSHAPE.NE.0) THEN
              Do iri = 1, meshn(icell)
                ir = iri + imt1
                r(ir, ih) = scale(icell)*alat*xrn(iri, icell)
                drdi(ir, ih) = scale(icell)*alat*drn(iri, icell)
                dror(ir, ih) = drdi(ir, ih)/r(ir, ih)
              End Do
            End If

            rws(ih) = r(irws1, ih)
!----------------------------------------------------------------------
! Kshape.eq.0 : calculate new rmt adapted to exp. mesh
!----------------------------------------------------------------------
            Call calrmt(ipf, ipfe, ipe, imt(ih), zat(ih), rmt(ih), rws(ih), &
              rmtnew(ih), alat, drdi(1,ih), a(ih), b(ih), irws1, r(1,ih), io, &
              ins)
! p.z. +                  R(1,IH),IO,KSHAPE)
!
            If (ins>0) Then
! p.z.            IF (KSHAPE.GT.0) THEN
              ircut(1, ih) = imt(ih)
              isum = imt(ih)
              Do ipan1 = 2, ipan(ih)
                isum = isum + nm(ipan1, icell)
                ircut(ipan1, ih) = isum
              End Do
              nr = isum
              If (irt1p/=nr) Then
                Write (*, *) 'STARTB1: Error: IRT1P.NE.NR', irt1p, nr, &
                  ' for atom', ih
                Stop 'STARTB1: IRT1P.NE.NR'
              End If
            Else ! (KSHAPE.GT.0)
              nr = irws(ih)
              If (kws>=1) Then
                ircut(1, ih) = irws1
              Else
                ircut(1, ih) = imt(ih)
              End If
            End If ! (KSHAPE.GT.0)
!
            irc(ih) = ircut(ipan(ih), ih)
!----------------------------------------------------------------------
! Fill array irmin in case of full potential
!----------------------------------------------------------------------
            If (ins/=0) Then
              If (fpradius(ih)>=0.E0_dp) Then
                irmin(ih) = min(floor(log(fpradius(ih)/b(ih)+ &
                  1.0E0_dp)/a(ih))+1, imt(ih))
                irns(ih) = nr - irmin(ih)
              Else If (irns(ih)>=meshn(icell)) Then
                irmin(ih) = nr - irns(ih)
              Else
                irns(ih) = irns1p
                irmin(ih) = nr - irns(ih)
              End If
            End If
!----------------------------------------------------------------------
! Generate arrays for the calculation of the wave functions
!----------------------------------------------------------------------
            z1 = zat(ih)
            Do l = 0, lmax
              If (krel>=1) Then
                s1 = sqrt(real(l*l+l+1,kind=dp)-4.0E0_dp*z1*z1/(cvlight* &
                  cvlight))
                If (z1==0.0E0_dp) s1 = real(l, kind=dp)
              Else
                s1 = real(l, kind=dp)
              End If

              s(l, ih) = s1
              rs(1, l, ih) = 0.0E0_dp
              Do ir = 2, nr
                rs(ir, l, ih) = r(ir, ih)**s1
              End Do
            End Do ! L = 0,LMAX
!-------------------------------------------------------------------
! Cut input potential at rmt if given only at exponential mesh
!-------------------------------------------------------------------
            If (kshape==1) Then
              imt1 = imt(ih)
              irc1 = ircut(ipan(ih), ih)
              Call potcut(imt1, irc1, ins, lmpot, r(1,ih), vm2z(1,i), vspsme, &
                vins(irmind,1,i), zat(ih), irm, irmind)
            End If
!-------------------------------------------------------------------
! First iteration : shift all potentials (only for test purpose)
! in case of test option 'atptshft' shift only potential of atom at position ivshift
!-------------------------------------------------------------------
            If (test('atptshft') .And. (ih==ivshift)) Then
              Write (1337, *) 'atptshft', ih, ivshift, vconst, nr, irmin(ih)
              Do j = 1, irmin(ih)
                vm2z(j, i) = vm2z(j, i) + vconst
              End Do
            Else If (.Not. test('atptshft')) Then
              Do j = 1, nr
                vm2z(j, i) = vm2z(j, i) + vconst
              End Do
            End If
          End If ! (ifile.ne.0)
!
          If (kshape==0 .And. kws==0) Then
!-------------------------------------------------------------------
! In case of a mt calculation cut potential at mt radius
!-------------------------------------------------------------------
            imt1 = imt(ih)
            irws1 = irws(ih)
            Call potcut(imt1, irws1, ins, lmpot, r(1,ih), vm2z(1,i), vspsme, &
              vins(irmind,1,i), zat(ih), irm, irmind)
          End If ! KSHAPE.EQ.0 .AND. KWS.EQ.0
        End Do ! ISPIN = 1,NSPIN
      End Do ! IH = NBEG,NEND

      If (ins/=0) Then
        i = 0
        Do ih = nbeg, nend
          If (irmin(ih)<irmind) Then
            Write (*, *) 'IRMIN < IRMIND for atom', ih
            Write (*, *) irmin(ih), irmind
            Write (*, *) 'Increase dimension IRNSD'
            i = 1
          End If
        End Do
        If (i/=0) Stop 'stop startb1 IRNS IRNSD'
      End If

      Return

100   Format (16I5)
110   Format (4D20.12)
120   Format (20A4)
130   Format (3F12.8)
140   Format (F10.5, /, F10.5, 2F15.10)
150   Format (I3, /, 2D15.8, /, 2I2)
160   Format (1P, 2D15.6, 1P, D15.8)
170   Format (I5, 1P, D20.11)
! 9080 format (10x,20a4)
180   Format (' < ', 20A4)
190   Format (' <#', 20A4)
200   Format (10I5)
210   Format (1P, 4D20.13)
    End Subroutine
