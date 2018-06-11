!-------------------------------------------------------------------------------
! SUBROUTINE: RITES
!> @brief this subroutine stores in 'ifile' the necessary results
!> (potentials etc.) to start self-consistency iterations
!
!> @ details Modified for the full potential case - if ins .gt. 0 there
!> is written a different potential card
!> if the sum of absolute values of an lm component of vins (non
!> spher. potential) is less than the given rms error qbound this
!> component will not be stored .
!
!> see to subroutine start , where most of the arrays are described)
!
!> @note modified by B. Drittler  aug. 1988
!> @note Jonathan Chico Apr. 2019: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
    Subroutine rites(ifile, natps, natyp, nspin, z, alat, rmt, rmtnew, rws, &
      ititle, r, drdi, vm2z, irws, a, b, txc, kxc, ins, irns, lpot, vins, &
      qbound, irc, kshape, efermi, vbc, ecore, lcore, ncore, ecorerel, nkcore, &
      kapcore, irm, irmind, lmpot)

      Use global_variables
      Use mod_datatypes, Only: dp

! .. Scalar Arguments
      Integer, Intent (In) :: ins !< 0 (MT), 1(ASA), 2(Full Potential)
      Integer, Intent (In) :: irm !< Maximum number of radial points
      Integer, Intent (In) :: kxc !< Type of xc-potential 0=vBH 1=MJW 2=VWN 3=PW91
      Integer, Intent (In) :: lpot !< Maximum l component in potential expansion
      Integer, Intent (In) :: lmpot !< (LPOT+1)**2
      Integer, Intent (In) :: ifile !< Unit specifier for potential card
      Integer, Intent (In) :: natps
      Integer, Intent (In) :: natyp !< Number of kinds of atoms in unit cell
      Integer, Intent (In) :: nspin !< Counter for spin directions
      Integer, Intent (In) :: kshape !< Exact treatment of WS cell
      Integer, Intent (In) :: irmind !< IRM-IRNSD
      Real (Kind=dp), Intent (In) :: alat !< Lattice constant in a.u.
      Real (Kind=dp), Intent (In) :: qbound !< Convergence parameter for the potential
      Real (Kind=dp), Intent (In) :: efermi !< Fermi energy
! .. Array Arguments
      Real (Kind=dp), Dimension (*), Intent (In) :: a !< Constants for exponential R mesh
      Real (Kind=dp), Dimension (*), Intent (In) :: b !< Constants for exponential R mesh
      Real (Kind=dp), Dimension (*), Intent (In) :: z
      Real (Kind=dp), Dimension (*), Intent (In) :: rws !< Wigner Seitz radius
      Real (Kind=dp), Dimension (2), Intent (In) :: vbc !< Potential constants
      Real (Kind=dp), Dimension (*), Intent (In) :: rmt !< Muffin-tin radius of true system
      Real (Kind=dp), Dimension (*), Intent (In) :: rmtnew !< Adapted muffin-tin radius
      Real (Kind=dp), Dimension (irm, *), Intent (In) :: r !< Radial mesh ( in units a Bohr)
      Real (Kind=dp), Dimension (irm, *), Intent (In) :: vm2z
      Real (Kind=dp), Dimension (irm, *), Intent (In) :: drdi !< Derivative dr/di
      Real (Kind=dp), Dimension (20, *), Intent (In) :: ecore !< Core energies !(2), 22.5,2000
      Real (Kind=dp), Dimension (irmind:irm, lmpot, *), Intent (In) :: vins !< Non-spherical part of the potential
!----------------------------------------------------------------------------
      Integer, Dimension (20, natyp), Intent (In) :: nkcore
      Integer, Dimension (20, 2*natyp), Intent (In) :: kapcore
      Real (Kind=dp), Dimension (krel*20+(1-krel), 2*natyp), &
        Intent (In) :: ecorerel ! relativistic core energies
!----------------------------------------------------------------------------
      Integer, Dimension (*), Intent (In) :: irc !< R point for potential cutting
      Integer, Dimension (*), Intent (In) :: irns !< Position of atoms in the unit cell in units of bravais vectors
      Integer, Dimension (*), Intent (In) :: irws !< R point at WS radius
      Integer, Dimension (*), Intent (In) :: ncore !< Number of core states
      Integer, Dimension (20, *), Intent (In) :: ititle
      Integer, Dimension (20, *), Intent (In) :: lcore !< Angular momentum of core states
      Character (Len=124), Dimension (*), Intent (In) :: txc
! .. Local Scalars
      Integer :: i, icore, ih, inew, ip, ir, irmin, irns1, is, isave, j, lm, &
        lmnr, ncore1, nr
      Real (Kind=dp) :: a1, b1, rmax, rmt1, rmtnw1, rv, sign, sum, z1
! .. Local Arrays
      Integer, Dimension (20) :: lcore1
      Real (Kind=dp), Dimension (20) :: ecore1
      Real (Kind=dp), Dimension (irm) :: dradi
      Real (Kind=dp), Dimension (irm) :: ra
      Real (Kind=dp), Dimension (irm) :: vm2za
      Real (Kind=dp), Dimension (20, 2) :: ecore2
      Character (Len=3), Dimension (4) :: txtk
      Character (Len=1), Dimension (0:3) :: txtl
!     ..
      Data txtl/'s', 'p', 'd', 'f'/
      Data txtk/'1/2', '3/2', '5/2', '7/2'/
!     ..
! -------------------------------------------------------------------
      isave = 1
      inew = 1
!
!
      Do ih = 1, natyp
        Do is = 1, nspin
          If (is==nspin) Then
            sign = 1.0E0_dp
          Else
            sign = -1.0E0_dp
          End If
          ip = nspin*(ih-1) + is

          rmt1 = rmt(ih)
          rmtnw1 = rmtnew(ih)
          z1 = z(ih)
          rmax = rws(ih)
          If (kshape==0) Then
            nr = irws(ih)
          Else
            nr = irc(ih)
          End If

          irns1 = irns(ih)
          irmin = nr - irns1
          a1 = a(ih)
          b1 = b(ih)
          ncore1 = ncore(ip)
!
          Do j = 1, nr
            ra(j) = r(j, ih)
            dradi(j) = drdi(j, ih)
!-------------------------------------------------------------------
! Store only lm=1 component of the potential
!-------------------------------------------------------------------
            vm2za(j) = vm2z(j, ip)
          End Do ! J
!
          Open (ifile, File='out_potential', Form='formatted')
          Write (ifile, Fmt=100)(ititle(i,ip), i=1, 7), txc(kxc+1)
          Write (ifile, Fmt=110) rmt1, alat, rmtnw1
          Write (ifile, Fmt=120) z1, rmax, efermi, vbc(is)
          Write (ifile, Fmt=130) nr, a1, b1, ncore1, inew
!
          If (ncore1>=1) Then
!
            If (krel==0) Then
              Do j = 1, ncore1
                lcore1(j) = lcore(j, ip)
                ecore1(j) = ecore(j, ip)
              End Do
              Write (ifile, Fmt=140)(lcore1(icore), ecore1(icore), icore=1, &
                ncore1)
            Else
              Do j = 1, ncore1
                lcore1(j) = lcore(j, ip)
                ecore2(j, 1) = ecorerel(j, 2*ih-1)
                ecore2(j, 2) = ecorerel(j, 2*ih)
              End Do
!----------------------------------------------------------------
! independent of spin, the \mu-averaged relativistic core energies
! are written out for \kappa = -l-1,l
! format compatible with the non-(scalar) relativistic mode
! however, the next read in has no meaning for the REL core-solver
! a detailed output of the core energies is supplied by < CORE >
!----------------------------------------------------------------
              Do icore = 1, ncore1
                Write (ifile, Fmt=150) lcore1(icore), &
                  (ecore2(icore,i+1), txtl(lcore1(icore)), txtk(abs( &
                  kapcore(icore,2*ih-1+i))), i=0, nkcore(icore,ih)-1)
              End Do
            End If
!
          End If
!
          If (ins==0 .Or. (ih<natps .And. ins<=2)) Then
!-------------------------------------------------------------------
! store only the spherically averaged potential
! (in mt or as - case)
! this is done always for the host
!-------------------------------------------------------------------
            If (inew==0) Then
              Write (ifile, Fmt=160)(ra(ir), dradi(ir), vm2za(ir), ir=1, nr)
            Else
              Write (ifile, Fmt=170)(vm2za(ir), ir=1, nr)
            End If
          Else
!-------------------------------------------------------------------
! store the full potential , but the non spherical contribution
! only from irns1 up to irws1 ;
! remember that the lm = 1 contribution is multiplied
! by a factor 1/sqrt(4 pi)
!-------------------------------------------------------------------
            Write (ifile, Fmt=180) nr, irns1, lmpot, isave
            Write (ifile, Fmt=190)(vm2za(ir), ir=1, nr)
            If (lpot>0) Then
              lmnr = 1
              Do lm = 2, lmpot
                sum = 0.0E0_dp
                Do ir = irmin, nr
                  rv = vins(ir, lm, ip)*ra(ir)
                  sum = sum + rv*rv*dradi(ir)
                End Do ! IR
                If (sqrt(sum)>qbound) Then
                  lmnr = lmnr + 1
                  Write (ifile, Fmt=180) lm
                  Write (ifile, Fmt=190)(vins(ir,lm,ip), ir=irmin, nr)
                End If
              End Do !LM
!----------------------------------------------------------------
! Write a one to mark the end
!----------------------------------------------------------------
              If (lmnr<lmpot) Write (ifile, Fmt=180) isave
            End If
          End If
        End Do !IS
      End Do !IH

      Close (ifile)

100   Format (7A4, 6X, '  exc:', A124, 3X, A10)
110   Format (3F12.8)
! 9010 FORMAT (3F) !f12.8) maybe change to higher accuracy in writeout?
120   Format (F10.5, /, F10.5, 2F20.15)
130   Format (I3, /, 2D15.8, /, 2I2)
140   Format (I5, 1P, D20.11)
150   Format (I5, 2(1P,D20.11,2X,A1,A3))
160   Format (1P, 2D15.6, 1P, D15.8)
170   Format (1P, 4D20.12)
180   Format (10I5)
190   Format (1P, 4D20.13)
    End Subroutine
