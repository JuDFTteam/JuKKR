!-------------------------------------------------------------------------------
! SUBROUTINE: TESTDIM
!> @brief Testing the dimension of several arrays
!> @note Jonathan Chico: Some of these tests seem unnecessary with the changes done to the
!> inc.p
!-------------------------------------------------------------------------------
    Subroutine testdim(nspin, naez, nemb, natyp, lmax, irm, ins, insref, nref, &
      irns, ncls, nlayer, krel, nspind, nclsd, nprincd, knosph, irnsd, korbit)

      Use mod_datatypes, Only: dp
      Implicit None

      Integer, Intent (In) :: ins !< 0 (MT), 1(ASA), 2(Full Potential)
      Integer, Intent (In) :: irm !< Maximum number of radial points
      Integer, Intent (In) :: lmax !< Maximum l component in wave function expansion
      Integer, Intent (In) :: naez !< Number of atoms in unit cell
      Integer, Intent (In) :: nref !< Number of diff. ref. potentials
      Integer, Intent (In) :: krel !< Switch for non-relativistic/relativistic (0/1) program. Attention: several other parameters depend explicitly on KREL, they are set automatically Used for Dirac solver in ASA
      Integer, Intent (In) :: ncls !< Number of reference clusters
      Integer, Intent (In) :: nspin !< Counter for spin directions
      Integer, Intent (In) :: natyp !< Number of kinds of atoms in unit cell
      Integer, Intent (In) :: nclsd !< Maximum number of different TB-clusters
      Integer, Intent (In) :: irnsd
      Integer, Intent (In) :: insref !< INS for reference pot. (usual 0)
      Integer, Intent (In) :: knosph !< switch for spherical/non-spherical (0/1) program.
      Integer, Intent (In) :: korbit !< Spin-orbit/non-spin-orbit (1/0) added to the Schroedinger or SRA equations. Works with FP. KREL and KORBIT cannot be both non-zero.
      Integer, Intent (In) :: nspind !< KREL+(1-KREL)*(NSPIN+1)
      Integer, Intent (In) :: nprincd !< Number of principle layers, set to a number >= NRPINC in output of main0
! .. In/Out variables
      Integer, Intent (Inout) :: nemb !< Number of 'embedding' positions
      Integer, Intent (Inout) :: nlayer !< Number of principal layer
      Integer, Dimension (natyp), Intent (Inout) :: irns !< Position of atoms in the unit cell in units of bravais vectors
!
      Integer :: stop_mark
      Integer :: i, j
      Logical :: test, opt
      External :: test, opt
!
! ---> dimension tests
!
      Write (1337, 120)
!
      stop_mark = 0
      If ((nspin>nspind) .And. (krel==0)) Then
        Write (6, *) 'There is an inconsistenciy between spin polarised &
          &calculation and relativistic options'
        stop_mark = 1
      End If
      If (max(ins,insref)>knosph) Then
        Write (6, *) 'Please, change the parameter insd in', &
          ' the inputcard to', max(ins, insref)
        stop_mark = 1
      End If
      j = 1
      Do i = 1, natyp
        j = max(j, irns(i))
      End Do
      If (ins==0 .And. j>1) Then
        Write (6, *) 'IRNS(*) is set to 1 in case of ', &
          'spherical potential treatment.'
        Do i = 1, natyp
          irns(i) = 1
        End Do
        j = 1
      End If
!
      If (j>irnsd) Then
        Write (6, *) 'Please, change the parameter irnsd in', &
          ' the inputcard to', j
        stop_mark = 1
      End If
!
      If (.Not. opt('VIRATOMS')) Then
        If (nref>natyp) Then
          Write (6, *) 'There are some inconsistencies in the input file./', &
            ' nref(=', nref, ') is greater than natyp (=', natyp, ').'
          stop_mark = 1
        End If
      End If

      If ((krel==1) .And. (korbit==1)) Then
        Write (6, *) 'Full relativistic for ASA and new SO solver', 'KREL', &
          krel, 'KORBIT', korbit
        stop_mark = 1
      End If

      If (.Not. opt('NEWSOSOL') .And. korbit==1) Then
        Write (6, *) &
          'Option NEWSOSOL not found, change KORBIT in the inputcard from', &
          korbit, 'to 0'
        stop_mark = 1
      End If

      If (opt('NEWSOSOL') .And. korbit==0) Then
        Write (6, *) &
          'Using option NEWSOSOL, change KORBIT in the inputcard from', &
          korbit, 'to 1'
        stop_mark = 1
      End If
!----------------------------------------------------------------------------
!  OPT 'WIRE' is only useful with OPT 'full inv' or
!  OPT 'SPARSE  ' because of sparsity of
!  the KKR matrix ( not tridiagonal like for 2D and 3D systems)
!----------------------------------------------------------------------------
      If (opt('WIRE    ') .And. .Not. (opt('full inv') .Or. opt('SPARSE  '))) &
        Then
        Write (6, *) 'Use option ''full inv'' or ''SPARSE  '' ', &
          'for WIRE calculation.'
        stop_mark = 1
      End If
!
      If (opt('COMPLEX ') .And. .Not. (opt('EigenV  ') .Or. opt( &
        'wfct    ') .Or. opt('iso surf'))) Then
        Write (6, *) &
          'Use option ''COMPLEX '' only for eigenvalue determination.'
        stop_mark = 1
      End If
!
      If (test('CONT    ')) Then
        nemb = 0
        Write (6, *) 'No usage of embedding points. NEMB is set to ', nemb, &
          '.'
      End If
!
      If (.Not. opt('full inv') .And. .Not. opt('SPARSE  ')) Then
!-------------------------------------------------------------------------
! Constants for O(N) algorithm for matrix inversion
!-------------------------------------------------------------------------
        nlayer = naez/nprincd
        Write (1337, 100) nprincd, nlayer
        Write (1337, 110)
! ignore this test if full inversion is done
        If (.Not. opt('full inv')) Then
          If (nlayer*nprincd/=naez) Then
            Write (6, *) 'NLAYER*NPRINCD ( = ', nlayer*nprincd, &
              ').NE.NAEZ ( = ', naez, ')'
            stop_mark = 1
          End If
        End If

      End If
!----------------------------------------------------------------------------
! STOP IF A DIMENSION ERROR OCCURED
!----------------------------------------------------------------------------
      If (stop_mark>0) Stop 'STOP : Dimension Error.'
100   Format (' NPRINCD  NLAYER', /, 2I8)
110   Format (2(7('-'),'+'), 63('-'))
120   Format (' Dimension and Input Data CHECK')
      Return
    End Subroutine
