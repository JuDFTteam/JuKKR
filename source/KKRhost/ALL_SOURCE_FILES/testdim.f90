! -------------------------------------------------------------------------------
! SUBROUTINE: TESTDIM
! > @brief Testing the dimension of several arrays
! > @note Jonathan Chico: Some of these tests seem unnecessary with the
! changes done to the
! > inc.p
! -------------------------------------------------------------------------------
subroutine testdim(nspin, naez, nemb, natyp, ins, insref, nref, &
  irns, nlayer, krel, nspind, nprincd, knosph, irnsd, korbit)

  implicit none

  integer, intent (in) :: ins      ! < 0 (MT), 1(ASA), 2(Full Potential)
  integer, intent (in) :: naez     ! < Number of atoms in unit cell
  integer, intent (in) :: nref     ! < Number of diff. ref. potentials
  integer, intent (in) :: krel     ! < Switch for
                                   ! non-relativistic/relativistic (0/1)
                                   ! program. Attention: several other
                                   ! parameters depend explicitly on KREL,
                                   ! they are set automatically Used for Dirac
                                   ! solver in ASA
  integer, intent (in) :: nspin    ! < Counter for spin directions
  integer, intent (in) :: natyp    ! < Number of kinds of atoms in unit cell
  integer, intent (in) :: irnsd
  integer, intent (in) :: insref   ! < INS for reference pot. (usual 0)
  integer, intent (in) :: knosph   ! < switch for spherical/non-spherical
                                   ! (0/1) program.
  integer, intent (in) :: korbit   ! < Spin-orbit/non-spin-orbit (1/0) added
                                   ! to the Schroedinger or SRA equations.
                                   ! Works with FP. KREL and KORBIT cannot be
                                   ! both non-zero.
  integer, intent (in) :: nspind   ! < KREL+(1-KREL)*(NSPIN+1)
  integer, intent (in) :: nprincd  ! < Number of principle layers, set to a
                                   ! number >= NRPINC in output of main0
  ! .. In/Out variables
  integer, intent (inout) :: nemb  ! < Number of 'embedding' positions
  integer, intent (inout) :: nlayer ! < Number of principal layer
  integer, dimension (natyp), intent (inout) :: irns ! < Position of atoms in
                                                     ! the unit cell in units
                                                     ! of bravais vectors

  integer :: stop_mark
  integer :: i, j
  logical :: test, opt
  external :: test, opt

  ! ---> dimension tests

  write (1337, 120)

  stop_mark = 0
  if ((nspin>nspind) .and. (krel==0)) then
    write (6, *) 'There is an inconsistenciy between spin polarised &
      &calculation and relativistic options'
    stop_mark = 1
  end if
  if (max(ins,insref)>knosph) then
    write (6, *) 'Please, change the parameter insd in', ' the inputcard to', &
      max(ins, insref)
    stop_mark = 1
  end if
  j = 1
  do i = 1, natyp
    j = max(j, irns(i))
  end do
  if (ins==0 .and. j>1) then
    write (6, *) 'IRNS(*) is set to 1 in case of ', &
      'spherical potential treatment.'
    do i = 1, natyp
      irns(i) = 1
    end do
    j = 1
  end if

  if (j>irnsd) then
    write (6, *) 'Please, change the parameter irnsd in', ' the inputcard to', &
      j
    stop_mark = 1
  end if

  if (.not. opt('VIRATOMS')) then
    if (nref>natyp) then
      write (6, *) 'There are some inconsistencies in the input file./', &
        ' nref(=', nref, ') is greater than natyp (=', natyp, ').'
      stop_mark = 1
    end if
  end if

  if ((krel==1) .and. (korbit==1)) then
    write (6, *) 'Full relativistic for ASA and new SO solver', 'KREL', krel, &
      'KORBIT', korbit
    stop_mark = 1
  end if

  if (.not. opt('NEWSOSOL') .and. korbit==1) then
    write (6, *) &
      'Option NEWSOSOL not found, change KORBIT in the inputcard from', &
      korbit, 'to 0'
    stop_mark = 1
  end if

  if (opt('NEWSOSOL') .and. korbit==0) then
    write (6, *) 'Using option NEWSOSOL, change KORBIT in the inputcard from', &
      korbit, 'to 1'
    stop_mark = 1
  end if
  ! ----------------------------------------------------------------------------
  ! OPT 'WIRE' is only useful with OPT 'full inv' or
  ! OPT 'SPARSE  ' because of sparsity of
  ! the KKR matrix ( not tridiagonal like for 2D and 3D systems)
  ! ----------------------------------------------------------------------------
  if (opt('WIRE    ') .and. .not. (opt('full inv') .or. opt('SPARSE  '))) then
    write (6, *) 'Use option ''full inv'' or ''SPARSE  '' ', &
      'for WIRE calculation.'
    stop_mark = 1
  end if

  if (opt('COMPLEX ') .and. .not. (opt('EigenV  ') .or. opt( &
    'wfct    ') .or. opt('iso surf'))) then
    write (6, *) 'Use option ''COMPLEX '' only for eigenvalue determination.'
    stop_mark = 1
  end if

  if (test('CONT    ')) then
    nemb = 0
    write (6, *) 'No usage of embedding points. NEMB is set to ', nemb, '.'
  end if

  if (.not. opt('full inv') .and. .not. opt('SPARSE  ')) then
    ! -------------------------------------------------------------------------
    ! Constants for O(N) algorithm for matrix inversion
    ! -------------------------------------------------------------------------
    nlayer = naez/nprincd
    write (1337, 100) nprincd, nlayer
    write (1337, 110)
    ! ignore this test if full inversion is done
    if (.not. opt('full inv')) then
      if (nlayer*nprincd/=naez) then
        write (6, *) 'NLAYER*NPRINCD ( = ', nlayer*nprincd, ').NE.NAEZ ( = ', &
          naez, ')'
        stop_mark = 1
      end if
    end if

  end if
  ! ----------------------------------------------------------------------------
  ! STOP IF A DIMENSION ERROR OCCURED
  ! ----------------------------------------------------------------------------
  if (stop_mark>0) stop 'STOP : Dimension Error.'
100 format (' NPRINCD  NLAYER', /, 2i8)
110 format (2(7('-'),'+'), 63('-'))
120 format (' Dimension and Input Data CHECK')
  return
end subroutine testdim
