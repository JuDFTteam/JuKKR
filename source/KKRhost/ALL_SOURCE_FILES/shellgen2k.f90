module mod_shellgen2k

contains

! 23.2.2000/ 27.9.2004 *************************************************
subroutine shellgen2k(icc, natom, rcls, atom, nofgij, iofgij, jofgij, nrot, &
  rsymat, isymindex, rotname, nshell, ratom, nsh1, nsh2, ish, jsh, ijtabsym, &
  ijtabsh, ijtabcalc, iprint, nsheld)
  use :: mod_datatypes, only: dp
  ! **********************************************************************
  ! *    Determines the number of different atomic pairs in a cluster by *
  ! * symmetry considerations, assigning a "shell" pointer (used to set  *
  ! * up the GF matrix), to each representative pair.                    *
  ! *                                                                    *
  ! * NATOM       number of atoms in the cluster                         *
  ! * RCLS(3,*)   atom positions in the cluster                          *
  ! * ATOM(*)     corresponding site index in the unit cell              *
  ! * NROT        actual number of symmetry operations                   *
  ! * RSYMAT      symmetry operation matrices                            *
  ! * ISYMINDEX   symmetry operation pointer                             *
  ! * NSHELD      dimension parameter ( max number of different shells)  *
  ! * IJTABCALC   flag to calculate the pair (I,J) - 1/0 for YES/NO      *
  ! *             (e.g. for impurity calc IJTABCALC(I,J) = 1 - delta_ij) *
  ! * NOFGIJ      total number of ij pairs (equals number of non-zero    *
  ! *             IJTABCALC elements                                     *
  ! * IOFGIJ      cluster indices i for pair ij                          *
  ! * JOFGIJ                      j for pair ij                          *
  ! *                                                                    *
  ! * NSHELL(0)   number of different shells (ij pairs)                  *
  ! * NSHELL(NS)  number of equivalent pairs in shell NS                 *
  ! * NSH1(NS),                                                          *
  ! * NSH2(NS)    site indices i,j of shell (representative pair) NS     *
  ! * ISH/JSH     cluster indices i,j of all NSHELL(NS) equivalent pairs *
  ! described by shell NS                                  *
  ! * IJTABSH     the index of the representative shell NS for G_ij      *
  ! * IJTABSYM    the index of the symmetry operation which brings G(NS) *
  ! *             into G_ij                                              *
  ! * RATOM(3,NS) diference vector R_i(NS) - R_j(NS)                     *
  ! *                                                                    *
  ! **********************************************************************
  implicit none
  ! ..
  ! .. Parameters
  integer :: nshell0
  parameter (nshell0=10000)
  ! ..
  ! .. Scalar arguments
  integer :: icc, nofgij, natom, nrot, iprint, nsheld
  ! ..
  ! .. Array arguments
  integer :: atom(*), isymindex(*), ijtabsym(*), ijtabsh(*), ijtabcalc(*)
  integer :: nshell(0:nsheld), nsh1(*), nsh2(*)
  integer :: ish(nsheld, *), jsh(nsheld, *)
  integer :: iofgij(*), jofgij(*)
  real (kind=dp) :: rcls(3, *), rsymat(64, 3, *)
  real (kind=dp) :: ratom(3, *)
  character (len=10) :: rotname(*)
  ! ..
  ! .. Local scalars
  integer :: ai, aj, i, j, k, ns, nsnew, nsgen, id, isym, ii, ij, igij
  real (kind=dp) :: r1, small
  logical :: lfound
  ! ..
  ! .. Local arrays
  real (kind=dp) :: ri(3), rj(3)
  integer :: nsh1i(:), nsh2i(:), nshelli(:)
  real (kind=dp) :: ratomi(:, :)
  allocatable :: nsh1i, nsh2i, nshelli, ratomi
  ! ..
  ! .. Data statements
  data small/1.0e-10_dp/
  ! ..

  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  write (1337, 120)
  if (iprint>1) call printijtab(natom, ijtabcalc)
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

  if (nsheld>=nshell0) then
    write (6, 110) 'local', 'NSHELL0', nsheld
    stop
  end if

  if (nofgij<=0) then
    write (6, '(A)') '      ICC set to 0'
    write (6, '(A)') '         maybe you should check your input?'
    icc = 0                        ! Bauer Long 2011-10-11
    return
  end if
  allocate (nsh1i(nshell0), nsh2i(nshell0), nshelli(nshell0), stat=ns)
  if (ns/=0) stop '   < shellgen2k > allocate NSHELLI arrays'
  allocate (ratomi(3,nshell0), stat=ns)
  if (ns/=0) stop '   < shellgen2k > allocate RATOMI array'
  ! ======================================================================

  ! --> initialise number of shells found for this cluster, setup the
  ! working arrays NSH1I,NSH2I,NSHELLI,RATOMI and set the number of
  ! new found shells (NSNEW) to zero

  do i = 1, nshell(0)
    nsh1i(i) = nsh1(i)
    nsh2i(i) = nsh2(i)
    nshelli(i) = nshell(i)
    do j = 1, 3
      ratomi(j, i) = ratom(j, i)
    end do
  end do
  nsnew = 0

  ! **********************************************************************
  ! loop over I,J-pairs in cluster
  do igij = 1, nofgij

    ! --> search for a symmetric equivalent pair of atoms, LFOUND takes
    ! on the value false/true if this equivalent pair is found

    i = iofgij(igij)
    j = jofgij(igij)
    ! take only those shells that have atom in them
    if (i/=0 .or. j/=0) then
      ai = atom(i)
      aj = atom(j)
      
      lfound = .false.
      nsgen = nshell(0) + nsnew
      
      ! RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
      do id = 1, nrot
        isym = isymindex(id)
        ! ----------------------------------------------------------------------
        do ii = 1, 3
          ri(ii) = rsymat(isym, ii, 1)*rcls(1, i) + rsymat(isym, ii, 2)*rcls(2, i) + rsymat(isym, ii, 3)*rcls(3, i)
      
          rj(ii) = rsymat(isym, ii, 1)*rcls(1, j) + rsymat(isym, ii, 2)*rcls(2, j) + rsymat(isym, ii, 3)*rcls(3, j)
        end do
        ! ----------------------------------------------------------------------
      
        ! --> search for an equivalent pair within the already generated
        ! shells (1..NSHELL(0)+NSNEW)
      
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ns = 0
        do while ((.not. lfound) .and. (ns<nsgen))
          ns = ns + 1
          ! ----------------------------------------------------------------------
          ! IF ( ( AI.EQ.NSH1I(NS) .AND. AJ.EQ.NSH2I(NS) ).OR.
          ! &              ( AI.EQ.NSH2I(NS) .AND. AJ.EQ.NSH1I(NS) )  ) THEN
          ! Commented out by Phivos Mavropoulos 31 Oct 2008. The problem is that
          ! if (I,J) and (J,I)
          ! are assigned to the same shell, then G(I,J) should be transposed to
          ! obtain G(J,I).
          ! However, this transposition is not performed in account in kkr1b
          ! (subr. tbxccpljij).
          ! There, only the real-space rotations (DSYMLL) are performed to
          ! generate each pair GF from the
          ! representative pair, but the transposition is forgotten. Thus there
          ! are two ways to resolve this:
          ! Either flag the pairs to be transposed, which is is a little faster
          ! but complicated
          ! to program, or do not consider the (I,J) and (J,I) pairs as
          ! belonging to the same shell,
          ! which is done now:
          if ((ai==nsh1i(ns) .and. aj==nsh2i(ns))) then
      
            r1 = (ri(1)-rj(1)+ratomi(1,ns))**2 + (ri(2)-rj(2)+ratomi(2,ns))**2 + &
              (ri(3)-rj(3)+ratomi(3,ns))**2
      
            if (r1<small) then
              lfound = .true.
              nshelli(ns) = nshelli(ns) + 1
              if (ns<=nshell(0)) write (1337, 130) ai, (rcls(ii,i), ii=1, 3), &
                aj, (rcls(ii,j), ii=1, 3), ns
              ish(ns, nshelli(ns)) = i
              jsh(ns, nshelli(ns)) = j
            end if
      
          end if
          ! ----------------------------------------------------------------------
        end do ! while loop ns = 1..nsgen while .not.lfound
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      end do ! id=1, nrot
      ! RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
      
      ! --> if the rotation and the representative pair (shell) that
      ! identify a pair of atoms was found LFOUND=.TRUE. and
      ! the search for a different pair of atoms starts; otherwise
      ! the pair (I,J) requires a new shell
      
      if (.not. lfound) then
        nsnew = nsnew + 1
        if (nsnew+nshell(0)>nshell0) then
          write (6, 110) 'local', 'NSHELL0', nsnew + nshell(0)
          stop
        end if
        if (nsnew+nshell(0)>nsheld) then
          write (6, 110) 'global', 'NSHELD', nsnew + nshell(0)
          stop
        end if
      
        nsh1i(nshell(0)+nsnew) = ai
        nsh2i(nshell(0)+nsnew) = aj
        nshelli(nshell(0)+nsnew) = 1
        ish(nshell(0)+nsnew, 1) = i
        jsh(nshell(0)+nsnew, 1) = j
        do ii = 1, 3
          ratomi(ii, nshell(0)+nsnew) = rcls(ii, j) - rcls(ii, i)
        end do
      end if ! .not. lfound

    end if !(i/=0 .or. j/=0) then

  end do ! igij = 1, nofgij
  ! **********************************************************************

  ! --> test number of shells

  if (nsnew+nshell(0)>nsheld) then
    write (6, 110) 'global', 'NSHELD', nsnew + nshell(0)
    stop
  end if

  ! --> update the argument arrays

  do i = 1, nshell(0) + nsnew
    nsh1(i) = nsh1i(i)
    nsh2(i) = nsh2i(i)
    nshell(i) = nshelli(i)
    do j = 1, 3
      ratom(j, i) = ratomi(j, i)
    end do
  end do

  nshell(0) = nshell(0) + nsnew
  deallocate (nsh1i, nsh2i, nshelli, ratomi, stat=ns)
  if (ns/=0) stop '   < shellgen2k > deallocate arrays'

  ! **********************************************************************

  ! --> scan once again the shells to find the corresponding symmetry
  ! index bringing GS(1..NSHELL(0)) to Gij.
  ! Setup the tables IJTABSH  assigning (I,J) --> NS
  ! IJTABSYM assigning (I,J) --> ISYM
  ! G_ij = D^\dagger(ISYM) * G(NS) * D(ISYM)

  ! **********************************************************************
  do i = 1, natom
    ai = (i-1)*natom
    do j = 1, natom
      ij = ai + j
      ijtabsh(ij) = 0
      ijtabsym(ij) = 0
    end do
  end do
  ! **********************************************************************
  do i = 1, natom
    ai = atom(i)
    do j = 1, natom
      aj = atom(j)
      ! =======================================================================
      do ii = 1, nshell(0)
        ! -----------------------------------------------------------------------
        do id = 1, nrot
          isym = isymindex(id)

          do k = 1, 3
            ri(k) = rsymat(isym, k, 1)*ratom(1, ii) + &
              rsymat(isym, k, 2)*ratom(2, ii) + rsymat(isym, k, 3)*ratom(3, ii &
              )
          end do

          if ((ai==nsh1(ii) .and. aj==nsh2(ii)) .or. (ai==nsh2( &
            ii) .and. aj==nsh1(ii))) then

            r1 = (rcls(1,j)-rcls(1,i)-ri(1))**2 + (rcls(2,j)-rcls(2,i)-ri(2)) &
              **2 + (rcls(3,j)-rcls(3,i)-ri(3))**2

            if (r1<small) then
              ij = (i-1)*natom + j
              ijtabsh(ij) = ii
              ijtabsym(ij) = id
              go to 100
            end if
          end if
        end do
        ! -----------------------------------------------------------------------
100     continue
      end do
      ! =======================================================================
    end do
  end do
  ! ***********************************************************************
  if (iprint<=0) return

  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  write (1337, 140) 'assigned shells and symmetries'
  do i = 1, natom
    ai = (i-1)*natom + j
    do j = 1, natom
      ij = ai + j
      if (ijtabcalc(ij)>0) write (1337, 150) i, j, ijtabsh(ij), ijtabsym(ij), &
        rotname(ijtabsym(ij))
    end do
  end do
  write (1337, 160)
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

110 format (6x, 'Dimension ERROR: please increase the ', a, ' parameter', /, &
    6x, a, ' to a value >=', i5, /)
120 format (9x, '< SHELLGEN2K > : assigning representative pairs', &
    ' (shells) ', /, 26x, 'for the off-diagonal elements Gij', /)
130 format (9x, 'INFO: For the atomic sites   I=', i3, ' :', 3f10.6, /, 9x, &
    29x, 'J=', i3, ' :', 3f10.6, /, 9x, 6x, &
    'an already generated equivalent shell (', i3, ') was found', /)
140 format (13x, 30('-'), /, 13x, a, /, 13x, 30('-'), /, 13x, ' I ', 1x, &
    ' J ', ' | ', 'shell', 4x, 'isym', /, 13x, 30('-'))
150 format (13x, i3, 1x, i3, ' | ', 1x, i4, 4x, i2, 2x, a10)
160 format (13x, 30('-'), /)
end subroutine shellgen2k          ! SUBROUTINE SHELLGEN

! **********************************************************************

subroutine printijtab(natom, ijtab)
  implicit none
  ! ..
  integer :: natom
  integer :: ijtab(*)
  ! ..
  integer :: i, j, ij
  integer :: lgmax
  ! ..
  lgmax = 59
  write (1337, 100, advance='no') &
    '  searched for pairs marked with 1 in the table below'
  do j = 1, min(natom+3, lgmax)
    write (1337, '("-")', advance='no')
  end do
  write (1337, *)
  do i = 1, natom
    write (1337, '(14X,I3," | ")', advance='no') i
    ij = (i-1)*natom
    do j = 1, natom
      write (1337, '(I1)', advance='no') ijtab(ij+j)
    end do
    write (1337, *)
  end do
  write (1337, '(13X,6("-"))', advance='no')
  do j = 1, min(natom+3, lgmax)
    write (1337, '("-")', advance='no')
  end do
  write (1337, '(/)')
  ! ...........................................
100 format (13x, 65('-'), /, 18x, a, /, 13x, 65('-'), /, 13x, '   J |', /, &
    13x, 'I    | 1..NATCLUS', /, 13x, 6('-'))
end subroutine printijtab

end module mod_shellgen2k
