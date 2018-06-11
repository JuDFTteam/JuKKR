! 23.2.2000/ 27.9.2004 *************************************************
    Subroutine shellgen2k(icc, natom, rcls, atom, nofgij, iofgij, jofgij, &
      nrot, rsymat, isymindex, rotname, nshell, ratom, nsh1, nsh2, ish, jsh, &
      ijtabsym, ijtabsh, ijtabcalc, iprint, nsheld)
      Use mod_datatypes, Only: dp
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
!               described by shell NS                                  *
! * IJTABSH     the index of the representative shell NS for G_ij      *
! * IJTABSYM    the index of the symmetry operation which brings G(NS) *
! *             into G_ij                                              *
! * RATOM(3,NS) diference vector R_i(NS) - R_j(NS)                     *
! *                                                                    *
! **********************************************************************
      Implicit None
!..
!.. Parameters
      Integer :: nshell0
      Parameter (nshell0=10000)
!..
!.. Scalar arguments
      Integer :: icc, nofgij, natom, nrot, iprint, nsheld
!..
!.. Array arguments
      Integer :: atom(*), isymindex(*), ijtabsym(*), ijtabsh(*), ijtabcalc(*)
      Integer :: nshell(0:nsheld), nsh1(*), nsh2(*)
      Integer :: ish(nsheld, *), jsh(nsheld, *)
      Integer :: iofgij(*), jofgij(*)
      Real (Kind=dp) :: rcls(3, *), rsymat(64, 3, *)
      Real (Kind=dp) :: ratom(3, *)
      Character (Len=10) :: rotname(*)
!..
!.. Local scalars
      Integer :: ai, aj, i, j, k, ns, nsnew, nsgen, id, isym, ii, ij, igij
      Real (Kind=dp) :: r1, small
      Logical :: lfound
!..
!.. Local arrays
      Real (Kind=dp) :: ri(3), rj(3)
      Integer :: nsh1i(:), nsh2i(:), nshelli(:)
      Real (Kind=dp) :: ratomi(:, :)
      Allocatable :: nsh1i, nsh2i, nshelli, ratomi
!..
!.. Data statements
      Data small/1.0E-10_dp/
!..

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      Write (1337, 120)
      If (iprint>1) Call printijtab(natom, ijtabcalc)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

      If (nsheld>=nshell0) Then
        Write (6, 110) 'local', 'NSHELL0', nsheld
        Stop
      End If

      If (nofgij<=0) Then
        Write (6, '(A)') '      ICC set to 0'
        Write (6, '(A)') '         maybe you should check your input?'
        icc = 0 ! Bauer Long 2011-10-11
        Return
      End If
      Allocate (nsh1i(nshell0), nsh2i(nshell0), nshelli(nshell0), Stat=ns)
      If (ns/=0) Stop '   < shellgen2k > allocate NSHELLI arrays'
      Allocate (ratomi(3,nshell0), Stat=ns)
      If (ns/=0) Stop '   < shellgen2k > allocate RATOMI array'
! ======================================================================

! --> initialise number of shells found for this cluster, setup the
!     working arrays NSH1I,NSH2I,NSHELLI,RATOMI and set the number of
!     new found shells (NSNEW) to zero

      Do i = 1, nshell(0)
        nsh1i(i) = nsh1(i)
        nsh2i(i) = nsh2(i)
        nshelli(i) = nshell(i)
        Do j = 1, 3
          ratomi(j, i) = ratom(j, i)
        End Do
      End Do
      nsnew = 0

! **********************************************************************
!                                         loop over I,J-pairs in cluster
      Do igij = 1, nofgij

! --> search for a symmetric equivalent pair of atoms, LFOUND takes
!     on the value false/true if this equivalent pair is found

        i = iofgij(igij)
        j = jofgij(igij)
        ai = atom(i)
        aj = atom(j)

        lfound = .False.
        nsgen = nshell(0) + nsnew

! RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
        Do id = 1, nrot
          isym = isymindex(id)
! ----------------------------------------------------------------------
          Do ii = 1, 3
            ri(ii) = rsymat(isym, ii, 1)*rcls(1, i) + &
              rsymat(isym, ii, 2)*rcls(2, i) + rsymat(isym, ii, 3)*rcls(3, i)

            rj(ii) = rsymat(isym, ii, 1)*rcls(1, j) + &
              rsymat(isym, ii, 2)*rcls(2, j) + rsymat(isym, ii, 3)*rcls(3, j)
          End Do
! ----------------------------------------------------------------------

! --> search for an equivalent pair within the already generated
!     shells (1..NSHELL(0)+NSNEW)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ns = 0
          Do While ((.Not. lfound) .And. (ns<nsgen))
            ns = ns + 1
! ----------------------------------------------------------------------
!              IF ( ( AI.EQ.NSH1I(NS) .AND. AJ.EQ.NSH2I(NS) ).OR.
!    &              ( AI.EQ.NSH2I(NS) .AND. AJ.EQ.NSH1I(NS) )  ) THEN
! Commented out by Phivos Mavropoulos 31 Oct 2008. The problem is that if (I,J) and (J,I)
! are assigned to the same shell, then G(I,J) should be transposed to obtain G(J,I).
! However, this transposition is not performed in account in kkr1b (subr. tbxccpljij).
! There, only the real-space rotations (DSYMLL) are performed to generate each pair GF from the
! representative pair, but the transposition is forgotten. Thus there are two ways to resolve this:
! Either flag the pairs to be transposed, which is is a little faster but complicated
! to program, or do not consider the (I,J) and (J,I) pairs as belonging to the same shell,
! which is done now:
            If ((ai==nsh1i(ns) .And. aj==nsh2i(ns))) Then

              r1 = (ri(1)-rj(1)+ratomi(1,ns))**2 + (ri(2)-rj(2)+ratomi(2,ns)) &
                **2 + (ri(3)-rj(3)+ratomi(3,ns))**2

              If (r1<small) Then
                lfound = .True.
                nshelli(ns) = nshelli(ns) + 1
                If (ns<=nshell(0)) Write (1337, 130) ai, &
                  (rcls(ii,i), ii=1, 3), aj, (rcls(ii,j), ii=1, 3), ns
                ish(ns, nshelli(ns)) = i
                jsh(ns, nshelli(ns)) = j
              End If

            End If
! ----------------------------------------------------------------------
          End Do ! NS = 1..NSGEN while .NOT.LFOUND
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        End Do
! RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR

! --> if the rotation and the representative pair (shell) that
!     identify a pair of atoms was found LFOUND=.TRUE. and
!     the search for a different pair of atoms starts; otherwise
!     the pair (I,J) requires a new shell

        If (.Not. lfound) Then
          nsnew = nsnew + 1
          If (nsnew+nshell(0)>nshell0) Then
            Write (6, 110) 'local', 'NSHELL0', nsnew + nshell(0)
            Stop
          End If
          If (nsnew+nshell(0)>nsheld) Then
            Write (6, 110) 'global', 'NSHELD', nsnew + nshell(0)
            Stop
          End If

          nsh1i(nshell(0)+nsnew) = ai
          nsh2i(nshell(0)+nsnew) = aj
          nshelli(nshell(0)+nsnew) = 1
          ish(nshell(0)+nsnew, 1) = i
          jsh(nshell(0)+nsnew, 1) = j
          Do ii = 1, 3
            ratomi(ii, nshell(0)+nsnew) = rcls(ii, j) - rcls(ii, i)
          End Do
        End If

      End Do
! **********************************************************************

! --> test number of shells

      If (nsnew+nshell(0)>nsheld) Then
        Write (6, 110) 'global', 'NSHELD', nsnew + nshell(0)
        Stop
      End If

! --> update the argument arrays

      Do i = 1, nshell(0) + nsnew
        nsh1(i) = nsh1i(i)
        nsh2(i) = nsh2i(i)
        nshell(i) = nshelli(i)
        Do j = 1, 3
          ratom(j, i) = ratomi(j, i)
        End Do
      End Do

      nshell(0) = nshell(0) + nsnew
      Deallocate (nsh1i, nsh2i, nshelli, ratomi, Stat=ns)
      If (ns/=0) Stop '   < shellgen2k > deallocate arrays'

! **********************************************************************

! --> scan once again the shells to find the corresponding symmetry
!     index bringing GS(1..NSHELL(0)) to Gij.
!     Setup the tables IJTABSH  assigning (I,J) --> NS
!                      IJTABSYM assigning (I,J) --> ISYM
!     G_ij = D^\dagger(ISYM) * G(NS) * D(ISYM)

! **********************************************************************
      Do i = 1, natom
        ai = (i-1)*natom
        Do j = 1, natom
          ij = ai + j
          ijtabsh(ij) = 0
          ijtabsym(ij) = 0
        End Do
      End Do
! **********************************************************************
      Do i = 1, natom
        ai = atom(i)
        Do j = 1, natom
          aj = atom(j)
!=======================================================================
          Do ii = 1, nshell(0)
!-----------------------------------------------------------------------
            Do id = 1, nrot
              isym = isymindex(id)

              Do k = 1, 3
                ri(k) = rsymat(isym, k, 1)*ratom(1, ii) + &
                  rsymat(isym, k, 2)*ratom(2, ii) + rsymat(isym, k, 3)*ratom(3 &
                  , ii)
              End Do

              If ((ai==nsh1(ii) .And. aj==nsh2(ii)) .Or. (ai==nsh2( &
                ii) .And. aj==nsh1(ii))) Then

                r1 = (rcls(1,j)-rcls(1,i)-ri(1))**2 + &
                  (rcls(2,j)-rcls(2,i)-ri(2))**2 + (rcls(3,j)-rcls(3,i)-ri(3)) &
                  **2

                If (r1<small) Then
                  ij = (i-1)*natom + j
                  ijtabsh(ij) = ii
                  ijtabsym(ij) = id
                  Go To 100
                End If
              End If
            End Do
!-----------------------------------------------------------------------
100         Continue
          End Do
!=======================================================================
        End Do
      End Do
!***********************************************************************
      If (iprint<=0) Return

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      Write (1337, 140) 'assigned shells and symmetries'
      Do i = 1, natom
        ai = (i-1)*natom + j
        Do j = 1, natom
          ij = ai + j
          If (ijtabcalc(ij)>0) Write (1337, 150) i, j, ijtabsh(ij), &
            ijtabsym(ij), rotname(ijtabsym(ij))
        End Do
      End Do
      Write (1337, 160)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

110   Format (6X, 'Dimension ERROR: please increase the ', A, ' parameter', /, &
        6X, A, ' to a value >=', I5, /)
120   Format (9X, '< SHELLGEN2K > : assigning representative pairs', &
        ' (shells) ', /, 26X, 'for the off-diagonal elements Gij', /)
130   Format (9X, 'INFO: For the atomic sites   I=', I3, ' :', 3F10.6, /, 9X, &
        29X, 'J=', I3, ' :', 3F10.6, /, 9X, 6X, &
        'an already generated equivalent shell (', I3, ') was found', /)
140   Format (13X, 30('-'), /, 13X, A, /, 13X, 30('-'), /, 13X, ' I ', 1X, &
        ' J ', ' | ', 'shell', 4X, 'isym', /, 13X, 30('-'))
150   Format (13X, I3, 1X, I3, ' | ', 1X, I4, 4X, I2, 2X, A10)
160   Format (13X, 30('-'), /)
    End Subroutine ! SUBROUTINE SHELLGEN

! **********************************************************************

    Subroutine printijtab(natom, ijtab)
      Implicit None
!     ..
      Integer :: natom
      Integer :: ijtab(*)
!     ..
      Integer :: i, j, ij
      Integer :: lgmax
!     ..
      lgmax = 59
      Write (1337, 100, Advance='no') &
        '  searched for pairs marked with 1 in the table below'
      Do j = 1, min(natom+3, lgmax)
        Write (1337, '("-")', Advance='no')
      End Do
      Write (1337, *)
      Do i = 1, natom
        Write (1337, '(14X,I3," | ")', Advance='no') i
        ij = (i-1)*natom
        Do j = 1, natom
          Write (1337, '(I1)', Advance='no') ijtab(ij+j)
        End Do
        Write (1337, *)
      End Do
      Write (1337, '(13X,6("-"))', Advance='no')
      Do j = 1, min(natom+3, lgmax)
        Write (1337, '("-")', Advance='no')
      End Do
      Write (1337, '(/)')
!     ...........................................
100   Format (13X, 65('-'), /, 18X, A, /, 13X, 65('-'), /, 13X, '   J |', /, &
        13X, 'I    | 1..NATCLUS', /, 13X, 6('-'))
    End Subroutine
