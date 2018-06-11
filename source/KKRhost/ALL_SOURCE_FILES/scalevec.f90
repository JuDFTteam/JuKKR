    Subroutine scalevec(lcartesian, rbasis, abasis, bbasis, cbasis, nlbasis, &
      nrbasis, nleft, nright, zperleft, zperight, tleft, tright, linterface, &
      naez, nemb, bravais, kaoez, noq, naezd, natypd, nembd)
      Use mod_datatypes, Only: dp

      Implicit None
      real(kind=dp), parameter :: eps=1E-14_dp
!..
!.. Arguments ..
      Integer :: naezd, natypd, nembd
      Integer :: naez, nemb, nlbasis, nleft, nrbasis, nright
      Real (Kind=dp) :: abasis, bbasis, cbasis
      Logical :: linterface
      Integer :: kaoez(natypd, *), noq(naezd)
      Real (Kind=dp) :: bravais(3, 3), rbasis(3, *), tleft(3, *), &
        tright(3, *), zperight(3), zperleft(3)
!..
!.. Locals ..
      Integer :: i, i1, j
      Logical :: lcartesian
      Real (Kind=dp) :: rbasis1(3, naezd+nembd), temp(3), tx, ty, tz

      Write (1337, '(79("="))')
      Write (1337, '(23X,A)') 'SCALEVEC: scale site coordinates'
      Write (1337, '(23X,A)') '          bring all to CARTESIAN system'
      Write (1337, '(79("="))')
      Write (1337, *)

! -->   normalization of basis vectors
!       multiplication instead of division 04/2004

      Do i = 1, naez + nemb
        rbasis1(1, i) = rbasis(1, i)*abasis
        rbasis1(2, i) = rbasis(2, i)*bbasis
        rbasis1(3, i) = rbasis(3, i)*cbasis
      End Do

      If (linterface) Then
        Do i = 1, nlbasis
          tleft(1, i) = tleft(1, i)*abasis
          tleft(2, i) = tleft(2, i)*bbasis
          tleft(3, i) = tleft(3, i)*cbasis
        End Do
        zperleft(1) = zperleft(1)*abasis
        zperleft(2) = zperleft(2)*bbasis
        zperleft(3) = zperleft(3)*cbasis

        Do i = 1, nrbasis
          tright(1, i) = tright(1, i)*abasis
          tright(2, i) = tright(2, i)*bbasis
          tright(3, i) = tright(3, i)*cbasis
        End Do
        zperight(1) = zperight(1)*abasis
        zperight(2) = zperight(2)*bbasis
        zperight(3) = zperight(3)*cbasis
      End If

      If (abs(abasis)<eps .Or. abs(bbasis)<eps .Or. abs(cbasis)<eps) Then
        Write (1337, '(5X,A,2(/,34X,F12.8,A))') &
          'Scaling site coordinates with:', abasis, '  x', bbasis, '  y'
        If (.Not. linterface) Write (1337, '(34X,F12.8,A)') cbasis, '  z'
        Write (1337, '(5X,44("-"))')
      Else
        Write (1337, '(5X,A)') 'Site coordinates will not be scaled'
      End If

! ---> normalization of atomic positions in the unit cell

!      if lcartesian is true cartesian coordinates are used
!      else the basis atoms are in units of the lattice vectors

      If (lcartesian) Then
        Write (1337, '(A)') ' CARTESIAN coordinates'
      Else
        Write (1337, '(A)') ' LATTICE VECTOR coordinates will be', &
          ' changed to CARTESIAN coordinates'
      End If

!**********************************************************************
! Change to cartesian coordinates
      If (linterface) Then
!======================================================================
        If (.Not. lcartesian) Then
!----------------------------------------------------------------------
          Write (1337, *)
          Write (1337, '(12X,49("-"))')
          Write (1337, '(13X,A)') &
            'Input positions transformed to CARTESIAN system'
          Write (1337, '(12X,49("-"),/,13X,A,/,12X,49("-"))') &
            'IQ        x             y             z        IT'
          Do i = 1, naez + nemb
            Do j = 1, 2
              rbasis(j, i) = (rbasis1(1,i)*bravais(j,1)+rbasis1(2,i)*bravais(j &
                ,2))
            End Do
            rbasis(3, i) = rbasis1(3, i)

            If (i<=naez) Then
              Write (1337, 140) i, (rbasis(j,i), j=1, 3), &
                (kaoez(j,i), j=1, noq(i))
            Else
              Write (1337, 140) i, (rbasis(j,i), j=1, 3), kaoez(1, i)
            End If
            If (i==naez) Write (1337, '(12X,49("."))')
          End Do
          Write (1337, '(12X,49("-"),/)')
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! -->  Do the same for the boundary vectors

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      left side

          Do i = 1, nlbasis
            Do i1 = 1, 2
              temp(i1) = tleft(i1, i)
            End Do
            Do j = 1, 2
              tleft(j, i) = (temp(1)*bravais(j,1)+temp(2)*bravais(j,2))
            End Do
          End Do

          Do i1 = 1, 2
            temp(i1) = zperleft(i1)
          End Do
          Do j = 1, 2
            zperleft(j) = (temp(1)*bravais(j,1)+temp(2)*bravais(j,2))
          End Do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      right side

          Do i = 1, nrbasis
            Do i1 = 1, 2
              temp(i1) = tright(i1, i)
            End Do
            Do j = 1, 2
              tright(j, i) = (temp(1)*bravais(j,1)+temp(2)*bravais(j,2))
            End Do
          End Do

          Do i1 = 1, 2
            temp(i1) = zperight(i1)
          End Do
          Do j = 1, 2
            zperight(j) = (temp(1)*bravais(j,1)+temp(2)*bravais(j,2))
          End Do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!----------------------------------------------------------------------
        Else
!----------------------------------------------------------------------
          Write (1337, '(42X,A)') '---> No transformation required'
          Do i = 1, 3
            Do j = 1, naez + nemb
              rbasis(i, j) = rbasis1(i, j)
            End Do
          End Do
!----------------------------------------------------------------------
        End If ! IF (.NOT.LCARTESIAN)
!======================================================================

        Write (1337, 110)
        Do i = nleft, 1, -1
          Do i1 = nlbasis, 1, -1
            tx = tleft(1, i1) + (i-1)*zperleft(1)
            ty = tleft(2, i1) + (i-1)*zperleft(2)
            tz = tleft(3, i1) + (i-1)*zperleft(3)
            Write (1337, 100)(i-1)*nlbasis + i1, tx, ty, tz, kaoez(1, naez+i1)
          End Do
        End Do

        Write (1337, 120)
        Do i = 1, naez
          Write (1337, 100) i, (rbasis(i1,i), i1=1, 3), &
            (kaoez(i1,i), i1=1, noq(i))
        End Do

        Write (1337, 130)
        Do i = 1, nright
          Do i1 = 1, nrbasis
            tx = tright(1, i1) + (i-1)*zperight(1)
            ty = tright(2, i1) + (i-1)*zperight(2)
            tz = tright(3, i1) + (i-1)*zperight(3)
            Write (1337, 100)(i-1)*nrbasis + i1, tx, ty, tz, &
              kaoez(1, naez+nlbasis+i1)
          End Do
        End Do
        Write (1337, '(14X,45("-"),/)')
!======================================================================
      Else If (.Not. lcartesian) Then ! Rescale lattice
!----------------------------------------------------------------------
        Write (1337, *)
        Write (1337, '(12X,49("-"))')
        Write (1337, '(13X,A)') &
          'Input positions transformed to CARTESIAN system'
        Write (1337, '(12X,49("-"),/,13X,A,/,12X,49("-"))') &
          'IQ        x             y             z        IT'
        Do i = 1, naez + nemb
          Do j = 1, 3
            rbasis(j, i) = (rbasis1(1,i)*bravais(j,1)+rbasis1(2,i)*bravais(j,2 &
              )+rbasis1(3,i)*bravais(j,3))
          End Do

          If (i<=naez) Then
            Write (1337, 140) i, (rbasis(j,i), j=1, 3), &
              (kaoez(j,i), j=1, noq(i))
          Else
            Write (1337, 140) i, (rbasis(j,i), j=1, 3), kaoez(1, i)
          End If
          If (i==naez .And. nemb>0) Write (1337, '(12X,49("."))')
        End Do
        Write (1337, '(12X,49("-"),/)')
!----------------------------------------------------------------------
      Else
!----------------------------------------------------------------------
        Write (1337, '(42X,A,/)') '---> No transformation required'
        Write (1337, 150)
!     changed by v.Bellini 21/10/99
        Do j = 1, naez + nemb
          Do i = 1, 3
            rbasis(i, j) = rbasis1(i, j)
          End Do

          If (j<=naez) Then
            Write (1337, 100) j, (rbasis(i,j), i=1, 3), &
              (kaoez(i,j), i=1, noq(j))
          Else
            Write (1337, 100) j, (rbasis(i,j), i=1, 3), kaoez(1, j)
          End If
          If (i==naez .And. nemb>0) Write (1337, '(12X,51("."))')
        End Do
!     end of the change
        Write (1337, '(12X,51("-"),/)')
!----------------------------------------------------------------------
!======================================================================
      End If !  IF (.NOT.LINTERFACE )
!**********************************************************************

! FROM NOW ON after < SCALEVEC > RBASIS are the basis vectors
! in units of au/alat in (xyz) reference

!**********************************************************************
100   Format (13X, I5, 3F12.6, 10I3)
110   Format (14X, 45('-'), /, 15X, '     Positions of ALL generated sites ', &
        /, 15X, '   in CARTESIAN coordinates (ALAT units)', /, 14X, 45('-'), &
        /, 15X, 'IQ       x           y           z       IT', /, 15X, &
        '**************** Left  Host ***************')
120   Format (15X, '****************   S L A B  ***************')
130   Format (15X, '**************** Right Host ***************')
140   Format (12X, I3, 3F14.8, 10I3)
150   Format (12X, 51('-'), /, 16X, '    Positions of (ALL) generated sites', &
        /, 16X, '   in CARTESIAN coordinates (ALAT units)', /, 12X, 51('-'), &
        /, 15X, 'IQ       x           y           z       IT', /, 12X, &
        51('-'))
    End Subroutine
