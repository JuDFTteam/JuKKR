    Subroutine findgroup(bravais, recbv, rbasis, nbasis, rsymat, rotname, &
      isymindex, nsymat, para, qmtet, qmphi, symunitary, krel, naezd, nembd, &
      nsymaxd)
      Use mod_datatypes, Only: dp
! **********************************************************
! This subroutine finds the rotation matrices that leave the
! real lattice unchanged.
! input:  bravais(i,j)    true bravais lattice vectors
!                         i = x,y,z ; j = A, B, C (a.u.)
!         recbv(i,j)      reciprocal basis vectors
!         rbasis          coordinates of basis atoms
!         nbasis          number of basis atoms
!         rsymat          all 64 rotation matrices.
!         rotname         names for the rotation matrices
! output: nsymat          number of rotations that restore the lattice.
!         ISYMINDEX       index for the symmeties found

! This sub makes all 64 rotations in the basis vectors and bravais
! vectors and checks if the new rotated vectror belongs in the
! lattice. The proper rotation must bring all vectors to a lattice
! vector. Information about the rotations found is printed in the end.
! The array ISYMINDEX holds the numbers of the symmetry operations
! that are stored in array RSYMAT
!----------------------------------------------------------------
! in case of relativistic calculation: take account of
! direction of the magnetic moment specified by (QMTET,QMPHI)
! if the PARA(magnetic) flag is set to .FALSE.

! **********************************************************
      Implicit None
      Integer :: krel
      Integer :: naezd, nembd, nsymaxd
!..
      Integer :: nbasis, nsymat
      Integer :: isymindex(nsymaxd)
      Real (Kind=dp) :: bravais(3, 3), rbasis(3, naezd+nembd)
      Real (Kind=dp) :: rsymat(64, 3, 3), recbv(3, 3)
      Real (Kind=dp) :: qmtet(naezd), qmphi(naezd)
      Logical :: symunitary(nsymaxd), para
!..
!.. Local variables
      Real (Kind=dp) :: r(3, 4), rotrbas(3, naezd+nembd)
      Real (Kind=dp) :: bravais1(3, 3)
      Integer :: i, j, isym, nsym, i0, ia
      Real (Kind=dp) :: mdotmp, mvecq(3, naezd), mvecqp(3, naezd)
      Real (Kind=dp) :: mrotr(3, 3), symdet, summdotmp
      Real (Kind=dp) :: stet, ddot, ddet33, pi
      Character (Len=10) :: rotname(64)
      Character (Len=10) :: char(64)
      Logical :: llatbas, latvec, lbulk
!..................................................................

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      Write (1337, 100)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

      nsym = 0
      pi = 4.E0_dp*atan(1.E0_dp)
      Do isym = 1, nsymaxd
        symunitary(isym) = .True.
      End Do
!     - ---------------------------------
      Do i = 1, 3
        Do j = 1, 3
          bravais1(j, i) = bravais(j, i)
        End Do
      End Do
!     Check for surface mode. If so, set bravais1(3,3) very large, so
!     that only the in-plane symmetries are found.
!     not checked, be careful of z--> -z!

      lbulk = .True.
!     Now check the bravais vectors if they have a z component
      If ((bravais(1,3)==0.E0_dp) .And. (bravais(2,3)==0.E0_dp) .And. (bravais &
        (3,3)==0.E0_dp)) Then
        lbulk = .False.
      End If

      Do isym = 1, 64

!--------------------------------- store rotation matrix

        Do i = 1, 3
          Do j = 1, 3
            mrotr(i, j) = rsymat(isym, i, j)
          End Do
        End Do

        summdotmp = 0E0_dp

        symdet = ddet33(mrotr)

!     rotate bravais lattice vectors

!     In the case of slab/interface geometry look only for
!     symmetry opperations that preserve the z axis..

        If (lbulk .Or. (rsymat(isym,3,3)==1)) Then
!     do rotation only in case bulk or if slab and z axis is restored..


          Do i = 1, 3 ! Loop on bravais vectors
            Do j = 1, 3 ! Loop on coordinates
              r(j, i) = rsymat(isym, j, 1)*bravais1(1, i) + &
                rsymat(isym, j, 2)*bravais1(2, i) + rsymat(isym, j, 3)* &
                bravais1(3, i)
            End Do
          End Do

!     rotate the basis atoms p and take RSYMAT.p - p then
!     find if R = (RSYMAT.bravais + RSYMAT.p - p) belongs to the
!     lattice. This is done by function latvec by checking
!     if R.q = integer (q reciprocal lattice vector)

          llatbas = .True.
          Do ia = 1, nbasis ! Loop on basis atoms
            Do j = 1, 3 ! Loop on coordinates
              rotrbas(j, ia) = rsymat(isym, j, 1)*rbasis(1, ia) + &
                rsymat(isym, j, 2)*rbasis(2, ia) + rsymat(isym, j, 3)*rbasis(3 &
                , ia)

              rotrbas(j, ia) = rotrbas(j, ia) - rbasis(j, ia)
              r(j, 4) = rotrbas(j, ia)
            End Do

            If (.Not. latvec(4,recbv,r)) llatbas = .False.

            If ((krel==1) .And. (.Not. para)) Then
              stet = sin(qmtet(ia)*pi/180E0_dp)
              mvecq(1, ia) = stet*cos(qmphi(ia)*pi/180E0_dp)
              mvecq(2, ia) = stet*sin(qmphi(ia)*pi/180E0_dp)
              mvecq(3, ia) = cos(qmtet(ia)*pi/180E0_dp)

              Call dgemv('N', 3, 3, 1E0_dp, mrotr, 3, mvecq(1,ia), 1, 0E0_dp, &
                mvecqp(1,ia), 1)

              Call dscal(3, real(symdet,kind=dp), mvecqp(1,ia), 1)

              mdotmp = ddot(3, mvecq(1,ia), 1, mvecqp(1,ia), 1)
              summdotmp = summdotmp + mdotmp

            End If

          End Do ! ia=1,nbasis

          If ((krel==1) .And. (.Not. para)) Then
            If (abs(abs(summdotmp)-nbasis)>0.00001E0_dp) Then
              llatbas = .False.
            Else
              If (summdotmp>0.00001E0_dp) Then
                symunitary(nsym+1) = .True.
              Else
                symunitary(nsym+1) = .False.
              End If
            End If
          End If

!     if llatbas=.true. the rotation does not change the lattice

          If (llatbas) Then
            nsym = nsym + 1
            isymindex(nsym) = isym
          End If
        End If ! (LBULK .OR. (RSYMAT(ISYM,3,3).EQ.1) )
      End Do ! isym=1,nmatd
!     nsym symmetries were found
!     the ISYMINDEX array has the numbers of the symmetries found


      nsymat = nsym

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      Write (1337, '(8X,60(1H-))')
      If (lbulk) Then
        Write (1337, 110)
      Else
        Write (1337, 120)
      End If
      Write (1337, 130) nsymat
      Do i = 1, nsymat
        i0 = isymindex(i)
        char(i) = rotname(i0)
      End Do
      nsym = nsymat/5
      Do i = 1, nsym + 1
        i0 = (i-1)*5
        isym = min(5, nsymat-i0)
        Write (1337, 140)(char(j), j=i0+1, i0+isym)
      End Do
      Write (1337, 150)
100   Format (5X, '< FINDGROUP > : Finding symmetry operations', /)
110   Format (8X, '3D symmetries:')
120   Format (8X, 'surface symmetries:')
130   Format (' found for this lattice: ', I2, /, 8X, 60('-'))
140   Format (8X, 5(A10,2X))
150   Format (8X, 60('-'), /)
    End Subroutine
