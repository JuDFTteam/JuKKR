    Subroutine impcheck(atomimp, natomimp, naez, rclsimp, rbasis, bravais, &
      ndim)
      Use mod_datatypes, Only: dp
! **********************************************************************
! * Checking the coordinates and site-index assignments of an impurity *
! * cluster read in from a file                                        *
! * The size of the Bravais lattice to be generated is determined      *
! * dynamically using the radius of the input cluster as reference     *
! *                                                                    *
! * For an input site not belonging to the Bravais lattice the program *
! * stops.                                                             *
! * A wrong site-index (unit cell) assignment is corrected and the     *
! * execution continues.                                               *
! **********************************************************************

      Implicit None
!..
!.. Scalar arguments
      Integer :: naez, natomimp, ndim
!.. Array arguments
      Integer :: atomimp(*)
      Real (Kind=dp) :: bravais(3, 3), rbasis(3, *), rclsimp(3, *)
!..
!.. Local scalars
      Integer :: i, iatok, iposok, iq, j, n1, n2, n3, nmax, nmaxz
      Real (Kind=dp) :: diff, rmaxclus, rmaxgen
      Character (Len=6) :: strat, strpos
!..
!.. Local arrays 
      Integer :: ain(natomimp), nbr(3)
      Real (Kind=dp) :: rclsnew(3, natomimp), vec1(3), vec2(3)
      Logical :: latom(natomimp), lpos(natomimp), labscord
      Real (Kind=dp) :: diffmin(natomimp)
!..
!.. External subroutine
      External :: getclusnxyz
!..

! ----------------------------------------------------------------------
!     initialize diffmin array with high value
      diffmin(:) = 1E+5_dp


!     LABSCORD - cluster coordinates are absolute atomic positions

      labscord = .False.
      j = 0
      Do While (j<3 .And. .Not. labscord)
        j = j + 1
        If (abs(rclsimp(j,1))>1E-8_dp) labscord = .True.
      End Do

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! --> set cluster coordinates in absolute atomic positions
!     RCLSNEW(3,*)

      Do i = 1, natomimp
        Call dcopy(3, rclsimp(1,i), 1, rclsnew(1,i), 1)
      End Do
      If (.Not. labscord) Then
        iq = atomimp(1)
        Do i = 1, natomimp
          Call daxpy(3, 1E0_dp, rbasis(1,iq), 1, rclsnew(1,i), 1)
        End Do
      End If

! --> determine the maximum radius of the input cluster
!     this will be then compared to the maximum generated radius
!     when testing the positions -- setting NMAX for generating the
!     lattice

      rmaxclus = 0E0_dp
      Do i = 2, natomimp
        diff = 0E0_dp
        Do j = 1, 3
          diff = diff + (rclsnew(j,i)-rclsnew(j,1))**2
        End Do
        diff = sqrt(diff)
        rmaxclus = max(rmaxclus, diff)
      End Do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      Do i = 1, 3
        nbr(i) = 0
      End Do
      rmaxclus = 3.0_dp*rmaxclus
!       RMAXCLUS = 1.5*RMAXCLUS
      Call getclusnxyz(rmaxclus, bravais, ndim, diff, nbr)
      nmax = max(nbr(1), nbr(2), nbr(3))
      nmaxz = nmax
      If (ndim==2) nmaxz = 0
      rmaxgen = 0E0_dp
      Write (1337, *) nmax, nmaxz

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Write (1337, *) 'rbasis of impurity (in cartesian coordinate)'
      Do i = 1, natomimp
        Write (1337, *)(rclsnew(j,i), j=1, 3)
        ain(i) = atomimp(i)
        lpos(i) = .False.
        latom(i) = .True.
!=======================================================================
        Do n1 = -nmax, nmax
          Do n2 = -nmax, nmax
            Do n3 = -nmaxz, nmaxz

              Do j = 1, 3
                vec1(j) = real(n1, kind=dp)*bravais(j, 1) + &
                  real(n2, kind=dp)*bravais(j, 2) + real(n3, kind=dp)*bravais( &
                  j, 3)
              End Do

!-----------------------------------------------------------------------
              Do iq = 1, naez
                diff = 0E0_dp
                Do j = 1, 3
                  vec2(j) = vec1(j) + rbasis(j, iq)
                  diff = diff + vec2(j)**2
                End Do
                rmaxgen = max(rmaxgen, sqrt(diff))

                diff = sqrt((rclsnew(1,i)-vec2(1))**2+(rclsnew(2, &
                  i)-vec2(2))**2+(rclsnew(3,i)-vec2(3))**2)

                If (diff<=(1E-5_dp)) Then
                  atomimp(i) = iq
                  If (ain(i)/=iq) latom(i) = .False.
                  lpos(i) = .True.
                  Go To 100
                End If

                If (diff<=diffmin(i)) diffmin(i) = diff

              End Do
!-----------------------------------------------------------------------
            End Do
          End Do
        End Do
!=======================================================================
100   End Do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Write (1337, 110) rmaxclus, rmaxgen
      Write (1337, 120) 'Input data for impurity sites - consistency check'


! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      iatok = 0
      iposok = 0
      Open (58, File='rimp.dat')
      Write (58, *) natomimp, ' NATOMIMP'
      Do i = 1, natomimp
        If (lpos(i)) Then
          Write (strpos, '(A6)') 'OK'
        Else
          Write (strpos, '(A6)') 'neq BL'
          iposok = iposok + 1
          Write (*, *) 'minimal difference for atom', i, '=', diffmin(i)
        End If
        Write (strat, '(I3,A3)') atomimp(i), ' <?'
        If (.Not. latom(i)) Then
          iatok = iatok + 1
        Else If (lpos(i)) Then
          Write (strat, '(I3)') atomimp(i)
        End If
        Write (1337, 130) i, (rclsimp(j,i), j=1, 3), ain(i), strpos, strat
        Write (58, Fmt='(I6,3E16.8,I6)') i, (rclsimp(j,i), j=1, 3), ain(i)
      End Do
      Close (58)
      Write (1337, 140)

      If (iposok/=0) Then
        Write (6, 150)
        Stop
      End If

      Do i = 1, natomimp
        If ((atomimp(i)>naez) .Or. (atomimp(i)==0)) Then
          Write (6, 180) i, atomimp(i)
          Stop
        End If
      End Do

      If (iatok/=0) Then
        Write (1337, 160)
      Else
        Write (1337, 170)
      End If
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

110   Format (12X, 'input-cluster R  : ', F11.6, /, 12X, &
        'test-cluster  R  : ', F11.6, /)
120   Format (13X, 63('-'), /, 15X, A, /, 13X, 63('-'), /, 13X, &
        ' imp |               READ IN DATA         host  |', ' CHECKED DATA ', &
        /, 13X, 'index|       x           y           z    site  |', &
        '  pos.   site ', /, 13X, 63('-'))
130   Format (13X, I3, 2X, '|', 3(F12.6), 1X, I4, 1X, '|', A6, 2X, A6)
140   Format (13X, 63('-'))
150   Format (/, 6X, 'ERROR: At least one of your input sites does not', &
        ' belong to the Bravais ', /, 13X, 'lattice (neq BL). ', &
        'Please check your input file', /)
160   Format (13X, 'WARNING: At least one inconsistent assignment of ', &
        'site indices', /, 13X, &
        '         was found in your input. The program will', ' override the', &
        /, 13X, '         input data. Crosscheck?  ', /, 13X, 63('-'), /)
170   Format (13X, 'Your cluster data is consistent', /, 13X, 63('-'), /)
180   Format (/, 6X, 'ERROR: Wrong assignment of impurity site ', I3, &
        ' to the unit-cell site ', I3, /, 13X, 'Please check your input file', &
        /)
    End Subroutine
