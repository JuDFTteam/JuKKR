    Subroutine gijxcpl(ido, naez, rbasis, bravais, linterface, niqcalc, &
      iqcalc, natomimp, rclsimp, atomimp, ijtabcalc, ijtabcalc_i, natomimpd)
      Use mod_datatypes, Only: dp
! **********************************************************************
! *                                                                    *
! * In case of tasks requiring Gij blocks calculation, set variables:  *
! *                                                                    *
! * NATOMIMP, RCLSIMP(3,1..NATOMIMP), ATOMIMP(1..NATOMIMP)             *
! * IJTABCALC flag to which pair is needed: I,J --> (I-1)*NATOMIMP + J *
! *           indexing refers to the generated cluster                 *
! * IJTABCALC_I same as IJTABCALC, but only flag to pair <IJ> is  set  *
! *             (in contrast to the pairs <IJ> and <JI> for NEWSOSOL)  *
! *           indexing refers to the generated cluster                 *
! * NIQCALC   number of sites in the first unit cell contained in the  *
! *           cluster                                                  *
! * IQCALC()  correspondence index to the 1..NAEZ sites                *
! * IDO takes on the value 1 or 0 if setting up process was OK or not  *
! *                                                                    *
! * EXCHANGE COUPLING CONSTANTS calculation case                       *
! *                                                                    *
! **********************************************************************
      Implicit None

! Arguments
      Integer :: ido, naez, natomimp, niqcalc, natomimpd
      Integer :: atomimp(*), ijtabcalc(*), ijtabcalc_i(*), iqcalc(*)
      Real (Kind=dp) :: bravais(3, 3), rbasis(3, *), rclsimp(3, *)
      Logical :: linterface

! Locals
      Integer :: i, i1, i2, i3, ibr(3), ieqvec, iq, iqs, j, jq, nbr(3), ndim, &
        nn, nout, jqs
!.......................................................................
!     uniquely identify a vector R_j - r_i = (r_j + T_n) - r_i by a
!     quintuple integer in IVECI2J(5,*):
!     index  meaning
!        1    unit-cell site of atom i (always i - first unit cell)
!        2    unit-cell site of atom j
!        3..5 the three coeficients of the translation vector
!             T_n = n1*a1 + n2*a2 + n3*a3
!     NVECI2J(I) gives the number of identical vectors Rj (initially 1)
!                in the array IVECI2J(2/3/4/5,I)
!     IREF(1..NVECI2J(I),I) points to the NVECI2J(I) identical vectors
!     in the array IVECI2J (one site is kept only once)
!.......................................................................
      Integer :: nb3max
      Integer :: njqcalc
      Integer :: iveci2j(:, :), nveci2j(:), iref(:, :), jqcalc(:)
      Allocatable :: iveci2j, nveci2j, iref, jqcalc
      Real (Kind=dp) :: clurad, cluradsq, dq(3), dr(3), drsq, tol, tolsq
      Real (Kind=dp) :: cluradxy, cluradxysq, drxysq
      Logical :: lspher
      Character (Len=256) :: uio ! NCOLIO=256

      Logical :: opt
!..
!.. Externals
      External :: getclusnxyz, ioinput, opt
!..
      ido = 0

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      Write (1337, 110)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

      tol = 1.0E-4_dp
      tolsq = tol*tol
      ndim = 3
      If (linterface) ndim = 2

      iq = 0
      clurad = 0E0_dp
      lspher = .True.
      Call ioinput('JIJRAD          ', uio, 0, 7, iq)
      If (iq==0) Then
        Read (Unit=uio, Fmt=*) clurad
        If (clurad<0E0_dp) clurad = 0E0_dp
        cluradxy = clurad
        If (clurad>0E0_dp) Then
          iq = 0
          Call ioinput('JIJRADXY        ', uio, 0, 7, iq)
          If (iq==0) Then
            Read (Unit=uio, Fmt=*) cluradxy
            If (cluradxy<=0E0_dp) cluradxy = clurad
          End If
        End If
        lspher = (abs(clurad-cluradxy)<tol)
      Else
        Write (1337, 120)
      End If

      Do i = 1, 3
        nbr(i) = 0
      End Do
      cluradxysq = max(clurad, cluradxy)
      Call getclusnxyz(cluradxysq, bravais, ndim, cluradsq, nbr)

      cluradsq = (clurad+tol)*(clurad+tol)
      cluradxysq = (cluradxy+tol)*(cluradxy+tol)

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      If (clurad>0E0_dp) Then
        If (lspher) Then
          Write (1337, 130) clurad, ((2*nbr(i)+1), i=1, 3)
        Else
          Write (1337, 140) cluradxy, clurad, ((2*nbr(i)+1), i=1, 3)
        End If
      Else
        Write (1337, 150)
      End If
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

! --> set the reference I and J sites for J_IJ within the unit cell
!     tokens JIJSITEI and JIJSITEJ allows one to specify which sites
!     in the unit cell should be considered in the calculation
!     (default ALL)
!     Usage:     JIJSITEx  NUMBER_OF_SITES SITE_INDICES_1_TO_NOS
!     with x=I/J
!     Examples:  assume a unit cell with 4 atomic sites
!           JIJSITEI= 1 1
!        will take only 1 'I' site having the index 1, as sites 'J'
!        all the others in the u.c.
!           JIJSITEI= 2 1 3  JIJSITEJ= 1 4
!        will take 2 sites 'I', 1st and 3rd, and 1 site 'J', the 4th
!        hence the calculated pairs will be (1,4),(3,4)

      Allocate (jqcalc(naez), Stat=iq)
      If (iq/=0) Stop '    Allocate JQCALC'
      niqcalc = naez
      njqcalc = niqcalc
      Do iq = 1, naez
        iqcalc(iq) = iq
        jqcalc(iq) = iq
      End Do
      nn = 0
      Call ioinput('JIJSITEI        ', uio, 0, 7, nn)
      If (nn==0) Then
        Read (Unit=uio, Fmt=*) niqcalc, (iqcalc(iq), iq=1, niqcalc)
        If (niqcalc<=0) Return
        Do iq = 1, niqcalc
          If ((iqcalc(iq)<=0) .Or. (iqcalc(iq)>naez)) Return
        End Do
      End If
      Call ioinput('JIJSITEJ        ', uio, 0, 7, nn)
      If (nn==0) Then
        Read (Unit=uio, Fmt=*) njqcalc, (jqcalc(iq), iq=1, njqcalc)
        If (njqcalc<=0) Return
        Do iq = 1, njqcalc
          If ((jqcalc(iq)<=0) .Or. (jqcalc(iq)>naez)) Return
        End Do
      End If

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      If (niqcalc==naez) Then
        Write (1337, 160) naez
      Else
        Write (1337, 170) niqcalc, naez
      End If
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

! ======================================================================

! --> determine the size of the arrays IVECI2J,NVECI2J,IREF

! ++++++++++++++++++++++++++++++++++++++++ (selected) sites in unit cell
      nn = 0
      Do iqs = 1, niqcalc
        nn = nn + 1
        iq = iqcalc(iqs)
! **************************************** (selected) sites in unit cell
        Do jqs = 1, njqcalc
          jq = jqcalc(jqs)
          Do i = 1, 3
            dq(i) = rbasis(i, jq) - rbasis(i, iq)
          End Do
! -------------------------------------------------- translation vectors
          Do i1 = -nbr(1), nbr(1)
            Do i2 = -nbr(2), nbr(2)
              Do i3 = -nbr(3), nbr(3)
                ibr(1) = i1
                ibr(2) = i2
                ibr(3) = i3
                Do i = 1, 3
                  dr(i) = dq(i)
                  Do j = 1, 3
                    dr(i) = dr(i) + real(ibr(j), kind=dp)*bravais(i, j)
                  End Do
                End Do
                drxysq = 0E0_dp
                Do i = 1, 2
                  drxysq = drxysq + dr(i)*dr(i)
                End Do
                drsq = dr(3)*dr(3)

                If ((drxysq<=cluradxysq) .And. (drsq<=cluradsq)) nn = nn + 1
              End Do
            End Do
          End Do
! ----------------------------------------------------------------------
        End Do
! **********************************************************************
      End Do
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      nb3max = nn
      Allocate (iveci2j(5,nb3max), nveci2j(nb3max), iref(nb3max,nb3max), &
        Stat=iq)
      If (iq/=0) Stop '    Allocate IVECI2J/NVECI2J/IREF'


! --> set the first NAEZ vectors (inside the unit-cell at [0,0,0])

      nn = 0
      Do iq = 1, niqcalc
        nn = nn + 1
        nveci2j(nn) = 1
        iref(1, nn) = nn
        iveci2j(1, nn) = iqcalc(iq)
        iveci2j(2, nn) = iqcalc(iq)
        Do j = 3, 5
          iveci2j(j, nn) = 0
        End Do
      End Do
! ++++++++++++++++++++++++++++++++++++++++ (selected) sites in unit cell
      Do iqs = 1, niqcalc
        iq = iqcalc(iqs)
! *************************************** (selected)  sites in unit cell
        Do jqs = 1, njqcalc
          jq = jqcalc(jqs)

! --> set up vector Rij = R_j - r_i = (r_j+T_n) - r_i = (r_j - r_i) + T_n
!                   DR  =                               DQ          + T_n

          Do i = 1, 3
            dq(i) = rbasis(i, jq) - rbasis(i, iq)
          End Do
! ================================================== translation vectors
          Do i1 = -nbr(1), nbr(1)
            Do i2 = -nbr(2), nbr(2)
              Do i3 = -nbr(3), nbr(3)
                ibr(1) = i1
                ibr(2) = i2
                ibr(3) = i3
                Do i = 1, 3
                  dr(i) = dq(i)
                  Do j = 1, 3
                    dr(i) = dr(i) + real(ibr(j), kind=dp)*bravais(i, j)
                  End Do
                End Do

! --> calculate Rij(xy)**2 -- DRXYSQ
!     and       Rij(z)**2  -- DRSQ

                drxysq = 0E0_dp
                Do i = 1, 2
                  drxysq = drxysq + dr(i)*dr(i)
                End Do
                drsq = dr(3)*dr(3)
                If (lspher) drsq = drsq + drxysq

! --> TOL <= Rij(xy)**2 <= CLURADXY**2 and
!     TOL <= Rij(z)**2  <= CLURADZ**2  --> keep the vector Rij by its
!     beginning and end points and the translation indices IBR(1..3)

! ------------------------------------------- TOL <= Rij**2 <= CLURAD**2
                If ((drxysq<=cluradxysq) .And. (drsq<=cluradsq)) Then
                  nn = nn + 1
                  nveci2j(nn) = 1
                  iref(1, nn) = nn
                  iveci2j(1, nn) = iq
                  iveci2j(2, nn) = jq
                  Do j = 3, 5
                    iveci2j(j, nn) = ibr(j-2)
                  End Do
                End If
! ----------------------------------------------------------------------
              End Do
            End Do
          End Do
! ======================================================================
        End Do
! **********************************************************************
      End Do
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! --> array IVECI2J contains now the positions of all NIQCALC clusters
!     next step eliminates common positions but keeps the information
!     about the connections I,J in IREF(*,*) to be used in setting up
!     IJTABCALC array
!     NOUT is the number of eliminated (repeated) sites

      nout = 0
! **********************************************************************
      Do i = 1, nn
! ======================================================================
        If (nveci2j(i)==1) Then

! --> check vector IVECI2J(I) only if NVECI2J(I).EQ.1
!     finding a J for which the R_j is the same, increment NVECI2J(I)
!     and set NVEVI2(J) to 0
!     same R_j means same site j, same translation vector (indices 2..5)

! ----------------------------------------------------------------------
          Do j = i + 1, nn

            If (nveci2j(j)==1) Then
              ieqvec = 0
              Do iq = 2, 5
                ieqvec = ieqvec + abs(iveci2j(iq,j)-iveci2j(iq,i))
              End Do
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
              If (ieqvec==0) Then
                nveci2j(j) = 0
                nveci2j(i) = nveci2j(i) + 1
                iref(nveci2j(i), i) = j
                nout = nout + 1
              End If
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            End If

          End Do
! ----------------------------------------------------------------------
        End If
! ======================================================================
      End Do
! **********************************************************************

! --> get now the actual NATOMIMP cluster positions R_j to be scanned
!     R_j is obtained from IVECI2J

      natomimp = 0
! **********************************************************************
      Do i = 1, nn
! ======================================================================
        If (nveci2j(i)/=0) Then

          Do j = 1, 3
            dr(j) = rbasis(j, iveci2j(2,i))
            Do iq = 3, 5
              dr(j) = dr(j) + iveci2j(iq, i)*bravais(j, iq-2)
            End Do
          End Do

          natomimp = natomimp + 1
          If (natomimp>natomimpd) Then
            Write (6, 180) 'global', 'NATOMIMPD', natomimp
            Stop
          End If
          Do j = 1, 3
            rclsimp(j, natomimp) = dr(j)
          End Do

          atomimp(natomimp) = iveci2j(2, i)
          nveci2j(natomimp) = nveci2j(i)
          Do j = 1, nveci2j(natomimp)
            iref(j, natomimp) = iref(j, i)
          End Do
        End If
! ======================================================================
      End Do
! **********************************************************************
! --> crosscheck -- if something went wrong return with IDO=0

      If ((nn-nout/=natomimp) .Or. (natomimp<=0)) Go To 100
      If ((naez==1) .And. (natomimp==naez)) Go To 100

! **********************************************************************
! --> set table IJTABCALC(I,J)
!     IJTABCALC(I,J) = 1 for each I = 1,NIQCALC
!                    and for each J which was previously obtained as
!                    connected to I ( IVECI2J(1,IREF)=I )

      Do iq = 1, 2*natomimp
        ijtabcalc(iq) = 0
      End Do
! ======================================================================
      Write (1337, 190)
      Do iq = 1, niqcalc
        nn = (iq-1)*natomimp
        nout = 0
        Do jq = 1, natomimp
! ----------------------------------------------------------------------
          If (jq/=iq) Then
            Do i = 1, nveci2j(jq)
              If (iveci2j(1,iref(i,jq))==atomimp(iq)) Then
!                  IF ( IVECI2J(1,IREF(I,JQ)).EQ.1.AND.
!     +               IVECI2J(2,IREF(I,JQ)).EQ.1 ) THEN
                ijtabcalc(nn+jq) = 1
                If (opt('NEWSOSOL')) Then !Jijtensor
                  ijtabcalc((jq-1)*natomimp+iq) = 1 !Jijtensor
                  ijtabcalc_i(nn+jq) = 1 !Jijtensor
                End If !Jijtensor
                nout = nout + 1
!                  END IF
              End If
            End Do
          End If
! ----------------------------------------------------------------------
        End Do
        Write (1337, 200) iq, nout
      End Do
      Write (1337, 210)
! ======================================================================
      ido = 1
100   Continue
      Deallocate (iveci2j, nveci2j, iref, Stat=iq)
      If (iq/=0) Stop '    Deallocate IVECI2J/NVECI2J/IREF'
      Deallocate (jqcalc, Stat=iq)
      If (iq/=0) Stop '    Deallocate JQCALC'
! ..
110   Format (5X, '< GIJXCPL > : Exchange coupling constants calculation', /)
120   Format (6X, 'WARNING: Calculation range JIJRAD missing from your input', &
        /, 6X, '         Default value JIJRAD = 0.0 will be assumed', /)
130   Format (6X, 'Range of calculating Jij around each atom', /, 6X, &
        '      spherical cluster of radius :', F7.4, ' (ALAT)', /, /, 6X, &
        'Sites j sought within a parallelipiped', /, 6X, 'of size  (', I3, &
        '*a_1) X (', I3, '*a_2) X (', I3, '*a_3)')
140   Format (6X, 'Range of calculating Jij around each atom', /, 6X, &
        '    cylindrical cluster of radius :', F7.4, ' (ALAT)', /, 6X, &
        '                           height :', F7.4, ' (ALAT)', /, /, 6X, &
        'Sites j sought within a parallelipiped', /, 6X, 'of size  (', I3, &
        '*a_1) X (', I3, '*a_2) X (', I3, '*a_3)')
150   Format (6X, 'Calculations restricted within the unit cell')
160   Format (6X, ' - all of the', I3, ' atoms of the u.c. ', &
        'will be taken into account', /)
170   Format (6X, ' - only', I3, ' atom(s) (out of', I3, ') in the u.c.', &
        'will be calculated', /)
180   Format (6X, 'Dimension ERROR: please increase the ', A, ' parameter', /, &
        6X, A, ' to a value >=', I5, /)
190   Format (8X, 'Jij connections set', /, 10X, 20('-'), /, 11X, ' I ', 3X, &
        'no. of J''S', /, 10X, 20('-'))
200   Format (10X, I4, 8X, I6)
210   Format (10X, 20('-'), /)
    End Subroutine
