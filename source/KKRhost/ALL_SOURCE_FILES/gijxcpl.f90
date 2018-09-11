module mod_gijxcpl

contains

subroutine gijxcpl(ido, naez, rbasis, bravais, linterface, niqcalc, iqcalc, &
  natomimp, rclsimp, atomimp, ijtabcalc, ijtabcalc_i, natomimpd)
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
  use :: mod_datatypes, only: dp
   use mod_getclusnxyz
   use mod_ioinput
  implicit none

  ! Arguments
  integer :: ido, naez, natomimp, niqcalc, natomimpd
  integer :: atomimp(*), ijtabcalc(*), ijtabcalc_i(*), iqcalc(*)
  real (kind=dp) :: bravais(3, 3), rbasis(3, *), rclsimp(3, *)
  logical :: linterface

  ! Locals
  integer :: i, i1, i2, i3, ibr(3), ieqvec, iq, iqs, j, jq, nbr(3), ndim, nn, &
    nout, jqs
  ! .......................................................................
  ! uniquely identify a vector R_j - r_i = (r_j + T_n) - r_i by a
  ! quintuple integer in IVECI2J(5,*):
  ! index  meaning
  ! 1    unit-cell site of atom i (always i - first unit cell)
  ! 2    unit-cell site of atom j
  ! 3..5 the three coeficients of the translation vector
  ! T_n = n1*a1 + n2*a2 + n3*a3
  ! NVECI2J(I) gives the number of identical vectors Rj (initially 1)
  ! in the array IVECI2J(2/3/4/5,I)
  ! IREF(1..NVECI2J(I),I) points to the NVECI2J(I) identical vectors
  ! in the array IVECI2J (one site is kept only once)
  ! .......................................................................
  integer :: nb3max
  integer :: njqcalc
  integer :: iveci2j(:, :), nveci2j(:), iref(:, :), jqcalc(:)
  allocatable :: iveci2j, nveci2j, iref, jqcalc
  real (kind=dp) :: clurad, cluradsq, dq(3), dr(3), drsq, tol
  real (kind=dp) :: cluradxy, cluradxysq, drxysq
  logical :: lspher
  character (len=256) :: uio                               ! NCOLIO=256
  ! ..
  ! .. Externals
  logical :: opt
  external :: opt
  ! ..
  ido = 0

  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  write (1337, 110)
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

  tol = 1.0e-4_dp
  ndim = 3
  if (linterface) ndim = 2

  iq = 0
  clurad = 0e0_dp
  lspher = .true.
  call ioinput('JIJRAD          ', uio, 0, 7, iq)
  if (iq==0) then
    read (unit=uio, fmt=*) clurad
    if (clurad<0e0_dp) clurad = 0e0_dp
    cluradxy = clurad
    if (clurad>0e0_dp) then
      iq = 0
      call ioinput('JIJRADXY        ', uio, 0, 7, iq)
      if (iq==0) then
        read (unit=uio, fmt=*) cluradxy
        if (cluradxy<=0e0_dp) cluradxy = clurad
      end if
    end if
    lspher = (abs(clurad-cluradxy)<tol)
  else
    write (1337, 120)
  end if

  do i = 1, 3
    nbr(i) = 0
  end do
  cluradxysq = max(clurad, cluradxy)
  call getclusnxyz(cluradxysq, bravais, ndim, cluradsq, nbr)

  cluradsq = (clurad+tol)*(clurad+tol)
  cluradxysq = (cluradxy+tol)*(cluradxy+tol)

  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  if (clurad>0e0_dp) then
    if (lspher) then
      write (1337, 130) clurad, ((2*nbr(i)+1), i=1, 3)
    else
      write (1337, 140) cluradxy, clurad, ((2*nbr(i)+1), i=1, 3)
    end if
  else
    write (1337, 150)
  end if
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

  ! --> set the reference I and J sites for J_IJ within the unit cell
  ! tokens JIJSITEI and JIJSITEJ allows one to specify which sites
  ! in the unit cell should be considered in the calculation
  ! (default ALL)
  ! Usage:     JIJSITEx  NUMBER_OF_SITES SITE_INDICES_1_TO_NOS
  ! with x=I/J
  ! Examples:  assume a unit cell with 4 atomic sites
  ! JIJSITEI= 1 1
  ! will take only 1 'I' site having the index 1, as sites 'J'
  ! all the others in the u.c.
  ! JIJSITEI= 2 1 3  JIJSITEJ= 1 4
  ! will take 2 sites 'I', 1st and 3rd, and 1 site 'J', the 4th
  ! hence the calculated pairs will be (1,4),(3,4)

  allocate (jqcalc(naez), stat=iq)
  if (iq/=0) stop '    Allocate JQCALC'
  niqcalc = naez
  njqcalc = niqcalc
  do iq = 1, naez
    iqcalc(iq) = iq
    jqcalc(iq) = iq
  end do
  nn = 0
  call ioinput('JIJSITEI        ', uio, 0, 7, nn)
  if (nn==0) then
    read (unit=uio, fmt=*) niqcalc, (iqcalc(iq), iq=1, niqcalc)
    if (niqcalc<=0) return
    do iq = 1, niqcalc
      if ((iqcalc(iq)<=0) .or. (iqcalc(iq)>naez)) return
    end do
  end if
  call ioinput('JIJSITEJ        ', uio, 0, 7, nn)
  if (nn==0) then
    read (unit=uio, fmt=*) njqcalc, (jqcalc(iq), iq=1, njqcalc)
    if (njqcalc<=0) return
    do iq = 1, njqcalc
      if ((jqcalc(iq)<=0) .or. (jqcalc(iq)>naez)) return
    end do
  end if

  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  if (niqcalc==naez) then
    write (1337, 160) naez
  else
    write (1337, 170) niqcalc, naez
  end if
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

  ! ======================================================================

  ! --> determine the size of the arrays IVECI2J,NVECI2J,IREF

  ! ++++++++++++++++++++++++++++++++++++++++ (selected) sites in unit cell
  nn = 0
  do iqs = 1, niqcalc
    nn = nn + 1
    iq = iqcalc(iqs)
    ! **************************************** (selected) sites in unit cell
    do jqs = 1, njqcalc
      jq = jqcalc(jqs)
      do i = 1, 3
        dq(i) = rbasis(i, jq) - rbasis(i, iq)
      end do
      ! -------------------------------------------------- translation vectors
      do i1 = -nbr(1), nbr(1)
        do i2 = -nbr(2), nbr(2)
          do i3 = -nbr(3), nbr(3)
            ibr(1) = i1
            ibr(2) = i2
            ibr(3) = i3
            do i = 1, 3
              dr(i) = dq(i)
              do j = 1, 3
                dr(i) = dr(i) + real(ibr(j), kind=dp)*bravais(i, j)
              end do
            end do
            drxysq = 0e0_dp
            do i = 1, 2
              drxysq = drxysq + dr(i)*dr(i)
            end do
            drsq = dr(3)*dr(3)

            if ((drxysq<=cluradxysq) .and. (drsq<=cluradsq)) nn = nn + 1
          end do
        end do
      end do
      ! ----------------------------------------------------------------------
    end do
    ! **********************************************************************
  end do
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  nb3max = nn
  allocate (iveci2j(5,nb3max), nveci2j(nb3max), iref(nb3max,nb3max), stat=iq)
  if (iq/=0) stop '    Allocate IVECI2J/NVECI2J/IREF'


  ! --> set the first NAEZ vectors (inside the unit-cell at [0,0,0])

  nn = 0
  do iq = 1, niqcalc
    nn = nn + 1
    nveci2j(nn) = 1
    iref(1, nn) = nn
    iveci2j(1, nn) = iqcalc(iq)
    iveci2j(2, nn) = iqcalc(iq)
    do j = 3, 5
      iveci2j(j, nn) = 0
    end do
  end do
  ! ++++++++++++++++++++++++++++++++++++++++ (selected) sites in unit cell
  do iqs = 1, niqcalc
    iq = iqcalc(iqs)
    ! *************************************** (selected)  sites in unit cell
    do jqs = 1, njqcalc
      jq = jqcalc(jqs)

      ! --> set up vector Rij = R_j - r_i = (r_j+T_n) - r_i = (r_j - r_i) +
      ! T_n
      ! DR  =                               DQ          + T_n

      do i = 1, 3
        dq(i) = rbasis(i, jq) - rbasis(i, iq)
      end do
      ! ================================================== translation vectors
      do i1 = -nbr(1), nbr(1)
        do i2 = -nbr(2), nbr(2)
          do i3 = -nbr(3), nbr(3)
            ibr(1) = i1
            ibr(2) = i2
            ibr(3) = i3
            do i = 1, 3
              dr(i) = dq(i)
              do j = 1, 3
                dr(i) = dr(i) + real(ibr(j), kind=dp)*bravais(i, j)
              end do
            end do

            ! --> calculate Rij(xy)**2 -- DRXYSQ
            ! and       Rij(z)**2  -- DRSQ

            drxysq = 0e0_dp
            do i = 1, 2
              drxysq = drxysq + dr(i)*dr(i)
            end do
            drsq = dr(3)*dr(3)
            if (lspher) drsq = drsq + drxysq

            ! --> TOL <= Rij(xy)**2 <= CLURADXY**2 and
            ! TOL <= Rij(z)**2  <= CLURADZ**2  --> keep the vector Rij by its
            ! beginning and end points and the translation indices IBR(1..3)

            ! ------------------------------------------- TOL <= Rij**2 <=
            ! CLURAD**2
            if ((drxysq<=cluradxysq) .and. (drsq<=cluradsq)) then
              nn = nn + 1
              nveci2j(nn) = 1
              iref(1, nn) = nn
              iveci2j(1, nn) = iq
              iveci2j(2, nn) = jq
              do j = 3, 5
                iveci2j(j, nn) = ibr(j-2)
              end do
            end if
            ! ----------------------------------------------------------------------
          end do
        end do
      end do
      ! ======================================================================
    end do
    ! **********************************************************************
  end do
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! --> array IVECI2J contains now the positions of all NIQCALC clusters
  ! next step eliminates common positions but keeps the information
  ! about the connections I,J in IREF(*,*) to be used in setting up
  ! IJTABCALC array
  ! NOUT is the number of eliminated (repeated) sites

  nout = 0
  ! **********************************************************************
  do i = 1, nn
    ! ======================================================================
    if (nveci2j(i)==1) then

      ! --> check vector IVECI2J(I) only if NVECI2J(I).EQ.1
      ! finding a J for which the R_j is the same, increment NVECI2J(I)
      ! and set NVEVI2(J) to 0
      ! same R_j means same site j, same translation vector (indices 2..5)

      ! ----------------------------------------------------------------------
      do j = i + 1, nn

        if (nveci2j(j)==1) then
          ieqvec = 0
          do iq = 2, 5
            ieqvec = ieqvec + abs(iveci2j(iq,j)-iveci2j(iq,i))
          end do
          ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
          if (ieqvec==0) then
            nveci2j(j) = 0
            nveci2j(i) = nveci2j(i) + 1
            iref(nveci2j(i), i) = j
            nout = nout + 1
          end if
          ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        end if

      end do
      ! ----------------------------------------------------------------------
    end if
    ! ======================================================================
  end do
  ! **********************************************************************

  ! --> get now the actual NATOMIMP cluster positions R_j to be scanned
  ! R_j is obtained from IVECI2J

  natomimp = 0
  ! **********************************************************************
  do i = 1, nn
    ! ======================================================================
    if (nveci2j(i)/=0) then

      do j = 1, 3
        dr(j) = rbasis(j, iveci2j(2,i))
        do iq = 3, 5
          dr(j) = dr(j) + iveci2j(iq, i)*bravais(j, iq-2)
        end do
      end do

      natomimp = natomimp + 1
      if (natomimp>natomimpd) then
        write (6, 180) 'global', 'NATOMIMPD', natomimp
        stop
      end if
      do j = 1, 3
        rclsimp(j, natomimp) = dr(j)
      end do

      atomimp(natomimp) = iveci2j(2, i)
      nveci2j(natomimp) = nveci2j(i)
      do j = 1, nveci2j(natomimp)
        iref(j, natomimp) = iref(j, i)
      end do
    end if
    ! ======================================================================
  end do
  ! **********************************************************************
  ! --> crosscheck -- if something went wrong return with IDO=0

  if ((nn-nout/=natomimp) .or. (natomimp<=0)) go to 100
  if ((naez==1) .and. (natomimp==naez)) go to 100

  ! **********************************************************************
  ! --> set table IJTABCALC(I,J)
  ! IJTABCALC(I,J) = 1 for each I = 1,NIQCALC
  ! and for each J which was previously obtained as
  ! connected to I ( IVECI2J(1,IREF)=I )

  do iq = 1, 2*natomimp
    ijtabcalc(iq) = 0
  end do
  ! ======================================================================
  write (1337, 190)
  do iq = 1, niqcalc
    nn = (iq-1)*natomimp
    nout = 0
    do jq = 1, natomimp
      ! ----------------------------------------------------------------------
      if (jq/=iq) then
        do i = 1, nveci2j(jq)
          if (iveci2j(1,iref(i,jq))==atomimp(iq)) then
            ! IF ( IVECI2J(1,IREF(I,JQ)).EQ.1.AND.
            ! +               IVECI2J(2,IREF(I,JQ)).EQ.1 ) THEN
            ijtabcalc(nn+jq) = 1
            if (opt('NEWSOSOL')) then ! Jijtensor
              ijtabcalc((jq-1)*natomimp+iq) = 1 ! Jijtensor
              ijtabcalc_i(nn+jq) = 1 ! Jijtensor
            end if                 ! Jijtensor
            nout = nout + 1
            ! END IF
          end if
        end do
      end if
      ! ----------------------------------------------------------------------
    end do
    write (1337, 200) iq, nout
  end do
  write (1337, 210)
  ! ======================================================================
  ido = 1
100 continue
  deallocate (iveci2j, nveci2j, iref, stat=iq)
  if (iq/=0) stop '    Deallocate IVECI2J/NVECI2J/IREF'
  deallocate (jqcalc, stat=iq)
  if (iq/=0) stop '    Deallocate JQCALC'
  ! ..
110 format (5x, '< GIJXCPL > : Exchange coupling constants calculation', /)
120 format (6x, 'WARNING: Calculation range JIJRAD missing from your input', &
    /, 6x, '         Default value JIJRAD = 0.0 will be assumed', /)
130 format (6x, 'Range of calculating Jij around each atom', /, 6x, &
    '      spherical cluster of radius :', f7.4, ' (ALAT)', /, /, 6x, &
    'Sites j sought within a parallelipiped', /, 6x, 'of size  (', i3, &
    '*a_1) X (', i3, '*a_2) X (', i3, '*a_3)')
140 format (6x, 'Range of calculating Jij around each atom', /, 6x, &
    '    cylindrical cluster of radius :', f7.4, ' (ALAT)', /, 6x, &
    '                           height :', f7.4, ' (ALAT)', /, /, 6x, &
    'Sites j sought within a parallelipiped', /, 6x, 'of size  (', i3, &
    '*a_1) X (', i3, '*a_2) X (', i3, '*a_3)')
150 format (6x, 'Calculations restricted within the unit cell')
160 format (6x, ' - all of the', i3, ' atoms of the u.c. ', &
    'will be taken into account', /)
170 format (6x, ' - only', i3, ' atom(s) (out of', i3, ') in the u.c.', &
    'will be calculated', /)
180 format (6x, 'Dimension ERROR: please increase the ', a, ' parameter', /, &
    6x, a, ' to a value >=', i5, /)
190 format (8x, 'Jij connections set', /, 10x, 20('-'), /, 11x, ' I ', 3x, &
    'no. of J''S', /, 10x, 20('-'))
200 format (10x, i4, 8x, i6)
210 format (10x, 20('-'), /)
end subroutine gijxcpl

end module mod_gijxcpl
