module mod_gfshells

contains

subroutine gfshells(icc, natomimp, nsh1, nsh2, ijtabsym, ijtabsh, ijtabcalc, &
  iofgij, jofgij, nofgij, ish, jsh, nshell, naez, natyp, noq, rbasis, bravais, &
  ifilimp, ratom, rclsimp, nsymat, isymindex, rsymat, kaoez, atomimp, rotname, &
  hostimp, lmaxd, lmmaxd, naezd, natypd, natomimpd, nembd, nsheld)
  ! **********************************************************************
  ! *                                                                    *
  ! * This subroutine constructs mainly the index arrays                 *
  ! * NSHELL, NSH1, NSH2 -- NSHELL(0) number of different GF blocks that *
  ! * have to be calculated, NSH1(I),NSH2(I) the sites connected for     *
  ! * the block I, I = 1,NSHELL(0)                                       *
  ! *                                                                    *
  ! **********************************************************************
  use :: mod_types, only: t_imp
  use :: mod_datatypes, only: dp
   use mod_impcheck
   use mod_impcoefs
   use mod_shellgen2k
  implicit none
  real (kind=dp), parameter :: eps = 1e-14_dp
  integer :: lmaxd, lmmaxd, naezd, natypd, natomimpd, nembd, nsheld
  ! ..
  ! .. Scalar arguments
  integer :: icc, naez, natomimp, natyp, nsymat, nofgij
  character (len=40) :: ifilimp
  ! ..
  ! .. Array arguments
  character (len=10) :: rotname(64)
  ! ..
  integer :: atomimp(natomimpd), hostimp(0:natypd)
  integer :: isymindex(*), kaoez(natypd, naezd+nembd)
  integer :: noq(naezd), nsh1(*), nsh2(*), nshell(0:nsheld)
  integer :: ish(nsheld, *), jsh(nsheld, *)
  integer :: ijtabsym(*), ijtabsh(*), ijtabcalc(*), iofgij(*), jofgij(*)
  ! ..
  real (kind=dp) :: bravais(3, 3), ratom(3, nsheld)
  real (kind=dp) :: rbasis(3, *), rclsimp(3, natomimpd)
  real (kind=dp) :: rsymat(64, 3, 3)
  ! ..
  ! .. Local scalars
  integer :: nb, i, j, pos, ii, io, ns, in, ndim, nsize, ihost, ierr
  character (len=9) :: str9
  logical :: lsurf
  integer :: nofgij_with_diag
  ! ..
  ! .. External subroutines
  logical :: opt
  external :: opt

  write (1337, 100)

  nsize = natomimpd*lmmaxd

  ! **********************************************************************

  ! --> construction of ratom, nsh1 and nsh2 for a self-consistent
  ! calculation

  if (.not. opt('VIRATOMS')) then
    nshell(0) = natyp
  else
    nshell(0) = naez
  end if                           ! ( .not. OPT('VIRATOMS') ) THEN

  if (nshell(0)>nsheld) then
    write (6, 110) 'NSHELD', nshell(0)
    stop
  end if

  do i = 1, nshell(0)
    ratom(1, i) = 0.0d0
    ratom(2, i) = 0.0d0
    ratom(3, i) = 0.0d0
    nshell(i) = 0

    do j = 1, naez
      do io = 1, noq(j)
        if (kaoez(io,j)==i) then
          nshell(i) = nshell(i) + 1
          if (nshell(i)==1) then
            nsh1(i) = j
            nsh2(i) = j
          end if
        end if
      end do
    end do
    if (opt('VIRATOMS')) then
      nshell(i) = 1
      nsh1(i) = i
      nsh2(i) = i
    end if



    if (nshell(i)==0) then
      write (6, 120)
      stop
    end if
  end do

  if (icc==0) then
    write (1337, 130) nshell(0)
    return
  end if

  ! end of simple SCF-calculation part.
  ! **********************************************************************

  ! heck if we are in surface mode

  lsurf = .false.
  if (abs(bravais(1,3))<eps .and. abs(bravais(2,3))<eps .and. abs(bravais(3, &
    3))<eps) lsurf = .true.
  ndim = 3
  if (lsurf) ndim = 2

  ! **********************************************************************
  ! NATOMIMP=0   ! BUG: This initialization breaks the shell generation for
  ! ! ICC=-1, which is set by option XCPL.  B. Zimmermann

  if (icc<0) then

    ! --->  ICC.LT.1 all shells are (should be) prepared

    write (1337, 210) natomimp
  else

    ! --> read-in the cluster coordinates from an external file

    rewind 25
    read (25, fmt=*) natomimp

    if (natomimp>natomimpd) then
      write (6, 110) 'NATOMIMPD', natomimp
      stop
    end if
    write (1337, 140) ifilimp, natomimp

    do i = 1, natomimp
      read (25, fmt=*)(rclsimp(j,i), j=1, 3), atomimp(i)
      atomimp(i) = atomimp(i) + icc - 1
    end do

    if (opt('GREENIMP') .or. opt('OPERATOR')) then
      ihost = 0
      outer: do i = 1, natypd
        inner: do j = 1, natomimp
          if (atomimp(j)==i) then
            ihost = ihost + 1
            hostimp(ihost) = atomimp(j)
            cycle outer
          end if
        end do inner
      end do outer

      ! save stuff to t_imp for later use
      t_imp%ihost = ihost
      t_imp%natomimp = natomimp
      allocate (t_imp%hostimp(ihost), stat=ierr)
      if (ierr/=0) stop 'Error allocating t_imp%HOSTIMP'
      t_imp%hostimp(1:ihost) = hostimp(1:ihost)
      allocate (t_imp%atomimp(natomimp), stat=ierr)
      if (ierr/=0) stop 'Error allocating t_imp%ATOMIMP'
      t_imp%atomimp(1:natomimp) = atomimp(1:natomimp)

    end if                         ! GREENIMP

  end if                           ! ICC>=0
  ! **********************************************************************

  call impcheck(atomimp, natomimp, naez, rclsimp, rbasis, bravais, ndim)

  ! **********************************************************************
  if (icc>0) then
    write (1337, 150)

    ! --> set up the number of all (I,J)-pairs to be looked for,
    ! avoid considering again the diagonal elements

    nofgij = 0
    nofgij_with_diag = 0
    do i = 1, natomimp
      nb = (i-1)*natomimp
      do j = 1, natomimp
        ijtabcalc(nb+j) = 0
      end do
      if (atomimp(i)>=0) then
        do j = 1, natomimp
          if ((atomimp(j)>=0) .and. (i/=j)) then
            nofgij = nofgij + 1
            if (nofgij>natomimpd*natomimpd) then
              write (6, 110) 'NATOMIMPD', nofgij/natomimp
              stop
            end if
            iofgij(nofgij) = i
            jofgij(nofgij) = j
            ijtabcalc(nb+j) = 1
          end if
          ! increment counter that includes diagonal elements
          if ((atomimp(j)>=0)) then
            nofgij_with_diag = nofgij_with_diag + 1
          end if
        end do
      end if
    end do
  end if
  ! **********************************************************************

  call shellgen2k(icc, natomimp, rclsimp(1,1), atomimp(1), nofgij, iofgij, &
    jofgij, nsymat, rsymat, isymindex, rotname, nshell, ratom(1,1), nsh1, &
    nsh2, ish, jsh, ijtabsym, ijtabsh, ijtabcalc, 2, nsheld)

  ! after shells have been created reset nofgij to nofgij_with_diag.
  ! Otherwise a segmentation fault occurs in kkrflex (in rotgll: ijtabsh etc too small)
  if (icc>0) nofgij = nofgij_with_diag

  ! **********************************************************************

  ! --> now write out the impurity.coefs file for the impurity calculation
  ! n.papanikolaou

  if (icc>0 .or. opt('KKRFLEX ')) call impcoefs(natomimp, naez, atomimp, &
    rclsimp, nshell, nsh1, nsh2, ratom, nsymat, isymindex, rotname, hostimp, &
    natypd, lmaxd, nsheld, nsize)

  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  write (1337, 130) nshell(0)
  write (1337, 160)
  nb = max(natyp, naez)
  do ns = 1, nshell(0)
    if (ns==nb+1) write (1337, 220)
    if (ns<=nb) then
      call setpairstr(nsh1(ns), nsh2(ns), str9)
      write (1337, 170) ns, nsh1(ns), nsh2(ns), (ratom(ii,ns), ii=1, 3), &
        sqrt(ratom(1,ns)**2+ratom(2,ns)**2+ratom(3,ns)**2), str9
    else
      write (1337, 180, advance='no') ns, nsh1(ns), nsh2(ns), &
        (ratom(ii,ns), ii=1, 3), sqrt(ratom(1,ns)**2+ratom(2,ns)**2+ratom(3,ns &
        )**2)
      io = min(2, nshell(ns))
      do i = 1, io
        call setpairstr(ish(ns,i), jsh(ns,i), str9)
        write (1337, '(A9)', advance='no') str9
      end do
      write (1337, *)
      pos = (nshell(ns)+1)/2
      do i = 2, pos
        io = (i-1)*2
        in = min(2, nshell(ns)-io)
        write (1337, 190, advance='no')
        do j = 1, in
          call setpairstr(ish(ns,io+j), jsh(ns,io+j), str9)
          write (1337, '(A9)', advance='no') str9
        end do
        write (1337, *)
      end do
    end if
  end do
  write (1337, '(6X,72("-"))')
  nb = 0
  do ns = 1, nshell(0)
    nb = nb + nshell(ns)
  end do
  write (1337, 200) nb
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

  ! ----------------------------------------------------------------------
100 format (5x, '< GFSHELLS > : setting up indices of the GF blocks', /)
110 format (6x, 'Dimension ERROR: please increase the global parameter', /, &
    6x, a, ' to a value >=', i5, /)
120 format (6x, 'ERROR: there are some inconsistencies in your input', /, 13x, &
    'not all atoms defined by NATYP have been found', /)
130 format (8x, 'Different shells for GF calculation : ', i3, /)
140 format (8x, 'Reading in cluster impurity sites from file', /, 12x, &
    'file name        : ', a, /, 12x, 'atoms in cluster : ', i3)
150 format (8x, 'Preparing indexing for impurity GF', /, 11x, &
    '- unsymmetrised GF is written out (v. 20.09.2001)', /, 11x, &
    '- files that will be created: impurity.coefs', /, 41x, 'intercell_ref', &
    /, 41x, 'green', /)
160 format (6x, 72('-'), /, 6x, 'shell|', ' IQ ', ' JQ', ' | ', 10x, &
    'vec R_IJ ', 11x, 'R_IJ   | equiv. pairs', /, 6x, 72('-'))
170 format (5x, i5, ' |', i3, 1x, i3, ' | ', 3f9.4, f9.5, 1x, '|', a9)
180 format (5x, i5, ' |', i3, 1x, i3, ' | ', 3f9.4, f9.5, 1x, '|')
190 format (5x, 5x, ' |', 7x, ' | ', 27x, 9x, 1x, '|')
200 format (8x, 'Number of block elements to be calculated : ', i3, /)
210 format (8x, 'Setting pairs for task-defined cluster sites ', &
    'and connections', /, 12x, 'atoms in cluster : ', i3)
220 format (6x, 72(':'), /, 22x, '(impurity) cluster related data/indexing', &
    /, 6x, 72(':'))

end subroutine gfshells
! **********************************************************************

subroutine setpairstr(i, j, str9)
  implicit none
  character (len=9) :: str9, strd
  integer :: i, j, l, lstr
  character (len=20) :: fmt1
  ! ..
  fmt1 = '("(",I'
  fmt1 = fmt1(1:6) // '1'
  lstr = 4
  if (i>=10) then
    fmt1 = fmt1(1:6) // '2'
    lstr = lstr + 1
    if (i>=100) then
      fmt1 = fmt1(1:6) // '3'
      lstr = lstr + 1
    end if
  end if
  fmt1 = fmt1(1:7) // ',",",I'
  fmt1 = fmt1(1:13) // '1'
  lstr = lstr + 1
  if (j>=10) then
    fmt1 = fmt1(1:13) // '2'
    lstr = lstr + 1
    if (j>=100) then
      fmt1 = fmt1(1:13) // '3'
      lstr = lstr + 1
    end if
  end if
  fmt1 = fmt1(1:14) // ',")")'
  write (strd, fmt1) i, j
  do l = 1, 9 - lstr
    str9(l:l) = ' '
  end do
  str9 = str9(1:9-lstr) // strd(1:lstr)
end subroutine setpairstr
! **********************************************************************

end module mod_gfshells
