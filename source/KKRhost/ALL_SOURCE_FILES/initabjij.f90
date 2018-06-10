subroutine initabjij(iprint, naez, natyp, natomimp, nofgij, nqcalc, nsmax, &
  nshell, iqcalc, atomimp, ish, jsh, ijtabcalc, ijtabsh, ijtabsym, nijcalc, &
  kijsh, nijmax, nshell0, nsheld)
!   ********************************************************************
!   *  subroutine called by < TBXCCPLJIJ > to set up some auxiliary    *
!   *  arrays allowing the indexing of shells, sites, atomic types     *
!   ********************************************************************

  implicit none

! Arguments
  integer :: iprint, naez, natomimp, natyp, nijmax, nofgij, nqcalc, nsheld, &
    nshell0, nsmax
  integer :: atomimp(*), ijtabcalc(*), ijtabsh(*), ijtabsym(*), iqcalc(*), &
    ish(nsheld, *), jsh(nsheld, *), kijsh(nijmax, nshell0), nijcalc(nshell0), &
    nshell(0:nsheld)

! Locals
  integer :: i1, ia, idone(naez), iqtojq(nijmax), j1, ja, lm1, lm2, ns
  integer :: nidone

! ======================================================================
  do ns = nsmax + 1, nshell(0)
    do i1 = 1, nijmax
      iqtojq(i1) = 0
    end do
! ----------------------------------------------------------------------
    do i1 = 1, nshell(ns)
      ia = atomimp(ish(ns,i1))
      ja = 0
      do j1 = 1, nijcalc(ns)
        if (ia==iqtojq(j1)) then
          ja = 1
          go to 100
        end if
      end do
100   continue
      if (ja==0) then
        nijcalc(ns) = nijcalc(ns) + 1
        if (nijcalc(ns)>nijmax) then
          write (6, 140) 'local', 'NIJMAX', nijcalc(ns)
          stop '       in < TBXCCPLJIJ > '
        end if
        iqtojq(nijcalc(ns)) = ia
        kijsh(nijcalc(ns), ns) = i1
      end if
    end do
  end do
! ======================================================================
  if (iprint<=0) return
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  nidone = 0
  do ia = 1, naez
    idone(ia) = 0
    do i1 = 1, nqcalc
      if (iqcalc(i1)==ia) then
        nidone = nidone + 1
        idone(nidone) = ia
      end if
    end do
  end do

  lm2 = min(25, natomimp)
  write (1337, 150) naez, natyp, natomimp, nofgij, nshell(0), lm2
  do i1 = 1, nidone
    do ia = 1, natomimp
      if (atomimp(ia)==idone(i1)) then
        lm1 = (ia-1)*natomimp
        write (1337, 160) ia, (ijtabcalc(lm1+ja), ja=1, lm2)
        go to 110
      end if
    end do
110 continue
  end do
  write (1337, 170) lm2
  do i1 = 1, nidone
    do ia = 1, natomimp
      if (atomimp(ia)==idone(i1)) then
        lm1 = (ia-1)*natomimp
        write (1337, 160) ia, (ijtabsh(lm1+ja), ja=1, lm2)
        go to 120
      end if
    end do
120 continue
  end do
  write (1337, 180) lm2
  do i1 = 1, nidone
    do ia = 1, natomimp
      if (atomimp(ia)==idone(i1)) then
        lm1 = (ia-1)*natomimp
        write (1337, 160) ia, (ijtabsym(lm1+ja), ja=1, lm2)
        go to 130
      end if
    end do
130 continue
  end do
  lm2 = 0
  do ns = nsmax + 1, nshell(0)
    lm2 = max(lm2, nijcalc(ns))
  end do
  lm2 = min(5, lm2)
  write (1337, 190)
  do ns = nsmax + 1, nshell(0)
    write (1337, 200) ns, (ish(ns,kijsh(i1,ns)), jsh(ns,kijsh(i1,ns)), i1=1, &
      min(nijcalc(ns),lm2))
    write (1337, 210)(atomimp(ish(ns,kijsh(i1,ns))), atomimp(jsh(ns, &
      kijsh(i1,ns))), ijtabsym((ish(ns,kijsh(i1,ns))-1)*natomimp+jsh(ns, &
      kijsh(i1,ns))), i1=1, min(nijcalc(ns),lm2))
  end do
!ccc      WRITE (6,99009) MIN(NATYP,25)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

140 format (6x, 'Dimension ERROR: please increase the ', a, ' parameter', /, &
    6x, a, ' to a value >=', i5, /)
150 format (8x, 60('-'), /, 8x, 'Data used for J_ij calculation:', /, /, 10x, &
    'Number of sites/types i        (NAEZ/NATYP) :', 2(1x,i3), /, 10x, &
    'Number of atoms in the cluster   (NATOMIMP) :', 1x, i3, /, 10x, &
    'Number of ij pairs                 (NOFGIJ) :', 1x, i3, /, 10x, &
    'Number of representative pairs     (NSHELL) :', 1x, i3, /, /, 10x, &
    'ij-pairs calculation table ( 1 = calculated )', /, 10x, 'IA   JA = 1 ..', &
    i3)
160 format (10x, i3, 3x, 25(i3))
170 format (/, 10x, 'ij-shells table ', /, 10x, 'IA   JA = 1 ..', i3)
180 format (/, 10x, 'ij-symmetries table ', /, 10x, 'IA   JA = 1 ..', i3)
190 format (/, 10x, 'effectively calculated pairs/shells', /, 10x, &
    'SHELL   (IAT,JAT) ', /, 10x, 'SHELL   (IQ,JQ - ISYM) ')
200 format (10x, i4, 3x, 5(i3,',',i3,5x))
210 format (10x, 4x, 3x, 5(i3,',',i3,' - ',i2))
220 format (/, 10x, 'effectively calculated type-type pairs (shells)', /, 10x, &
    'IT   JT = 1 ..', i3)
end subroutine
