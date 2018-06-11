    Subroutine initabjij(iprint, naez, natyp, natomimp, nofgij, nqcalc, nsmax, &
      nshell, iqcalc, atomimp, ish, jsh, ijtabcalc, ijtabsh, ijtabsym, &
      nijcalc, kijsh, nijmax, nshell0, nsheld)
!   ********************************************************************
!   *  subroutine called by < TBXCCPLJIJ > to set up some auxiliary    *
!   *  arrays allowing the indexing of shells, sites, atomic types     *
!   ********************************************************************

      Implicit None

! Arguments
      Integer :: iprint, naez, natomimp, natyp, nijmax, nofgij, nqcalc, &
        nsheld, nshell0, nsmax
      Integer :: atomimp(*), ijtabcalc(*), ijtabsh(*), ijtabsym(*), iqcalc(*), &
        ish(nsheld, *), jsh(nsheld, *), kijsh(nijmax, nshell0), &
        nijcalc(nshell0), nshell(0:nsheld)

! Locals
      Integer :: i1, ia, idone(naez), iqtojq(nijmax), j1, ja, lm1, lm2, ns
      Integer :: nidone

! ======================================================================
      Do ns = nsmax + 1, nshell(0)
        Do i1 = 1, nijmax
          iqtojq(i1) = 0
        End Do
! ----------------------------------------------------------------------
        Do i1 = 1, nshell(ns)
          ia = atomimp(ish(ns,i1))
          ja = 0
          Do j1 = 1, nijcalc(ns)
            If (ia==iqtojq(j1)) Then
              ja = 1
              Go To 100
            End If
          End Do
100       Continue
          If (ja==0) Then
            nijcalc(ns) = nijcalc(ns) + 1
            If (nijcalc(ns)>nijmax) Then
              Write (6, 140) 'local', 'NIJMAX', nijcalc(ns)
              Stop '       in < TBXCCPLJIJ > '
            End If
            iqtojq(nijcalc(ns)) = ia
            kijsh(nijcalc(ns), ns) = i1
          End If
        End Do
      End Do
! ======================================================================
      If (iprint<=0) Return
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      nidone = 0
      Do ia = 1, naez
        idone(ia) = 0
        Do i1 = 1, nqcalc
          If (iqcalc(i1)==ia) Then
            nidone = nidone + 1
            idone(nidone) = ia
          End If
        End Do
      End Do

      lm2 = min(25, natomimp)
      Write (1337, 150) naez, natyp, natomimp, nofgij, nshell(0), lm2
      Do i1 = 1, nidone
        Do ia = 1, natomimp
          If (atomimp(ia)==idone(i1)) Then
            lm1 = (ia-1)*natomimp
            Write (1337, 160) ia, (ijtabcalc(lm1+ja), ja=1, lm2)
            Go To 110
          End If
        End Do
110     Continue
      End Do
      Write (1337, 170) lm2
      Do i1 = 1, nidone
        Do ia = 1, natomimp
          If (atomimp(ia)==idone(i1)) Then
            lm1 = (ia-1)*natomimp
            Write (1337, 160) ia, (ijtabsh(lm1+ja), ja=1, lm2)
            Go To 120
          End If
        End Do
120     Continue
      End Do
      Write (1337, 180) lm2
      Do i1 = 1, nidone
        Do ia = 1, natomimp
          If (atomimp(ia)==idone(i1)) Then
            lm1 = (ia-1)*natomimp
            Write (1337, 160) ia, (ijtabsym(lm1+ja), ja=1, lm2)
            Go To 130
          End If
        End Do
130     Continue
      End Do
      lm2 = 0
      Do ns = nsmax + 1, nshell(0)
        lm2 = max(lm2, nijcalc(ns))
      End Do
      lm2 = min(5, lm2)
      Write (1337, 190)
      Do ns = nsmax + 1, nshell(0)
        Write (1337, 200) ns, (ish(ns,kijsh(i1,ns)), jsh(ns,kijsh(i1,ns)), i1= &
          1, min(nijcalc(ns),lm2))
        Write (1337, 210)(atomimp(ish(ns,kijsh(i1,ns))), atomimp(jsh(ns, &
          kijsh(i1,ns))), ijtabsym((ish(ns,kijsh(i1,ns))-1)*natomimp+jsh(ns, &
          kijsh(i1,ns))), i1=1, min(nijcalc(ns),lm2))
      End Do
!ccc      WRITE (6,99009) MIN(NATYP,25)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

140   Format (6X, 'Dimension ERROR: please increase the ', A, ' parameter', /, &
        6X, A, ' to a value >=', I5, /)
150   Format (8X, 60('-'), /, 8X, 'Data used for J_ij calculation:', /, /, &
        10X, 'Number of sites/types i        (NAEZ/NATYP) :', 2(1X,I3), /, &
        10X, 'Number of atoms in the cluster   (NATOMIMP) :', 1X, I3, /, 10X, &
        'Number of ij pairs                 (NOFGIJ) :', 1X, I3, /, 10X, &
        'Number of representative pairs     (NSHELL) :', 1X, I3, /, /, 10X, &
        'ij-pairs calculation table ( 1 = calculated )', /, 10X, &
        'IA   JA = 1 ..', I3)
160   Format (10X, I3, 3X, 25(I3))
170   Format (/, 10X, 'ij-shells table ', /, 10X, 'IA   JA = 1 ..', I3)
180   Format (/, 10X, 'ij-symmetries table ', /, 10X, 'IA   JA = 1 ..', I3)
190   Format (/, 10X, 'effectively calculated pairs/shells', /, 10X, &
        'SHELL   (IAT,JAT) ', /, 10X, 'SHELL   (IQ,JQ - ISYM) ')
200   Format (10X, I4, 3X, 5(I3,',',I3,5X))
210   Format (10X, 4X, 3X, 5(I3,',',I3,' - ',I2))
220   Format (/, 10X, 'effectively calculated type-type pairs (shells)', /, &
        10X, 'IT   JT = 1 ..', I3)
    End Subroutine
