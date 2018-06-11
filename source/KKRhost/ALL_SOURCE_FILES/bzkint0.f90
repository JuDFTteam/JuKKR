    Subroutine bzkint0(nshell, naez, natyp, noq, rbasis, kaoez, icc, bravais, &
      recbv, atomimp, rsymat, isymindex, nsymat, ifilimp, natomimp, nsh1, &
      nsh2, rclsimp, ratom, ijtabsym, ijtabsh, ijtabcalc, iofgij, jofgij, &
      nofgij, ish, jsh, rrot, dsymll, para, qmtet, qmphi, symunitary, hostimp, &
      intervx, intervy, intervz, ielast, ez, kmesh, maxmesh, maxmshd, nsymaxd, &
      krel, lmaxd, lmmaxd, kpoibz, naezd, natypd, natomimpd, nsheld, nembd)
      Use mod_datatypes, Only: dp
      Implicit None
!.. Parameters ..
      Integer :: nsymaxd, krel, lmaxd, lmmaxd
      Integer :: kpoibz, naezd, natypd, natomimpd, nsheld, nembd
!..
!.. Scalar Arguments ..
      Integer :: icc, naez, natomimp, natyp, nsymat, nofgij
      Integer :: intervx, intervy, intervz, maxmesh, maxmshd, ielast
      Character (Len=40) :: ifilimp
!..
!.. Array Arguments ..
      Complex (Kind=dp) :: dsymll(lmmaxd, lmmaxd, nsymaxd), ez(*)
      Real (Kind=dp) :: bravais(3, 3), ratom(3, nsheld), &
        rbasis(3, naezd+nembd), rclsimp(3, natomimpd), recbv(3, 3), &
        rrot(48, 3, nsheld), rsymat(64, 3, 3)
      Integer :: atomimp(natomimpd), isymindex(nsymaxd), &
        kaoez(natypd, naezd+nembd), noq(naezd), kmesh(*), nsh1(*), nsh2(*), &
        nshell(0:nsheld), ijtabsym(*), ijtabsh(*), ijtabcalc(*), iofgij(*), &
        jofgij(*), ish(nsheld, *), jsh(nsheld, *)

!..anges for impurity 20/02/2004 -- v.popescu according to 
!..                                 n.papanikolaou 

      Integer :: hostimp(0:natypd)
!..
!.. Local Scalars ..
      Integer :: i, ishell, iu, iprint
      Logical :: lirr
!..
!.. Local Arrays ..
      Character (Len=10) :: rotname(64)
!.. magnetisation angles ..
      Real (Kind=dp) :: qmtet(naezd), qmphi(naezd)
!.. unitary/antiunitary symmetry flag
      Logical :: symunitary(nsymaxd), para
!..
!.. External Functions ..
      Logical :: test, opt
      External :: test, opt
!..
!.. External Subroutines ..
      External :: bzkmesh, crtstar, findgroup, gfshells, pointgrp, symtaumat
!..

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      Write (1337, '(79(1H=),/,15X,A)') &
        'BZKINT0: finding symmetry, setting BZ integration'
      Write (1337, '(79(1H=),/)')
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

      Call pointgrp(rsymat, rotname)
      Call findgroup(bravais, recbv, rbasis, naez, rsymat, rotname, isymindex, &
        nsymat, para, qmtet, qmphi, symunitary, krel, naezd, nembd, nsymaxd)

      lirr = .True.
      iprint = 0
      If (test('TAUSTRUC')) iprint = 2

! --> test: full BZ integration

      If (test('fullBZ  ') .Or. opt('NEWSOSOL')) Then
        nsymat = 1
        lirr = .False.
        Write (1337, '(8X,2A,/)') &
          'Test option < fullBZ > or Run option < NEWSOSOL >: ', &
          ' overriding NSYMAT, generate full BZ k-mesh'
      End If

! --> generate BZ k-mesh

      Call bzkmesh(intervx, intervy, intervz, maxmesh, lirr, bravais, recbv, &
        nsymat, rsymat, isymindex, symunitary, ielast, ez, kmesh, iprint, &
        krel, kpoibz, maxmshd)

      Call symtaumat(rotname, rsymat, dsymll, nsymat, isymindex, symunitary, &
        naezd, lmmaxd, naez, lmaxd+1, krel, iprint, nsymaxd)

! Now DSYMLL hold NSYMAT symmetrization matrices

! 20.02.2004
      Call gfshells(icc, natomimp, nsh1, nsh2, ijtabsym, ijtabsh, ijtabcalc, &
        iofgij, jofgij, nofgij, ish, jsh, nshell, naez, natyp, noq, rbasis, &
        bravais, ifilimp, ratom, rclsimp, nsymat, isymindex, rsymat, kaoez, &
        atomimp, rotname, hostimp, lmaxd, lmmaxd, naezd, natypd, natomimpd, &
        nembd, nsheld)

! -->  creates difference vectors RROT for BZ integration in KKRMAT01

      Call crtstar(ratom, nshell(0), rsymat, nsymat, isymindex, rrot)
! ----------------------------------------------------------------------
      If (iprint>2) Then
        Do ishell = 1, nshell(0)
          Write (1337, Fmt='(I4)') ishell
          Write (1337, Fmt='((I4,3F10.1))')(iu, (rrot(iu,i, &
            ishell),i=1,3), iu=1, nsymat)
        End Do
      End If
! ----------------------------------------------------------------------
    End Subroutine
