    Subroutine drvreltmat(eryd, tmatll, vt, bt, r, drdi, r2drdi, zat, jws, &
      solver, soctl, ctl, lmmaxd, lmaxd, irmd)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   * driving routine to call relativistic < SSITE > routine           *
!   * only to calculate the single-site t matrix                       *
!   * v.popescu, munich, may 2004                                      *
!   *                                                                  *
!   ********************************************************************

      Implicit None

! PARAMETER definitions
      Integer :: nrmax
      Parameter (nrmax=900)
      Integer :: nlamax, nqmax, ntmax, nmmax
      Parameter (nlamax=1, nqmax=1, ntmax=1, nmmax=1)
      Integer :: nlmax, nkmmax, nmuemax, nkmpmax, nkmax, linmax
      Parameter (nlmax=5) ! this should be >= LMAXD + 1
      Parameter (nkmmax=2*nlmax**2, nkmax=2*nlmax-1)
      Parameter (nkmpmax=nkmmax+2*nlmax, nmuemax=2*nlmax)
      Parameter (linmax=2*nlmax*(2*nlmax-1))

! Dummy arguments
      Integer :: lmaxd, lmmaxd, irmd
      Integer :: zat(ntmax), jws(nmmax)
      Complex (Kind=dp) :: tmatll(lmmaxd, lmmaxd)
      Real (Kind=dp) :: soctl(nlmax)
      Real (Kind=dp) :: ctl(nlmax)
      Real (Kind=dp) :: vt(nrmax), bt(nrmax)
      Real (Kind=dp) :: r(nrmax, nmmax), r2drdi(nrmax, nmmax)
      Real (Kind=dp) :: drdi(nrmax, nmmax)

! Local variables
      Real (Kind=dp) :: ameopo(nkmmax, nkmmax, nlamax, 3), &
        at(nrmax, nlamax, 3, ntmax)
      Complex (Kind=dp) :: bzj(linmax, ntmax), bzz(linmax, ntmax), &
        dzj(linmax, ntmax), dzz(linmax, ntmax), eryd, &
        msst(nkmmax, nkmmax, ntmax), ozj(linmax, ntmax), ozz(linmax, ntmax), &
        p, qzj(linmax, ntmax), qzz(linmax, ntmax), szj(linmax, ntmax), &
        szz(linmax, ntmax)
      Complex (Kind=dp) :: tsst(nkmmax, nkmmax, ntmax), &
        tsstlin(linmax, ntmax), tzj(linmax, ntmax), tzz(linmax, ntmax)
      Complex (Kind=dp) :: ozzs(linmax, ntmax, 2), ozjs(linmax, ntmax, 2)
      Logical :: calcint, getirrsol
      Real (Kind=dp) :: cgc(nkmpmax, 2)
      Integer :: i, ihyper, ikm1lin(linmax), ikm2lin(linmax), il, imt(ntmax), &
        imue, iprint, iq, iqat(nqmax, ntmax), it, iwrirrwf, iwrregwf, j, &
        lopt(ntmax), mmax, nkm, nkmq(nqmax), nl, nlinq(nqmax), nlq(nqmax), nt, &
        nucleus
      Integer :: nfilcbwf
      Integer :: nsollm(nlmax, nmuemax), ltab(nmuemax), kaptab(nmuemax), &
        nmuetab(nmuemax)
      Character (Len=10) :: solver
      Integer :: icall

      Data icall/0/

      Save :: icall, ikm1lin, ikm2lin, lopt, nlq, nkmq, iqat, imt, nkm, &
        ihyper, iprint, it, nt, nucleus, iwrregwf, iwrirrwf, calcint, &
        getirrsol, nfilcbwf

      icall = icall + 1

!=======================================================================
!       initialise relativistic and dummy variables and SAVE them
!=======================================================================
      If (icall==1) Then

        If (lmaxd>nlmax-1) Then
          Write (6, *) ' LMAXD = ', lmaxd, ' > NLMAX-1 = ', nlmax - 1
          Stop ' Increase NLMAX in < DRVRELTMAT > '
        End If

        If (irmd>nrmax) Then
          Write (6, *) ' IRMD = ', irmd, ' > NRMAX = ', nrmax
          Write (6, *) ' Increase NRMAX in < sprkkr_rmesh.dim > '
          Stop ' and in < DRVRELTMAT > '
        End If

        nl = lmaxd + 1 ! no need to save, used only here once
        iprint = 0

        Do i = 1, nmuemax
          ltab(i) = i/2
          If (2*ltab(i)==i) Then
            kaptab(i) = ltab(i)
          Else
            kaptab(i) = -ltab(i) - 1
          End If
          nmuetab(i) = 2*abs(kaptab(i))
        End Do

        Do il = 1, nlmax
          mmax = 2*il
          Do imue = 1, mmax
            If ((imue==1) .Or. (imue==mmax)) Then
              nsollm(il, imue) = 1
            Else
              nsollm(il, imue) = 2
            End If
          End Do
        End Do

        Call ikmlin(iprint, nsollm, ikm1lin, ikm2lin, nlmax, nmuemax, linmax, &
          nlmax)

        Call calccgc(ltab, kaptab, nmuetab, cgc, nkmax, nmuemax, nkmpmax)

        Do it = 1, ntmax
          imt(it) = 1
          lopt(it) = -1 ! this should change for Brooks' OP
        End Do

        Do iq = 1, nqmax
          nlq(iq) = nl
          nkmq(iq) = lmmaxd
          nlinq(iq) = 2*nlq(iq)*(2*nlq(iq)-1)
          iqat(iq, 1) = 1
        End Do

        nkm = lmmaxd
        ihyper = 0
        nt = 1
        it = 1
        nucleus = 0

        iwrregwf = 0
        iwrirrwf = 0
        calcint = .False.
        getirrsol = .False.
        nfilcbwf = 0
      End If ! ICALL.EQ.1
!=======================================================================

      Call ssite(iwrregwf, iwrirrwf, nfilcbwf, calcint, getirrsol, soctl, ctl, &
        eryd, p, ihyper, iprint, ikm1lin, ikm2lin, nlq, nkmq, nlinq, nt, nkm, &
        iqat, tsst, msst, tsstlin, dzz, dzj, szz, szj, ozz, ozj, bzz, bzj, &
        qzz, qzj, tzz, tzj, vt, bt, at, zat, nucleus, r, drdi, r2drdi, jws, &
        imt, ameopo, lopt, solver, cgc, ozzs, ozjs, nlmax, nqmax, linmax, &
        nrmax, nmmax, ntmax, nkmmax, nkmpmax, nlamax)

      Do j = 1, nkm
        Call zcopy(nkm, tsst(1,j,it), 1, tmatll(1,j), 1)
      End Do

      Return
    End Subroutine
