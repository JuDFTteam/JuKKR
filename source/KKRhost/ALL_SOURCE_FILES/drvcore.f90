    Subroutine drvcore(iprint, itprt, lcore, ncore, cscl, vtin, btin, rin, a, &
      b, drdiin, r2drdiin, zat, jws, ishift, rhoc, ecorerel, nkcore, kapcore, &
      ecore, lmaxd, irmd)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   * driving routine to call relativistic < CORE > routine            *
!   * counterpart of < COREL > of the non- or scalar-relativistic mode *
!   *                                                                  *
!   * ATTENTION: all the variables connected with hyperfine fields are *
!   *            OFF                                                   *
!   *                                                                  *
!   * The non-relativistic variables NCORE and LCORE are used as read  *
!   * in from the potential file; they are not modified and are again  *
!   * written out for the next iteration ( routine <RITES>)            *
!   *                                                                  *
!   * Relativistic variables passed outside this routine:              *
!   *    NKCORE(1..NCORE)       = number of KAPPA values for a given   *
!   *                             (n,l) core state                     *
!   *    KAPCORE(1..NCORE,1..NCORE) = the (maximum 2) values of KAPPA  *
!   *
!   *    ECOREREL(1..NCORE,1..NCORE) = for a given (n,l) state the core*
!   *                              energies corresponding first/second *
!   *                              KAPPA value, AVERAGED over \mu's    *
!   *                              These values are written out to the *
!   *                              potential file (routine <RITES>),   *
!   *                              but the read in (routine <STARTB1>) *
!   *                              updates the ECORE array             *
!   *     Please note that ALL the core energies (also \mu-resolved)   *
!   *     are output by <CORE> routine but not passed out of here      *
!   *                                                                  *
!   *    ECORE(1..NCORE,1..2) = this (non- or scalar relativistic)     *
!   *                           variable is updated here to be used in *
!   *                           calculating/printing out the spin-     *
!   *                           resolved energies (see <ESPCB> )       *
!   *                           A SUMMATION is done here:              *
!   *        ECORE(L,1/2) = SUM_{\kappa=-L-1,L} E(\kappa,\mu)          *
!   *                       /(2*L+1)                                   *
!   *                           with negative \mu's for E(*,1) and     *
!   *                           positive \mu's for E(*,2)              *
!   *                                                                  *
!   *                           v.popescu July/2002                    *
!   ********************************************************************
      Implicit None

! PARAMETER definitions

      Integer :: ntmax, nmmax, ncstmax, nmemax, nlmax, nkmmax
      Parameter (ntmax=1, nmmax=1, ncstmax=6, nmemax=5, nlmax=5, &
        nkmmax=2*nlmax**2) ! NLMAX should be >= LCOREMAX + 1
      Integer :: nrmax
      Parameter (nrmax=900)
      Real (Kind=dp) :: dzero
      Parameter (dzero=0.0E0_dp)

! Dummy arguments
      Real (Kind=dp) :: a, b
      Integer :: iprint, irmd, itprt, lmaxd, ncore, ishift

!obs: here, in contrast to DRVRHO, one works with local
!arrays, since they have to be set up as far as NRMAX (see CORE)

      Real (Kind=dp) :: vtin(irmd), btin(irmd)
      Real (Kind=dp) :: drdiin(irmd), r2drdiin(irmd)
      Real (Kind=dp) :: rin(irmd), cscl(lmaxd+1)

      Real (Kind=dp) :: ecore(20, 2), ecorerel(20*2)
      Integer :: kapcore(20*2), lcore(20, 2), nkcore(20)
      Integer :: zat(ntmax), jws(nmmax)
      Real (Kind=dp) :: rhoc(irmd, 2)

! Local variables
      Real (Kind=dp) :: bcor(ntmax), bcors(ntmax), ea, ecor(ncstmax), &
        ecortab(120, ntmax), fcor(nrmax, 2, ncstmax), gcor(nrmax, 2, ncstmax), &
        qdia(nkmmax), qmdia(nkmmax), qmoff(nkmmax), qoff(nkmmax), &
        rhochr(nrmax, ntmax), rhospn(nrmax, ntmax), sdia(nkmmax), &
        smdia(nkmmax), smoff(nkmmax), soff(nkmmax), szcor(ncstmax)
      Real (Kind=dp) :: vt(nrmax, ntmax), bt(nrmax, ntmax), ctl(ntmax, nlmax)
      Real (Kind=dp) :: drdi(nrmax, nmmax), r2drdi(nrmax, nmmax)
      Real (Kind=dp) :: r(nrmax, nmmax)
      Integer :: i, icall, ikmcor(ncstmax, 2), imt(ntmax), ip, ismqhfi, it, &
        itxray, izero(ncstmax), j, kapcor(ncstmax), lcxray(ntmax), &
        mm05cor(ncstmax), ncort(ntmax), ncxray(ntmax), nkpcor(ncstmax), nt, &
        nucleus
      Integer :: lcoremax

      Save :: imt, itxray, nt, nucleus, qdia, qmdia, qmoff, qoff, sdia, smdia, &
        smoff, soff

      Data ncxray/ntmax*0/, lcxray/ntmax*0/, ismqhfi/0/

      Data icall/0/

      icall = icall + 1

!=======================================================================
!       initialise relativistic and dummy variables and SAVE them
!=======================================================================
      If (icall==1) Then

        If (lmaxd>nlmax-1) Then
          Write (6, *) ' LMAXD = ', lmaxd, ' > NLMAX-1 = ', nlmax - 1
          Stop ' Increase NLMAX in < DRVCORE > '
        End If

        If (irmd>nrmax) Then
          Write (6, *) ' IRMD = ', irmd, ' > NRMAX = ', nrmax
          Write (6, *) ' Increase NRMAX in < sprkkr_rmesh.dim > '
          Stop ' In < DRVCORE > '
        End If

        itxray = 0

        Do it = 1, ntmax
          imt(it) = 1
        End Do

        nt = 1
        nucleus = 0

        Do it = 1, nkmmax
          sdia(it) = dzero
          smdia(it) = dzero
          soff(it) = dzero
          smoff(it) = dzero

          qdia(it) = dzero
          qmdia(it) = dzero
          qoff(it) = dzero
          qmoff(it) = dzero
        End Do

      End If ! ICALL.EQ.1
!=======================================================================

! --> fill up CTL array for the case of core states with higher L values
!     than those used in the valence band

      lcoremax = 0
      Do it = 1, ncore
        j = lcore(it, 1)
        lcoremax = max(lcoremax, j)
      End Do
      If (lcoremax>nlmax-1) Then
        Write (6, *) ' LCOREMAX = ', lcoremax, ' > NLMAX-1 = ', nlmax - 1
        Stop ' Increase NLMAX in < DRVCORE > '
      End If
      Do j = 1, lmaxd + 1
        ctl(1, j) = cscl(j)
      End Do
      If (lcoremax>0) Then
        Do j = lcoremax + 1, nlmax
          ctl(1, j) = cscl(lcoremax)
        End Do
      End If

      Call dcopy(jws(1), vtin, 1, vt(1,1), 1)
      Call dcopy(jws(1), btin, 1, bt(1,1), 1)
      Call dcopy(jws(1), rin, 1, r(1,1), 1)
      Call dcopy(jws(1), drdiin, 1, drdi(1,1), 1)
      Call dcopy(jws(1), r2drdiin, 1, r2drdi(1,1), 1)

      Do j = jws(1) + 1, nrmax
        ea = exp(a*real(j+ishift-1,kind=dp)) ! corrected from (J-1) 07.05.2004
        r(j, 1) = b*(ea-1E0_dp)
        drdi(j, 1) = a*b*ea
        r2drdi(j, 1) = r(j, 1)*r(j, 1)*drdi(j, 1)
        vt(j, 1) = 0E0_dp
        bt(j, 1) = 0E0_dp
      End Do

      ncort(1) = 0 ! no. of core electrons = no. of diff. core states

      Do it = 1, ncore
        ncort(1) = ncort(1) + 2*(2*lcore(it,1)+1)
      End Do

      Call core(iprint, itprt, nt, ncort, ctl, vt, bt, zat, nucleus, r, &
        r2drdi, drdi, jws, imt, rhochr, rhospn, ecortab, gcor, fcor, ecor, &
        szcor, kapcor, mm05cor, nkpcor, ikmcor, izero, ncxray, lcxray, itxray, &
        bcor, bcors, sdia, smdia, soff, smoff, qdia, qoff, qmdia, qmoff, &
        nkmmax, nmemax, ismqhfi, ntmax, nrmax, nmmax, ncstmax, nlmax)

      Call rinit(2*irmd, rhoc(1,1))

      Do i = 1, jws(1)
        ip = i + ishift
        rhoc(ip, 2) = (rhochr(i,1)+rhospn(i,1))*0.5E0_dp*(r(i,1)**2)
        rhoc(ip, 1) = (rhochr(i,1)-rhospn(i,1))*0.5E0_dp*(r(i,1)**2)
      End Do

      Call sumecore(ncore, lcore(1,1), ecortab(1,1), nkcore, ecorerel, ecore, &
        kapcore)

    End Subroutine

    Subroutine sumecore(ncore, lcore, ecortab, nkcore, ecorerel, ecore, &
      kapcore)
      Use mod_datatypes, Only: dp
      Implicit None

! Dummy arguments
      Integer :: ncore
      Real (Kind=dp) :: ecore(20, 2), ecorerel(20*2), ecortab(*)
      Integer :: kapcore(20*2), lcore(*), nkcore(*)

! Local variables
      Real (Kind=dp) :: dble
      Integer :: i, ic, icrel, jrel, kfg(4), l, lmp1, lmxc, lp1, muem05, nc, &
        nmax, nn, nsol, wgt(2)
      Real (Kind=dp) :: mj
      Intrinsic :: abs

! --> find the principal quantum numbers

      Do ic = 1, 4
        kfg(ic) = 0
      End Do
      Do ic = 1, ncore
        If (lcore(ic)==0) kfg(1) = kfg(1) + 1
        If (lcore(ic)==1) kfg(2) = kfg(2) + 1
        If (lcore(ic)==2) kfg(3) = kfg(3) + 1
        If (lcore(ic)==3) kfg(4) = kfg(4) + 1
      End Do

      If (kfg(2)/=0) kfg(2) = kfg(2) + 1
      If (kfg(3)/=0) kfg(3) = kfg(3) + 2
      If (kfg(4)/=0) kfg(4) = kfg(4) + 3

      lmxc = 0
      If (kfg(2)/=0) lmxc = 1
      If (kfg(3)/=0) lmxc = 2
      If (kfg(4)/=0) lmxc = 3

      lmp1 = lmxc + 1
      nc = 0
      Do lp1 = 1, lmp1
        l = lp1 - 1
        nmax = kfg(lp1)
        Do nn = lp1, nmax
          nc = nc + 1
          ecorerel(nc) = 0.0E0_dp
          ecorerel(nc+20) = 0.0E0_dp

!  ECOREREL(NC..NC+20) = 1st/2nd value of \kappa for current l


! --> icrel is pointing a core-state (n,l) = (nn,l) from the
!     relativistic-routine sequence 1s,2s,2p,3s,3p,... to the
!     ECOREREL array  1s,2s,3s,2p,3p,...

!     (nn,l) -> nn*(nn-1)*(2*nn-1)/3 + 2*l**2

          icrel = nn*(nn-1)*(2*nn-1)/3 + 2*l**2
          jrel = 0
          wgt(1) = 0
          wgt(2) = 0
          Do muem05 = -l - 1, +l
            mj = muem05 + 0.5E0_dp
            If (abs(mj)>l) Then
              nsol = 1
            Else
              nsol = 2
            End If
            Do i = 1, nsol
              wgt(i) = wgt(i) + 1
              jrel = jrel + 1
              ecorerel((i-1)*20+nc) = ecorerel((i-1)*20+nc) + &
                ecortab(icrel+jrel)
            End Do
          End Do
          nkcore(nc) = 1
          If (l/=0) nkcore(nc) = 2
          kapcore(nc) = -l - 1
          kapcore(nc+20) = l

          Do i = 1, nkcore(nc)
            ecorerel((i-1)*20+nc) = ecorerel((i-1)*20+nc)/dble(wgt(i))
          End Do

! --> update the array ECORE(1..NCORE,UP/DOWN) as

!         ECORE(L,SIGMA) = 1/(2*L+1) *
!              SUM_KAPPA(L) SUM_(MUE,SIGN(MUE)=SIGN(SIGMA)) ECORTAB(KAP,MUE)
!     i.e., states with negative MUE are added to SPIN DOWN, those with
!     positive MUE to SPIN UP.

!     ECORE is used later in calculating the total energy

          ecore(nc, 1) = 0.0E0_dp
          ecore(nc, 2) = 0.0E0_dp
          Do i = 1, 2*l + 1
            ecore(nc, 1) = ecore(nc, 1) + ecortab(icrel+i)
            ecore(nc, 2) = ecore(nc, 2) + ecortab(icrel+2*l+1+i)
          End Do
          Do i = 1, 2
            ecore(nc, i) = ecore(nc, i)/dble(2*l+1)
          End Do
        End Do
      End Do
    End Subroutine
