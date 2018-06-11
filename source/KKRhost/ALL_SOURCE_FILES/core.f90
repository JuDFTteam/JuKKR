subroutine core(iprint, itprt, nt, ncort, ctl, vt, bt, z, nucleus, r, r2drdi, &
  drdi, jws, imt, rhochr, rhospn, ecortab, gcor, fcor, ecor, szcor, kapcor, &
  mm05cor, nkpcor, ikmcor, izero, ncxray, lcxray, itxray, bcor, bcors, sdia, &
  smdia, soff, smoff, qdia, qoff, qmdia, qmoff, nkmmax, nmemax, ismqhfi, &
  ntmax, nrmax, nmmax, ncstmax, nlmax)
!   ********************************************************************
!   *                                                                  *
!   *   SUBROUTINE TO CALCULATE THE RELATIVISTIC CORE WAVE             *
!   *   FUNCTIONS FOR A SPIN-DEPENDENT POTENTIAL                       *
!   *                                                                  *
!   *   FOR A GIVEN POTENTIAL THE NUMBER OF CORE AND VALENCE           *
!   *   ELECTRONS IS DETERMINED AND ALL CORE STATES THEN CALCULATED    *
!   *   > THE ROUTINE IS ORGANIZED AS DESCLAUX'S ROUTINE <RESLD>       *
!   *     BUT FINDS THE CORRECTION TO THE E-EIGENVALUE AND THE         *
!   *     MATCHING PARAMETERS BY A NEWTON RAPHSON ALGORITHM            *
!   *     THIS IS IN VARIANCE TO THE METHOD SUGGESTED BY CORTONA       *
!   *   > SET THE SWITCH 'CHECK'  TO COPARE E-EIGENVALUES WITH         *
!   *     RESULTS OBTAINED WITH THE CONVENTIONAL E-CORRECTION          *
!   *     ALGORITHM, WHICH WORKS ONLY IF NO COUPLING IS PRESENT !      *
!   *   > THE FUNCTIONS  {GC,FC}(I,J) J=1,NSOL ARE THE LINEAR          *
!   *     INDEPENDENT SOLUTIONS TO THE DIFFERENTIAL EQUATIONS WITH     *
!   *     KAPPA-CHARACTER I=1,NSOL;   FOR OUTWARD AND INWARD           *
!   *     INTEGRATION THE SAME ARRAYS ARE USED !                       *
!   *   > THE PROPER SOLUTIONS SATISFYING THE BOUNDARY CONDITIONS      *
!   *     AT R=0 AND(!) R=INFINITY ARE STORED IN {GCK,FCK}(K,S)        *
!   *     S,K=1,NSOL   SOLUTION S=1 FOR  KAPPA = - L - 1               *
!   *                           S=2 FOR  KAPPA = + L (IF EXISTENT)     *
!   *   > THE SWITCH NUCLEUS SELECTS WHETHER A FINITE NUCLEUS          *
!   *     SHOULD BE USED                                               *
!   *                                                                  *
!   *   ADAPTED FOR FINITE NUCLEUS       MB MAR. 1995                  *
!   *   HYPERFINE FIELD SPLITTING introduced if icore=1 MB JUN. 1995   *
!   *                                                                  *
!   *   SCALEB:                                                        *
!   *   if the B-field is quite high it might happen that the routine  *
!   *   fails to find both 'spin-orbit-split' solutions.               *
!   *   in that case the whole l-shell is rerun with the B-field       *
!   *   gradually switched on, i.e. scaled with a parameter that       *
!   *   increases from 0 to 1 during the iteration loop  HE Nov. 95    *
!   *                                                                  *
!   *   ITXRAY =  0  run over all core states to get charge density    *
!   *   ITXRAY >  0  deal exclusively with state  NCXRAY,LCXRAY        *
!   *   ITXRAY <  0  state  NCXRAY,LCXRAY  is checked to be a          *
!   *                bound state or not. on return:                    *
!   *                ITXRAY = |ITXRAY| indicates bound state           *
!   *                ITXRAY = -1       indicates NO bound state found  *
!   *                                                                  *
!   *                                                                  *
!   *   few changes in the TB-KKR implementation as compared to SPR    *
!   *        IPRINT values between 0 and 2                             *
!   *        ITPRT  correct value of the atom-type index               *
!   ********************************************************************
  use :: mod_types, only: t_inc
      Use mod_datatypes, Only: dp
  implicit none


! PARAMETER definitions
  real (kind=dp) :: unend, tolvar, trymix, dvstep
  parameter (unend=600.0d0, tolvar=1.0d-6, trymix=0.01d0, dvstep=0.01d0)
  integer :: itermax, nlshellmax
  parameter (itermax=200, nlshellmax=15)

! Dummy arguments
  integer :: iprint, ismqhfi, itprt, itxray, ncstmax, nkmmax, nlmax, nmemax, &
    nmmax, nrmax, nt, ntmax, nucleus
  real (kind=dp) :: bcor(ntmax), bcors(ntmax), bt(nrmax, ntmax), &
    ctl(ntmax, nlmax), drdi(nrmax, nmmax), ecor(ncstmax), ecortab(120, ntmax), &
    fcor(nrmax, 2, ncstmax), gcor(nrmax, 2, ncstmax), qdia(nkmmax), &
    qmdia(nkmmax), qmoff(nkmmax), qoff(nkmmax), r(nrmax, nmmax), &
    r2drdi(nrmax, nmmax), rhochr(nrmax, ntmax), rhospn(nrmax, ntmax), &
    sdia(nkmmax), smdia(nkmmax), smoff(nkmmax), soff(nkmmax), szcor(ncstmax), &
    vt(nrmax, ntmax)
  integer :: ikmcor(ncstmax, 2), imt(ntmax), izero(ncstmax), jws(nmmax), &
    kapcor(ncstmax), lcxray(ntmax), mm05cor(ncstmax), ncxray(ntmax), &
    nkpcor(ncstmax), z(ntmax), ncort(ntmax)

! Local variables
  real (kind=dp) :: aux, bb(nrmax*2), bhf(2, 2), bhf1(2, 2), bhf2(2, 2), &
    bsh, bsol, bsum, cgd(2), cgmd(2), cgo, dec, dedv(4, 4), dovrc(nrmax*2), &
    d_p(2, 2, nrmax*2), dq(2, 2, nrmax*2), drdic(nrmax*2), drovrn(2*nrmax), &
    dv(4), dvde(4, 4), ec, ecc, elim, err(4), errnew(4), fc(2, 2, nrmax*2), &
    fck(2, 2, nrmax*2), gc(2, 2, nrmax*2), gck(2, 2, nrmax*2), mj, niw(2), &
    norm, now(2), piw(2, 2), pow(2, 2), qiw(2, 2), qow(2, 2), &
    r2drdic(nrmax*2), rat, ratt, rc(nrmax*2), rint(nrmax), rnuc, rr, scale, &
    shf(2, 2, nmemax), simp, split(nmemax), split1(nmemax), &
    split2(nmemax, ntmax), split3(nmemax, ntmax), sz, val, var(4), varnew(4), &
    vartab(4, 20), vv(nrmax*2), vz, w, wp(2, 2, nrmax*2), wq(2, 2, nrmax*2)
  logical :: check, ferro, scaleb, suppressb, bndstachk
  real (kind=dp) :: dble, dsign
  integer :: i, ic, ic1, ic2, icst, ie, iflag, ii, il, ilc, ilshell, im, imin, &
    in, info, ipiv(4), ish, istart, it, iter, iv, j, jlim, jtop, jv, k, &
    kap(2), kap1, kap2, kc, l, lcp1, lll, loop, lqntab(nlshellmax), muem05, n, &
    nlshell, nmatch, nn, node, nqn, nqntab(nlshellmax), nrc, nsh, nsol, nvar, &
    nzero, s, t
  integer :: iabs
  integer :: ikapmue
  real (kind=dp) :: rnuctab
  character (len=10) :: txtb(1:5)
  character (len=3) :: txtk(4)
  character (len=1) :: txtl(0:3)

  data txtb/'B_nses', 'B_nseo', 'B_ssc ', 'B_sum ', 'B_tot '/
  data txtl/'s', 'p', 'd', 'f'/
  data txtk/'1/2', '3/2', '5/2', '7/2'/
  data nqntab/1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 4, 5, 6, 6/
  data lqntab/0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 3, 2, 0, 1/
  data check/.false./

  nrc = 2*nrmax

  call rinit(120*ntmax, ecortab)
  call rinit(4*20, vartab)

  if (itxray<0) then
    itxray = abs(itxray)
    bndstachk = .true.
  else
    bndstachk = .false.
  end if

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
  do it = 1, nt
    suppressb = .false.
    scaleb = .false.
    scale = 1d0

    im = imt(it)
    jtop = jws(im)

    rat = r(nrmax, im)/r(nrmax-1, im)
    do n = 1, nrmax
      rc(n) = r(n, im)
      drdic(n) = drdi(n, im)
      r2drdic(n) = r2drdi(n, im)
      dovrc(n) = drdic(n)/rc(n)
    end do
    do n = (nrmax+1), nrc
      rc(n) = rat*rc(n-1)
      drdic(n) = (rat-1.0d0)*rc(n-1)
      r2drdic(n) = rc(n)*rc(n)*drdic(n)
      dovrc(n) = drdic(n)/rc(n)
    end do
    if (nucleus/=0) then
      rnuc = rnuctab(z(it))
      in = 1
      do while (rc(in)<=rnuc)
        in = in + 1
      end do
!     INTEGRATION BOUNDARY FOR HYPERFINE FIELDS FOR FINITE NUCLEUS
!     2 MESH POINTS MORE FOR EXECUTING APPROPRIATE INTERPOLATION
!     TO REAL NUCLEAR RADIUS RNUC
      jlim = in + 2
      if (mod(jlim,2)==0) jlim = jlim - 1
    end if
    do i = 1, nrc
      if (nucleus/=0) drovrn(i) = (rc(i)/rnuc)**3*drdic(i)
    end do

    loop = 1
100 continue
    bcor(it) = 0.0d0
    bcors(it) = 0.0d0
    do i = 1, nmemax
      split2(i, it) = 0.0d0
      split3(i, it) = 0.0d0
    end do
    bsum = 0.0d0
    do n = 1, nrmax
      rhochr(n, it) = 0.0d0
      rhospn(n, it) = 0.0d0
    end do
    do n = 1, jws(im)
      vv(n) = vt(n, it)
      bb(n) = bt(n, it)
      bsum = bsum + abs(bb(n))
    end do

    if (suppressb) then
      do n = 1, jws(im)
        bb(n) = 0.0d0
      end do
      bsum = 0.0d0
    end if

    do n = (jws(im)+1), nrc
      vv(n) = 0.0d0
      bb(n) = 0.0d0
    end do

    nlshell = 0
    if (z(it)>2) nlshell = 1
    if (z(it)>10) nlshell = 3
    if (z(it)>18) nlshell = 5
    if (z(it)>30) nlshell = 6
    if (z(it)>36) nlshell = 8
    if (z(it)>48) nlshell = 9
    if (z(it)>54) nlshell = 11
    if (z(it)>70) nlshell = 12
    if (z(it)>80) nlshell = 13
    if (z(it)>86) nlshell = 15

    if (ncort(it)/=0) then
      nlshell = 0
      n = 0
      do ilshell = 1, nlshellmax
        l = lqntab(ilshell)
        n = n + 2*(2*l+1)
        if (n==ncort(it)) nlshell = ilshell
      end do
      if (nlshell==0) then
        write (*, *) 'NLSHELL not found for IT=', it, ' NCORT=', ncort(it)
        stop ' in <CORE>'
      end if
    end if

    if (bsum>1.0d-8) then
      ferro = .true.
    else
      ferro = .false.
      if (iprint>=1 .and. (t_inc%i_write>0)) write (1337, 170)
    end if

    if (itxray==0 .and. iprint>0 .and. (t_inc%i_write>0)) write (1337, 180) &
      itprt, z(it)


! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
!                   ---------------------------------------
!                   INITIALIZE QUANTUM NUMBERS  NQN  AND  L
!                   ---------------------------------------
    ic = 0
    do ilshell = 1, nlshell
      nqn = nqntab(ilshell)
      l = lqntab(ilshell)
      il = l + 1
      ilc = min(nlmax, il)
      nsh = 2*(2*l+1)

! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!     SKIP SHELL IF NOT NEEDED IN A  XRAY - CALCULATION
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      if (itxray/=0) then
        if (it/=itxray) go to 160
        if ((nqn/=ncxray(it)) .or. (l/=lcxray(it))) go to 160
        do icst = 1, ncstmax
          do kc = 1, 2
            do n = 1, nrmax
              gcor(n, kc, icst) = 0.0d0
              fcor(n, kc, icst) = 0.0d0
            end do
          end do
        end do
      end if
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      ish = 0
      bsh = 0.0d0
      do i = 1, nmemax
        split1(i) = 0.0d0
      end do
! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
      do muem05 = -l - 1, +l
        mj = muem05 + 0.5d0


        kap1 = -l - 1
        kap2 = l
        kap(1) = kap1
        kap(2) = kap2

        lll = l*(l+1)
        if (abs(mj)>l) then
          nsol = 1
        else
          nsol = 2
        end if

        if (ferro) then
          nvar = 2*nsol
        else
          nvar = 2
        end if

! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
        do s = 1, nsol
          ic = ic + 1
          ish = ish + 1
          t = 3 - s
!                  ----------------------------------------
!                   USE EC OF PREVIOUS RUNS AS START-VALUE
!                   TAKE SPIN-ORBIT SPLITTING INTO ACCOUNT
!                  ----------------------------------------
          if (ish>1) then
            ec = ecortab(ic-1, it)
            if (s==2) ec = ecortab(ic-1, it)*1.1d0
            if (ish>=4) ec = ecortab(ic-2, it)
            go to 130
          end if


!                                      --------------------
!                                         FIND  E-LIMIT
!                                      --------------------
          if (lll==0) then
            elim = -2*dble(z(it)**2)/(1.5d0*nqn*nqn)
          else
            elim = vv(1) + lll/rc(1)**2
            do n = 2, nrc
              val = vv(n) + lll/rc(n)**2
              if (val<=elim) elim = val
            end do
          end if

          ec = -dble(z(it)**2)/(2.0d0*nqn*nqn)

          istart = 1
110       continue
          if (ec<=elim) ec = elim*0.7d0

!                                      --------------------
!                                         FIND    NZERO
!                                      --------------------
          do n = 1, (nrc-1)
            if ((vv(n)-ec)*rc(n)**2>unend) then
              if (mod(n,2)==0) then
                nzero = n + 1
              else
                nzero = n
              end if
              go to 120
            end if
          end do
          nzero = nrc - 1
          write (6, 190) itprt, nqn, l, (nrc-1)
          stop
!                                      --------------------
!                                         FIND    NMATCH
!                                      --------------------
120       continue
          n = nzero + 1
          do nn = 1, nzero
            n = n - 1
            if ((vv(n)+lll/rc(n)**2-ec)<0.0) then
              nmatch = n
              go to 130
            end if
          end do
          if (t_inc%i_write>0) write (1337, 200) itprt, nqn, l, ec

130       continue
          call coredir(it, ctl(it,ilc), ec, l, mj, 'OUT', vv, bb, rc, drdic, &
            dovrc, nmatch, nzero, gc, fc, d_p, dq, wp, wq, pow, qow, piw, qiw, &
            cgd, cgmd, cgo, nrc, z(it), nucleus)

          node = 0
          do n = 2, nmatch
            if (gc(s,s,n)*gc(s,s,n-1)<0.0) node = node + 1
          end do


          if (iprint>=2 .and. (t_inc%i_write>0)) write (1337, 320) itprt, nqn, &
            l, kap(s), (2*muem05+1), ic, ish, 0, ec, nmatch, rc(nmatch), &
            nzero, rc(nzero), node, (gc(s,s,nmatch)/gc(s,s,nmatch-1))

          if (node/=(nqn-l-1)) then
            if (node>(nqn-l-1)) then
              ec = 1.2d0*ec
            else
              ec = 0.8d0*ec
            end if
            go to 110
          else if ((gc(s,s,nmatch)/gc(s,s,nmatch-1)<=0.0) .or. (gc(s,s, &
              nmatch)/gc(s,s,nmatch-1)>=1.0)) then
            ec = 0.9d0*ec
            go to 110
          end if


          call coredir(it, ctl(it,ilc), ec, l, mj, 'INW', vv, bb, rc, drdic, &
            dovrc, nmatch, nzero, gc, fc, d_p, dq, wp, wq, pow, qow, piw, qiw, &
            cgd, cgmd, cgo, nrc, z(it), nucleus)

!                                      --------------------
!                                       START VALUES FOR
!                                           PARAMETERS
!                                      --------------------

          var(1) = ec
          var(2) = pow(s, s)/piw(s, s)

          if (nsol/=2 .or. nvar==2) then
            do iv = 3, 4
              err(iv) = 0.0d0
              errnew(iv) = 0.0d0
              var(iv) = 0.0d0
              varnew(iv) = 0.0d0
              dv(iv) = 0.0d0
            end do
          else if (ish>=4) then
            do iv = 1, 4
              var(iv) = vartab(iv, ish-2)
            end do
          else

            do j = 1, nsol
              now(j) = 0.0d0
            end do
            do n = 1, nmatch - 1
              rr = rc(n)**3
              do j = 1, nsol
                now(j) = now(j) + gc(j, j, n)**2*rr
              end do
            end do

            do j = 1, nsol
              niw(j) = 0.0d0
            end do
            do n = nmatch, nzero - 1
              rr = rc(n)**3
              do j = 1, nsol
                niw(j) = niw(j) + gc(j, j, n)**2*rr
              end do
            end do

            ratt = pow(t, t)/piw(t, t)
            var(3) = trymix*(now(s)+niw(s)*var(2))/(now(t)+niw(t)*ratt)
            var(4) = ratt*var(3)/var(2)
          end if


          call coreerr(err, var, s, nsol, pow, qow, piw, qiw)

          do iv = 1, nvar
            dv(iv) = var(iv)
          end do

          iter = 0
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
140       continue
          iter = iter + 1

          if (scaleb) then
            scale = min(1.0d0, 0.1d0*iter)
            if (nsol==1) scale = 1.0d0
            do n = 1, jws(im)
              bb(n) = bt(n, it)*scale
            end do
          end if
!                         ----------------------------------
!                         CHECK WHETHER NUMBER OF NODES O.K.
!                         ----------------------------------
          if (iter>1) then
            node = 0
            do n = 2, (nmatch-1)
              if (gc(s,s,n)*gc(s,s,n-1)<0.0) node = node + 1
            end do
            if (iprint>=2 .and. (t_inc%i_write>0)) write (1337, 320) itprt, &
              nqn, l, kap(s), (2*muem05+1), ic, ish, iter, ec, nmatch, &
              rc(nmatch), nzero, rc(nzero), node, (gc(s,s,nmatch)/gc(s,s, &
              nmatch-1))

            if (node/=(nqn-l-1)) then
              if (node>(nqn-l-1)) then
                ec = 1.2d0*ec
              else
                ec = 0.8d0*ec
              end if
              istart = istart + 1
              if (istart<20) go to 110
            end if
          end if

          do iv = 2, nvar
            do jv = 1, nvar
              varnew(jv) = var(jv)
            end do
            varnew(iv) = var(iv) + dv(iv)*dvstep

            if (abs(var(iv))>1d-16) then
              if (abs(dv(iv)/var(iv))<tolvar) varnew(iv) = var(iv)* &
                (1.0d0+dsign(dvstep*tolvar,dv(iv)))
            else if (ferro) then
              if (t_inc%i_write>0) write (1337, 270) ' VAR(', iv, &
                ') = 0 for (T,N,L,K,M;S,NSOL) ', itprt, nqn, l, kap(s), &
                (2*muem05+1), '/2  ', s, nsol, '  --- suppress B'
              loop = 2
              suppressb = .true.
              go to 100
            else if (suppressb) then
              write (6, *) 'suppressing B did not help !!'
              stop 'in <CORE>'
            end if

            call coreerr(errnew, varnew, s, nsol, pow, qow, piw, qiw)

            do ie = 1, nvar
              if (abs(errnew(ie)-err(ie))<1d-16) then
                dedv(ie, iv) = 0.0d0
                if ((ie==iv) .and. .not. ferro) dedv(ie, iv) = 1.0d0
              else
                dedv(ie, iv) = (errnew(ie)-err(ie))/(varnew(iv)-var(iv))
              end if
            end do
          end do

          do jv = 1, nvar
            varnew(jv) = var(jv)
          end do
          varnew(1) = var(1) + dv(1)*dvstep
          if (abs(dv(1)/var(1))<tolvar) varnew(1) = var(1)* &
            (1.0d0+dsign(dvstep*tolvar,dv(1)))
          call coredir(it, ctl(it,ilc), varnew(1), l, mj, 'OUT', vv, bb, rc, &
            drdic, dovrc, nmatch, nzero, gc, fc, d_p, dq, wp, wq, pow, qow, &
            piw, qiw, cgd, cgmd, cgo, nrc, z(it), nucleus)
          call coredir(it, ctl(it,ilc), varnew(1), l, mj, 'INW', vv, bb, rc, &
            drdic, dovrc, nmatch, nzero, gc, fc, d_p, dq, wp, wq, pow, qow, &
            piw, qiw, cgd, cgmd, cgo, nrc, z(it), nucleus)

          call coreerr(errnew, varnew, s, nsol, pow, qow, piw, qiw)

          do ie = 1, nvar
            dedv(ie, 1) = (errnew(ie)-err(ie))/(varnew(1)-var(1))
          end do

          do ie = 1, nvar
            call dcopy(nvar, dedv(1,ie), 1, dvde(1,ie), 1)
          end do
          call dgetrf(nvar, nvar, dvde, 4, ipiv, info)
          call dgetri(nvar, dvde, 4, ipiv, dedv, 4*4, info)

          do iv = 1, nvar
            dv(iv) = 0.0d0
            do ie = 1, nvar
              dv(iv) = dv(iv) + dvde(iv, ie)*err(ie)
            end do
            var(iv) = var(iv) - dv(iv)
          end do

          if (var(1)>0.0d0) then
            if (iprint>=1 .and. (t_inc%i_write>0)) write (1337, *) &
              ' warning from <CORE> E=', var(1), it, nqn, l
            var(1) = -0.2d0
          end if

          call coredir(it, ctl(it,ilc), var(1), l, mj, 'OUT', vv, bb, rc, &
            drdic, dovrc, nmatch, nzero, gc, fc, d_p, dq, wp, wq, pow, qow, &
            piw, qiw, cgd, cgmd, cgo, nrc, z(it), nucleus)
          call coredir(it, ctl(it,ilc), var(1), l, mj, 'INW', vv, bb, rc, &
            drdic, dovrc, nmatch, nzero, gc, fc, d_p, dq, wp, wq, pow, qow, &
            piw, qiw, cgd, cgmd, cgo, nrc, z(it), nucleus)

          call coreerr(err, var, s, nsol, pow, qow, piw, qiw)

          ec = var(1)

          if (iprint>=2 .and. (t_inc%i_write>0)) write (1337, 210) loop, &
            scale, var(1), (var(iv), iv=1, 4), (dv(iv), iv=1, 4), &
            (err(ie), ie=1, 4)

!----------------------------------  check relative change in parameters
! ----------------------- parameters 3 and 4 = 0 for paramagnetic case !
          if (iter<itermax) then
            do iv = 1, nvar
              vartab(iv, ish) = var(iv)
!           IF( ABS(VAR(IV)) .EQ. 0.0D0 ) THEN
              if ((abs(var(iv))+abs(var(iv)))<1.0d-30) then
                if (ferro .and. (t_inc%i_write>0)) write (1337, '(A,I3,A)') &
                  ' VAR ', iv, ' = 0 ??????!!!!!'
              else if (abs(dv(iv)/var(iv))>tolvar) then
                go to 140
              end if
            end do
          else
            if (bndstachk) then
              itxray = -1
              return
            end if
            if (t_inc%i_write>0) write (1337, 220) itermax, &
              (var(iv), iv=1, 4), (dv(iv), iv=1, 4), (err(ie), ie=1, 4)
          end if
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

!                         ---------------------------------
!                         NORMALIZE WAVEFUNCTIONS ACCORDING
!                               TO MATCHING CONDITIONS
!                         ---------------------------------

!                                    INWARD - SOLUTION
          do n = nmatch, nzero
            do j = 1, nsol
              do i = 1, nsol
                gc(i, j, n) = gc(i, j, n)*var(2)
                fc(i, j, n) = fc(i, j, n)*var(2)
              end do
            end do
          end do

          if (nsol==2) then
!                                   OUTWARD - SOLUTION
            do n = 1, (nmatch-1)
              do i = 1, nsol
                gc(i, t, n) = gc(i, t, n)*var(3)
                fc(i, t, n) = fc(i, t, n)*var(3)
              end do
            end do
!                                    INWARD - SOLUTION
            do n = nmatch, nzero
              do i = 1, nsol
                gc(i, t, n) = gc(i, t, n)*var(4)
                fc(i, t, n) = fc(i, t, n)*var(4)
              end do
            end do
          end if

!                                    SUM FOR EACH KAPPA
          do n = 1, nzero
            do k = 1, nsol
              gck(k, s, n) = 0.0d0
              fck(k, s, n) = 0.0d0
              do j = 1, nsol
                gck(k, s, n) = gck(k, s, n) + gc(k, j, n)
                fck(k, s, n) = fck(k, s, n) + fc(k, j, n)
              end do
            end do
          end do

!                       -----------------------------------
!                       CALCULATE  NORM  AND NORMALIZE TO 1
!                       -----------------------------------
          do k = 1, nsol
            norm = r2drdic(1)*(gck(k,s,1)**2+fck(k,s,1)**2)
          end do

          simp = -1.0d0
          do n = 2, nzero
            simp = -simp
            w = (3.0d0+simp)*r2drdic(n)
            do k = 1, nsol
              norm = norm + w*(gck(k,s,n)**2+fck(k,s,n)**2)
            end do
          end do

          n = nzero
          do k = 1, nsol
            norm = norm - r2drdic(n)*(gck(k,s,n)**2+fck(k,s,n)**2)
          end do
          norm = norm/3.0d0

          do k = 1, nsol
            norm = norm + 0.5d0*r2drdic(1)*(gck(k,s,1)**2+fck(k,s,1)**2)
          end do
          norm = 1.0d0/sqrt(norm)

          do n = 1, nzero
            do k = 1, nsol
              gck(k, s, n) = gck(k, s, n)*norm
              fck(k, s, n) = fck(k, s, n)*norm
            end do
          end do
          if (nzero<jtop) then
            do n = (nzero+1), jtop
              do k = 1, nsol
                gck(k, s, n) = 0.0d0
                fck(k, s, n) = 0.0d0
              end do
            end do
          end if

          call rinit(nrmax, rint)

          do n = 1, jtop
            do k = 1, nsol
              rint(n) = rint(n) + r2drdi(n, im)*(gck(k,s,n)*gck(k,s,n)+fck(k,s &
                ,n)*fck(k,s,n))
            end do
          end do

          call rintsimp(rint, jtop, aux)

! ------------------------------ omit normalization for XRAY calculation
! ------------------------------ to recover old (errounous data) -------
          if (itxray>0) then
            norm = 1.0d0
          else
            norm = 1.0d0/sqrt(aux)
          end if

          do n = 1, max(nzero, jtop)
            do k = 1, nsol
              gck(k, s, n) = gck(k, s, n)*norm
              fck(k, s, n) = fck(k, s, n)*norm
            end do
          end do

!                       -----------------------------------
!                       CALCULATE  CHARGE AND SPIN DENSITY
!                       -----------------------------------

          do n = 1, jws(im)
            do k = 1, nsol
              rhochr(n, it) = rhochr(n, it) + (gck(k,s,n)**2+fck(k,s,n)**2)
              rhospn(n, it) = rhospn(n, it) + (gck(k,s,n)**2*cgd(k)-fck(k,s,n) &
                **2*cgmd(k))
            end do
          end do

          if (nsol>1) then
            do n = 1, jws(im)
              rhospn(n, it) = rhospn(n, it) + gck(1, s, n)*gck(2, s, n)*cgo*2
            end do
          end if


!                       -----------------------------------
!                            CALCULATE  SPIN CHARACTER
!                       -----------------------------------

          w = r2drdic(1)
          sz = 0.0d0
          do k = 1, nsol
            sz = sz + w*(gck(k,s,1)**2*cgd(k)+fck(k,s,1)**2*cgmd(k))
          end do

          simp = -1.0d0
          do n = 2, nzero
            simp = -simp
            w = (3.0d0+simp)*r2drdic(n)
            do k = 1, nsol
              sz = sz + w*(gck(k,s,n)**2*cgd(k)+fck(k,s,n)**2*cgmd(k))
            end do
          end do

          n = nzero
          w = r2drdic(n)
          do k = 1, nsol
            sz = sz + w*(gck(k,s,n)**2*cgd(k)+fck(k,s,n)**2*cgmd(k))
          end do


          if (nsol>1) then

            w = r2drdic(1)
            sz = sz + w*gck(1, s, 1)*gck(2, s, 1)*cgo*2

            simp = -1.0d0
            do n = 2, nzero
              simp = -simp
              w = (3.0d0+simp)*r2drdic(n)
              do k = 1, nsol
                sz = sz + w*gck(1, s, n)*gck(2, s, n)*cgo*2
              end do
            end do

            n = nzero
            w = r2drdic(n)
            sz = sz - w*gck(1, s, n)*gck(2, s, n)*cgo*2

          end if

          sz = sz/3.0d0


!                         ------------------------------
!                         CALCULATE   HYPERFINE - FIELD
!                         ------------------------------

          call corehff(kap1, kap2, mj, s, nsol, bhf, gck, fck, rc, drdic, &
            0.0d0, nzero, nrc)
          if (nucleus/=0) then
            call corehff(kap1, kap2, mj, s, nsol, bhf1, gck, fck, rc, drdic, &
              rnuc, jlim, nrc)
            call corehff(kap1, kap2, mj, s, nsol, bhf2, gck, fck, rc, drovrn, &
              rnuc, jlim, nrc)
            do i = 1, nsol
              do j = 1, nsol
                bhf(i, j) = bhf(i, j) - bhf1(i, j) + bhf2(i, j)
              end do
            end do
          end if !end of nucleus.eq.0

          bsol = 0.0d0
          do j = 1, nsol
            do i = 1, nsol
              bsol = bsol + bhf(i, j)
              bsh = bsh + bhf(i, j)
              bcor(it) = bcor(it) + bhf(i, j)
            end do
          end do
          if (kap1==-1) bcors(it) = bcors(it) + bhf(1, 1)

          ecortab(ic, it) = ec

!     ------------------
!     SPLIT HFF-FIELD
!     ------------------
          if (ismqhfi==1) then
            call hffcore(rnuc, nzero, kap1, kap2, nsol, mj, gck, fck, nrc, &
              shf, s, nmemax, nkmmax, rc, drdic, sdia, smdia, soff, smoff, &
              qdia, qoff, qmdia, qmoff, nucleus, jlim)

            do k = 1, nmemax
              split(k) = 0.0d0
              do j = 1, nsol
                do i = 1, nsol
                  split(k) = split(k) + shf(i, j, k)
                  split1(k) = split1(k) + shf(i, j, k)
                  split2(k, it) = split2(k, it) + shf(i, j, k)
                end do
              end do
            end do
            do k = 1, nmemax
              if (kap1==-1) split3(k, it) = split3(k, it) + shf(1, 1, k)
            end do
          end if
!MBE

!---------------------------------------------------- l-shell UNCOMPLETE
          if (ish>=nsh) then
!----------------------------------------------------- l-shell completed
            if (itxray==0 .and. iprint>0 .and. (t_inc%i_write>0)) then
              write (1337, 280) itprt, nqn, txtl(l), txtk(iabs(kap(s))), &
                (2*muem05+1), kap(s), iter, ec, bsol*.001d0, bsh*.001d0
              if (ismqhfi==1) then
                do k = 1, nmemax
                  write (1337, 290) txtb(k), split(k)*.001d0, split1(k)*.001d0
                end do
                write (1337, 300) 'total error in %', &
                  100.0d0*(1.0d0-split(4)/split(5))
              end if
            end if
!                              ----------------------------
!                              CHECK CONSISTENCY OF RESULTS
!                              ----------------------------
            if (l/=0) then
              ic1 = ic - nsh + 1
              ic2 = ic
              if (ecortab(ic2,it)>=ecortab(ic1,it)) then
                imin = ic1
                vz = +1.0d0
              else
                imin = ic2
                vz = -1.0d0
              end if
              iflag = 0
              ii = 1
              do i = ic1 + 1, ic2, 2
                if (vz*(ecortab(i,it)-ecortab(i-ii,it))<0.0) iflag = 1
                ii = 2
              end do
              if (ecortab(ic1+2,it)>ecortab(imin,it)) iflag = 1
              do i = ic1 + 4, ic2 - 1, 2
                if (ecortab(i,it)>ecortab(imin,it)) iflag = 1
                if (vz*(ecortab(i,it)-ecortab(i-ii,it))>0.0) iflag = 1
              end do

              if (ferro .and. (iflag==1)) then
                if (t_inc%i_write>0) write (1337, 230)
                scaleb = .true.
                if (loop==1) then
                  loop = 2
                  if (t_inc%i_write>0) write (1337, 240) itprt
                  go to 100
                end if
              end if

            end if
          else if (itxray==0 .and. iprint>0) then
            if (t_inc%i_write>0) write (1337, 280) itprt, nqn, txtl(l), &
              txtk(iabs(kap(s))), (2*muem05+1), kap(s), iter, ec, bsol*.001d0
            if (ismqhfi==1) then
              do k = 1, nmemax
                if (t_inc%i_write>0) write (1337, 290) txtb(k), &
                  split(k)*.001d0
              end do
              if (t_inc%i_write>0) write (1337, 300) 'total error in %', &
                100.0d0*(1.0d0-split(4)/split(5))
            end if
          end if
!-----------------------------------------------------------------------

          if (iprint>=1 .and. (t_inc%i_write>0)) write (1337, 310)((bhf(i, &
            j)*.001d0,i=1,nsol), j=1, nsol)


!                          --------------------------------
!                            IF THE SWITCH CHECK IS SET:
!                            RECALCULATE THE EIGENVALUE
!                          USING THE CONVENTIONAL ALGORITHM
!                          --------------------------------
          if (check) then
            ecc = 0.95d0*ec
150         continue
            call coredir(it, ctl(it,ilc), ecc, l, mj, 'OUT', vv, bb, rc, &
              drdic, dovrc, nmatch, nzero, gc, fc, d_p, dq, wp, wq, pow, qow, &
              piw, qiw, cgd, cgmd, cgo, nrc, z(it), nucleus)
            call coredir(it, ctl(it,ilc), ecc, l, mj, 'INW', vv, bb, rc, &
              drdic, dovrc, nmatch, nzero, gc, fc, d_p, dq, wp, wq, pow, qow, &
              piw, qiw, cgd, cgmd, cgo, nrc, z(it), nucleus)

            norm = pow(s, s)/piw(s, s)
            do n = nmatch, nzero
              gc(s, s, n) = gc(s, s, n)*norm
              fc(s, s, n) = fc(s, s, n)*norm
            end do

            norm = 0.0d0
            do n = 3, nzero, 2
              norm = norm + r2drdic(n)*(gc(s,s,n)**2+fc(s,s,n)**2) + &
                4.d0*r2drdic(n-1)*(gc(s,s,n-1)**2+fc(s,s,n-1)**2) + &
                r2drdic(n-2)*(gc(s,s,n-2)**2+fc(s,s,n-2)**2)
            end do
            norm = norm/3.0d0

            lcp1 = min(nlmax, l+1)
            dec = pow(s, s)*(qow(s,s)-rc(nmatch)*ctl(it,lcp1)*fc(s,s,nmatch))/ &
              norm
            ecc = ecc + dec
            if (abs(dec/ecc)>tolvar) go to 150
            if (t_inc%i_write>0) write (1337, '(7X,''CHECK-E:'',10X,F12.5,/)') &
              ecc
          end if


! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!           STORE CORE WAVE FUNCTIONS IF REQUIRED
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
          if (itxray/=0) then
            if (nsol==2) then
              if (s==2) then
                icst = l + 1 + muem05
              else
                icst = l + 1 + muem05 + 2*l + 1
              end if
            else if (mj<0) then
              icst = 2*l + 1
            else
              icst = 4*l + 2
            end if

            mm05cor(icst) = muem05
            nkpcor(icst) = nsol
            kapcor(icst) = kap(s)
            ikmcor(icst, 1) = ikapmue(kap(s), muem05)
            izero(icst) = min(nzero, jws(im))
            szcor(icst) = sz
            ecor(icst) = ec

            do n = 1, izero(icst)
              gcor(n, 1, icst) = gck(s, s, n)
              fcor(n, 1, icst) = fck(s, s, n)
            end do
            if (nsol==2) then
              do n = 1, izero(icst)
                gcor(n, 2, icst) = gck(t, s, n)
                fcor(n, 2, icst) = fck(t, s, n)
              end do
              ikmcor(icst, 2) = ikapmue(kap(t), muem05)
            end if
          end if
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


        end do
! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

      end do
! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
160 end do
! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL


    if (itxray==0 .and. iprint>0 .and. (t_inc%i_write>0)) write (1337, 250) &
      bcor(it)*.001d0, bcors(it)*.001d0
    if (ismqhfi==1) then
      do n = 1, nmemax
        if (t_inc%i_write>0) write (1337, 260) split2(n, it)*.001d0, &
          split3(n, it)*.001d0
      end do
      if (t_inc%i_write>0) write (1337, 300) 'total error', &
        100.0d0*(1.0d0-split2(4,it)/split2(5,it))
      if (t_inc%i_write>0) write (1337, 300) 'total error', &
        100.0d0*(1.0d0-split3(4,it)/split3(5,it))
    end if

    if ((itxray==0) .or. (it==itxray)) then
      do n = 1, jtop
        rint(n) = rhochr(n, it)*r2drdi(n, im)
      end do
      call rintsimp(rint, jtop, aux)
      if (iprint>-2 .and. (t_inc%i_write>0)) write (1337, 330) 'charge', &
        itprt, aux
      do n = 1, jtop
        rint(n) = rhospn(n, it)*r2drdi(n, im)
      end do
      call rintsimp(rint, jtop, aux)
      if (iprint>-2 .and. (t_inc%i_write>0)) write (1337, 330) ' spin ', &
        itprt, aux
    end if

  end do
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

170 format (/, 10x, 'potential is not exchange split ')
180 format (/, '  ATOM   IT      : ', i5, /, '  ATOMIC NUMBER  : ', i5, /, /, &
    '  IT', 10x, 'MUE  KAP ITER    ENERGY       B  (k-Gauss)  ')
190 format ('  IT=', i2, ' NQN=', i2, ' L=', i2, '  NZERO set to  (NRC-1) =', &
    i4)
200 format (/, /, '  STOP IN <<CORE>>', /, '  IT=', i2, ' NQN=', i2, ' L=', &
    i2, /, '  no matching-radius found for  EC=', f10.3)
210 format (' LOOP    =  ', i3, ' BSCL=', f10.5, /, ' E=', f25.16, ' VAR  ', &
    4e11.4, /, 17x, ' CORR ', 4e11.4, /, 17x, ' ERR  ', 4e11.4)
220 format (' iteration not converged after', i3, ' steps !', /, &
    ' parameters:', 4e18.10, /, ' last corr.:', 4e18.10, /, ' last error:', &
    4e18.10)
230 format (' >>> check data E(KAP,MJ) should be monotonous ', &
    ' and  E(+L,MJ) < E(-L-1,MJ) ', /, /)
240 format (' all core states for atom type ', i2, ' will be recalculated ', &
    /, ' with the spin dependent exchange field gradually switched on')
250 format (2x, 57('-'), /, 42x, f17.3, /, 38x, '(S) ', f17.3, /, 2x, 57('*'), &
    /, /)
260 format (2x, 57('-'), /, 37x, f17.3, /, 33x, '(S) ', f17.3, /, 2x, 57('*'), &
    /, /)
270 format (a, i1, a, 5i3, a, 2i2, a)
280 format (2i4, a1, a3, i3, '/2', 2i4, 2x, f15.8, f17.3, :, f17.3, /)
290 format (a, :, 32x, f17.3, :, f17.3, /)
300 format (a, :, 37x, f6.3)
310 format (37x, f17.3)
320 format (/, ' IT=', i2, '  NQN=', i2, '  L=', i2, '  KAP=', i2, '  MJ=', &
    i2, '/2    IC=', i3, '  ISH=', i2, /, ' E(', i2, ')   =', f15.5, /, &
    ' NMATCH  =', i5, '    R=', f10.5, /, ' NZERO   =', i5, '    R=', f10.5, &
    /, ' NODES   =', i5, '  RAT=', e11.4)
330 format (' integrated core ', a, ' density for atom type ', i4, ':', f12.8)
end subroutine
