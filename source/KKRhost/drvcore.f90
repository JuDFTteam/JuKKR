!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_drvcore

private
public :: drvcore

contains

  !-------------------------------------------------------------------------------
  !> Summary: Driver relativistic core routine
  !> Author: 
  !> Category: KKRhost, core-electrons
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Driving routine to call relativistic < CORE > routine           
  !> counterpart of < COREL > of the non- or scalar-relativistic mode
  !>
  !> ATTENTION: all the variables connected with hyperfine fields are
  !>            OFF                                                  
  !>
  !> The non-relativistic variables NCORE and LCORE are used as read 
  !> in from the potential file; they are not modified and are again 
  !> written out for the next iteration ( routine <RITES>)           
  !>
  !> Relativistic variables passed outside this routine:             
  !>    NKCORE(1..NCORE)       = number of KAPPA values for a given  
  !>                             (n,l) core state                    
  !>    KAPCORE(1..NCORE,1..NCORE) = the (maximum 2) values of KAPPA 
  !>    ECOREREL(1..NCORE,1..NCORE) = for a given (n,l) state the core
  !>                              energies corresponding first/second
  !>                              KAPPA value, AVERAGED over \mu's   
  !>                              These values are written out to the
  !>                              potential file (routine <RITES>),  
  !>                              but the read in (routine <STARTB1>)
  !>                              updates the ECORE array            
  !>     Please note that ALL the core energies (also \mu-resolved)  
  !>     are output by <CORE> routine but not passed out of here     
  !>
  !>    ECORE(1..NCORE,1..2) = this (non- or scalar relativistic)    
  !>                           variable is updated here to be used in
  !>                           calculating/printing out the spin-    
  !>                           resolved energies (see <ESPCB> )      
  !>                           A SUMMATION is done here:             
  !>        ECORE(L,1/2) = SUM_{\kappa=-L-1,L} E(\kappa,\mu)         
  !>                       /(2*L+1)                                  
  !>                           with negative \mu's for E(*,1) and    
  !>                           positive \mu's for E(*,2)             
  !>                                                                 
  !>                           V. Popescu July/2002
  !-------------------------------------------------------------------------------
  subroutine drvcore(iprint, itprt, lcore, ncore, cscl, vtin, btin, rin, a, b, drdiin, r2drdiin, zat_in, jws_in, ishift, rhoc, ecorerel, nkcore, kapcore, ecore, lmaxd, irmd)

    use :: mod_datatypes, only: dp
    use :: mod_core, only: core
    use :: mod_rinit, only: rinit
    implicit none

    ! PARAMETER definitions
    integer :: ntmax, nmmax, ncstmax, nmemax, nlmax, nkmmax
    parameter (ntmax=1, nmmax=1, ncstmax=6, nmemax=5, nlmax=5, nkmmax=2*nlmax**2) ! NLMAX should be >= LCOREMAX + 1
    integer :: nrmax
    parameter (nrmax=900)
    real (kind=dp) :: dzero
    parameter (dzero=0.0e0_dp)

    ! Dummy arguments
    real (kind=dp) :: a, b
    integer :: iprint, irmd, itprt, lmaxd, ncore, ishift

    ! obs: here, in contrast to DRVRHO, one works with local
    ! arrays, since they have to be set up as far as NRMAX (see CORE)

    real (kind=dp) :: vtin(irmd), btin(irmd)
    real (kind=dp) :: drdiin(irmd), r2drdiin(irmd)
    real (kind=dp) :: rin(irmd), cscl(lmaxd+1)

    real (kind=dp) :: ecore(20, 2), ecorerel(20*2)
    integer :: kapcore(20*2), lcore(20, 2), nkcore(20)
    integer :: zat_in, zat(ntmax), jws_in, jws(nmmax)
    real (kind=dp) :: rhoc(irmd, 2)

    ! Local variables
    real (kind=dp) :: bcor(ntmax), bcors(ntmax), ea, ecor(ncstmax), ecortab(120, ntmax), fcor(nrmax, 2, ncstmax), gcor(nrmax, 2, ncstmax), qdia(nkmmax), qmdia(nkmmax), &
      qmoff(nkmmax), qoff(nkmmax), rhochr(nrmax, ntmax), rhospn(nrmax, ntmax), sdia(nkmmax), smdia(nkmmax), smoff(nkmmax), soff(nkmmax), szcor(ncstmax)
    real (kind=dp) :: vt(nrmax, ntmax), bt(nrmax, ntmax), ctl(ntmax, nlmax)
    real (kind=dp) :: drdi(nrmax, nmmax), r2drdi(nrmax, nmmax)
    real (kind=dp) :: r(nrmax, nmmax)
    integer :: i, icall, ikmcor(ncstmax, 2), imt(ntmax), ip, ismqhfi, it, itxray, izero(ncstmax), j, kapcor(ncstmax), lcxray(ntmax), mm05cor(ncstmax), ncort(ntmax), ncxray(ntmax), &
      nkpcor(ncstmax), nt, nucleus
    integer :: lcoremax

    save :: imt, itxray, nt, nucleus, qdia, qmdia, qmoff, qoff, sdia, smdia, smoff, soff

    data ncxray/ntmax*0/, lcxray/ntmax*0/, ismqhfi/0/

    data icall/0/

    icall = icall + 1
    zat(1) = zat_in
    jws(1) = jws_in

    ! =======================================================================
    ! initialise relativistic and dummy variables and SAVE them
    ! =======================================================================
    if (icall==1) then

      if (lmaxd>nlmax-1) then
        write (6, *) ' LMAXD = ', lmaxd, ' > NLMAX-1 = ', nlmax - 1
        stop ' Increase NLMAX in < DRVCORE > '
      end if

      if (irmd>nrmax) then
        write (6, *) ' IRMD = ', irmd, ' > NRMAX = ', nrmax
        write (6, *) ' Increase NRMAX in < sprkkr_rmesh.dim > '
        stop ' In < DRVCORE > '
      end if

      itxray = 0

      do it = 1, ntmax
        imt(it) = 1
      end do

      nt = 1
      nucleus = 0

      do it = 1, nkmmax
        sdia(it) = dzero
        smdia(it) = dzero
        soff(it) = dzero
        smoff(it) = dzero

        qdia(it) = dzero
        qmdia(it) = dzero
        qoff(it) = dzero
        qmoff(it) = dzero
      end do

    end if                         ! ICALL.EQ.1
    ! =======================================================================

    ! --> fill up CTL array for the case of core states with higher L values
    ! than those used in the valence band

    lcoremax = 0
    do it = 1, ncore
      j = lcore(it, 1)
      lcoremax = max(lcoremax, j)
    end do
    if (lcoremax>nlmax-1) then
      write (6, *) ' LCOREMAX = ', lcoremax, ' > NLMAX-1 = ', nlmax - 1
      stop ' Increase NLMAX in < DRVCORE > '
    end if
    do j = 1, lmaxd + 1
      ctl(1, j) = cscl(j)
    end do
    if (lcoremax>0) then
      do j = lcoremax + 1, nlmax
        ctl(1, j) = cscl(lcoremax)
      end do
    end if

    call dcopy(jws(1), vtin, 1, vt(1,1), 1)
    call dcopy(jws(1), btin, 1, bt(1,1), 1)
    call dcopy(jws(1), rin, 1, r(1,1), 1)
    call dcopy(jws(1), drdiin, 1, drdi(1,1), 1)
    call dcopy(jws(1), r2drdiin, 1, r2drdi(1,1), 1)

    do j = jws(1) + 1, nrmax
      ea = exp(a*real(j+ishift-1,kind=dp)) ! corrected from (J-1) 07.05.2004
      r(j, 1) = b*(ea-1e0_dp)
      drdi(j, 1) = a*b*ea
      r2drdi(j, 1) = r(j, 1)*r(j, 1)*drdi(j, 1)
      vt(j, 1) = 0e0_dp
      bt(j, 1) = 0e0_dp
    end do

    ncort(1) = 0                   ! no. of core electrons = no. of diff. core
    ! states

    do it = 1, ncore
      ncort(1) = ncort(1) + 2*(2*lcore(it,1)+1)
    end do

    call core(iprint, itprt, nt, ncort, ctl, vt, bt, zat, nucleus, r, r2drdi, drdi, jws, imt, rhochr, rhospn, ecortab, gcor, fcor, ecor, szcor, kapcor, mm05cor, nkpcor, ikmcor, &
      izero, ncxray, lcxray, itxray, bcor, bcors, sdia, smdia, soff, smoff, qdia, qoff, qmdia, qmoff, nkmmax, nmemax, ismqhfi, ntmax, nrmax, nmmax, ncstmax, nlmax)

    call rinit(2*irmd, rhoc(1,1))

    do i = 1, jws(1)
      ip = i + ishift
      rhoc(ip, 2) = (rhochr(i,1)+rhospn(i,1))*0.5e0_dp*(r(i,1)**2)
      rhoc(ip, 1) = (rhochr(i,1)-rhospn(i,1))*0.5e0_dp*(r(i,1)**2)
    end do

    call sumecore(ncore, lcore(1,1), ecortab(1,1), nkcore, ecorerel, ecore, kapcore)

  end subroutine drvcore

  !-------------------------------------------------------------------------------
  !> Summary: Summation of core-electron energies
  !> Author: 
  !> Category: KKRhost, core-electrons
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> 
  !-------------------------------------------------------------------------------
  subroutine sumecore(ncore, lcore, ecortab, nkcore, ecorerel, ecore, kapcore)
    use :: mod_datatypes, only: dp
    implicit none

    ! Dummy arguments
    integer :: ncore
    real (kind=dp) :: ecore(20, 2), ecorerel(20*2), ecortab(*)
    integer :: kapcore(20*2), lcore(*), nkcore(*)

    ! Local variables
    real (kind=dp) :: dble
    integer :: i, ic, icrel, jrel, kfg(4), l, lmp1, lmxc, lp1, muem05, nc, nmax, nn, nsol, wgt(2)
    real (kind=dp) :: mj
    intrinsic :: abs

    ! --> find the principal quantum numbers

    do ic = 1, 4
      kfg(ic) = 0
    end do
    do ic = 1, ncore
      if (lcore(ic)==0) kfg(1) = kfg(1) + 1
      if (lcore(ic)==1) kfg(2) = kfg(2) + 1
      if (lcore(ic)==2) kfg(3) = kfg(3) + 1
      if (lcore(ic)==3) kfg(4) = kfg(4) + 1
    end do

    if (kfg(2)/=0) kfg(2) = kfg(2) + 1
    if (kfg(3)/=0) kfg(3) = kfg(3) + 2
    if (kfg(4)/=0) kfg(4) = kfg(4) + 3

    lmxc = 0
    if (kfg(2)/=0) lmxc = 1
    if (kfg(3)/=0) lmxc = 2
    if (kfg(4)/=0) lmxc = 3

    lmp1 = lmxc + 1
    nc = 0
    do lp1 = 1, lmp1
      l = lp1 - 1
      nmax = kfg(lp1)
      do nn = lp1, nmax
        nc = nc + 1
        ecorerel(nc) = 0.0e0_dp
        ecorerel(nc+20) = 0.0e0_dp

        ! ECOREREL(NC..NC+20) = 1st/2nd value of \kappa for current l


        ! --> icrel is pointing a core-state (n,l) = (nn,l) from the
        ! relativistic-routine sequence 1s,2s,2p,3s,3p,... to the
        ! ECOREREL array  1s,2s,3s,2p,3p,...

        ! (nn,l) -> nn*(nn-1)*(2*nn-1)/3 + 2*l**2

        icrel = nn*(nn-1)*(2*nn-1)/3 + 2*l**2
        jrel = 0
        wgt(1) = 0
        wgt(2) = 0
        do muem05 = -l - 1, +l
          mj = muem05 + 0.5e0_dp
          if (abs(mj)>l) then
            nsol = 1
          else
            nsol = 2
          end if
          do i = 1, nsol
            wgt(i) = wgt(i) + 1
            jrel = jrel + 1
            ecorerel((i-1)*20+nc) = ecorerel((i-1)*20+nc) + ecortab(icrel+jrel)
          end do
        end do
        nkcore(nc) = 1
        if (l/=0) nkcore(nc) = 2
        kapcore(nc) = -l - 1
        kapcore(nc+20) = l

        do i = 1, nkcore(nc)
          ecorerel((i-1)*20+nc) = ecorerel((i-1)*20+nc)/dble(wgt(i))
        end do

        ! --> update the array ECORE(1..NCORE,UP/DOWN) as

        ! ECORE(L,SIGMA) = 1/(2*L+1) *
        ! SUM_KAPPA(L) SUM_(MUE,SIGN(MUE)=SIGN(SIGMA)) ECORTAB(KAP,MUE)
        ! i.e., states with negative MUE are added to SPIN DOWN, those with
        ! positive MUE to SPIN UP.

        ! ECORE is used later in calculating the total energy

        ecore(nc, 1) = 0.0e0_dp
        ecore(nc, 2) = 0.0e0_dp
        do i = 1, 2*l + 1
          ecore(nc, 1) = ecore(nc, 1) + ecortab(icrel+i)
          ecore(nc, 2) = ecore(nc, 2) + ecortab(icrel+2*l+1+i)
        end do
        do i = 1, 2
          ecore(nc, i) = ecore(nc, i)/dble(2*l+1)
        end do
      end do
    end do
  end subroutine sumecore

end module mod_drvcore
