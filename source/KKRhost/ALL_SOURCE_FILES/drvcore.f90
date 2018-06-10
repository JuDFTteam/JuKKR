SUBROUTINE drvcore(iprint,itprt,lcore,ncore,cscl,vtin,btin,rin,  &
    a,b,drdiin,r2drdiin,zat,jws,ishift,rhoc,  &
    ecorerel,nkcore,kapcore,ecore,lmaxd,irmd)
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
IMPLICIT NONE

! PARAMETER definitions

INTEGER NTMAX,NMMAX,NCSTMAX,NMEMAX,NLMAX,NKMMAX
PARAMETER (NTMAX=1,NMMAX=1,NCSTMAX=6,NMEMAX=5,NLMAX=5, &
           NKMMAX=2*NLMAX**2) ! NLMAX should be >= LCOREMAX + 1
INTEGER NRMAX
PARAMETER (NRMAX = 900)
DOUBLE PRECISION DZERO
PARAMETER (DZERO=0.0D0)

! Dummy arguments
DOUBLE PRECISION A,B
INTEGER IPRINT,IRMD,ITPRT,LMAXD,NCORE,ISHIFT

!obs: here, in contrast to DRVRHO, one works with local
!arrays, since they have to be set up as far as NRMAX (see CORE)

DOUBLE PRECISION VTIN(IRMD),BTIN(IRMD)
DOUBLE PRECISION DRDIIN(IRMD),R2DRDIIN(IRMD)
DOUBLE PRECISION RIN(IRMD),CSCL(LMAXD+1)

DOUBLE PRECISION ECORE(20,2),ECOREREL(20*2)
INTEGER KAPCORE(20*2),LCORE(20,2),NKCORE(20)
INTEGER ZAT(NTMAX),JWS(NMMAX)
DOUBLE PRECISION RHOC(IRMD,2)

! Local variables
DOUBLE PRECISION BCOR(NTMAX),BCORS(NTMAX), &
       EA,ECOR(NCSTMAX),ECORTAB(120,NTMAX), &
       FCOR(NRMAX,2,NCSTMAX),GCOR(NRMAX,2,NCSTMAX),QDIA(NKMMAX), &
       QMDIA(NKMMAX),QMOFF(NKMMAX),QOFF(NKMMAX), &
       RHOCHR(NRMAX,NTMAX),RHOSPN(NRMAX,NTMAX), &
       SDIA(NKMMAX),SMDIA(NKMMAX),SMOFF(NKMMAX),SOFF(NKMMAX), &
       SZCOR(NCSTMAX)
DOUBLE PRECISION VT(NRMAX,NTMAX),BT(NRMAX,NTMAX),CTL(NTMAX,NLMAX)
DOUBLE PRECISION DRDI(NRMAX,NMMAX),R2DRDI(NRMAX,NMMAX)
DOUBLE PRECISION R(NRMAX,NMMAX)
INTEGER I,ICALL,IKMCOR(NCSTMAX,2),IMT(NTMAX),IP,ISMQHFI, &
        IT,ITXRAY,IZERO(NCSTMAX),J,KAPCOR(NCSTMAX), &
        LCXRAY(NTMAX),MM05COR(NCSTMAX),NCORT(NTMAX),NCXRAY(NTMAX), &
        NKPCOR(NCSTMAX),NT,NUCLEUS
INTEGER LCOREMAX

SAVE IMT,ITXRAY,NT,NUCLEUS,QDIA,QMDIA,QMOFF,QOFF,SDIA,SMDIA, &
     SMOFF,SOFF

DATA NCXRAY/NTMAX*0/,LCXRAY/NTMAX*0/,ISMQHFI/0/

DATA ICALL/0/

icall = icall + 1

!=======================================================================
!       initialise relativistic and dummy variables and SAVE them
!=======================================================================
IF ( icall == 1 ) THEN
  
  IF ( lmaxd > nlmax-1 ) THEN
    WRITE (6,*) ' LMAXD = ',lmaxd,' > NLMAX-1 = ',nlmax - 1
    STOP ' Increase NLMAX in < DRVCORE > '
  END IF
  
  IF ( irmd > nrmax ) THEN
    WRITE (6,*) ' IRMD = ',irmd,' > NRMAX = ',nrmax
    WRITE (6,*) ' Increase NRMAX in < sprkkr_rmesh.dim > '
    STOP ' In < DRVCORE > '
  END IF
  
  itxray = 0
  
  DO it = 1,ntmax
    imt(it) = 1
  END DO
  
  nt = 1
  nucleus = 0
  
  DO it = 1,nkmmax
    sdia(it) = dzero
    smdia(it) = dzero
    soff(it) = dzero
    smoff(it) = dzero
    
    qdia(it) = dzero
    qmdia(it) = dzero
    qoff(it) = dzero
    qmoff(it) = dzero
  END DO
  
END IF                    ! ICALL.EQ.1
!=======================================================================

! --> fill up CTL array for the case of core states with higher L values
!     than those used in the valence band

lcoremax = 0
DO it = 1,ncore
  j = lcore(it,1)
  lcoremax = MAX(lcoremax,j)
END DO
IF ( lcoremax > nlmax-1 ) THEN
  WRITE (6,*) ' LCOREMAX = ',lcoremax,' > NLMAX-1 = ',nlmax - 1
  STOP ' Increase NLMAX in < DRVCORE > '
END IF
DO j = 1,lmaxd+1
  ctl(1,j) = cscl(j)
END DO
IF ( lcoremax > 0 ) THEN
  DO j = lcoremax+1,nlmax
    ctl(1,j) = cscl(lcoremax)
  END DO
END IF

CALL dcopy(jws(1),vtin,1,vt(1,1),1)
CALL dcopy(jws(1),btin,1,bt(1,1),1)
CALL dcopy(jws(1),rin,1,r(1,1),1)
CALL dcopy(jws(1),drdiin,1,drdi(1,1),1)
CALL dcopy(jws(1),r2drdiin,1,r2drdi(1,1),1)

DO j = jws(1) + 1,nrmax
  ea = DEXP(a*DBLE(j+ishift-1)) ! corrected from (J-1) 07.05.2004
  r(j,1) = b*(ea-1D0)
  drdi(j,1) = a*b*ea
  r2drdi(j,1) = r(j,1)*r(j,1)*drdi(j,1)
  vt(j,1) = 0D0
  bt(j,1) = 0D0
END DO

ncort(1) = 0  ! no. of core electrons = no. of diff. core states

DO it = 1,ncore
  ncort(1) = ncort(1) + 2*(2*lcore(it,1)+1)
END DO

CALL core(iprint,itprt,nt,ncort,ctl,vt,bt,zat,nucleus,r,r2drdi,  &
    drdi,jws,imt,rhochr,rhospn,ecortab,gcor,fcor,ecor,szcor,  &
    kapcor,mm05cor,nkpcor,ikmcor,izero,ncxray,lcxray,itxray,  &
    bcor,bcors,sdia,smdia,soff,smoff,qdia,qoff,qmdia,qmoff,  &
    nkmmax,nmemax,ismqhfi,ntmax,nrmax,nmmax,ncstmax,nlmax)

CALL rinit(2*irmd,rhoc(1,1))

DO i = 1,jws(1)
  ip = i + ishift
  rhoc(ip,2) = (rhochr(i,1)+rhospn(i,1))*0.5D0*(r(i,1)**2)
  rhoc(ip,1) = (rhochr(i,1)-rhospn(i,1))*0.5D0*(r(i,1)**2)
END DO

CALL sumecore(ncore,lcore(1,1),ecortab(1,1),nkcore,ecorerel,ecore, kapcore)

END SUBROUTINE drvcore
!*==sumecore.f    processed by SPAG 6.05Rc at 11:35 on 10 May 2004

SUBROUTINE sumecore(ncore,lcore,ecortab,nkcore,ecorerel,ecore, kapcore)
IMPLICIT NONE

! Dummy arguments
INTEGER NCORE
DOUBLE PRECISION ECORE(20,2),ECOREREL(20*2),ECORTAB(*)
INTEGER KAPCORE(20*2),LCORE(*),NKCORE(*)

! Local variables
DOUBLE PRECISION DBLE
INTEGER I,IC,ICREL,JREL,KFG(4),L,LMP1,LMXC,LP1,MUEM05,NC,NMAX,NN, &
        NSOL,WGT(2)
DOUBLE PRECISION MJ
INTRINSIC ABS

! --> find the principal quantum numbers

DO ic = 1,4
  kfg(ic) = 0
END DO
DO ic = 1,ncore
  IF ( lcore(ic) == 0 ) kfg(1) = kfg(1) + 1
  IF ( lcore(ic) == 1 ) kfg(2) = kfg(2) + 1
  IF ( lcore(ic) == 2 ) kfg(3) = kfg(3) + 1
  IF ( lcore(ic) == 3 ) kfg(4) = kfg(4) + 1
END DO

IF ( kfg(2) /= 0 ) kfg(2) = kfg(2) + 1
IF ( kfg(3) /= 0 ) kfg(3) = kfg(3) + 2
IF ( kfg(4) /= 0 ) kfg(4) = kfg(4) + 3

lmxc = 0
IF ( kfg(2) /= 0 ) lmxc = 1
IF ( kfg(3) /= 0 ) lmxc = 2
IF ( kfg(4) /= 0 ) lmxc = 3

lmp1 = lmxc + 1
nc = 0
DO lp1 = 1,lmp1
  l = lp1 - 1
  nmax = kfg(lp1)
  DO nn = lp1,nmax
    nc = nc + 1
    ecorerel(nc) = 0.0D0
    ecorerel(nc+20) = 0.0D0
    
!  ECOREREL(NC..NC+20) = 1st/2nd value of \kappa for current l
    
    
! --> icrel is pointing a core-state (n,l) = (nn,l) from the
!     relativistic-routine sequence 1s,2s,2p,3s,3p,... to the
!     ECOREREL array  1s,2s,3s,2p,3p,...
    
!     (nn,l) -> nn*(nn-1)*(2*nn-1)/3 + 2*l**2
    
    icrel = nn*(nn-1)*(2*nn-1)/3 + 2*l**2
    jrel = 0
    wgt(1) = 0
    wgt(2) = 0
    DO muem05 = -l - 1, + l
      mj = muem05 + 0.5D0
      IF ( ABS(mj) > l ) THEN
        nsol = 1
      ELSE
        nsol = 2
      END IF
      DO i = 1,nsol
        wgt(i) = wgt(i) + 1
        jrel = jrel + 1
        ecorerel((i-1)*20+nc) = ecorerel((i-1)*20+nc) + ecortab(icrel+jrel)
      END DO
    END DO
    nkcore(nc) = 1
    IF ( l /= 0 ) nkcore(nc) = 2
    kapcore(nc) = -l - 1
    kapcore(nc+20) = l
    
    DO i = 1,nkcore(nc)
      ecorerel((i-1)*20+nc) = ecorerel((i-1)*20+nc) /DBLE(wgt(i))
    END DO
    
! --> update the array ECORE(1..NCORE,UP/DOWN) as
    
!         ECORE(L,SIGMA) = 1/(2*L+1) *
!              SUM_KAPPA(L) SUM_(MUE,SIGN(MUE)=SIGN(SIGMA)) ECORTAB(KAP,MUE)
!     i.e., states with negative MUE are added to SPIN DOWN, those with
!     positive MUE to SPIN UP.
    
!     ECORE is used later in calculating the total energy
    
    ecore(nc,1) = 0.0D0
    ecore(nc,2) = 0.0D0
    DO i = 1,2*l + 1
      ecore(nc,1) = ecore(nc,1) + ecortab(icrel+i)
      ecore(nc,2) = ecore(nc,2) + ecortab(icrel+2*l+1+i)
    END DO
    DO i=1,2
      ecore(nc,i) = ecore(nc,i)/DBLE(2*l+1)
    END DO
  END DO
END DO
END SUBROUTINE sumecore
