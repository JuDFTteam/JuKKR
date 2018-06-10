SUBROUTINE core(iprint,itprt,nt,ncort,ctl,vt,bt,z,nucleus,  &
    r,r2drdi,drdi,jws,imt,rhochr,rhospn, ecortab,gcor,fcor,ecor,szcor,  &
    kapcor,mm05cor,nkpcor,ikmcor,izero,ncxray,lcxray,  &
    itxray,bcor,bcors,sdia,smdia,soff,smoff,qdia,qoff,  &
    qmdia,qmoff,nkmmax,nmemax,ismqhfi,ntmax,nrmax, nmmax,ncstmax,nlmax)
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
use mod_types, only: t_inc
implicit none


! PARAMETER definitions
DOUBLE PRECISION UNEND,TOLVAR,TRYMIX,DVSTEP
PARAMETER (UNEND=600.0D0,TOLVAR=1.0D-6,TRYMIX=0.01D0, &
           DVSTEP=0.01D0)
INTEGER ITERMAX,NLSHELLMAX
PARAMETER (ITERMAX=200,NLSHELLMAX=15)

! Dummy arguments
INTEGER IPRINT,ISMQHFI,ITPRT,ITXRAY,NCSTMAX,NKMMAX,NLMAX,NMEMAX, &
        NMMAX,NRMAX,NT,NTMAX,NUCLEUS
DOUBLE PRECISION BCOR(NTMAX),BCORS(NTMAX),BT(NRMAX,NTMAX), &
       CTL(NTMAX,NLMAX), &
       DRDI(NRMAX,NMMAX),ECOR(NCSTMAX),ECORTAB(120,NTMAX), &
       FCOR(NRMAX,2,NCSTMAX),GCOR(NRMAX,2,NCSTMAX),QDIA(NKMMAX), &
       QMDIA(NKMMAX),QMOFF(NKMMAX),QOFF(NKMMAX),R(NRMAX,NMMAX), &
       R2DRDI(NRMAX,NMMAX),RHOCHR(NRMAX,NTMAX),RHOSPN(NRMAX,NTMAX) &
       ,SDIA(NKMMAX),SMDIA(NKMMAX),SMOFF(NKMMAX),SOFF(NKMMAX), &
       SZCOR(NCSTMAX),VT(NRMAX,NTMAX)
INTEGER IKMCOR(NCSTMAX,2),IMT(NTMAX),IZERO(NCSTMAX),JWS(NMMAX), &
        KAPCOR(NCSTMAX),LCXRAY(NTMAX),MM05COR(NCSTMAX), &
        NCXRAY(NTMAX),NKPCOR(NCSTMAX),Z(NTMAX),NCORT(NTMAX)

! Local variables
DOUBLE PRECISION AUX,BB(NRMAX*2),BHF(2,2),BHF1(2,2),BHF2(2,2), &
       BSH,BSOL,BSUM, &
       CGD(2),CGMD(2),CGO,DEC,DEDV(4,4),DOVRC(NRMAX*2), &
       DP(2,2,NRMAX*2),DQ(2,2,NRMAX*2),DRDIC(NRMAX*2), &
       DROVRN(2*NRMAX),DV(4),DVDE(4,4),EC,ECC,ELIM,ERR(4), &
       ERRNEW(4),FC(2,2,NRMAX*2),FCK(2,2,NRMAX*2),GC(2,2,NRMAX*2), &
       GCK(2,2,NRMAX*2),MJ,NIW(2),NORM,NOW(2),PIW(2,2),POW(2,2), &
       QIW(2,2),QOW(2,2),R2DRDIC(NRMAX*2),RAT,RATT,RC(NRMAX*2), &
       RINT(NRMAX),RNUC,RR,SCALE,SHF(2,2,NMEMAX),SIMP, &
       SPLIT(NMEMAX),SPLIT1(NMEMAX),SPLIT2(NMEMAX,NTMAX), &
       SPLIT3(NMEMAX,NTMAX),SZ,VAL,VAR(4),VARNEW(4),VARTAB(4,20), &
       VV(NRMAX*2),VZ,W,WP(2,2,NRMAX*2),WQ(2,2,NRMAX*2)
LOGICAL CHECK,FERRO,SCALEB,SUPPRESSB, BNDSTACHK
DOUBLE PRECISION DBLE,DSIGN
INTEGER I,IC,IC1,IC2,ICST,IE,IFLAG,II,IL,ILC,ILSHELL,IM,IMIN,IN, &
        INFO,IPIV(4),ISH,ISTART,IT,ITER,IV,J,JLIM,JTOP,JV,K,KAP(2) &
        ,KAP1,KAP2,KC,L,LCP1,LLL,LOOP,LQNTAB(NLSHELLMAX),MUEM05,N, &
        NLSHELL,NMATCH,NN,NODE,NQN,NQNTAB(NLSHELLMAX), &
        NRC,NSH,NSOL,NVAR,NZERO,S,T
INTEGER IABS
INTEGER IKAPMUE
DOUBLE PRECISION RNUCTAB
CHARACTER (len=10) :: TXTB(1:5)
CHARACTER (len=3) :: TXTK(4)
CHARACTER (len=1) :: TXTL(0:3)

DATA TXTB/'B_nses','B_nseo','B_ssc ','B_sum ','B_tot '/
DATA TXTL/'s','p','d','f'/
DATA TXTK/'1/2','3/2','5/2','7/2'/
DATA NQNTAB/1,2,2,3,3,3,4,4,4,5,5,4,5,6,6/
DATA LQNTAB/0,0,1,0,1,2,0,1,2,0,1,3,2,0,1/
DATA CHECK/.FALSE./

nrc = 2*nrmax

CALL rinit(120*ntmax,ecortab)
CALL rinit(4*20,vartab)

IF( itxray < 0 ) THEN
  itxray = ABS(itxray)
  bndstachk = .true.
ELSE
  bndstachk = .false.
END IF

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
DO it = 1,nt
  suppressb = .false.
  scaleb = .false.
  scale = 1D0
  
  im = imt(it)
  jtop = jws(im)
  
  rat = r(nrmax,im)/r(nrmax-1,im)
  DO n = 1,nrmax
    rc(n) = r(n,im)
    drdic(n) = drdi(n,im)
    r2drdic(n) = r2drdi(n,im)
    dovrc(n) = drdic(n)/rc(n)
  END DO
  DO n = (nrmax+1),nrc
    rc(n) = rat*rc(n-1)
    drdic(n) = (rat-1.0D0)*rc(n-1)
    r2drdic(n) = rc(n)*rc(n)*drdic(n)
    dovrc(n) = drdic(n)/rc(n)
  END DO
  IF ( nucleus /= 0 ) THEN
    rnuc = rnuctab(z(it))
    in = 1
    DO WHILE ( rc(in) <= rnuc )
      in = in + 1
    END DO
!     INTEGRATION BOUNDARY FOR HYPERFINE FIELDS FOR FINITE NUCLEUS
!     2 MESH POINTS MORE FOR EXECUTING APPROPRIATE INTERPOLATION
!     TO REAL NUCLEAR RADIUS RNUC
    jlim = in + 2
    IF ( MOD(jlim,2) == 0 ) jlim = jlim - 1
  END IF
  DO i = 1,nrc
    IF ( nucleus /= 0 ) drovrn(i) = (rc(i)/rnuc)**3*drdic(i)
  END DO
  
  loop = 1
  50      CONTINUE
  bcor(it) = 0.0D0
  bcors(it) = 0.0D0
  DO i = 1,nmemax
    split2(i,it) = 0.0D0
    split3(i,it) = 0.0D0
  END DO
  bsum = 0.0D0
  DO n = 1,nrmax
    rhochr(n,it) = 0.0D0
    rhospn(n,it) = 0.0D0
  END DO
  DO n = 1,jws(im)
    vv(n) = vt(n,it)
    bb(n) = bt(n,it)
    bsum = bsum + ABS(bb(n))
  END DO
  
  IF ( suppressb ) THEN
    DO n = 1,jws(im)
      bb(n) = 0.0D0
    END DO
    bsum = 0.0D0
  END IF
  
  DO n = (jws(im)+1),nrc
    vv(n) = 0.0D0
    bb(n) = 0.0D0
  END DO
  
  nlshell = 0
  IF ( z(it) > 2 ) nlshell = 1
  IF ( z(it) > 10 ) nlshell = 3
  IF ( z(it) > 18 ) nlshell = 5
  IF ( z(it) > 30 ) nlshell = 6
  IF ( z(it) > 36 ) nlshell = 8
  IF ( z(it) > 48 ) nlshell = 9
  IF ( z(it) > 54 ) nlshell = 11
  IF ( z(it) > 70 ) nlshell = 12
  IF ( z(it) > 80 ) nlshell = 13
  IF ( z(it) > 86 ) nlshell = 15
  
  IF( ncort(it) /= 0 ) THEN
    nlshell = 0
    n = 0
    DO ilshell = 1,nlshellmax
      l = lqntab(ilshell)
      n = n + 2*(2*l+1)
      IF( n == ncort(it) ) nlshell = ilshell
    END DO
    IF( nlshell == 0 ) THEN
      WRITE(*,*) 'NLSHELL not found for IT=',it,' NCORT=', ncort(it)
      STOP ' in <CORE>'
    END IF
  END IF
  
  IF ( bsum > 1.0D-8 ) THEN
    ferro = .true.
  ELSE
    ferro = .false.
    IF ( iprint >= 1 .AND. (t_inc%i_write>0)) WRITE (1337,99001)
  END IF
  
  IF ( itxray == 0 .AND. iprint > 0 .AND. (t_inc%i_write>0))  &
      WRITE (1337,99002) itprt,z(it)
  
  
! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
!                   ---------------------------------------
!                   INITIALIZE QUANTUM NUMBERS  NQN  AND  L
!                   ---------------------------------------
  ic = 0
  DO ilshell = 1,nlshell
    nqn = nqntab(ilshell)
    l = lqntab(ilshell)
    il = l + 1
    ilc = MIN(nlmax,il)
    nsh = 2*(2*l+1)
    
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!     SKIP SHELL IF NOT NEEDED IN A  XRAY - CALCULATION
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    IF ( itxray /= 0 ) THEN
      IF ( it /= itxray ) GO TO 100
      IF ( (nqn /= ncxray(it)) .OR. (l /= lcxray(it)) ) GO TO 100
      DO icst = 1,ncstmax
        DO kc = 1,2
          DO n = 1,nrmax
            gcor(n,kc,icst) = 0.0D0
            fcor(n,kc,icst) = 0.0D0
          END DO
        END DO
      END DO
    END IF
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    
    ish = 0
    bsh = 0.0D0
    DO i = 1,nmemax
      split1(i) = 0.0D0
    END DO
! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    DO muem05 = -l - 1, + l
      mj = muem05 + 0.5D0
      
      
      kap1 = -l - 1
      kap2 = l
      kap(1) = kap1
      kap(2) = kap2
      
      lll = l*(l+1)
      IF ( ABS(mj) > l ) THEN
        nsol = 1
      ELSE
        nsol = 2
      END IF
      
      IF ( ferro ) THEN
        nvar = 2*nsol
      ELSE
        nvar = 2
      END IF
      
! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      DO s = 1,nsol
        ic = ic + 1
        ish = ish + 1
        t = 3 - s
!                  ----------------------------------------
!                   USE EC OF PREVIOUS RUNS AS START-VALUE
!                   TAKE SPIN-ORBIT SPLITTING INTO ACCOUNT
!                  ----------------------------------------
        IF ( ish > 1 ) THEN
          ec = ecortab(ic-1,it)
          IF ( s == 2 ) ec = ecortab(ic-1,it)*1.1D0
          IF ( ish >= 4 ) ec = ecortab(ic-2,it)
          GO TO 65
        END IF
        
        
!                                      --------------------
!                                         FIND  E-LIMIT
!                                      --------------------
        IF ( lll == 0 ) THEN
          elim = -2*DBLE(z(it)**2)/(1.5D0*nqn*nqn)
        ELSE
          elim = vv(1) + lll/rc(1)**2
          DO n = 2,nrc
            val = vv(n) + lll/rc(n)**2
            IF ( val <= elim ) elim = val
          END DO
        END IF
        
        ec = -DBLE(z(it)**2)/(2.0D0*nqn*nqn)
        
        istart = 1
        55               CONTINUE
        IF ( ec <= elim ) ec = elim*0.7D0
        
!                                      --------------------
!                                         FIND    NZERO
!                                      --------------------
        DO n = 1,(nrc-1)
          IF ( (vv(n)-ec)*rc(n)**2 > unend ) THEN
            IF ( MOD(n,2) == 0 ) THEN
              nzero = n + 1
            ELSE
              nzero = n
            END IF
            GO TO 60
          END IF
        END DO
        nzero = nrc - 1
        WRITE (6,99003) itprt,nqn,l,(nrc-1)
        STOP
!                                      --------------------
!                                         FIND    NMATCH
!                                      --------------------
        60               CONTINUE
        n = nzero + 1
        DO nn = 1,nzero
          n = n - 1
          IF ( (vv(n)+lll/rc(n)**2-ec) < 0.0 ) THEN
            nmatch = n
            GO TO 65
          END IF
        END DO
        IF(t_inc%i_write>0) WRITE (1337,99004) itprt,nqn,l,ec
        
        65               CONTINUE
        CALL coredir(it,ctl(it,ilc),ec,l,mj,'OUT',vv,bb,rc,  &
            drdic,dovrc,nmatch,nzero,gc,fc,dp,dq,wp,  &
            wq,pow,qow,piw,qiw,cgd,cgmd,cgo,nrc, z(it),nucleus)
        
        node = 0
        DO n = 2,nmatch
          IF ( gc(s,s,n)*gc(s,s,n-1) < 0.0 ) node = node + 1
        END DO
        
        
        IF ( iprint >= 2 .AND. (t_inc%i_write>0))  &
            WRITE (1337,99016) itprt,nqn,l,  &
            kap(s),(2*muem05+1),ic,ish,0,ec,nmatch,rc(nmatch),nzero,  &
            rc(nzero),node,(gc(s,s,nmatch)/gc(s,s,nmatch-1))
        
        IF ( node /= (nqn-l-1) ) THEN
          IF ( node > (nqn-l-1) ) THEN
            ec = 1.2D0*ec
          ELSE
            ec = 0.8D0*ec
          END IF
          GO TO 55
        ELSE IF ( (gc(s,s,nmatch)/gc(s,s,nmatch-1) <= 0.0)  &
              .OR. (gc(s,s,nmatch)/gc(s,s,nmatch-1) >= 1.0) )  &
              THEN
          ec = 0.9D0*ec
          GO TO 55
        END IF
        
        
        CALL coredir(it,ctl(it,ilc),ec,l,mj,'INW',vv,bb,rc,  &
            drdic,dovrc,nmatch,nzero,gc,fc,dp,dq,wp,  &
            wq,pow,qow,piw,qiw,cgd,cgmd,cgo,nrc, z(it),nucleus)
        
!                                      --------------------
!                                       START VALUES FOR
!                                           PARAMETERS
!                                      --------------------
        
        var(1) = ec
        var(2) = pow(s,s)/piw(s,s)
        
        IF ( nsol /= 2 .OR. nvar == 2 ) THEN
          DO iv = 3,4
            ERR(iv) = 0.0D0
            errnew(iv) = 0.0D0
            var(iv) = 0.0D0
            varnew(iv) = 0.0D0
            dv(iv) = 0.0D0
          END DO
        ELSE IF ( ish >= 4 ) THEN
          DO iv = 1,4
            var(iv) = vartab(iv,ish-2)
          END DO
        ELSE
          
          DO j = 1,nsol
            now(j) = 0.0D0
          END DO
          DO n = 1,nmatch - 1
            rr = rc(n)**3
            DO j = 1,nsol
              now(j) = now(j) + gc(j,j,n)**2*rr
            END DO
          END DO
          
          DO j = 1,nsol
            niw(j) = 0.0D0
          END DO
          DO n = nmatch,nzero - 1
            rr = rc(n)**3
            DO j = 1,nsol
              niw(j) = niw(j) + gc(j,j,n)**2*rr
            END DO
          END DO
          
          ratt = pow(t,t)/piw(t,t)
          var(3) = trymix*(now(s)+niw(s)*var(2)) /(now(t)+niw(t)*ratt)
          var(4) = ratt*var(3)/var(2)
        END IF
        
        
        CALL coreerr(ERR,var,s,nsol,pow,qow,piw,qiw)
        
        DO iv = 1,nvar
          dv(iv) = var(iv)
        END DO
        
        iter = 0
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
        70               CONTINUE
        iter = iter + 1
        
        IF ( scaleb ) THEN
          scale = MIN(1.0D0,0.1D0*iter)
          IF ( nsol == 1 ) scale = 1.0D0
          DO n = 1,jws(im)
            bb(n) = bt(n,it)*scale
          END DO
        END IF
!                         ----------------------------------
!                         CHECK WHETHER NUMBER OF NODES O.K.
!                         ----------------------------------
        IF ( iter > 1 ) THEN
          node = 0
          DO n = 2,(nmatch-1)
            IF ( gc(s,s,n)*gc(s,s,n-1) < 0.0 ) node = node + 1
          END DO
          IF ( iprint >= 2 .AND. (t_inc%i_write>0))  &
              WRITE (1337,99016) itprt,nqn,l,  &
              kap(s),(2*muem05+1),ic,ish,iter,ec,nmatch,  &
              rc(nmatch),nzero,rc(nzero),node, (gc(s,s,nmatch)/gc(s,s,nmatch-1))
          
          IF ( node /= (nqn-l-1) ) THEN
            IF ( node > (nqn-l-1) ) THEN
              ec = 1.2D0*ec
            ELSE
              ec = 0.8D0*ec
            END IF
            istart = istart + 1
            IF ( istart < 20 ) GO TO 55
          END IF
        END IF
        
        DO iv = 2,nvar
          DO jv = 1,nvar
            varnew(jv) = var(jv)
          END DO
          varnew(iv) = var(iv) + dv(iv)*dvstep
          
          IF ( ABS(var(iv)) > 1D-16 ) THEN
            IF ( ABS(dv(iv)/var(iv)) < tolvar ) varnew(iv) = var(iv)  &
                *(1.0D0+DSIGN(dvstep*tolvar,dv(iv)))
          ELSE IF ( ferro ) THEN
            IF(t_inc%i_write>0) WRITE(1337,99011) ' VAR(',iv,  &
                ') = 0 for (T,N,L,K,M;S,NSOL) ',  &
                itprt,nqn,l,kap(s),(2*muem05+1), '/2  ',s,nsol,'  --- suppress B'
            loop = 2
            suppressb = .true.
            GO TO 50
          ELSE IF ( suppressb ) THEN
            WRITE (6,*) 'suppressing B did not help !!'
            STOP 'in <CORE>'
          END IF
          
          CALL coreerr(errnew,varnew,s,nsol,pow,qow,piw,qiw)
          
          DO ie = 1,nvar
            IF ( ABS(errnew(ie)-ERR(ie)) < 1D-16 ) THEN
              dedv(ie,iv) = 0.0D0
              IF ( (ie == iv) .AND. .NOT.ferro ) dedv(ie,iv) = 1.0D0
            ELSE
              dedv(ie,iv) = (errnew(ie)-ERR(ie)) /(varnew(iv)-var(iv))
            END IF
          END DO
        END DO
        
        DO jv = 1,nvar
          varnew(jv) = var(jv)
        END DO
        varnew(1) = var(1) + dv(1)*dvstep
        IF ( ABS(dv(1)/var(1)) < tolvar ) varnew(1) = var(1)  &
            *(1.0D0+DSIGN(dvstep*tolvar,dv(1)))
        CALL coredir(it,ctl(it,ilc),varnew(1),l,mj,'OUT',vv,  &
            bb,rc,drdic,dovrc,nmatch,nzero,gc,fc,dp,  &
            dq,wp,wq,pow,qow,piw,qiw,cgd,cgmd,cgo, nrc,z(it),nucleus)
        CALL coredir(it,ctl(it,ilc),varnew(1),l,mj,'INW',vv,  &
            bb,rc,drdic,dovrc,nmatch,nzero,gc,fc,dp,  &
            dq,wp,wq,pow,qow,piw,qiw,cgd,cgmd,cgo, nrc,z(it),nucleus)
        
        CALL coreerr(errnew,varnew,s,nsol,pow,qow,piw,qiw)
        
        DO ie = 1,nvar
          dedv(ie,1) = (errnew(ie)-ERR(ie)) /(varnew(1)-var(1))
        END DO
        
        DO ie = 1,nvar
          CALL dcopy(nvar,dedv(1,ie),1,dvde(1,ie),1)
        END DO
        CALL dgetrf(nvar,nvar,dvde,4,ipiv,info)
        CALL dgetri(nvar,dvde,4,ipiv,dedv,4*4,info)
        
        DO iv = 1,nvar
          dv(iv) = 0.0D0
          DO ie = 1,nvar
            dv(iv) = dv(iv) + dvde(iv,ie)*ERR(ie)
          END DO
          var(iv) = var(iv) - dv(iv)
        END DO
        
        IF ( var(1) > 0.0D0 ) THEN
          IF ( iprint >= 1 .AND. (t_inc%i_write>0)) WRITE (1337,*)  &
              ' warning from <CORE> E=',var(1),it,nqn,l
          var(1) = -0.2D0
        END IF
        
        CALL coredir(it,ctl(it,ilc),var(1),l,mj,'OUT',vv,bb,  &
            rc,drdic,dovrc,nmatch,nzero,gc,fc,dp,dq,  &
            wp,wq,pow,qow,piw,qiw,cgd,cgmd,cgo,nrc, z(it),nucleus)
        CALL coredir(it,ctl(it,ilc),var(1),l,mj,'INW',vv,bb,  &
            rc,drdic,dovrc,nmatch,nzero,gc,fc,dp,dq,  &
            wp,wq,pow,qow,piw,qiw,cgd,cgmd,cgo,nrc, z(it),nucleus)
        
        CALL coreerr(ERR,var,s,nsol,pow,qow,piw,qiw)
        
        ec = var(1)
        
        IF ( iprint >= 2 .AND. (t_inc%i_write>0))  &
            WRITE (1337,99005) loop,scale,  &
            var(1),(var(iv),iv=1,4),(dv(iv),iv=1,4),(ERR(ie),ie=1,4)
        
!----------------------------------  check relative change in parameters
! ----------------------- parameters 3 and 4 = 0 for paramagnetic case !
        IF ( iter < itermax ) THEN
          DO iv = 1,nvar
            vartab(iv,ish) = var(iv)
!           IF( ABS(VAR(IV)) .EQ. 0.0D0 ) THEN
            IF ( (ABS(var(iv))+ABS(var(iv))) < 1.0D-30 ) THEN
              IF ( ferro .AND. (t_inc%i_write>0))  &
                  WRITE (1337,'(A,I3,A)') ' VAR ',iv,' = 0 ??????!!!!!'
            ELSE IF ( ABS(dv(iv)/var(iv)) > tolvar ) THEN
              GO TO 70
            END IF
          END DO
        ELSE
          IF( bndstachk ) THEN
            itxray = -1
            RETURN
          END IF
          IF(t_inc%i_write>0) WRITE (1337,99006) itermax,(var(iv),iv=1,4),  &
              (dv(iv),iv=1,4),(ERR(ie),ie=1,4)
        END IF
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
        
!                         ---------------------------------
!                         NORMALIZE WAVEFUNCTIONS ACCORDING
!                               TO MATCHING CONDITIONS
!                         ---------------------------------
        
!                                    INWARD - SOLUTION
        DO n = nmatch,nzero
          DO j = 1,nsol
            DO i = 1,nsol
              gc(i,j,n) = gc(i,j,n)*var(2)
              fc(i,j,n) = fc(i,j,n)*var(2)
            END DO
          END DO
        END DO
        
        IF ( nsol == 2 ) THEN
!                                   OUTWARD - SOLUTION
          DO n = 1,(nmatch-1)
            DO i = 1,nsol
              gc(i,t,n) = gc(i,t,n)*var(3)
              fc(i,t,n) = fc(i,t,n)*var(3)
            END DO
          END DO
!                                    INWARD - SOLUTION
          DO n = nmatch,nzero
            DO i = 1,nsol
              gc(i,t,n) = gc(i,t,n)*var(4)
              fc(i,t,n) = fc(i,t,n)*var(4)
            END DO
          END DO
        END IF
        
!                                    SUM FOR EACH KAPPA
        DO n = 1,nzero
          DO k = 1,nsol
            gck(k,s,n) = 0.0D0
            fck(k,s,n) = 0.0D0
            DO j = 1,nsol
              gck(k,s,n) = gck(k,s,n) + gc(k,j,n)
              fck(k,s,n) = fck(k,s,n) + fc(k,j,n)
            END DO
          END DO
        END DO
        
!                       -----------------------------------
!                       CALCULATE  NORM  AND NORMALIZE TO 1
!                       -----------------------------------
        DO k = 1,nsol
          norm = r2drdic(1)*(gck(k,s,1)**2+fck(k,s,1)**2)
        END DO
        
        simp = -1.0D0
        DO n = 2,nzero
          simp = -simp
          w = (3.0D0+simp)*r2drdic(n)
          DO k = 1,nsol
            norm = norm + w*(gck(k,s,n)**2+fck(k,s,n)**2)
          END DO
        END DO
        
        n = nzero
        DO k = 1,nsol
          norm = norm - r2drdic(n) *(gck(k,s,n)**2+fck(k,s,n)**2)
        END DO
        norm = norm/3.0D0
        
        DO k = 1,nsol
          norm = norm + 0.5D0*r2drdic(1) *(gck(k,s,1)**2+fck(k,s,1)**2)
        END DO
        norm = 1.0D0/SQRT(norm)
        
        DO n = 1,nzero
          DO k = 1,nsol
            gck(k,s,n) = gck(k,s,n)*norm
            fck(k,s,n) = fck(k,s,n)*norm
          END DO
        END DO
        IF ( nzero < jtop ) THEN
          DO n = (nzero+1),jtop
            DO k = 1,nsol
              gck(k,s,n) = 0.0D0
              fck(k,s,n) = 0.0D0
            END DO
          END DO
        END IF
        
        CALL rinit(nrmax,rint)
        
        DO n = 1,jtop
          DO k = 1,nsol
            rint(n) = rint(n) + r2drdi(n,im)  &
                *(gck(k,s,n)*gck(k,s,n)+fck(k,s,n) *fck(k,s,n))
          END DO
        END DO
        
        CALL rintsimp(rint,jtop,aux)
        
! ------------------------------ omit normalization for XRAY calculation
! ------------------------------ to recover old (errounous data) -------
        IF ( itxray > 0 ) THEN
          norm = 1.0D0
        ELSE
          norm = 1.0D0/SQRT(aux)
        END IF
        
        DO n = 1,MAX(nzero,jtop)
          DO k = 1,nsol
            gck(k,s,n) = gck(k,s,n)*norm
            fck(k,s,n) = fck(k,s,n)*norm
          END DO
        END DO
        
!                       -----------------------------------
!                       CALCULATE  CHARGE AND SPIN DENSITY
!                       -----------------------------------
        
        DO n = 1,jws(im)
          DO k = 1,nsol
            rhochr(n,it) = rhochr(n,it) + (gck(k,s,n)**2+fck(k,s,n)**2)
            rhospn(n,it) = rhospn(n,it) + (gck(k,s,n)**2*cgd(k)  &
                -fck(k,s,n)**2*cgmd(k))
          END DO
        END DO
        
        IF ( nsol > 1 ) THEN
          DO n = 1,jws(im)
            rhospn(n,it) = rhospn(n,it) + gck(1,s,n) *gck(2,s,n)*cgo*2
          END DO
        END IF
        
        
!                       -----------------------------------
!                            CALCULATE  SPIN CHARACTER
!                       -----------------------------------
        
        w = r2drdic(1)
        sz = 0.0D0
        DO k = 1,nsol
          sz = sz + w*(gck(k,s,1)**2*cgd(k)+fck(k,s,1) **2*cgmd(k))
        END DO
        
        simp = -1.0D0
        DO n = 2,nzero
          simp = -simp
          w = (3.0D0+simp)*r2drdic(n)
          DO k = 1,nsol
            sz = sz + w*(gck(k,s,n)**2*cgd(k)+fck(k,s,n) **2*cgmd(k))
          END DO
        END DO
        
        n = nzero
        w = r2drdic(n)
        DO k = 1,nsol
          sz = sz + w*(gck(k,s,n)**2*cgd(k)+fck(k,s,n) **2*cgmd(k))
        END DO
        
        
        IF ( nsol > 1 ) THEN
          
          w = r2drdic(1)
          sz = sz + w*gck(1,s,1)*gck(2,s,1)*cgo*2
          
          simp = -1.0D0
          DO n = 2,nzero
            simp = -simp
            w = (3.0D0+simp)*r2drdic(n)
            DO k = 1,nsol
              sz = sz + w*gck(1,s,n)*gck(2,s,n)*cgo*2
            END DO
          END DO
          
          n = nzero
          w = r2drdic(n)
          sz = sz - w*gck(1,s,n)*gck(2,s,n)*cgo*2
          
        END IF
        
        sz = sz/3.0D0
        
        
!                         ------------------------------
!                         CALCULATE   HYPERFINE - FIELD
!                         ------------------------------
        
        CALL corehff(kap1,kap2,mj,s,nsol,bhf,gck,fck,rc,drdic,  &
            0.0D0,nzero,nrc)
        IF ( nucleus /= 0 ) THEN
          CALL corehff(kap1,kap2,mj,s,nsol,bhf1,gck,fck,rc,  &
              drdic,rnuc,jlim,nrc)
          CALL corehff(kap1,kap2,mj,s,nsol,bhf2,gck,fck,rc,  &
              drovrn,rnuc,jlim,nrc)
          DO i = 1,nsol
            DO j = 1,nsol
              bhf(i,j) = bhf(i,j) - bhf1(i,j) + bhf2(i,j)
            END DO
          END DO
        END IF        !end of nucleus.eq.0
        
        bsol = 0.0D0
        DO j = 1,nsol
          DO i = 1,nsol
            bsol = bsol + bhf(i,j)
            bsh = bsh + bhf(i,j)
            bcor(it) = bcor(it) + bhf(i,j)
          END DO
        END DO
        IF ( kap1 == -1 ) bcors(it) = bcors(it) + bhf(1,1)
        
        ecortab(ic,it) = ec
        
!     ------------------
!     SPLIT HFF-FIELD
!     ------------------
        IF ( ismqhfi == 1 ) THEN
          CALL hffcore(rnuc,nzero,kap1,kap2,nsol,mj,gck,fck,  &
              nrc,shf,s,nmemax,nkmmax,rc,drdic,sdia,  &
              smdia,soff,smoff,qdia,qoff,qmdia, qmoff,nucleus,jlim)
          
          DO k = 1,nmemax
            split(k) = 0.0D0
            DO j = 1,nsol
              DO i = 1,nsol
                split(k) = split(k) + shf(i,j,k)
                split1(k) = split1(k) + shf(i,j,k)
                split2(k,it) = split2(k,it) + shf(i,j,k)
              END DO
            END DO
          END DO
          DO k = 1,nmemax
            IF ( kap1 == -1 ) split3(k,it) = split3(k,it) + shf(1,1,k)
          END DO
        END IF
!MBE
        
!---------------------------------------------------- l-shell UNCOMPLETE
        IF ( ish >= nsh ) THEN
!----------------------------------------------------- l-shell completed
          IF ( itxray == 0 .AND. iprint > 0 .AND. (t_inc%i_write>0)) THEN
            WRITE (1337,99012) itprt,nqn,txtl(l),  &
                txtk(IABS(kap(s))),(2*muem05+1), kap(s),iter,ec,bsol*.001D0,  &
                bsh*.001D0
            IF ( ismqhfi == 1 ) THEN
              DO k = 1,nmemax
                WRITE (1337,99013) txtb(k),split(k)*.001D0 ,split1(k)*.001D0
              END DO
              WRITE (1337,99014) 'total error in %',  &
                  100.0D0*(1.0D0-split(4)/split(5))
            END IF
          END IF
!                              ----------------------------
!                              CHECK CONSISTENCY OF RESULTS
!                              ----------------------------
          IF ( l /= 0 ) THEN
            ic1 = ic - nsh + 1
            ic2 = ic
            IF ( ecortab(ic2,it) >= ecortab(ic1,it) ) THEN
              imin = ic1
              vz = +1.0D0
            ELSE
              imin = ic2
              vz = -1.0D0
            END IF
            iflag = 0
            ii = 1
            DO i = ic1 + 1,ic2,2
              IF ( vz*(ecortab(i,it)-ecortab(i-ii,it))  < 0.0 ) iflag = 1
              ii = 2
            END DO
            IF ( ecortab(ic1+2,it) > ecortab(imin,it) ) iflag = 1
            DO i = ic1 + 4,ic2 - 1,2
              IF ( ecortab(i,it) > ecortab(imin,it) ) iflag = 1
              IF ( vz*(ecortab(i,it)-ecortab(i-ii,it))  > 0.0 ) iflag = 1
            END DO
            
            IF ( ferro .AND. (iflag == 1) ) THEN
              IF(t_inc%i_write>0) WRITE (1337,99007)
              scaleb = .true.
              IF ( loop == 1 ) THEN
                loop = 2
                IF(t_inc%i_write>0) WRITE(1337,99008) itprt
                GO TO 50
              END IF
            END IF
            
          END IF
        ELSE IF ( itxray == 0 .AND. iprint > 0 ) THEN
          IF(t_inc%i_write>0) WRITE (1337,99012) itprt,nqn,  &
              txtl(l),txtk(IABS(kap(s))), (2*muem05+1),kap(s),iter,ec,  &
              bsol*.001D0
          IF ( ismqhfi == 1 ) THEN
            DO k = 1,nmemax
              IF(t_inc%i_write>0) WRITE (1337,99013) txtb(k),split(k)*.001D0
            END DO
            IF(t_inc%i_write>0) WRITE (1337,99014) 'total error in %',  &
                100.0D0*(1.0D0-split(4)/split(5) )
          END IF
        END IF
!-----------------------------------------------------------------------
        
        IF ( iprint >= 1 .AND. (t_inc%i_write>0)) WRITE (1337,99015)  &
            ((bhf(i,j)*.001D0,i=1,nsol),j=1,nsol)
        
        
!                          --------------------------------
!                            IF THE SWITCH CHECK IS SET:
!                            RECALCULATE THE EIGENVALUE
!                          USING THE CONVENTIONAL ALGORITHM
!                          --------------------------------
        IF ( check ) THEN
          ecc = 0.95D0*ec
          72                  CONTINUE
          CALL coredir(it,ctl(it,ilc),ecc,l,mj,'OUT',vv,bb,  &
              rc,drdic,dovrc,nmatch,nzero,gc,fc,dp,  &
              dq,wp,wq,pow,qow,piw,qiw,cgd,cgmd,cgo, nrc,z(it),nucleus)
          CALL coredir(it,ctl(it,ilc),ecc,l,mj,'INW',vv,bb,  &
              rc,drdic,dovrc,nmatch,nzero,gc,fc,dp,  &
              dq,wp,wq,pow,qow,piw,qiw,cgd,cgmd,cgo, nrc,z(it),nucleus)
          
          norm = pow(s,s)/piw(s,s)
          DO n = nmatch,nzero
            gc(s,s,n) = gc(s,s,n)*norm
            fc(s,s,n) = fc(s,s,n)*norm
          END DO
          
          norm = 0.0D0
          DO n = 3,nzero,2
            norm = norm + r2drdic(n) *(gc(s,s,n)**2+fc(s,s,n)**2)  &
                + 4.d0*r2drdic(n-1) *(gc(s,s,n-1)**2+fc(s,s,n-1)**2)  &
                + r2drdic(n-2) *(gc(s,s,n-2)**2+fc(s,s,n-2)**2)
          END DO
          norm = norm/3.0D0
          
          lcp1 = MIN(nlmax,l+1)
          dec = pow(s,s) *(qow(s,s)-rc(nmatch)*ctl(it,lcp1)*fc(s,s,  &
              nmatch))/norm
          ecc = ecc + dec
          IF ( ABS(dec/ecc) > tolvar ) GO TO 72
          IF(t_inc%i_write>0)  &
              WRITE (1337,'(7X,''CHECK-E:'',10X,F12.5,/)') ecc
        END IF
        
        
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!           STORE CORE WAVE FUNCTIONS IF REQUIRED
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        IF ( itxray /= 0 ) THEN
          IF ( nsol == 2 ) THEN
            IF ( s == 2 ) THEN
              icst = l + 1 + muem05
            ELSE
              icst = l + 1 + muem05 + 2*l + 1
            END IF
          ELSE IF ( mj < 0 ) THEN
            icst = 2*l + 1
          ELSE
            icst = 4*l + 2
          END IF
          
          mm05cor(icst) = muem05
          nkpcor(icst) = nsol
          kapcor(icst) = kap(s)
          ikmcor(icst,1) = ikapmue(kap(s),muem05)
          izero(icst) = MIN(nzero,jws(im))
          szcor(icst) = sz
          ecor(icst) = ec
          
          DO n = 1,izero(icst)
            gcor(n,1,icst) = gck(s,s,n)
            fcor(n,1,icst) = fck(s,s,n)
          END DO
          IF ( nsol == 2 ) THEN
            DO n = 1,izero(icst)
              gcor(n,2,icst) = gck(t,s,n)
              fcor(n,2,icst) = fck(t,s,n)
            END DO
            ikmcor(icst,2) = ikapmue(kap(t),muem05)
          END IF
        END IF
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        
        
      END DO
! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      
    END DO
! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
  100     END DO
! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
  
  
  IF ( itxray == 0 .AND. iprint > 0 .AND. (t_inc%i_write>0))  &
      WRITE (1337,99009) bcor(it)*.001D0,bcors(it)*.001D0
  IF ( ismqhfi == 1 ) THEN
    DO n = 1,nmemax
      IF(t_inc%i_write>0) WRITE (1337,99010) split2(n,it)*.001D0,  &
          split3(n,it)*.001D0
    END DO
    IF(t_inc%i_write>0) WRITE (1337,99014) 'total error',  &
        100.0D0*(1.0D0-split2(4,it)/split2(5,it))
    IF(t_inc%i_write>0) WRITE (1337,99014) 'total error',  &
        100.0D0*(1.0D0-split3(4,it)/split3(5,it))
  END IF
  
  IF ( (itxray == 0) .OR. (it == itxray) ) THEN
    DO n = 1,jtop
      rint(n) = rhochr(n,it)*r2drdi(n,im)
    END DO
    CALL rintsimp(rint,jtop,aux)
    IF( iprint > -2 .AND. (t_inc%i_write>0))  &
        WRITE (1337,99017) 'charge',itprt,aux
    DO n = 1,jtop
      rint(n) = rhospn(n,it)*r2drdi(n,im)
    END DO
    CALL rintsimp(rint,jtop,aux)
    IF( iprint > -2 .AND. (t_inc%i_write>0))  &
        WRITE (1337,99017) ' spin ',itprt,aux
  END IF
  
END DO
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

99001 FORMAT (/,10X,'potential is not exchange split ')
99002 FORMAT (/,'  ATOM   IT      : ',i5,/,'  ATOMIC NUMBER  : ',i5,//,  &
    '  IT',10X,'MUE  KAP ITER    ENERGY       B  (k-Gauss)  ')
99003 FORMAT ('  IT=',i2,' NQN=',i2,' L=',i2, '  NZERO set to  (NRC-1) =',i4)
99004 FORMAT (//,'  STOP IN <<CORE>>',/,'  IT=',i2,' NQN=',i2,' L=',i2,  &
    /,'  no matching-radius found for  EC=',f10.3)
99005 FORMAT (' LOOP    =  ',i3,' BSCL=',f10.5,/,' E=',f25.16,' VAR  ',  &
    4E11.4,/,17X,' CORR ',4E11.4,/,17X,' ERR  ',4E11.4)
99006 FORMAT (' iteration not converged after',i3,' steps !',/,  &
    ' parameters:',4E18.10,/,' last corr.:',4E18.10,/, ' last error:',4E18.10)
99007 FORMAT (' >>> check data E(KAP,MJ) should be monotonous ',  &
    ' and  E(+L,MJ) < E(-L-1,MJ) ',//)
99008 FORMAT (' all core states for atom type ',i2,  &
    ' will be recalculated ',/,  &
    ' with the spin dependent exchange field gradually switched on' )
99009 FORMAT (2X,57('-'),/,42X,f17.3,/,38X,'(S) ',f17.3,/,2X,57('*'),//)
99010 FORMAT (2X,57('-'),/,37X,f17.3,/,33X,'(S) ',f17.3,/,2X,57('*'),//)
99011 FORMAT (a,i1,a,5I3,a,2I2,a)
99012 FORMAT (2I4,a1,a3,i3,'/2',2I4,2X,f15.8,f17.3,:,f17.3,/)
99013 FORMAT (a,:,32X,f17.3,:,f17.3,/)
99014 FORMAT (a,:,37X,f6.3)
99015 FORMAT (37X,f17.3)
99016 FORMAT (/,' IT=',i2,'  NQN=',i2,'  L=',i2,'  KAP=',i2,'  MJ=',i2,  &
    '/2    IC=',i3,'  ISH=',i2,/,' E(',i2,')   =',f15.5,/,  &
    ' NMATCH  =',i5,'    R=',f10.5,/,' NZERO   =',i5,'    R=',  &
    f10.5,/,' NODES   =',i5,'  RAT=',e11.4)
99017 FORMAT (' integrated core ',a,' density for atom type ',i4,':', f12.8)
END SUBROUTINE core
