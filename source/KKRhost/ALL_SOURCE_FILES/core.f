C*==core.f    processed by SPAG 6.05Rc at 12:36 on 29 Apr 2001
      SUBROUTINE CORE(IPRINT,ITPRT,NT,NCORT,CTL,VT,BT,Z,NUCLEUS,
     &                R,R2DRDI,DRDI,JWS,IMT,RHOCHR,RHOSPN,
     &                ECORTAB,GCOR,FCOR,ECOR,SZCOR,
     &                KAPCOR,MM05COR,NKPCOR,IKMCOR,IZERO,NCXRAY,LCXRAY,
     &                ITXRAY,BCOR,BCORS,SDIA,SMDIA,SOFF,SMOFF,QDIA,QOFF,
     &                QMDIA,QMOFF,NKMMAX,NMEMAX,ISMQHFI,NTMAX,NRMAX,
     &                NMMAX,NCSTMAX,NLMAX)
C   ********************************************************************
C   *                                                                  *
C   *   SUBROUTINE TO CALCULATE THE RELATIVISTIC CORE WAVE             *
C   *   FUNCTIONS FOR A SPIN-DEPENDENT POTENTIAL                       *
C   *                                                                  *
C   *   FOR A GIVEN POTENTIAL THE NUMBER OF CORE AND VALENCE           *
C   *   ELECTRONS IS DETERMINED AND ALL CORE STATES THEN CALCULATED    *
C   *   > THE ROUTINE IS ORGANIZED AS DESCLAUX'S ROUTINE <RESLD>       *
C   *     BUT FINDS THE CORRECTION TO THE E-EIGENVALUE AND THE         *
C   *     MATCHING PARAMETERS BY A NEWTON RAPHSON ALGORITHM            *
C   *     THIS IS IN VARIANCE TO THE METHOD SUGGESTED BY CORTONA       *
C   *   > SET THE SWITCH 'CHECK'  TO COPARE E-EIGENVALUES WITH         *
C   *     RESULTS OBTAINED WITH THE CONVENTIONAL E-CORRECTION          *
C   *     ALGORITHM, WHICH WORKS ONLY IF NO COUPLING IS PRESENT !      *
C   *   > THE FUNCTIONS  {GC,FC}(I,J) J=1,NSOL ARE THE LINEAR          *
C   *     INDEPENDENT SOLUTIONS TO THE DIFFERENTIAL EQUATIONS WITH     *
C   *     KAPPA-CHARACTER I=1,NSOL;   FOR OUTWARD AND INWARD           *
C   *     INTEGRATION THE SAME ARRAYS ARE USED !                       *
C   *   > THE PROPER SOLUTIONS SATISFYING THE BOUNDARY CONDITIONS      *
C   *     AT R=0 AND(!) R=INFINITY ARE STORED IN {GCK,FCK}(K,S)        *
C   *     S,K=1,NSOL   SOLUTION S=1 FOR  KAPPA = - L - 1               *
C   *                           S=2 FOR  KAPPA = + L (IF EXISTENT)     *
C   *   > THE SWITCH NUCLEUS SELECTS WHETHER A FINITE NUCLEUS          *
C   *     SHOULD BE USED                                               *
C   *                                                                  *
C   *   ADAPTED FOR FINITE NUCLEUS       MB MAR. 1995                  *
C   *   HYPERFINE FIELD SPLITTING introduced if icore=1 MB JUN. 1995   *
C   *                                                                  *
C   *   SCALEB:                                                        *
C   *   if the B-field is quite high it might happen that the routine  *
C   *   fails to find both 'spin-orbit-split' solutions.               *
C   *   in that case the whole l-shell is rerun with the B-field       *
C   *   gradually switched on, i.e. scaled with a parameter that       *
C   *   increases from 0 to 1 during the iteration loop  HE Nov. 95    *
C   *                                                                  *
C   *   ITXRAY =  0  run over all core states to get charge density    *
C   *   ITXRAY >  0  deal exclusively with state  NCXRAY,LCXRAY        *
C   *   ITXRAY <  0  state  NCXRAY,LCXRAY  is checked to be a          *
C   *                bound state or not. on return:                    *
C   *                ITXRAY = |ITXRAY| indicates bound state           *
C   *                ITXRAY = -1       indicates NO bound state found  *
C   *                                                                  *
C   *                                                                  *
C   *   few changes in the TB-KKR implementation as compared to SPR    *
C   *        IPRINT values between 0 and 2                             *
C   *        ITPRT  correct value of the atom-type index               *
C   ********************************************************************
      use mod_types, only: t_inc
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      REAL*8 UNEND,TOLVAR,TRYMIX,DVSTEP
      PARAMETER (UNEND=600.0D0,TOLVAR=1.0D-6,TRYMIX=0.01D0,
     &           DVSTEP=0.01D0)
      INTEGER ITERMAX,NLSHELLMAX
      PARAMETER (ITERMAX=200,NLSHELLMAX=15)
C
C Dummy arguments
C
      INTEGER IPRINT,ISMQHFI,ITPRT,ITXRAY,NCSTMAX,NKMMAX,NLMAX,NMEMAX,
     &        NMMAX,NRMAX,NT,NTMAX,NUCLEUS
      REAL*8 BCOR(NTMAX),BCORS(NTMAX),BT(NRMAX,NTMAX),CTL(NTMAX,NLMAX),
     &       DRDI(NRMAX,NMMAX),ECOR(NCSTMAX),ECORTAB(120,NTMAX),
     &       FCOR(NRMAX,2,NCSTMAX),GCOR(NRMAX,2,NCSTMAX),QDIA(NKMMAX),
     &       QMDIA(NKMMAX),QMOFF(NKMMAX),QOFF(NKMMAX),R(NRMAX,NMMAX),
     &       R2DRDI(NRMAX,NMMAX),RHOCHR(NRMAX,NTMAX),RHOSPN(NRMAX,NTMAX)
     &       ,SDIA(NKMMAX),SMDIA(NKMMAX),SMOFF(NKMMAX),SOFF(NKMMAX),
     &       SZCOR(NCSTMAX),VT(NRMAX,NTMAX)
      INTEGER IKMCOR(NCSTMAX,2),IMT(NTMAX),IZERO(NCSTMAX),JWS(NMMAX),
     &        KAPCOR(NCSTMAX),LCXRAY(NTMAX),MM05COR(NCSTMAX),
     &        NCXRAY(NTMAX),NKPCOR(NCSTMAX),Z(NTMAX),NCORT(NTMAX)
C
C Local variables
C
      REAL*8 AUX,BB(NRMAX*2),BHF(2,2),BHF1(2,2),BHF2(2,2),BSH,BSOL,BSUM,
     &       CGD(2),CGMD(2),CGO,DEC,DEDV(4,4),DOVRC(NRMAX*2),
     &       DP(2,2,NRMAX*2),DQ(2,2,NRMAX*2),DRDIC(NRMAX*2),
     &       DROVRN(2*NRMAX),DV(4),DVDE(4,4),EC,ECC,ELIM,ERR(4),
     &       ERRNEW(4),FC(2,2,NRMAX*2),FCK(2,2,NRMAX*2),GC(2,2,NRMAX*2),
     &       GCK(2,2,NRMAX*2),MJ,NIW(2),NORM,NOW(2),PIW(2,2),POW(2,2),
     &       QIW(2,2),QOW(2,2),R2DRDIC(NRMAX*2),RAT,RATT,RC(NRMAX*2),
     &       RINT(NRMAX),RNUC,RR,SCALE,SHF(2,2,NMEMAX),SIMP,
     &       SPLIT(NMEMAX),SPLIT1(NMEMAX),SPLIT2(NMEMAX,NTMAX),
     &       SPLIT3(NMEMAX,NTMAX),SZ,VAL,VAR(4),VARNEW(4),VARTAB(4,20),
     &       VV(NRMAX*2),VZ,W,WP(2,2,NRMAX*2),WQ(2,2,NRMAX*2)
      LOGICAL CHECK,FERRO,SCALEB,SUPPRESSB, BNDSTACHK
      DOUBLE PRECISION DBLE,DSIGN
      INTEGER I,IC,IC1,IC2,ICST,IE,IFLAG,II,IL,ILC,ILSHELL,IM,IMIN,IN,
     &        INFO,IPIV(4),ISH,ISTART,IT,ITER,IV,J,JLIM,JTOP,JV,K,KAP(2)
     &        ,KAP1,KAP2,KC,L,LCP1,LLL,LOOP,LQNTAB(NLSHELLMAX),MUEM05,N,
     &        NLSHELL,NMATCH,NN,NODE,NQN,NQNTAB(NLSHELLMAX),
     &        NRC,NSH,NSOL,NVAR,NZERO,S,T
      INTEGER IABS
      INTEGER IKAPMUE
      REAL*8 RNUCTAB
      CHARACTER*10 TXTB(1:5)
      CHARACTER*3 TXTK(4)
      CHARACTER*1 TXTL(0:3)
C
C*** End of declarations rewritten by SPAG
C
      DATA TXTB/'B_nses','B_nseo','B_ssc ','B_sum ','B_tot '/
      DATA TXTL/'s','p','d','f'/
      DATA TXTK/'1/2','3/2','5/2','7/2'/
      DATA NQNTAB/1,2,2,3,3,3,4,4,4,5,5,4,5,6,6/
      DATA LQNTAB/0,0,1,0,1,2,0,1,2,0,1,3,2,0,1/
      DATA CHECK/.FALSE./
C
      NRC = 2*NRMAX
C
      CALL RINIT(120*NTMAX,ECORTAB)
      CALL RINIT(4*20,VARTAB)
C
      IF( ITXRAY .LT. 0 ) THEN 
         ITXRAY = ABS(ITXRAY)
         BNDSTACHK = .TRUE.
      ELSE
         BNDSTACHK = .FALSE.
      END IF
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = 1,NT
         SUPPRESSB = .FALSE.
         SCALEB = .FALSE.
         SCALE = 1D0
C
         IM = IMT(IT)
         JTOP = JWS(IM)
C
         RAT = R(NRMAX,IM)/R(NRMAX-1,IM)
         DO N = 1,NRMAX
            RC(N) = R(N,IM)
            DRDIC(N) = DRDI(N,IM)
            R2DRDIC(N) = R2DRDI(N,IM)
            DOVRC(N) = DRDIC(N)/RC(N)
         END DO
         DO N = (NRMAX+1),NRC
            RC(N) = RAT*RC(N-1)
            DRDIC(N) = (RAT-1.0D0)*RC(N-1)
            R2DRDIC(N) = RC(N)*RC(N)*DRDIC(N)
            DOVRC(N) = DRDIC(N)/RC(N)
         END DO
         IF ( NUCLEUS.NE.0 ) THEN
            RNUC = RNUCTAB(Z(IT))
            IN = 1
            DO WHILE ( RC(IN).LE.RNUC )
               IN = IN + 1
            END DO
C     INTEGRATION BOUNDARY FOR HYPERFINE FIELDS FOR FINITE NUCLEUS
C     2 MESH POINTS MORE FOR EXECUTING APPROPRIATE INTERPOLATION
C     TO REAL NUCLEAR RADIUS RNUC
            JLIM = IN + 2
            IF ( MOD(JLIM,2).EQ.0 ) JLIM = JLIM - 1
         END IF
         DO I = 1,NRC
            IF ( NUCLEUS.NE.0 ) DROVRN(I) = (RC(I)/RNUC)**3*DRDIC(I)
         END DO
C
         LOOP = 1
 50      CONTINUE
         BCOR(IT) = 0.0D0
         BCORS(IT) = 0.0D0
         DO I = 1,NMEMAX
            SPLIT2(I,IT) = 0.0D0
            SPLIT3(I,IT) = 0.0D0
         END DO
         BSUM = 0.0D0
         DO N = 1,NRMAX
            RHOCHR(N,IT) = 0.0D0
            RHOSPN(N,IT) = 0.0D0
         END DO
         DO N = 1,JWS(IM)
            VV(N) = VT(N,IT)
            BB(N) = BT(N,IT)
            BSUM = BSUM + ABS(BB(N))
         END DO
C
         IF ( SUPPRESSB ) THEN
            DO N = 1,JWS(IM)
               BB(N) = 0.0D0
            END DO
            BSUM = 0.0D0
         END IF
C
         DO N = (JWS(IM)+1),NRC
            VV(N) = 0.0D0
            BB(N) = 0.0D0
         END DO
C
         NLSHELL = 0
         IF ( Z(IT).GT.2 ) NLSHELL = 1
         IF ( Z(IT).GT.10 ) NLSHELL = 3
         IF ( Z(IT).GT.18 ) NLSHELL = 5
         IF ( Z(IT).GT.30 ) NLSHELL = 6
         IF ( Z(IT).GT.36 ) NLSHELL = 8
         IF ( Z(IT).GT.48 ) NLSHELL = 9
         IF ( Z(IT).GT.54 ) NLSHELL = 11
         IF ( Z(IT).GT.70 ) NLSHELL = 12
         IF ( Z(IT).GT.80 ) NLSHELL = 13
         IF ( Z(IT).GT.86 ) NLSHELL = 15
C
         IF( NCORT(IT) .NE. 0 ) THEN
            NLSHELL = 0
            N = 0 
            DO ILSHELL = 1,NLSHELLMAX
               L = LQNTAB(ILSHELL)
               N = N + 2*(2*L+1)
               IF( N .EQ. NCORT(IT) ) NLSHELL = ILSHELL
            END DO
            IF( NLSHELL .EQ. 0 ) THEN 
               WRITE(*,*) 'NLSHELL not found for IT=',IT,' NCORT=',
     &                     NCORT(IT)
               STOP ' in <CORE>'
            END IF   
         END IF
C
         IF ( BSUM.GT.1.0D-8 ) THEN
            FERRO = .TRUE.
         ELSE
            FERRO = .FALSE.
            IF ( IPRINT.GE.1 .and. (t_inc%i_write>0)) WRITE (1337,99001)
         END IF
C
         IF ( ITXRAY.EQ.0 .AND. IPRINT.GT.0 .and. (t_inc%i_write>0)) 
     &        WRITE (1337,99002) ITPRT,Z(IT)
C
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C                   ---------------------------------------
C                   INITIALIZE QUANTUM NUMBERS  NQN  AND  L
C                   ---------------------------------------
         IC = 0
         DO ILSHELL = 1,NLSHELL
            NQN = NQNTAB(ILSHELL)
            L = LQNTAB(ILSHELL)
            IL = L + 1
            ILC = MIN(NLMAX,IL)
            NSH = 2*(2*L+1)
C
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C     SKIP SHELL IF NOT NEEDED IN A  XRAY - CALCULATION
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            IF ( ITXRAY.NE.0 ) THEN
               IF ( IT.NE.ITXRAY ) GOTO 100
               IF ( (NQN.NE.NCXRAY(IT)) .OR. (L.NE.LCXRAY(IT)) )
     &              GOTO 100
               DO ICST = 1,NCSTMAX
                  DO KC = 1,2
                     DO N = 1,NRMAX
                        GCOR(N,KC,ICST) = 0.0D0
                        FCOR(N,KC,ICST) = 0.0D0
                     END DO
                  END DO
               END DO
            END IF
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C
            ISH = 0
            BSH = 0.0D0
            DO I = 1,NMEMAX
               SPLIT1(I) = 0.0D0
            END DO
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
            DO MUEM05 = -L - 1, + L
               MJ = MUEM05 + 0.5D0
C
C
               KAP1 = -L - 1
               KAP2 = L
               KAP(1) = KAP1
               KAP(2) = KAP2
C
               LLL = L*(L+1)
               IF ( ABS(MJ).GT.L ) THEN
                  NSOL = 1
               ELSE
                  NSOL = 2
               END IF
C
               IF ( FERRO ) THEN
                  NVAR = 2*NSOL
               ELSE
                  NVAR = 2
               END IF
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
               DO S = 1,NSOL
                  IC = IC + 1
                  ISH = ISH + 1
                  T = 3 - S
C                  ----------------------------------------
C                   USE EC OF PREVIOUS RUNS AS START-VALUE
C                   TAKE SPIN-ORBIT SPLITTING INTO ACCOUNT
C                  ----------------------------------------
                  IF ( ISH.GT.1 ) THEN
                     EC = ECORTAB(IC-1,IT)
                     IF ( S.EQ.2 ) EC = ECORTAB(IC-1,IT)*1.1D0
                     IF ( ISH.GE.4 ) EC = ECORTAB(IC-2,IT)
                     GOTO 65
                  END IF
C
C
C                                      --------------------
C                                         FIND  E-LIMIT
C                                      --------------------
                  IF ( LLL.EQ.0 ) THEN
                     ELIM = -2*DBLE(Z(IT)**2)/(1.5D0*NQN*NQN)
                  ELSE
                     ELIM = VV(1) + LLL/RC(1)**2
                     DO N = 2,NRC
                        VAL = VV(N) + LLL/RC(N)**2
                        IF ( VAL.LE.ELIM ) ELIM = VAL
                     END DO
                  END IF
C
                  EC = -DBLE(Z(IT)**2)/(2.0D0*NQN*NQN)
C
                  ISTART = 1
 55               CONTINUE
                  IF ( EC.LE.ELIM ) EC = ELIM*0.7D0
C
C                                      --------------------
C                                         FIND    NZERO
C                                      --------------------
                  DO N = 1,(NRC-1)
                     IF ( (VV(N)-EC)*RC(N)**2.GT.UNEND ) THEN
                        IF ( MOD(N,2).EQ.0 ) THEN
                           NZERO = N + 1
                        ELSE
                           NZERO = N
                        END IF
                        GOTO 60
                     END IF
                  END DO
                  NZERO = NRC - 1
                  WRITE (6,99003) ITPRT,NQN,L,(NRC-1)
                  STOP
C                                      --------------------
C                                         FIND    NMATCH
C                                      --------------------
 60               CONTINUE
                  N = NZERO + 1
                  DO NN = 1,NZERO
                     N = N - 1
                     IF ( (VV(N)+LLL/RC(N)**2-EC).LT.0.0 ) THEN
                        NMATCH = N
                        GOTO 65
                     END IF
                  END DO
                  if(t_inc%i_write>0) WRITE (1337,99004) ITPRT,NQN,L,EC
C
 65               CONTINUE
                  CALL COREDIR(IT,CTL(IT,ILC),EC,L,MJ,'OUT',VV,BB,RC,
     &                         DRDIC,DOVRC,NMATCH,NZERO,GC,FC,DP,DQ,WP,
     &                         WQ,POW,QOW,PIW,QIW,CGD,CGMD,CGO,NRC,
     &                         Z(IT),NUCLEUS)
C
                  NODE = 0
                  DO N = 2,NMATCH
                     IF ( GC(S,S,N)*GC(S,S,N-1).LT.0.0 ) NODE = NODE + 1
                  END DO
C
C
                  IF ( IPRINT.GE.2 .and. (t_inc%i_write>0)) 
     &         WRITE (1337,99016) ITPRT,NQN,L,
     &         KAP(S),(2*MUEM05+1),IC,ISH,0,EC,NMATCH,RC(NMATCH),NZERO,
     &                 RC(NZERO),NODE,(GC(S,S,NMATCH)/GC(S,S,NMATCH-1))
C
                  IF ( NODE.NE.(NQN-L-1) ) THEN
                     IF ( NODE.GT.(NQN-L-1) ) THEN
                        EC = 1.2D0*EC
                     ELSE
                        EC = 0.8D0*EC
                     END IF
                     GOTO 55
                  ELSE IF ( (GC(S,S,NMATCH)/GC(S,S,NMATCH-1).LE.0.0)
     &                      .OR. 
     &                      (GC(S,S,NMATCH)/GC(S,S,NMATCH-1).GE.1.0) )
     &                      THEN
                     EC = 0.9D0*EC
                     GOTO 55
                  END IF
C
C
                  CALL COREDIR(IT,CTL(IT,ILC),EC,L,MJ,'INW',VV,BB,RC,
     &                         DRDIC,DOVRC,NMATCH,NZERO,GC,FC,DP,DQ,WP,
     &                         WQ,POW,QOW,PIW,QIW,CGD,CGMD,CGO,NRC,
     &                         Z(IT),NUCLEUS)
C
C                                      --------------------
C                                       START VALUES FOR
C                                           PARAMETERS
C                                      --------------------
C
                  VAR(1) = EC
                  VAR(2) = POW(S,S)/PIW(S,S)
C
                  IF ( NSOL.NE.2 .OR. NVAR.EQ.2 ) THEN
                     DO IV = 3,4
                        ERR(IV) = 0.0D0
                        ERRNEW(IV) = 0.0D0
                        VAR(IV) = 0.0D0
                        VARNEW(IV) = 0.0D0
                        DV(IV) = 0.0D0
                     END DO
                  ELSE IF ( ISH.GE.4 ) THEN
                     DO IV = 1,4
                        VAR(IV) = VARTAB(IV,ISH-2)
                     END DO
                  ELSE
C
                     DO J = 1,NSOL
                        NOW(J) = 0.0D0
                     END DO
                     DO N = 1,NMATCH - 1
                        RR = RC(N)**3
                        DO J = 1,NSOL
                           NOW(J) = NOW(J) + GC(J,J,N)**2*RR
                        END DO
                     END DO
C
                     DO J = 1,NSOL
                        NIW(J) = 0.0D0
                     END DO
                     DO N = NMATCH,NZERO - 1
                        RR = RC(N)**3
                        DO J = 1,NSOL
                           NIW(J) = NIW(J) + GC(J,J,N)**2*RR
                        END DO
                     END DO
C
                     RATT = POW(T,T)/PIW(T,T)
                     VAR(3) = TRYMIX*(NOW(S)+NIW(S)*VAR(2))
     &                        /(NOW(T)+NIW(T)*RATT)
                     VAR(4) = RATT*VAR(3)/VAR(2)
                  END IF
C
C
                  CALL COREERR(ERR,VAR,S,NSOL,POW,QOW,PIW,QIW)
C
                  DO IV = 1,NVAR
                     DV(IV) = VAR(IV)
                  END DO
C
                  ITER = 0
C IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
 70               CONTINUE
                  ITER = ITER + 1
C
                  IF ( SCALEB ) THEN
                     SCALE = MIN(1.0D0,0.1D0*ITER)
                     IF ( NSOL.EQ.1 ) SCALE = 1.0D0
                     DO N = 1,JWS(IM)
                        BB(N) = BT(N,IT)*SCALE
                     END DO
                  END IF
C                         ----------------------------------
C                         CHECK WHETHER NUMBER OF NODES O.K.
C                         ----------------------------------
                  IF ( ITER.GT.1 ) THEN
                     NODE = 0
                     DO N = 2,(NMATCH-1)
                        IF ( GC(S,S,N)*GC(S,S,N-1).LT.0.0 )
     &                       NODE = NODE + 1
                     END DO
                     IF ( IPRINT.GE.2 .and. (t_inc%i_write>0)) 
     &                    WRITE (1337,99016) ITPRT,NQN,L,
     &                    KAP(S),(2*MUEM05+1),IC,ISH,ITER,EC,NMATCH,
     &                    RC(NMATCH),NZERO,RC(NZERO),NODE,
     &                    (GC(S,S,NMATCH)/GC(S,S,NMATCH-1))
C
                     IF ( NODE.NE.(NQN-L-1) ) THEN
                        IF ( NODE.GT.(NQN-L-1) ) THEN
                           EC = 1.2D0*EC
                        ELSE
                           EC = 0.8D0*EC
                        END IF
                        ISTART = ISTART + 1
                        IF ( ISTART.LT.20 ) GOTO 55
                     END IF
                  END IF
C
                  DO IV = 2,NVAR
                     DO JV = 1,NVAR
                        VARNEW(JV) = VAR(JV)
                     END DO
                     VARNEW(IV) = VAR(IV) + DV(IV)*DVSTEP
C
                     IF ( ABS(VAR(IV)).GT.1D-16 ) THEN
                        IF ( ABS(DV(IV)/VAR(IV)).LT.TOLVAR ) VARNEW(IV)
     &                       = VAR(IV)
     &                         *(1.0D0+DSIGN(DVSTEP*TOLVAR,DV(IV)))
                     ELSE IF ( FERRO ) THEN
                       if(t_inc%i_write>0) WRITE(1337,99011) ' VAR(',IV,
     &                                  ') = 0 for (T,N,L,K,M;S,NSOL) ',
     &                                  ITPRT,NQN,L,KAP(S),(2*MUEM05+1),
     &                                  '/2  ',S,NSOL,'  --- suppress B'
                        LOOP = 2
                        SUPPRESSB = .TRUE.
                        GOTO 50
                     ELSE IF ( SUPPRESSB ) THEN
                        WRITE (6,*) 'suppressing B did not help !!'
                        STOP 'in <CORE>'
                     END IF
C
                     CALL COREERR(ERRNEW,VARNEW,S,NSOL,POW,QOW,PIW,QIW)
C
                     DO IE = 1,NVAR
                        IF ( ABS(ERRNEW(IE)-ERR(IE)).LT.1D-16 ) THEN
                           DEDV(IE,IV) = 0.0D0
                           IF ( (IE.EQ.IV) .AND. .NOT.FERRO )
     &                          DEDV(IE,IV) = 1.0D0
                        ELSE
                           DEDV(IE,IV) = (ERRNEW(IE)-ERR(IE))
     &                        /(VARNEW(IV)-VAR(IV))
                        END IF
                     END DO
                  END DO
C
                  DO JV = 1,NVAR
                     VARNEW(JV) = VAR(JV)
                  END DO
                  VARNEW(1) = VAR(1) + DV(1)*DVSTEP
                  IF ( ABS(DV(1)/VAR(1)).LT.TOLVAR ) VARNEW(1) = VAR(1)
     &                 *(1.0D0+DSIGN(DVSTEP*TOLVAR,DV(1)))
                  CALL COREDIR(IT,CTL(IT,ILC),VARNEW(1),L,MJ,'OUT',VV,
     &                         BB,RC,DRDIC,DOVRC,NMATCH,NZERO,GC,FC,DP,
     &                         DQ,WP,WQ,POW,QOW,PIW,QIW,CGD,CGMD,CGO,
     &                         NRC,Z(IT),NUCLEUS)
                  CALL COREDIR(IT,CTL(IT,ILC),VARNEW(1),L,MJ,'INW',VV,
     &                         BB,RC,DRDIC,DOVRC,NMATCH,NZERO,GC,FC,DP,
     &                         DQ,WP,WQ,POW,QOW,PIW,QIW,CGD,CGMD,CGO,
     &                         NRC,Z(IT),NUCLEUS)
C
                  CALL COREERR(ERRNEW,VARNEW,S,NSOL,POW,QOW,PIW,QIW)
C
                  DO IE = 1,NVAR
                     DEDV(IE,1) = (ERRNEW(IE)-ERR(IE))
     &                            /(VARNEW(1)-VAR(1))
                  END DO
C
                  DO IE = 1,NVAR
                     CALL DCOPY(NVAR,DEDV(1,IE),1,DVDE(1,IE),1)
                  END DO
                  CALL DGETRF(NVAR,NVAR,DVDE,4,IPIV,INFO)
                  CALL DGETRI(NVAR,DVDE,4,IPIV,DEDV,4*4,INFO)
C
                  DO IV = 1,NVAR
                     DV(IV) = 0.0D0
                     DO IE = 1,NVAR
                        DV(IV) = DV(IV) + DVDE(IV,IE)*ERR(IE)
                     END DO
                     VAR(IV) = VAR(IV) - DV(IV)
                  END DO
C
                  IF ( VAR(1).GT.0.0D0 ) THEN
                     IF ( IPRINT.GE.1 .and. (t_inc%i_write>0)) 
     &                     WRITE (1337,*)
     &                     ' warning from <CORE> E=',VAR(1),IT,NQN,L
                     VAR(1) = -0.2D0
                  END IF
C
                  CALL COREDIR(IT,CTL(IT,ILC),VAR(1),L,MJ,'OUT',VV,BB,
     &                         RC,DRDIC,DOVRC,NMATCH,NZERO,GC,FC,DP,DQ,
     &                         WP,WQ,POW,QOW,PIW,QIW,CGD,CGMD,CGO,NRC,
     &                         Z(IT),NUCLEUS)
                  CALL COREDIR(IT,CTL(IT,ILC),VAR(1),L,MJ,'INW',VV,BB,
     &                         RC,DRDIC,DOVRC,NMATCH,NZERO,GC,FC,DP,DQ,
     &                         WP,WQ,POW,QOW,PIW,QIW,CGD,CGMD,CGO,NRC,
     &                         Z(IT),NUCLEUS)
C
                  CALL COREERR(ERR,VAR,S,NSOL,POW,QOW,PIW,QIW)
C
                  EC = VAR(1)
C
                  IF ( IPRINT.GE.2 .and. (t_inc%i_write>0)) 
     &        WRITE (1337,99005) LOOP,SCALE,
     &        VAR(1),(VAR(IV),IV=1,4),(DV(IV),IV=1,4),(ERR(IE),IE=1,4)
C
C----------------------------------  check relative change in parameters
C ----------------------- parameters 3 and 4 = 0 for paramagnetic case !
                  IF ( ITER.LT.ITERMAX ) THEN
                     DO IV = 1,NVAR
                        VARTAB(IV,ISH) = VAR(IV)
C           IF( ABS(VAR(IV)) .EQ. 0.0D0 ) THEN
                        IF ( (ABS(VAR(IV))+ABS(VAR(IV))).LT.1.0D-30 )
     &                       THEN
                           IF ( FERRO .and. (t_inc%i_write>0))
     &           WRITE (1337,'(A,I3,A)') ' VAR ',IV,' = 0 ??????!!!!!'
                        ELSE IF ( ABS(DV(IV)/VAR(IV)).GT.TOLVAR ) THEN
                           GOTO 70
                        END IF
                     END DO
                  ELSE
                     IF( BNDSTACHK ) THEN 
                        ITXRAY = -1
                        RETURN
                     END IF
                     if(t_inc%i_write>0) 
     &                    WRITE (1337,99006) ITERMAX,(VAR(IV),IV=1,4),
     &                               (DV(IV),IV=1,4),(ERR(IE),IE=1,4)
                  END IF
C IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
C
C                         ---------------------------------
C                         NORMALIZE WAVEFUNCTIONS ACCORDING
C                               TO MATCHING CONDITIONS
C                         ---------------------------------
C
C                                    INWARD - SOLUTION
                  DO N = NMATCH,NZERO
                     DO J = 1,NSOL
                        DO I = 1,NSOL
                           GC(I,J,N) = GC(I,J,N)*VAR(2)
                           FC(I,J,N) = FC(I,J,N)*VAR(2)
                        END DO
                     END DO
                  END DO
C
                  IF ( NSOL.EQ.2 ) THEN
C                                   OUTWARD - SOLUTION
                     DO N = 1,(NMATCH-1)
                        DO I = 1,NSOL
                           GC(I,T,N) = GC(I,T,N)*VAR(3)
                           FC(I,T,N) = FC(I,T,N)*VAR(3)
                        END DO
                     END DO
C                                    INWARD - SOLUTION
                     DO N = NMATCH,NZERO
                        DO I = 1,NSOL
                           GC(I,T,N) = GC(I,T,N)*VAR(4)
                           FC(I,T,N) = FC(I,T,N)*VAR(4)
                        END DO
                     END DO
                  END IF
C
C                                    SUM FOR EACH KAPPA
                  DO N = 1,NZERO
                     DO K = 1,NSOL
                        GCK(K,S,N) = 0.0D0
                        FCK(K,S,N) = 0.0D0
                        DO J = 1,NSOL
                           GCK(K,S,N) = GCK(K,S,N) + GC(K,J,N)
                           FCK(K,S,N) = FCK(K,S,N) + FC(K,J,N)
                        END DO
                     END DO
                  END DO
C
C                       -----------------------------------
C                       CALCULATE  NORM  AND NORMALIZE TO 1
C                       -----------------------------------
                  DO K = 1,NSOL
                     NORM = R2DRDIC(1)*(GCK(K,S,1)**2+FCK(K,S,1)**2)
                  END DO
C
                  SIMP = -1.0D0
                  DO N = 2,NZERO
                     SIMP = -SIMP
                     W = (3.0D0+SIMP)*R2DRDIC(N)
                     DO K = 1,NSOL
                        NORM = NORM + W*(GCK(K,S,N)**2+FCK(K,S,N)**2)
                     END DO
                  END DO
C
                  N = NZERO
                  DO K = 1,NSOL
                     NORM = NORM - R2DRDIC(N)
     &                      *(GCK(K,S,N)**2+FCK(K,S,N)**2)
                  END DO
                  NORM = NORM/3.0D0
C
                  DO K = 1,NSOL
                     NORM = NORM + 0.5D0*R2DRDIC(1)
     &                      *(GCK(K,S,1)**2+FCK(K,S,1)**2)
                  END DO
                  NORM = 1.0D0/SQRT(NORM)
C
                  DO N = 1,NZERO
                     DO K = 1,NSOL
                        GCK(K,S,N) = GCK(K,S,N)*NORM
                        FCK(K,S,N) = FCK(K,S,N)*NORM
                     END DO
                  END DO
                  IF ( NZERO.LT.JTOP ) THEN
                     DO N = (NZERO+1),JTOP
                        DO K = 1,NSOL
                           GCK(K,S,N) = 0.0D0
                           FCK(K,S,N) = 0.0D0
                        END DO
                     END DO
                  END IF
C
                  CALL RINIT(NRMAX,RINT)
C
                  DO N = 1,JTOP
                     DO K = 1,NSOL
                        RINT(N) = RINT(N) + R2DRDI(N,IM)
     &                            *(GCK(K,S,N)*GCK(K,S,N)+FCK(K,S,N)
     &                            *FCK(K,S,N))
                     END DO
                  END DO
C
                  CALL RINTSIMP(RINT,JTOP,AUX)
C
C ------------------------------ omit normalization for XRAY calculation
C ------------------------------ to recover old (errounous data) -------
                  IF ( ITXRAY.GT.0 ) THEN
                     NORM = 1.0D0
                  ELSE
                     NORM = 1.0D0/SQRT(AUX)
                  END IF
C
                  DO N = 1,MAX(NZERO,JTOP)
                     DO K = 1,NSOL
                        GCK(K,S,N) = GCK(K,S,N)*NORM
                        FCK(K,S,N) = FCK(K,S,N)*NORM
                     END DO
                  END DO
C
C                       -----------------------------------
C                       CALCULATE  CHARGE AND SPIN DENSITY
C                       -----------------------------------
C
                  DO N = 1,JWS(IM)
                     DO K = 1,NSOL
                        RHOCHR(N,IT) = RHOCHR(N,IT)
     &                                 + (GCK(K,S,N)**2+FCK(K,S,N)**2)
                        RHOSPN(N,IT) = RHOSPN(N,IT)
     &                                 + (GCK(K,S,N)**2*CGD(K)
     &                                 -FCK(K,S,N)**2*CGMD(K))
                     END DO
                  END DO
C
                  IF ( NSOL.GT.1 ) THEN
                     DO N = 1,JWS(IM)
                        RHOSPN(N,IT) = RHOSPN(N,IT) + GCK(1,S,N)
     &                                 *GCK(2,S,N)*CGO*2
                     END DO
                  END IF
C
C
C                       -----------------------------------
C                            CALCULATE  SPIN CHARACTER
C                       -----------------------------------
C
                  W = R2DRDIC(1)
                  SZ = 0.0D0
                  DO K = 1,NSOL
                     SZ = SZ + W*(GCK(K,S,1)**2*CGD(K)+FCK(K,S,1)
     &                    **2*CGMD(K))
                  END DO
C
                  SIMP = -1.0D0
                  DO N = 2,NZERO
                     SIMP = -SIMP
                     W = (3.0D0+SIMP)*R2DRDIC(N)
                     DO K = 1,NSOL
                        SZ = SZ + W*(GCK(K,S,N)**2*CGD(K)+FCK(K,S,N)
     &                       **2*CGMD(K))
                     END DO
                  END DO
C
                  N = NZERO
                  W = R2DRDIC(N)
                  DO K = 1,NSOL
                     SZ = SZ + W*(GCK(K,S,N)**2*CGD(K)+FCK(K,S,N)
     &                    **2*CGMD(K))
                  END DO
C
C
                  IF ( NSOL.GT.1 ) THEN
C
                     W = R2DRDIC(1)
                     SZ = SZ + W*GCK(1,S,1)*GCK(2,S,1)*CGO*2
C
                     SIMP = -1.0D0
                     DO N = 2,NZERO
                        SIMP = -SIMP
                        W = (3.0D0+SIMP)*R2DRDIC(N)
                        DO K = 1,NSOL
                           SZ = SZ + W*GCK(1,S,N)*GCK(2,S,N)*CGO*2
                        END DO
                     END DO
C
                     N = NZERO
                     W = R2DRDIC(N)
                     SZ = SZ - W*GCK(1,S,N)*GCK(2,S,N)*CGO*2
C
                  END IF
C
                  SZ = SZ/3.0D0
C
C
C                         ------------------------------
C                         CALCULATE   HYPERFINE - FIELD
C                         ------------------------------
C
                  CALL COREHFF(KAP1,KAP2,MJ,S,NSOL,BHF,GCK,FCK,RC,DRDIC,
     &                         0.0D0,NZERO,NRC)
                  IF ( NUCLEUS.NE.0 ) THEN
                     CALL COREHFF(KAP1,KAP2,MJ,S,NSOL,BHF1,GCK,FCK,RC,
     &                            DRDIC,RNUC,JLIM,NRC)
                     CALL COREHFF(KAP1,KAP2,MJ,S,NSOL,BHF2,GCK,FCK,RC,
     &                            DROVRN,RNUC,JLIM,NRC)
                     DO I = 1,NSOL
                        DO J = 1,NSOL
                           BHF(I,J) = BHF(I,J) - BHF1(I,J) + BHF2(I,J)
                        END DO
                     END DO
                  END IF        !end of nucleus.eq.0
C
                  BSOL = 0.0D0
                  DO J = 1,NSOL
                     DO I = 1,NSOL
                        BSOL = BSOL + BHF(I,J)
                        BSH = BSH + BHF(I,J)
                        BCOR(IT) = BCOR(IT) + BHF(I,J)
                     END DO
                  END DO
                  IF ( KAP1.EQ.-1 ) BCORS(IT) = BCORS(IT) + BHF(1,1)
C
                  ECORTAB(IC,IT) = EC
C
C     ------------------
C     SPLIT HFF-FIELD
C     ------------------
                  IF ( ISMQHFI.EQ.1 ) THEN
                     CALL HFFCORE(RNUC,NZERO,KAP1,KAP2,NSOL,MJ,GCK,FCK,
     &                            NRC,SHF,S,NMEMAX,NKMMAX,RC,DRDIC,SDIA,
     &                            SMDIA,SOFF,SMOFF,QDIA,QOFF,QMDIA,
     &                            QMOFF,NUCLEUS,JLIM)
C
                     DO K = 1,NMEMAX
                        SPLIT(K) = 0.0D0
                        DO J = 1,NSOL
                           DO I = 1,NSOL
                              SPLIT(K) = SPLIT(K) + SHF(I,J,K)
                              SPLIT1(K) = SPLIT1(K) + SHF(I,J,K)
                              SPLIT2(K,IT) = SPLIT2(K,IT) + SHF(I,J,K)
                           END DO
                        END DO
                     END DO
                     DO K = 1,NMEMAX
                        IF ( KAP1.EQ.-1 ) SPLIT3(K,IT) = SPLIT3(K,IT)
     &                       + SHF(1,1,K)
                     END DO
                  END IF
CMBE
C
C---------------------------------------------------- l-shell UNCOMPLETE
                  IF ( ISH.GE.NSH ) THEN
C----------------------------------------------------- l-shell completed
                     IF ( ITXRAY.EQ.0 .AND. IPRINT.GT.0 .and. 
     &                                         (t_inc%i_write>0)) THEN
                        WRITE (1337,99012) ITPRT,NQN,TXTL(L),
     &                                  TXTK(IABS(KAP(S))),(2*MUEM05+1),
     &                                  KAP(S),ITER,EC,BSOL*.001D0,
     &                                  BSH*.001D0
                        IF ( ISMQHFI.EQ.1 ) THEN
                           DO K = 1,NMEMAX
                              WRITE (1337,99013) TXTB(K),SPLIT(K)*.001D0
     &                               ,SPLIT1(K)*.001D0
                           END DO
                           WRITE (1337,99014) 'total error in %',
     &                            100.0D0*(1.0D0-SPLIT(4)/SPLIT(5))
                        END IF
                     END IF
C                              ----------------------------
C                              CHECK CONSISTENCY OF RESULTS
C                              ----------------------------
                     IF ( L.NE.0 ) THEN
                        IC1 = IC - NSH + 1
                        IC2 = IC
                        IF ( ECORTAB(IC2,IT).GE.ECORTAB(IC1,IT) ) THEN
                           IMIN = IC1
                           VZ = +1.0D0
                        ELSE
                           IMIN = IC2
                           VZ = -1.0D0
                        END IF
                        IFLAG = 0
                        II = 1
                        DO I = IC1 + 1,IC2,2
                           IF ( VZ*(ECORTAB(I,IT)-ECORTAB(I-II,IT))
     &                          .LT.0.0 ) IFLAG = 1
                           II = 2
                        END DO
                        IF ( ECORTAB(IC1+2,IT).GT.ECORTAB(IMIN,IT) )
     &                       IFLAG = 1
                        DO I = IC1 + 4,IC2 - 1,2
                           IF ( ECORTAB(I,IT).GT.ECORTAB(IMIN,IT) )
     &                          IFLAG = 1
                           IF ( VZ*(ECORTAB(I,IT)-ECORTAB(I-II,IT))
     &                          .GT.0.0 ) IFLAG = 1
                        END DO
C
                        IF ( FERRO .AND. (IFLAG.EQ.1) ) THEN
                           if(t_inc%i_write>0) WRITE (1337,99007)
                           SCALEB = .TRUE.
                           IF ( LOOP.EQ.1 ) THEN
                              LOOP = 2
                            if(t_inc%i_write>0) WRITE(1337,99008) ITPRT
                              GOTO 50
                           END IF
                        END IF
C
                     END IF
                  ELSE IF ( ITXRAY.EQ.0 .AND. IPRINT.GT.0 ) THEN
                     if(t_inc%i_write>0) WRITE (1337,99012) ITPRT,NQN,
     &                               TXTL(L),TXTK(IABS(KAP(S))),
     &                               (2*MUEM05+1),KAP(S),ITER,EC,
     &                               BSOL*.001D0
                     IF ( ISMQHFI.EQ.1 ) THEN
                        DO K = 1,NMEMAX
                           if(t_inc%i_write>0) WRITE (1337,99013) 
     &                                         TXTB(K),SPLIT(K)*.001D0
                        END DO
                        if(t_inc%i_write>0) WRITE (1337,99014) 
     &                                 'total error in %',
     &                                  100.0D0*(1.0D0-SPLIT(4)/SPLIT(5)
     &                                  )
                     END IF
                  END IF
C-----------------------------------------------------------------------
C
                  IF ( IPRINT.GE.1 .and. (t_inc%i_write>0)) 
     &                 WRITE (1337,99015)
     &                 ((BHF(I,J)*.001D0,I=1,NSOL),J=1,NSOL)
C
C
C                          --------------------------------
C                            IF THE SWITCH CHECK IS SET:
C                            RECALCULATE THE EIGENVALUE
C                          USING THE CONVENTIONAL ALGORITHM
C                          --------------------------------
                  IF ( CHECK ) THEN
                     ECC = 0.95D0*EC
 72                  CONTINUE
                     CALL COREDIR(IT,CTL(IT,ILC),ECC,L,MJ,'OUT',VV,BB,
     &                            RC,DRDIC,DOVRC,NMATCH,NZERO,GC,FC,DP,
     &                            DQ,WP,WQ,POW,QOW,PIW,QIW,CGD,CGMD,CGO,
     &                            NRC,Z(IT),NUCLEUS)
                     CALL COREDIR(IT,CTL(IT,ILC),ECC,L,MJ,'INW',VV,BB,
     &                            RC,DRDIC,DOVRC,NMATCH,NZERO,GC,FC,DP,
     &                            DQ,WP,WQ,POW,QOW,PIW,QIW,CGD,CGMD,CGO,
     &                            NRC,Z(IT),NUCLEUS)
C
                     NORM = POW(S,S)/PIW(S,S)
                     DO N = NMATCH,NZERO
                        GC(S,S,N) = GC(S,S,N)*NORM
                        FC(S,S,N) = FC(S,S,N)*NORM
                     END DO
C
                     NORM = 0.0D0
                     DO N = 3,NZERO,2
                        NORM = NORM + R2DRDIC(N)
     &                         *(GC(S,S,N)**2+FC(S,S,N)**2)
     &                         + 4.D0*R2DRDIC(N-1)
     &                         *(GC(S,S,N-1)**2+FC(S,S,N-1)**2)
     &                         + R2DRDIC(N-2)
     &                         *(GC(S,S,N-2)**2+FC(S,S,N-2)**2)
                     END DO
                     NORM = NORM/3.0D0
C
                     LCP1 = MIN(NLMAX,L+1)
                     DEC = POW(S,S)
     &                     *(QOW(S,S)-RC(NMATCH)*CTL(IT,LCP1)*FC(S,S,
     &                     NMATCH))/NORM
                     ECC = ECC + DEC
                     IF ( ABS(DEC/ECC).GT.TOLVAR ) GOTO 72
                     if(t_inc%i_write>0) 
     &                  WRITE (1337,'(7X,''CHECK-E:'',10X,F12.5,/)') ECC
                  END IF
C
C
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C           STORE CORE WAVE FUNCTIONS IF REQUIRED
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                  IF ( ITXRAY.NE.0 ) THEN
                     IF ( NSOL.EQ.2 ) THEN
                        IF ( S.EQ.2 ) THEN
                           ICST = L + 1 + MUEM05
                        ELSE
                           ICST = L + 1 + MUEM05 + 2*L + 1
                        END IF
                     ELSE IF ( MJ.LT.0 ) THEN
                        ICST = 2*L + 1
                     ELSE
                        ICST = 4*L + 2
                     END IF
C
                     MM05COR(ICST) = MUEM05
                     NKPCOR(ICST) = NSOL
                     KAPCOR(ICST) = KAP(S)
                     IKMCOR(ICST,1) = IKAPMUE(KAP(S),MUEM05)
                     IZERO(ICST) = MIN(NZERO,JWS(IM))
                     SZCOR(ICST) = SZ
                     ECOR(ICST) = EC
C
                     DO N = 1,IZERO(ICST)
                        GCOR(N,1,ICST) = GCK(S,S,N)
                        FCOR(N,1,ICST) = FCK(S,S,N)
                     END DO
                     IF ( NSOL.EQ.2 ) THEN
                        DO N = 1,IZERO(ICST)
                           GCOR(N,2,ICST) = GCK(T,S,N)
                           FCOR(N,2,ICST) = FCK(T,S,N)
                        END DO
                        IKMCOR(ICST,2) = IKAPMUE(KAP(T),MUEM05)
                     END IF
                  END IF
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C
C
               END DO
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
            END DO
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
 100     END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
C
         IF ( ITXRAY.EQ.0 .AND. IPRINT.GT.0 .and. (t_inc%i_write>0)) 
     &        WRITE (1337,99009) 
     &        BCOR(IT)*.001D0,BCORS(IT)*.001D0
         IF ( ISMQHFI.EQ.1 ) THEN
            DO N = 1,NMEMAX
              if(t_inc%i_write>0) WRITE (1337,99010) 
     &                            SPLIT2(N,IT)*.001D0,
     &                            SPLIT3(N,IT)*.001D0
            END DO
            if(t_inc%i_write>0) WRITE (1337,99014) 'total error',
     &                      100.0D0*(1.0D0-SPLIT2(4,IT)/SPLIT2(5,IT))
            if(t_inc%i_write>0) WRITE (1337,99014) 'total error',
     &                      100.0D0*(1.0D0-SPLIT3(4,IT)/SPLIT3(5,IT))
         END IF
C
         IF ( (ITXRAY.EQ.0) .OR. (IT.EQ.ITXRAY) ) THEN
            DO N = 1,JTOP
               RINT(N) = RHOCHR(N,IT)*R2DRDI(N,IM)
            END DO
            CALL RINTSIMP(RINT,JTOP,AUX)
            IF( IPRINT.GT. -2 .and. (t_inc%i_write>0)) 
     &                          WRITE (1337,99017) 'charge',ITPRT,AUX
            DO N = 1,JTOP
               RINT(N) = RHOSPN(N,IT)*R2DRDI(N,IM)
            END DO
            CALL RINTSIMP(RINT,JTOP,AUX)
            IF( IPRINT.GT. -2 .and. (t_inc%i_write>0)) 
     &                          WRITE (1337,99017) ' spin ',ITPRT,AUX
         END IF
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
99001 FORMAT (/,10X,'potential is not exchange split ')
99002 FORMAT (/,'  ATOM   IT      : ',I5,/,'  ATOMIC NUMBER  : ',I5,//,
     &        '  IT',10X,'MUE  KAP ITER    ENERGY       B  (k-Gauss)  ')
99003 FORMAT ('  IT=',I2,' NQN=',I2,' L=',I2,
     &        '  NZERO set to  (NRC-1) =',I4)
99004 FORMAT (//,'  STOP IN <<CORE>>',/,'  IT=',I2,' NQN=',I2,' L=',I2,
     &        /,'  no matching-radius found for  EC=',F10.3)
99005 FORMAT (' LOOP    =  ',I3,' BSCL=',F10.5,/,' E=',F14.7,' VAR  ',
     &        4E11.4,/,17X,' CORR ',4E11.4,/,17X,' ERR  ',4E11.4)
99006 FORMAT (' iteration not converged after',I3,' steps !',/,
     &        ' parameters:',4E18.10,/,' last corr.:',4E18.10,/,
     &        ' last error:',4E18.10)
99007 FORMAT (' >>> check data E(KAP,MJ) should be monotonous ',
     &        ' and  E(+L,MJ) < E(-L-1,MJ) ',//)
99008 FORMAT (' all core states for atom type ',I2,
     &        ' will be recalculated ',/,
     &   ' with the spin dependent exchange field gradually switched on'
     &   )
99009 FORMAT (2X,57('-'),/,42X,F17.3,/,38X,'(S) ',F17.3,/,2X,57('*'),//)
99010 FORMAT (2X,57('-'),/,37X,F17.3,/,33X,'(S) ',F17.3,/,2X,57('*'),//)
99011 FORMAT (A,I1,A,5I3,A,2I2,A)
99012 FORMAT (2I4,A1,A3,I3,'/2',2I4,2X,F15.8,F17.3,:,F17.3,/)
99013 FORMAT (A,:,32X,F17.3,:,F17.3,/)
99014 FORMAT (A,:,37X,F6.3)
99015 FORMAT (37X,F17.3)
99016 FORMAT (/,' IT=',I2,'  NQN=',I2,'  L=',I2,'  KAP=',I2,'  MJ=',I2,
     &        '/2    IC=',I3,'  ISH=',I2,/,' E(',I2,')   =',F15.5,/,
     &        ' NMATCH  =',I5,'    R=',F10.5,/,' NZERO   =',I5,'    R=',
     &        F10.5,/,' NODES   =',I5,'  RAT=',E11.4)
99017 FORMAT (' integrated core ',A,' density for atom type ',I4,':',
     &        F12.8)
      END
