      SUBROUTINE DIRABMSOC(GETIRRSOL,C,SOCSCL,IT,E,L,MJ,KAP1,KAP2,PIS,
     &                     CG1,CG2,CG4,CG5,CG8,V,B,Z,NUCLEUS,R,DRDI,
     &                     DOVR,NMESH,DXP,PR,QR,PI,QI,DP,DQ,NRMAX)
C   ********************************************************************
C   *                                                                  *
C   *   ROUTINE TO SOLVE THE SPIN-POLARISED RADIAL DIRAC EQUATIONS     *
C   *                                                                  *
C   *               scaling the SPIN-ORBIT-COUPLING                    *
C   *                                                                  *
C   *   the outward integration is started by a power expansion        *
C   *   and continued by ADAMS-BASHFORTH-MOULTON - pred./corr.-method  *
C   *   NABM = 4(5) selects the 4(5)-point formula                     *
C   *                                                                  *
C   *   the inward integration is started analytically                 *
C   *                                                                  *
C   *   returns the wave functions up to the mesh point NMESH          *
C   *   PR,QR and PI,QI  with   P=r*g and Q=r*c*f                      *
C   *   and    R/I standing for regular/irregular solution             *
C   *                                                                  *
C   *  19/12/94  HE                                                    *
C   *  28/06/95  HF: corrected init of inward integration              *
C   *  21/01/98  HE  finite nucelus                                    *
C   ********************************************************************
      IMPLICIT NONE
C
C
C PARAMETER definitions
C
      INTEGER MPSMAX,NPEMAX,NABM
      PARAMETER (MPSMAX=40,NPEMAX=4,NABM=4)
      COMPLEX*16 C0
      PARAMETER (C0=(0.0D0,0.0D0))
      REAL*8 TOL
      PARAMETER (TOL=1.0D-9)
      INTEGER ITMAX
      PARAMETER (ITMAX=50)
C
C Dummy arguments
C
      REAL*8 C,CG1,CG2,CG4,CG5,CG8,MJ,SOCSCL
      COMPLEX*16 E
      LOGICAL GETIRRSOL
      INTEGER IT,KAP1,KAP2,L,NMESH,NRMAX,NUCLEUS,Z
      COMPLEX*16 PIS
      REAL*8 B(NRMAX),DOVR(NRMAX),DRDI(NRMAX),R(NRMAX),V(NRMAX)
      COMPLEX*16 DP(2,2,NRMAX),DQ(2,2,NRMAX),DXP(2,2),PI(2,2,NRMAX),
     &           PR(2,2,NRMAX),QI(2,2,NRMAX),QR(2,2,NRMAX)
C
C Local variables
C
      COMPLEX*16 AA11,AA12,AA21,AA22,ARG,BB1,BB2,BPP,BQQ,CFAC,CGO,D14,
     &           DH,DIFFA,DIFFB,EMVPP,EMVQQ
      REAL*8 ACORR(0:NABM-1),ACORR0(0:NABM-1),APRED(NABM),APRED0(NABM),
     &       ASTEP,B14,BC(0:NPEMAX),BH,BHLP(NABM+4),CGD(2),CGMD(2),
     &       CM(NPEMAX,NPEMAX),CMI(NPEMAX,NPEMAX),CSQR,DHLP(NABM+4),
     &       GAM(2),GPM,HLP(NABM+4),HLP1,KAP(2),KPX(2),KPY(2),LMK(2),
     &       R14,RH,RHLP(NABM+4),RPWGPM,RR,SK(2),SK1,SK2,SO2,SO6,SRK,TZ,
     &       V14,VC(0:NPEMAX),VH,VHLP(NABM+4),X14,XH
      COMPLEX*16 CJLZ
      DOUBLE PRECISION DABS,DBLE,DSQRT
      COMPLEX*16 DETD,MP1(2,2),MP2(2,2),MP3(2,2),MP4(2,2),MQ1(2,2),
     &           MQ2(2,2),MQ3(2,2),MQ4(2,2),P1(2,2),P2(2,2),P3(2,2),
     &           P4(2,2),PC(2,2,-NPEMAX:MPSMAX),PNEW(2,2),POLD(2,2),
     &           Q1(2,2),Q2(2,2),Q3(2,2),Q4(2,2),QC(2,2,-NPEMAX:MPSMAX),
     &           QNEW(2,2),QOLD(2,2),S0,SOCPP(2),T0,ZZ
      INTEGER I,IC,IP,IRK,ISK1,ISK2,IV,J,JCORR,K,LB(2),LB1,LB2,M,MPS,N,
     &        NACORR,NDIV,NHLP,NM,NPE,NSOL,NTOP
      INTEGER INT,ISIGN,NINT
      REAL*8 YLAG
C
      DATA APRED0/55.0D0, - 59.0D0, + 37.0D0, - 9.0D0/
      DATA ACORR0/9.0D0, + 19.0D0, - 5.0D0, + 1.0D0/
      DATA ASTEP/24.0D0/
C

C
      CSQR = C*C
      CFAC = PIS*C/(E+CSQR)
C
C find   NPE  expansion coefficients for the potential and b-field
      NPE = 4
C
      TZ = DBLE(2*Z)
C
      DO IV = 1,NPE
         DO N = 1,NPE
            CM(N,IV) = R(N)**(IV-1)
         END DO
      END DO
C
      CALL RINVGJ(CMI,CM,NPEMAX,NPE)
C
      DO IV = 1,NPE
         VC(IV-1) = 0.0D0
         DO N = 1,NPE
            VC(IV-1) = VC(IV-1) + CMI(IV,N)*(V(N)+TZ/R(N))
         END DO
      END DO
C
      DO IV = 1,NPE
         BC(IV-1) = 0.0D0
         DO N = 1,NPE
            BC(IV-1) = BC(IV-1) + CMI(IV,N)*B(N)
         END DO
      END DO
C
C
C    calculate g-coefficients of b-field
C
      ISK1 = ISIGN(1,KAP1)
      ISK2 = ISIGN(1,KAP2)
      SK1 = DBLE(ISK1)
      SK2 = DBLE(ISK2)
      LB1 = L - ISK1
      LB2 = L - ISK2
C
      CG1 = -MJ/(KAP1+0.5D0)
      CG5 = -MJ/(-KAP1+0.5D0)
      CGD(1) = CG1
      CGMD(1) = CG5
      KAP(1) = DBLE(KAP1)
      GAM(1) = DSQRT(KAP(1)**2-(TZ/C)**2)
      LB(1) = LB1
      SK(1) = SK1
      IF ( DABS(MJ).GT.L ) THEN
         CG2 = 0.0D0
         CG4 = 0.0D0
         CG8 = 0.0D0
         NSOL = 1
         CGD(2) = 0.0D0
         CGO = 0.0D0
         CGMD(2) = 0.0D0
         GAM(2) = 0.0D0
         KAP(2) = 0.0D0
         LB(2) = 0
         SK(2) = 0.0D0
      ELSE
         CG2 = -DSQRT(1.0D0-(MJ/(KAP1+0.5D0))**2)
         CG4 = -MJ/(KAP2+0.5D0)
         CG8 = -MJ/(-KAP2+0.5D0)
         NSOL = 2
         CGD(2) = CG4
         CGO = CG2
         CGMD(2) = CG8
         KAP(2) = DBLE(KAP2)
         GAM(2) = DSQRT(KAP(2)**2-(TZ/C)**2)
         LB(2) = LB2
         SK(2) = SK2
      END IF
      DO I = 1,NSOL
         KPX(I) = -1.0D0 + SOCSCL*DBLE(1+KAP(I))
         LMK(I) = DBLE(L*(L+1)) - KPX(I)*(KPX(I)+1.0D0)
        
         CGMD(I)= 0.0D0
c-------------------------------------- causes numerical inconsistencies
c        GAM(I) = DSQRT( KPX(I)**2 - (TZ/C)**2 )
C        KPY(I) = KPX(I)
C
         GAM(I) = DSQRT( KAP(I)**2 - (TZ/C)**2 )
         KPY(I) = KAP(I)
      END DO
C
      DO IP = 1,NABM
         IC = IP - 1
         APRED(IP) = APRED0(IP)/ASTEP
         ACORR(IC) = ACORR0(IC)/ASTEP
      END DO
      NACORR = NABM - 1
C
      DO I = 1,2
         DO J = 1,2
            DO IP = -NPEMAX,MPSMAX
               PC(I,J,IP) = C0
               QC(I,J,IP) = C0
            END DO
         END DO
      END DO
C
C ======================================================================
      IF ( (TZ.GE.2) .AND. (NUCLEUS.EQ.0) ) THEN
C
         DO J = 1,NSOL
            I = 3 - J
            PC(J,J,0) = DSQRT(ABS(KPY(J))-GAM(J))
            QC(J,J,0) = (KPY(J)+GAM(J))*(CSQR/TZ)*PC(J,J,0)
            PC(I,J,0) = C0
            QC(I,J,0) = C0
         END DO
C
C  DETERMINE HIGHER EXPANSION COEFFICIENTS FOR THE WAVE FUNCTIONS
C
         MPS = 40
C
         AA12 = -TZ/CSQR
         AA21 = TZ
         EMVQQ = (E-VC(0)+CSQR)/CSQR
         EMVPP = -E + VC(0)
         BQQ = BC(0)/CSQR
         DO I = 1,NSOL
            SOCPP(I) = LMK(I)/(EMVQQ+BQQ*CGMD(I))
         END DO
C
         DO J = 1,NSOL
C
            DO M = 1,MPS
               DO I = 1,NSOL
                  K = 3 - I
                  BB1 = (EMVQQ+BQQ*CGMD(I))*QC(I,J,M-1)
                  BB2 = (EMVPP+BC(0)*CGD(I))*PC(I,J,M-1) + BC(0)
     &                  *CGO*PC(K,J,M-1)
                  DO IP = 1,NPE - 1
                    
                     BB1 = BB1 + (-VC(IP)+BC(IP)*CGMD(I))*QC(I,J,M-1-IP)
     &                     /CSQR
                     BB2 = BB2 + (+VC(IP)+BC(IP)*CGD(I))*PC(I,J,M-1-IP)
     &                     + BC(IP)*CGO*PC(K,J,M-1-IP)
                  END DO
C
                  AA11 = GAM(J) + M + KPY(I)
                  AA22 = GAM(J) + M - KPY(I)
                  DETD = AA11*AA22 - AA12*AA21
                  PC(I,J,M) = (BB1*AA22-AA12*BB2)/DETD
                  QC(I,J,M) = (AA11*BB2-BB1*AA21)/DETD
C
               END DO
            END DO
C
         END DO
C
C
C  PERFORM SUMMATION OVER WAVE FUNCTION - EXPANSION COEFFICIENTS
C  FOR THE FIRST   NABM   R - MESH - POINTS
C
         DO N = 1,NABM
            RR = R(N)
C
            DO J = 1,NSOL
               RPWGPM = RR**GAM(J)
C
               DO I = 1,NSOL
                  PR(I,J,N) = PC(I,J,0)*RPWGPM
                  QR(I,J,N) = QC(I,J,0)*RPWGPM
                  DP(I,J,N) = PC(I,J,0)*RPWGPM*GAM(J)*DOVR(N)
                  DQ(I,J,N) = QC(I,J,0)*RPWGPM*GAM(J)*DOVR(N)
               END DO
C
               DO M = 1,MPS
                  RPWGPM = RPWGPM*RR
                  GPM = GAM(J) + M
C
                  DO I = 1,NSOL
                     PR(I,J,N) = PR(I,J,N) + PC(I,J,M)*RPWGPM
                     QR(I,J,N) = QR(I,J,N) + QC(I,J,M)*RPWGPM
                     DP(I,J,N) = DP(I,J,N) + PC(I,J,M)
     &                           *RPWGPM*GPM*DOVR(N)
                     DQ(I,J,N) = DQ(I,J,N) + QC(I,J,M)
     &                           *RPWGPM*GPM*DOVR(N)
                  END DO
C
               END DO
            END DO
         END DO
C ======================================================================
C                                  == EMPTY SPHERE  or FINITE NUCLEUS ==
      ELSE
C
C        assume constant pot: V=V(1)   ignore coupling: B=0
C
         T0 = E - V(1)
         S0 = (E-V(1))/CSQR + 1
C
         DO N = 1,NABM
C
            DO J = 1,NSOL
               DO I = 1,NSOL
                  PR(I,J,N) = C0
                  QR(I,J,N) = C0
                  DP(I,J,N) = C0
                  DQ(I,J,N) = C0
               END DO
            END DO
C
            ZZ = CDSQRT(S0*T0)*R(N)
C
            DO J = 1,NSOL
               PR(J,J,N) = CJLZ(L,ZZ)*R(N)
               DP(J,J,N) = (DBLE(L+1)*CJLZ(L,ZZ)-ZZ*CJLZ(L+1,ZZ))
     &                     *DRDI(N)
C
               QR(J,J,N) = (DP(J,J,N)/DRDI(N)+PR(J,J,N)*(KPX(J)/R(N)))
     &                     /S0
               DQ(J,J,N) = QR(J,J,N)*(KPX(J)/R(N)) - PR(J,J,N)*T0
            END DO
         END DO
C
      END IF
C ===================================================================
C
C
C =============================================================== N ====
C     CALCULATE ALL NEXT POINTS BY PRE/CORR(ADAMS-BASHFORTH-MOULTON)
C
      DO N = NABM + 1,NMESH
C
C    EVALUATE PREDICTOR
C
         DO J = 1,NSOL
            DO I = 1,NSOL
               PNEW(I,J) = PR(I,J,N-1)
               QNEW(I,J) = QR(I,J,N-1)
C
               DO IP = 1,NABM
                  PNEW(I,J) = PNEW(I,J) + APRED(IP)*DP(I,J,N-IP)
                  QNEW(I,J) = QNEW(I,J) + APRED(IP)*DQ(I,J,N-IP)
               END DO
            END DO
         END DO
C
         EMVQQ = (E-V(N)+CSQR)*DRDI(N)/CSQR
         EMVPP = -(E-V(N))*DRDI(N)
         BQQ = B(N)*DRDI(N)/CSQR
         BPP = B(N)*DRDI(N)
         DO I = 1,NSOL
            SOCPP(I) = LMK(I)*DOVR(N)**2/(EMVQQ+BQQ*CGMD(I))
         END DO
C
C    EVALUATE CORRECTOR
C
C
         DO JCORR = 1,ITMAX
C
            DO J = 1,NSOL
               DO I = 1,NSOL
                  K = 3 - I
                  POLD(I,J) = PNEW(I,J)
                  QOLD(I,J) = QNEW(I,J)
                  DP(I,J,N) = -KPX(I)*PNEW(I,J)*DOVR(N)
     &                        + (EMVQQ+BQQ*CGMD(I))*QNEW(I,J)
                  DQ(I,J,N) = KPX(I)*QNEW(I,J)*DOVR(N)
     &                        + (EMVPP+BPP*CGD(I))*PNEW(I,J)
     &                        + BPP*CGO*PNEW(K,J) + SOCPP(I)*PNEW(I,J)
C
                  PNEW(I,J) = PR(I,J,N-1)
                  QNEW(I,J) = QR(I,J,N-1)
                  DO IC = 0,NACORR
                     PNEW(I,J) = PNEW(I,J) + ACORR(IC)*DP(I,J,N-IC)
                     QNEW(I,J) = QNEW(I,J) + ACORR(IC)*DQ(I,J,N-IC)
                  END DO
               END DO
            END DO
C
            DO J = 1,NSOL
               DO I = 1,NSOL
                  DIFFA = POLD(I,J) - PNEW(I,J)
                  IF ( ABS(DREAL(DIFFA)).GT.(TOL*ABS(DREAL(PNEW(I,J))))
     &                 ) GOTO 50
                  IF ( ABS(DIMAG(DIFFA)).GT.(TOL*ABS(DIMAG(PNEW(I,J))))
     &                 ) GOTO 50
C
                  DIFFB = QOLD(I,J) - QNEW(I,J)
                  IF ( ABS(DREAL(DIFFB)).GT.(TOL*ABS(DREAL(QNEW(I,J))))
     &                 ) GOTO 50
                  IF ( ABS(DIMAG(DIFFB)).GT.(TOL*ABS(DIMAG(QNEW(I,J))))
     &                 ) GOTO 50
               END DO
            END DO
            GOTO 100
C
 50      END DO
         WRITE (6,99001) KAP1,N,R(N),DIFFA,DIFFB,IT,L,INT(2*MJ),'REG'
C
C                   SORRY NOT CONVERGED IN  ITMAX  ITERATIONS
C
C
C
 100     CONTINUE
         DO J = 1,NSOL
            DO I = 1,NSOL
               K = 3 - I
               PR(I,J,N) = PNEW(I,J)
               QR(I,J,N) = QNEW(I,J)
               DP(I,J,N) = -KPX(I)*PNEW(I,J)*DOVR(N)
     &                     + (EMVQQ+BQQ*CGMD(I))*QNEW(I,J)
               DQ(I,J,N) = KPX(I)*QNEW(I,J)*DOVR(N) + (EMVPP+BPP*CGD(I))
     &                     *PNEW(I,J) + BPP*CGO*PNEW(K,J) + SOCPP(I)
     &                     *PNEW(I,J)
            END DO
         END DO
C
C
      END DO
C =============================================================== N ====
C
      DO I = 1,2
         DO J = 1,2
            DXP(I,J) = DP(I,J,NMESH)
         END DO
      END DO
C
      IF ( .NOT.GETIRRSOL ) RETURN
C
C #####################################################################
C #####################################################################
C #####################################################################
C
C             IRREGULAR SOLUTION IRREGULAR SOLUTION  IRREGULAR SOLUTION
C
C
C  CALCULATE THE INITIAL VALUES OF THE WAVEFUNCTION AT THE SPHERE
C  BOUNDARY
C
C
      DO N = NMESH,NMESH + NABM
         ARG = PIS*R(N)
C
         DO J = 1,NSOL
            I = 3 - J
            PI(J,J,N) = CJLZ(L,ARG)*R(N)
            QI(J,J,N) = CFAC*SK(J)*CJLZ(LB(J),ARG)*R(N)*C
            DP(J,J,N) = (DBLE(L+1)*CJLZ(L,ARG)-ARG*CJLZ(L+1,ARG))
     &                  *DRDI(N)
            M = LB(J)
            DQ(J,J,N) = CFAC*SK(J)
     &                  *(DBLE(M+1)*CJLZ(M,ARG)-ARG*CJLZ(M+1,ARG))
     &                  *DRDI(N)*C
C
            PI(I,J,N) = C0
            QI(I,J,N) = C0
            DP(I,J,N) = C0
            DQ(I,J,N) = C0
         END DO
      END DO
C        ------------------------------------------------------------
C              INITIALIZE INWARD INTEGRATION WITH RUNGE - KUTTA
C        ------------------------------------------------------------
      NDIV = 60
      IF ( NDIV.NE.0 ) THEN
C
         SRK = 1.0D0/DBLE(NDIV)
         SO2 = SRK/2.0D0
         SO6 = SRK/6.0D0
C
         N = NMESH
C
         EMVQQ = (E-V(N)+CSQR)*DRDI(N)/CSQR
         EMVPP = -(E-V(N))*DRDI(N)
         BQQ = B(N)*DRDI(N)/CSQR
         BPP = B(N)*DRDI(N)
         DO I = 1,NSOL
            SOCPP(I) = LMK(I)*DOVR(N)**2/(EMVQQ+BQQ*CGMD(I))
         END DO
C
C *** reinitialize Q using only DP and PI
         DO J = 1,NSOL
            I = 3 - J
            QI(J,J,N) = (DP(J,J,N)+KPX(J)*PI(J,J,N)*DOVR(N))
     &                  /(EMVQQ+BQQ*CGMD(J))
            QI(I,J,N) = C0
         END DO
C
         DO J = 1,NSOL
            DO I = 1,NSOL
               K = 3 - I
               DQ(I,J,N) = KPX(I)*QI(I,J,N)*DOVR(N) + (EMVPP+BPP*CGD(I))
     &                     *PI(I,J,N) + BPP*CGO*PI(K,J,N) + SOCPP(I)
     &                     *PI(I,J,N)
            END DO
         END DO
C
         DO J = 1,NSOL
            DO I = 1,NSOL
               P1(I,J) = PI(I,J,N)
               Q1(I,J) = QI(I,J,N)
               MP1(I,J) = DP(I,J,N)
               MQ1(I,J) = DQ(I,J,N)
            END DO
         END DO
C
         X14 = DBLE(N)
         NHLP = NABM + 4
         HLP1 = DBLE(NMESH-NHLP)
         DO I = 1,NHLP
            HLP(I) = DBLE(I)
            VHLP(I) = V(NMESH-NHLP+I)
            BHLP(I) = B(NMESH-NHLP+I)
            DHLP(I) = DRDI(NMESH-NHLP+I)
            RHLP(I) = R(NMESH-NHLP+I)
         END DO
C
         DO IRK = 1,(NABM-1)*NDIV
C
            XH = X14 - SO2
            VH = YLAG(XH-HLP1,HLP,VHLP,0,3,NHLP)
            BH = YLAG(XH-HLP1,HLP,BHLP,0,3,NHLP)
            RH = YLAG(XH-HLP1,HLP,RHLP,0,3,NHLP)
            DH = YLAG(XH-HLP1,HLP,DHLP,0,3,NHLP)
C
            EMVQQ = (E-VH+CSQR)*DH/CSQR
            EMVPP = -(E-VH)*DH
            BQQ = BH*DH/CSQR
            BPP = BH*DH
            DO I = 1,NSOL
               SOCPP(I) = LMK(I)*(DH/RH)**2/(EMVQQ+BQQ*CGMD(I))
            END DO
            N = NMESH - IRK/NDIV - N + NINT(XH)
C
            DO J = 1,NSOL
               DO I = 1,NSOL
                  P2(I,J) = P1(I,J) - SO2*MP1(I,J)
                  Q2(I,J) = Q1(I,J) - SO2*MQ1(I,J)
               END DO
            END DO
C
            DO J = 1,NSOL
               DO I = 1,NSOL
                  K = 3 - I
                  MP2(I,J) = -KPX(I)*P2(I,J)*DH/RH + (EMVQQ+BQQ*CGMD(I))
     &                       *Q2(I,J)
                  MQ2(I,J) = KPX(I)*Q2(I,J)*DH/RH + (EMVPP+BPP*CGD(I))
     &                       *P2(I,J) + BPP*CGO*P2(K,J) + SOCPP(I)
     &                       *P2(I,J)
               END DO
            END DO
C
            DO J = 1,NSOL
               DO I = 1,NSOL
                  P3(I,J) = P1(I,J) - SO2*MP2(I,J)
                  Q3(I,J) = Q1(I,J) - SO2*MQ2(I,J)
               END DO
            END DO
C
            DO J = 1,NSOL
               DO I = 1,NSOL
                  K = 3 - I
                  MP3(I,J) = -KPX(I)*P3(I,J)*DH/RH + (EMVQQ+BQQ*CGMD(I))
     &                       *Q3(I,J)
                  MQ3(I,J) = KPX(I)*Q3(I,J)*DH/RH + (EMVPP+BPP*CGD(I))
     &                       *P3(I,J) + BPP*CGO*P3(K,J) + SOCPP(I)
     &                       *P3(I,J)
               END DO
            END DO
C
            X14 = X14 - SRK
            V14 = YLAG(X14-HLP1,HLP,VHLP,0,3,NHLP)
            B14 = YLAG(X14-HLP1,HLP,BHLP,0,3,NHLP)
            R14 = YLAG(X14-HLP1,HLP,RHLP,0,3,NHLP)
            D14 = YLAG(X14-HLP1,HLP,DHLP,0,3,NHLP)
C
C
            EMVQQ = (E-V14+CSQR)*D14/CSQR
            EMVPP = -(E-V14)*D14
            BQQ = B14*D14/CSQR
            BPP = B14*D14
            DO I = 1,NSOL
               SOCPP(I) = LMK(I)*(D14/R14)**2/(EMVQQ+BQQ*CGMD(I))
            END DO
C
            DO J = 1,NSOL
               DO I = 1,NSOL
                  P4(I,J) = P1(I,J) - SRK*MP3(I,J)
                  Q4(I,J) = Q1(I,J) - SRK*MQ3(I,J)
               END DO
            END DO
C
            DO J = 1,NSOL
               DO I = 1,NSOL
                  K = 3 - I
                  MP4(I,J) = -KPX(I)*P4(I,J)*D14/R14 + 
     &                       (EMVQQ+BQQ*CGMD(I))*Q4(I,J)
                  MQ4(I,J) = KPX(I)*Q4(I,J)*D14/R14 + (EMVPP+BPP*CGD(I))
     &                       *P4(I,J) + BPP*CGO*P4(K,J) + SOCPP(I)
     &                       *P4(I,J)
               END DO
            END DO
C
            DO J = 1,NSOL
               DO I = 1,NSOL
                  P1(I,J) = P1(I,J)
     &                      - SO6*(MP1(I,J)+2*(MP2(I,J)+MP3(I,J))
     &                      +MP4(I,J))
                  Q1(I,J) = Q1(I,J)
     &                      - SO6*(MQ1(I,J)+2*(MQ2(I,J)+MQ3(I,J))
     &                      +MQ4(I,J))
               END DO
            END DO
C
            DO J = 1,NSOL
               DO I = 1,NSOL
                  K = 3 - I
                  MP1(I,J) = -KPX(I)*P1(I,J)*D14/R14 + 
     &                       (EMVQQ+BQQ*CGMD(I))*Q1(I,J)
                  MQ1(I,J) = KPX(I)*Q1(I,J)*D14/R14 + (EMVPP+BPP*CGD(I))
     &                       *P1(I,J) + BPP*CGO*P1(K,J) + SOCPP(I)
     &                       *P1(I,J)
               END DO
            END DO
C
            IF ( MOD(IRK,NDIV).EQ.0 ) THEN
               N = NMESH - IRK/NDIV
               IF ( ABS(X14-DBLE(N)).GT.1.0D-5 ) THEN
                  WRITE (*,*) ' <DIRAC> RUNGE-KUTTA: ',IRK,NDIV,N,X14
                  STOP
               END IF
               DO J = 1,NSOL
                  DO I = 1,NSOL
                     PI(I,J,N) = P1(I,J)
                     QI(I,J,N) = Q1(I,J)
                     DP(I,J,N) = MP1(I,J)
                     DQ(I,J,N) = MQ1(I,J)
                  END DO
               END DO
C
            END IF
C
         END DO
C
      END IF
C
C =============================================================== N ====
C
C     CALCULATE ALL NEXT POINTS BY PRE/CORR(ADAMS-BASHFORTH-MOULTON)
C
      IF ( NDIV.NE.0 ) THEN
         NTOP = NMESH - NABM
      ELSE
         NTOP = NMESH
      END IF
C
      DO NM = 1,NTOP
         N = 1 + NTOP - NM
C
C    EVALUATE PREDICTOR
C
         DO J = 1,NSOL
            DO I = 1,NSOL
               PNEW(I,J) = PI(I,J,N+1)
               QNEW(I,J) = QI(I,J,N+1)
C
               DO IP = 1,NABM
                  PNEW(I,J) = PNEW(I,J) - APRED(IP)*DP(I,J,N+IP)
                  QNEW(I,J) = QNEW(I,J) - APRED(IP)*DQ(I,J,N+IP)
               END DO
            END DO
         END DO
C
         EMVQQ = (E-V(N)+CSQR)*DRDI(N)/CSQR
         EMVPP = -(E-V(N))*DRDI(N)
         BQQ = B(N)*DRDI(N)/CSQR
         BPP = B(N)*DRDI(N)
         DO I = 1,NSOL
            SOCPP(I) = LMK(I)*DOVR(N)**2/(EMVQQ+BQQ*CGMD(I))
         END DO
C
C    EVALUATE CORRECTOR
C
         DO JCORR = 1,ITMAX
            DO J = 1,NSOL
               DO I = 1,NSOL
                  K = 3 - I
                  POLD(I,J) = PNEW(I,J)
                  QOLD(I,J) = QNEW(I,J)
                  DP(I,J,N) = -KPX(I)*PNEW(I,J)*DOVR(N)
     &                        + (EMVQQ+BQQ*CGMD(I))*QNEW(I,J)
                  DQ(I,J,N) = KPX(I)*QNEW(I,J)*DOVR(N)
     &                        + (EMVPP+BPP*CGD(I))*PNEW(I,J)
     &                        + BPP*CGO*PNEW(K,J) + SOCPP(I)*PNEW(I,J)
C
                  PNEW(I,J) = PI(I,J,N+1)
                  QNEW(I,J) = QI(I,J,N+1)
                  DO IC = 0,NACORR
                     PNEW(I,J) = PNEW(I,J) - ACORR(IC)*DP(I,J,N+IC)
                     QNEW(I,J) = QNEW(I,J) - ACORR(IC)*DQ(I,J,N+IC)
                  END DO
               END DO
            END DO
C
            DO J = 1,NSOL
               DO I = 1,NSOL
                  DIFFA = POLD(I,J) - PNEW(I,J)
                  IF ( ABS(DREAL(DIFFA)).GT.(TOL*ABS(DREAL(PNEW(I,J))))
     &                 ) GOTO 150
                  IF ( ABS(DIMAG(DIFFA)).GT.(TOL*ABS(DIMAG(PNEW(I,J))))
     &                 ) GOTO 150
C
                  DIFFB = QOLD(I,J) - QNEW(I,J)
                  IF ( ABS(DREAL(DIFFB)).GT.(TOL*ABS(DREAL(QNEW(I,J))))
     &                 ) GOTO 150
                  IF ( ABS(DIMAG(DIFFB)).GT.(TOL*ABS(DIMAG(QNEW(I,J))))
     &                 ) GOTO 150
               END DO
            END DO
            GOTO 200
C
 150     END DO
         WRITE (6,99001) KAP1,N,R(N),DIFFA,DIFFB,IT,L,INT(2*MJ),'IRR'
C
C                   SORRY NOT CONVERGED IN  ITMAX  ITERATIONS
C
C
C
C
 200     CONTINUE
         DO J = 1,NSOL
            DO I = 1,NSOL
               K = 3 - I
               PI(I,J,N) = PNEW(I,J)
               QI(I,J,N) = QNEW(I,J)
               DP(I,J,N) = -KPX(I)*PNEW(I,J)*DOVR(N)
     &                     + (EMVQQ+BQQ*CGMD(I))*QNEW(I,J)
               DQ(I,J,N) = KPX(I)*QNEW(I,J)*DOVR(N) + (EMVPP+BPP*CGD(I))
     &                     *PNEW(I,J) + BPP*CGO*PNEW(K,J) + SOCPP(I)
     &                     *PNEW(I,J)
            END DO
         END DO
C
      END DO

99001 FORMAT (' PRE/CORR NOT CONV. IN <DIRABMSOC> ',2I4,F10.7,2X,
     &        4E12.4,3I2,'/2 ',A3)
C =============================================================== N ====
C
c     the minor component for the soc-manipulated wf is meaningless 
c     =>  set it to zero
c
      CALL CINIT(2*2*NRMAX,QR)
      CALL CINIT(2*2*NRMAX,QI)

      RETURN
      END
