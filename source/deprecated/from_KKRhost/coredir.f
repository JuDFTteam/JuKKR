      SUBROUTINE COREDIR(IT,C,E,L,MJ,WAY,VV,BB,RC,DRDIC,DOVRC,NMATCH,
     &                   NZERO,GC,FC,DP,DQ,WP,WQ,POW,QOW,PIW,QIW,CGD,
     &                   CGMD,CGO,NRC,Z,NUCLEUS)
C   ********************************************************************
C   *                                                                  *
C   *   ROUTINE TO SOLVE RADIAL SPIN-POLARISED DIRAC EQUATIONS         *
C   *   FOR THE CORE WAVE FUNCTIONS                                    *
C   *                                                                  *
C   *   SIMILAR TO LOUCKS' METHOD TO SOLVE THE COUPLED SET OF          *
C   *   DIFFERENTIAL EQUATIONS                                         *
C   *                                    HE JAN. 1989                  *
C   *                                                                  *
C   *   ADOPTED FOR FINITE NUCLEUS       MB MAR. 1995                  *
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
C
      use mod_types, only: t_inc
      IMPLICIT NONE
C
C
C PARAMETER definitions
C
      INTEGER MPSMAX,NPEMAX,INVMAX
      PARAMETER (MPSMAX=20,NPEMAX=20,INVMAX=3)
      DOUBLE PRECISION TOL
      PARAMETER (TOL=1.0D-9)
      INTEGER ITMAX
      PARAMETER (ITMAX=50)
C
C Dummy arguments
C
      DOUBLE PRECISION C,CGO,MJ
      DOUBLE PRECISION E
      INTEGER IT,L,NMATCH,NRC,NUCLEUS,NZERO,Z
      CHARACTER*3 WAY
      DOUBLE PRECISION BB(NRC),CGD(2),CGMD(2),DOVRC(NRC),DP(2,2,NRC),
     &       DQ(2,2,NRC),
     &       DRDIC(NRC),FC(2,2,NRC),GC(2,2,NRC),PIW(2,2),POW(2,2),
     &       QIW(2,2),QOW(2,2),RC(NRC),VV(NRC),WP(2,2,NRC),WQ(2,2,NRC)
C
C Local variables
C
      DOUBLE PRECISION A11,A12,A21,A22,AA11,AA12,AA21,AA22,ALPHA,BB1,
     &       BB2,BETA,
     &       BOVA,BPP,BQQ,DIFFA,DIFFB,DMUE,EMVPP,EMVQQ,W1,W2,W3,W4,W5,
     &       W6,W7
      DOUBLE PRECISION BC(0:NPEMAX),CG1,CG2,CG4,CG5,CG8,CSQR,DET,DVC,
     &       GAM(2),GPM,
     &       H24,KAP(2),PC(2,2,0:MPSMAX),PNEW(2,2),POLD(2,2),
     &       QC(2,2,0:MPSMAX),QNEW(2,2),QOLD(2,2),RPWGPM,RR,TZ,
     &       VC(0:NPEMAX)
      DOUBLE PRECISION DABS,DBLE,DEXP,DSQRT
      INTEGER I,IV,J,JCORR,K,KAP1,KAP2,M,MPS,N,NN,NSOL
      INTEGER INT,NINT
      SAVE A11,A12,A21,A22,AA11,AA12,AA21,AA22,ALPHA,BB1,BB2,BC,BETA,
     &     BOVA,BPP,BQQ,CG1,CG2,CG4,CG5,CG8,CSQR,DET,DIFFA,DIFFB,DMUE,
     &     DVC,EMVPP,EMVQQ,GAM,GPM,H24,I,IV,J,JCORR,K,KAP,KAP1,KAP2,M,
     &     MPS,N,NN,NSOL,PC,PNEW,POLD,QC,QNEW,QOLD,RPWGPM,RR,TZ,VC,W1,
     &     W2,W3,W4,W5,W6,W7
C
C MB
C     DOUBLE PRECISION CM(INVMAX,INVMAX),CMI(INVMAX,INVMAX)
C MB
C
C
      H24 = 1.0D0/24.0D0
      DVC = C
      CSQR = DVC*DVC
C
C     EXPANSION COEFFICIENTS FOR THE POTENTIAL AND B-FIELD
C MB
      IF ( NUCLEUS.EQ.0 ) THEN
         TZ = DBLE(NINT(-VV(1)*RC(1)))
         VC(0) = VV(1) - (-TZ)/RC(1)
      ELSE
         TZ = 2.0D0*DBLE(Z)
         VC(0) = VV(1)
      END IF
      DO I = 1,2
         DO J = 1,2
            DO K = 1,NPEMAX
               PC(I,J,K) = 0.0D0
               QC(I,J,K) = 0.0D0
            END DO
         END DO
      END DO
C MB
C
      BC(0) = BB(1)
C
C
C    CALCULATE G-COEFFICIENTS OF B-FIELD
C
      KAP1 = -L - 1
      KAP2 = +L
C
      CG1 = -MJ/(KAP1+0.5D0)
      CG5 = -MJ/(-KAP1+0.5D0)
      CGD(1) = CG1
      CGMD(1) = CG5
      KAP(1) = DBLE(KAP1)
C MB
      IF ( NUCLEUS.EQ.0 ) THEN
         GAM(1) = DSQRT(KAP(1)**2-(TZ/DVC)**2)
      ELSE
         GAM(1) = DABS(KAP(1))
      END IF
C MB
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
      ELSE
         CG2 = -DSQRT(1.0D0-(MJ/(KAP1+0.5D0))**2)
         CG4 = -MJ/(KAP2+0.5D0)
         CG8 = -MJ/(-KAP2+0.5D0)
         NSOL = 2
         CGD(2) = CG4
         CGO = CG2
         CGMD(2) = CG8
         KAP(2) = DBLE(KAP2)
CMBA
         IF ( NUCLEUS.EQ.0 ) THEN
            GAM(2) = DSQRT(KAP(2)**2-(TZ/DVC)**2)
         ELSE
            GAM(2) = DABS(KAP(2))
         END IF
CMBE
      END IF
C
C
      IF ( WAY.EQ.'INW' ) THEN
C
C #####################################################################
C #####################################################################
C #####################################################################
C
C             INWARD INTEGRATION
C
C
C
C
         DMUE = SQRT(-E-E*E/CSQR)
         BOVA = -DMUE/(1.0D0+E/CSQR)
C
         DO N = (NZERO-3),NZERO
C
            RR = RC(N)
C
            DO J = 1,NSOL
               I = 3 - J
               WP(J,J,N) = DEXP(-DMUE*RR)
               DP(J,J,N) = -DMUE*DRDIC(N)*WP(J,J,N)
               WQ(J,J,N) = BOVA*WP(J,J,N)
               DQ(J,J,N) = BOVA*DP(J,J,N)
C
               WP(I,J,N) = 0.0D0
               WQ(I,J,N) = 0.0D0
               DP(I,J,N) = 0.0D0
               DQ(I,J,N) = 0.0D0
            END DO
         END DO
C
C =============================================================== N ====
C
C     CALCULATE ALL NEXT POINTS BY PRE/CORR(ADAMS-MOULTON-BASHFORTH)
C
         DO NN = 1,(NZERO-3-NMATCH)
            N = NZERO - 3 - NN
C
C    EVALUATE PREDICTOR
C
            DO J = 1,NSOL
               DO I = 1,NSOL
                  PNEW(I,J) = WP(I,J,N+1)
     &                        - H24*(55.0D0*DP(I,J,N+1)-59.0D0*DP(I,J,
     &                        N+2)+37.0D0*DP(I,J,N+3)-9.0D0*DP(I,J,N+4))
                  QNEW(I,J) = WQ(I,J,N+1)
     &                        - H24*(55.0D0*DQ(I,J,N+1)-59.0D0*DQ(I,J,
     &                        N+2)+37.0D0*DQ(I,J,N+3)-9.0D0*DQ(I,J,N+4))
               END DO
            END DO
C
            EMVQQ = (E-VV(N)+CSQR)*DRDIC(N)/CSQR
            EMVPP = -(E-VV(N))*DRDIC(N)
            BQQ = BB(N)*DRDIC(N)/CSQR
            BPP = BB(N)*DRDIC(N)
C
C    EVALUATE CORRECTOR
C
            DO JCORR = 1,ITMAX
C
               DO J = 1,NSOL
                  DO I = 1,NSOL
                     POLD(I,J) = PNEW(I,J)
                     QOLD(I,J) = QNEW(I,J)
                     DP(I,J,N) = -KAP(I)*PNEW(I,J)*DOVRC(N)
     &                           + (EMVQQ+BQQ*CGMD(I))*QNEW(I,J)
                     DQ(I,J,N) = KAP(I)*QNEW(I,J)*DOVRC(N)
     &                           + (EMVPP+BPP*CGD(I))*PNEW(I,J)
     &                           + BPP*CGO*PNEW(3-I,J)
C
                     PNEW(I,J) = WP(I,J,N+1)
     &                           - H24*(9.0D0*DP(I,J,N)+19.0D0*DP(I,J,
     &                           N+1)-5.0D0*DP(I,J,N+2)+DP(I,J,N+3))
                     QNEW(I,J) = WQ(I,J,N+1)
     &                           - H24*(9.0D0*DQ(I,J,N)+19.0D0*DQ(I,J,
     &                           N+1)-5.0D0*DQ(I,J,N+2)+DQ(I,J,N+3))
                  END DO
               END DO
C
               DO J = 1,NSOL
                  DO I = 1,NSOL
                     DIFFA = POLD(I,J) - PNEW(I,J)
                     IF ( ABS(DIFFA).GT.(TOL*ABS(PNEW(I,J))) ) GOTO 20
                     IF ( ABS(DIFFA).GT.(TOL*ABS(PNEW(I,J))) ) GOTO 20
C
                     DIFFB = QOLD(I,J) - QNEW(I,J)
                     IF ( ABS(DIFFB).GT.(TOL*ABS(QNEW(I,J))) ) GOTO 20
                     IF ( ABS(DIFFB).GT.(TOL*ABS(QNEW(I,J))) ) GOTO 20
                  END DO
               END DO
               GOTO 40
C
 20         END DO
            if(t_inc%i_write>0) WRITE (1337,99001) KAP1,N,RC(N),DIFFA,
     &                     DIFFB,IT,L,INT(2*MJ),' IN'
C
C                   SORRY NOT CONVERGED IN  ITMAX  ITERATIONS
C
C
C
 40         CONTINUE
            DO J = 1,NSOL
               DO I = 1,NSOL
                  WP(I,J,N) = PNEW(I,J)
                  WQ(I,J,N) = QNEW(I,J)
                  DP(I,J,N) = -KAP(I)*PNEW(I,J)*DOVRC(N)
     &                        + (EMVQQ+BQQ*CGMD(I))*QNEW(I,J)
                  DQ(I,J,N) = KAP(I)*QNEW(I,J)*DOVRC(N)
     &                        + (EMVPP+BPP*CGD(I))*PNEW(I,J)
     &                        + BPP*CGO*PNEW(3-I,J)
               END DO
            END DO
C
         END DO
C =============================================================== N ====
C
C
C     NOW TRANSFORM TO THE PROPER WAVEFUNCTIONS
C
         DO N = NMATCH,NZERO
            DO J = 1,NSOL
               DO I = 1,NSOL
                  GC(I,J,N) = WP(I,J,N)/RC(N)
                  FC(I,J,N) = WQ(I,J,N)/(RC(N)*DVC)
               END DO
            END DO
         END DO
C
         DO J = 1,NSOL
            DO I = 1,NSOL
               PIW(I,J) = WP(I,J,NMATCH)
               QIW(I,J) = WQ(I,J,NMATCH)
            END DO
         END DO
         GOTO 99999
      END IF
C
C #####################################################################
C #####################################################################
C #####################################################################
C
C             OUTWARD INTEGRATION
C
C
C
C
C
C  DETERMINE HIGHER EXPANSION COEFFICIENTS FOR THE WAVE FUNCTIONS
C
      MPS = 20
C
      AA12 = -TZ/CSQR
      AA21 = TZ
      EMVQQ = (E-VC(0)+CSQR)/CSQR
      EMVPP = -E + VC(0)
      BQQ = BC(0)/CSQR
CMBA
      IF ( NUCLEUS.EQ.0 ) THEN
CMBE
         DO J = 1,NSOL
            I = 3 - J
            PC(J,J,0) = DSQRT(ABS(KAP(J))-GAM(J))
            QC(J,J,0) = (KAP(J)+GAM(J))*(CSQR/TZ)*PC(J,J,0)
            PC(I,J,0) = 0.0D0
            QC(I,J,0) = 0.0D0
         END DO
C
         DO J = 1,NSOL
C
            DO M = 1,MPS
               DO I = 1,NSOL
                  BB1 = (EMVQQ+BQQ*CGMD(I))*QC(I,J,M-1)
                  BB2 = (EMVPP+BC(0)*CGD(I))*PC(I,J,M-1) + BC(0)
     &                  *CGO*PC(3-I,J,M-1)
                  AA11 = GAM(J) + M + KAP(I)
                  AA22 = GAM(J) + M - KAP(I)
                  DET = AA11*AA22 - AA12*AA21
                  PC(I,J,M) = (BB1*AA22-AA12*BB2)/DET
                  QC(I,J,M) = (AA11*BB2-BB1*AA21)/DET
               END DO
            END DO
C
         END DO
C MBA
      ELSE
C EXPANSION ADAPTED FOR POTENTIALS WITH FINITE NUCLEUS
C EXPANSION OF POTENTIAL actually UP TO zeroth ORDER
C
C       DO IV=1,INVMAX
C        DO N=1,INVMAX
C         CM(N,IV)=RC(N)**(IV-1)
C        ENDDO
C       ENDDO
C
C       CALL RINVGJ(CMI,CM,INVMAX,INVMAX)
         DO IV = 1,INVMAX
            VC(IV-1) = 0.0D0
C        DO N=1,INVMAX
C         VC(IV-1)=VC(IV-1)+CMI(IV,N)*VV(N)
C        ENDDO
         END DO
         DO J = 1,NSOL
            I = 3 - J
            IF ( KAP(J).GT.0 ) THEN
C ARBITRARY STARTING VALUES
               ALPHA = 0.0D0
               BETA = 0.174D0
            ELSE
               BETA = 0.0D0
               ALPHA = 0.174D0
            END IF
            PC(J,J,0) = ALPHA
            QC(J,J,0) = BETA
            PC(I,J,0) = 0.0D0
            QC(I,J,0) = 0.0D0
         END DO
C
         W4 = BC(0)*CGO
         W2 = VC(1)/CSQR
         W5 = VC(1)
         W6 = VC(2)/CSQR
         W7 = VC(2)
         DO J = 1,NSOL
            DO I = 1,NSOL
               W1 = EMVQQ + BQQ*CGMD(I)
               W3 = -EMVPP + BC(0)*CGD(I)
               A11 = GAM(J) + KAP(I) + 1D0
               A12 = GAM(J) - KAP(I) + 1D0
               IF ( A11.NE.0 ) PC(I,J,1) = W1/A11*QC(I,J,0)
               IF ( A12.NE.0 ) QC(I,J,1)
     &              = (-W3*PC(I,J,0)+W4*PC(3-I,J,0))/A12
C
            END DO
         END DO
         DO J = 1,NSOL
            DO I = 1,NSOL
               W1 = EMVQQ + BQQ*CGMD(I)
               W3 = -EMVPP + BC(0)*CGD(I)
               A11 = GAM(J) + KAP(I) + 2D0
               A12 = GAM(J) - KAP(I) + 2D0
               IF ( A11.NE.0 ) PC(I,J,2) = (W1*QC(I,J,1)-W2*QC(I,J,0))
     &              /A11
               IF ( A12.NE.0 ) QC(I,J,2)
     &              = (-W3*PC(I,J,1)+W4*PC(3-I,J,1)+W5*PC(I,J,0))/A12
            END DO
         END DO
         DO J = 1,NSOL
            DO M = 3,MPS
               DO I = 1,NSOL
                  W1 = EMVQQ + BQQ*CGMD(I)
                  W3 = -EMVPP + BC(0)*CGD(I)
                  A21 = GAM(J) + KAP(I) + DBLE(M)
                  A22 = GAM(J) - KAP(I) + DBLE(M)
                  IF ( A21.NE.0 ) PC(I,J,M)
     &                 = (W1*QC(I,J,M-1)-W2*QC(I,J,M-2)-W6*QC(I,J,M-3))
     &                 /A21
                  IF ( A22.NE.0 ) QC(I,J,M)
     &                 = (-W3*PC(I,J,M-1)+W4*PC(3-I,J,M-1)
     &                 +W5*PC(I,J,M-2)+W7*PC(I,J,M-3))/A22
               END DO
            END DO
         END DO
      END IF
CMBE
C
C  PERFORM SUMMATION OVER WAVE FUNCTION - EXPANSION COEFFICIENTS
C  FOR THE FIRST 4 R - MESH - POINTS
C
      DO N = 1,4
         RR = RC(N)
C
         DO J = 1,NSOL
            RPWGPM = RR**GAM(J)
C
            DO I = 1,NSOL
               WP(I,J,N) = PC(I,J,0)*RPWGPM
               WQ(I,J,N) = QC(I,J,0)*RPWGPM
C      print*,WP(i,j,N), WQ(I,J,N),0,i,j
               DP(I,J,N) = PC(I,J,0)*RPWGPM*GAM(J)*DOVRC(N)
               DQ(I,J,N) = QC(I,J,0)*RPWGPM*GAM(J)*DOVRC(N)
            END DO
C
            DO M = 1,MPS
               RPWGPM = RPWGPM*RR
               GPM = GAM(J) + M
C
               DO I = 1,NSOL
                  WP(I,J,N) = WP(I,J,N) + PC(I,J,M)*RPWGPM
                  WQ(I,J,N) = WQ(I,J,N) + QC(I,J,M)*RPWGPM
C       print*,WP(i,j,N),WQ(I,J,N) ,m,i,j
                  DP(I,J,N) = DP(I,J,N) + PC(I,J,M)*RPWGPM*GPM*DOVRC(N)
                  DQ(I,J,N) = DQ(I,J,N) + QC(I,J,M)*RPWGPM*GPM*DOVRC(N)
               END DO
C
            END DO
C      if((nsol.eq.2).and.(N.gt.1))stop
         END DO
      END DO
C
C
C
C =============================================================== N ====
C     CALCULATE ALL NEXT POINTS BY PRE/CORR(ADAMS-MOULTON-BASHFORTH)
C
      DO N = 5,NMATCH
C
C    EVALUATE PREDICTOR
C
         DO J = 1,NSOL
            DO I = 1,NSOL
               PNEW(I,J) = WP(I,J,N-1)
     &                     + H24*(55.0D0*DP(I,J,N-1)-59.0D0*DP(I,J,N-2)
     &                     +37.0D0*DP(I,J,N-3)-9.0D0*DP(I,J,N-4))
               QNEW(I,J) = WQ(I,J,N-1)
     &                     + H24*(55.0D0*DQ(I,J,N-1)-59.0D0*DQ(I,J,N-2)
     &                     +37.0D0*DQ(I,J,N-3)-9.0D0*DQ(I,J,N-4))
            END DO
         END DO
C
         EMVQQ = (E-VV(N)+CSQR)*DRDIC(N)/CSQR
         EMVPP = -(E-VV(N))*DRDIC(N)
         BQQ = BB(N)*DRDIC(N)/CSQR
         BPP = BB(N)*DRDIC(N)
C
C    EVALUATE CORRECTOR
C
C
         DO JCORR = 1,ITMAX
C
            DO J = 1,NSOL
               DO I = 1,NSOL
                  POLD(I,J) = PNEW(I,J)
                  QOLD(I,J) = QNEW(I,J)
                  DP(I,J,N) = -KAP(I)*PNEW(I,J)*DOVRC(N)
     &                        + (EMVQQ+BQQ*CGMD(I))*QNEW(I,J)
                  DQ(I,J,N) = KAP(I)*QNEW(I,J)*DOVRC(N)
     &                        + (EMVPP+BPP*CGD(I))*PNEW(I,J)
     &                        + BPP*CGO*PNEW(3-I,J)
C
                  PNEW(I,J) = WP(I,J,N-1)
     &                        + H24*(9.0D0*DP(I,J,N)+19.0D0*DP(I,J,N-1)
     &                        -5.0D0*DP(I,J,N-2)+DP(I,J,N-3))
                  QNEW(I,J) = WQ(I,J,N-1)
     &                        + H24*(9.0D0*DQ(I,J,N)+19.0D0*DQ(I,J,N-1)
     &                        -5.0D0*DQ(I,J,N-2)+DQ(I,J,N-3))
               END DO
            END DO
C
            DO J = 1,NSOL
               DO I = 1,NSOL
                  DIFFA = POLD(I,J) - PNEW(I,J)
                  IF ( ABS(DIFFA).GT.(TOL*ABS(PNEW(I,J))) ) GOTO 50
                  IF ( ABS(DIFFA).GT.(TOL*ABS(PNEW(I,J))) ) GOTO 50
C
                  DIFFB = QOLD(I,J) - QNEW(I,J)
                  IF ( ABS(DIFFB).GT.(TOL*ABS(QNEW(I,J))) ) GOTO 50
                  IF ( ABS(DIFFB).GT.(TOL*ABS(QNEW(I,J))) ) GOTO 50
               END DO
            END DO
            GOTO 100
C
 50      END DO
         if(t_inc%i_write>0) WRITE (1337,99001) KAP1,N,RC(N),DIFFA,
     &                      DIFFB,IT,L,INT(2*MJ),'OUT'
C
C                   SORRY NOT CONVERGED IN  ITMAX  ITERATIONS
C
C
 100     CONTINUE
         DO J = 1,NSOL
            DO I = 1,NSOL
               WP(I,J,N) = PNEW(I,J)
               WQ(I,J,N) = QNEW(I,J)
               DP(I,J,N) = -KAP(I)*PNEW(I,J)*DOVRC(N)
     &                     + (EMVQQ+BQQ*CGMD(I))*QNEW(I,J)
               DQ(I,J,N) = KAP(I)*QNEW(I,J)*DOVRC(N)
     &                     + (EMVPP+BPP*CGD(I))*PNEW(I,J)
     &                     + BPP*CGO*PNEW(3-I,J)
            END DO
         END DO
C
      END DO
C =============================================================== N ====
C
C     NOW TRANSFORM TO THE PROPER WAVEFUNCTIONS
C
C
      DO N = 1,NMATCH
         DO J = 1,NSOL
            DO I = 1,NSOL
               GC(I,J,N) = WP(I,J,N)/RC(N)
C
C
               FC(I,J,N) = WQ(I,J,N)/(RC(N)*DVC)
            END DO
         END DO
      END DO
C
      DO J = 1,NSOL
         DO I = 1,NSOL
            POW(I,J) = WP(I,J,NMATCH)
            QOW(I,J) = WQ(I,J,NMATCH)
         END DO
      END DO
C
      RETURN
C
99001 FORMAT (' P/C NOT CONV. IN <DIRAC> ',2I4,2X,F10.7,2X,2E12.4,3I2,
     &        '/2 ',A3)
99999 CONTINUE
      END
