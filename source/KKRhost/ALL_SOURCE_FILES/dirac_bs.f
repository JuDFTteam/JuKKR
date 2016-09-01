      SUBROUTINE DIRBS(GETIRRSOL,C,E,L,MJ,KAP1,KAP2,PIS,CG1,CG2,CG4,
     &                 CG5,CG8,V,B,Z,NUCLEUS,R,DRDI,DOVR,NMESH,PR,QR,PI,
     &                 QI,DP,DQ)
C   ********************************************************************
C   *                                                                  *
C   *   ROUTINE TO SOLVE THE SPIN-POLARISED RADIAL DIRAC EQUATIONS     *
C   *                                                                  *
C   *   the outward integration is started by a power expansion        *
C   *   and the inward integration is started analytically             *
C   *   the integration itself is done by the BURLISCH-STOER method    *
C   *   see: numerical recipes chapter 15.4                            *
C   *                                                                  *
C   *   returns the wave functions up to the mesh point NMESH          *
C   *   PR,QR and PI,QI  with   P=r*g and Q=r*c*f                      *
C   *   and    R/I standing for regular/irregular solution             *
C   *                                                                  *
C   *   bug fixed 93/11/24                                             *
C   *  31/10/94  HE  arg. list changed - return P,Q instead of g,f     *
C   *  06/12/94  HE  CM real                                           *
C   *  29/04/95  MB  Adopted for finite nucleus                        *
C   ********************************************************************
C

      IMPLICIT NONE
      INCLUDE 'sprkkr_rmesh.dim'
C
C PARAMETER definitions
C
      INTEGER MPSMAX,NPEMAX,NABM
      PARAMETER (MPSMAX=40,NPEMAX=4,NABM=5)
      COMPLEX*16 C0
      PARAMETER (C0=(0.0D0,0.0D0))
      REAL*8 EPSBS
      PARAMETER (EPSBS=2.0D-7)
C
C COMMON variables
C
      REAL*8 CGD(2),CGMD(2),CGO,CSQR,KAP(2)
      COMPLEX*16 EBS
      INTEGER NRADBS,NSOLBS
      COMMON /COMMBS/ EBS,CSQR,CGD,CGMD,CGO,KAP,NSOLBS,NRADBS
C
C Dummy arguments
C
      REAL*8 C,CG1,CG2,CG4,CG5,CG8,MJ
      COMPLEX*16 E
      LOGICAL GETIRRSOL
      INTEGER KAP1,KAP2,L,NMESH,NUCLEUS,Z
      COMPLEX*16 PIS
      REAL*8 B(NRMAX),DOVR(NRMAX),DRDI(NRMAX),R(NRMAX),V(NRMAX)
      COMPLEX*16 DP(2,2,NRMAX),DQ(2,2,NRMAX),PI(2,2,NRMAX),PR(2,2,NRMAX)
      COMPLEX*16 QI(2,2,NRMAX),QR(2,2,NRMAX)
C
C Local variables
C
      COMPLEX*16 A11,A12,A21,A22,AA11,AA12,AA21,AA22,ALPHA,BB1,BB2,BETA
      COMPLEX*16 BQQ,CFAC,EMVPP,EMVQQ,W1,W2,W3,W4,W5,W6,W7
      REAL*8 BC(0:NPEMAX),CM(NPEMAX,NPEMAX),CMI(NPEMAX,NPEMAX),DIX
      REAL*8 GAM(2),GPM,HBS,RPWGPM,RR,SK(2),SK1,SK2,TZ,VC(0:NPEMAX),X
      COMPLEX*16 CJLZ
C      
      COMPLEX*16 DETD,DY(NCFMAX),EFAC,FY(NCFMAX),PC(2,2,-NPEMAX:MPSMAX)
      COMPLEX*16 QC(2,2,-NPEMAX:MPSMAX),ZZ
C
      INTEGER I,IP,ISK1,ISK2,IV,J,K,LB(2),LB1,LB2,M,MPS,N,NFY,NPE,NSOL
      INTEGER ISIGN
C
      SAVE A11,A12,A21,A22,AA11,AA12,AA21,AA22,ALPHA,BB1,BB2,BC,BETA,
     &     BQQ,CFAC,CM,CMI,DETD,DIX,DY,EFAC,EMVPP,EMVQQ,FY,GAM,GPM,HBS,
     &     I,IP,ISK1,ISK2,IV,J,K,LB,LB1,LB2,M,MPS,N,NFY,NPE,NSOL,PC,QC,
     &     RPWGPM,RR,SK,SK1,SK2,TZ,VC,W1,W2,W3,W4,W5,W6,W7,X,ZZ
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
            IF ( NUCLEUS.EQ.0 ) THEN
               VC(IV-1) = VC(IV-1) + CMI(IV,N)*(V(N)+TZ/R(N))
            ELSE
               VC(IV-1) = VC(IV-1) + CMI(IV,N)*V(N)
            END IF
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
C MB
      IF ( NUCLEUS.EQ.0 ) THEN
         GAM(1) = DSQRT(KAP(1)**2-(TZ/C)**2)
      ELSE
         GAM(1) = DABS(KAP(1))
      END IF
C MB
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
C MB
         IF ( NUCLEUS.EQ.0 ) THEN
            GAM(2) = DSQRT(KAP(2)**2-(TZ/C)**2)
         ELSE
            GAM(2) = DABS(KAP(2))
         END IF
C MB
         LB(2) = LB2
         SK(2) = SK2
      END IF
C
      NSOLBS = NSOL
      EBS = E
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
      IF ( TZ.GE.2 ) THEN
C
         DO J = 1,NSOL
            I = 3 - J
            PC(J,J,0) = DSQRT(ABS(KAP(J))-GAM(J))
            QC(J,J,0) = (KAP(J)+GAM(J))*(CSQR/TZ)*PC(J,J,0)
            PC(I,J,0) = C0
            QC(I,J,0) = C0
         END DO
C
C  determine higher expansion coefficients for the wave functions
C
         MPS = 40
C
         AA12 = -TZ/CSQR
         AA21 = TZ
         EMVQQ = (E-VC(0)+CSQR)/CSQR
         EMVPP = -E + VC(0)
         BQQ = BC(0)/CSQR
CMBA
         IF ( NUCLEUS.EQ.0 ) THEN
CMBE
C
            DO J = 1,NSOL
C
               DO M = 1,MPS
                  DO I = 1,NSOL
                     K = 3 - I
                     BB1 = (EMVQQ+BQQ*CGMD(I))*QC(I,J,M-1)
                     BB2 = (EMVPP+BC(0)*CGD(I))*PC(I,J,M-1) + BC(0)
     &                     *CGO*PC(K,J,M-1)
                     DO IP = 1,NPE - 1
                        BB1 = BB1 + (-VC(IP)+BC(IP)*CGMD(I))
     &                        *QC(I,J,M-1-IP)/CSQR
                        BB2 = BB2 + (+VC(IP)+BC(IP)*CGD(I))
     &                        *PC(I,J,M-1-IP) + (+BC(IP)*CGO)
     &                        *PC(K,J,M-1-IP)
                     END DO
C
                     AA11 = GAM(J) + M + KAP(I)
                     AA22 = GAM(J) + M - KAP(I)
                     DETD = AA11*AA22 - AA12*AA21
                     PC(I,J,M) = (BB1*AA22-AA12*BB2)/DETD
                     QC(I,J,M) = (AA11*BB2-BB1*AA21)/DETD
                  END DO
               END DO
C
            END DO
C MBA
         ELSE
C EXPANSION ADAPTED FOR POTENTIALS WITH FINITE NUCLEUS
C EXPANSION OF POTENTIAL UP TO SECOND ORDER: V_O+V_1*R+V_2*R*R
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
     &                 = (-W3*PC(I,J,0)+W4*PC(3-I,J,0))/A12
               END DO
            END DO
            DO J = 1,NSOL
               DO I = 1,NSOL
                  W1 = EMVQQ + BQQ*CGMD(I)
                  W3 = -EMVPP + BC(0)*CGD(I)
                  A11 = GAM(J) + KAP(I) + 2D0
                  A12 = GAM(J) - KAP(I) + 2D0
                  IF ( A11.NE.0 ) PC(I,J,2)
     &                 = (W1*QC(I,J,1)-W2*QC(I,J,0))/A11
                  IF ( A12.NE.0 ) QC(I,J,2)
     &                 = (-W3*PC(I,J,1)+W4*PC(3-I,J,1)+W5*PC(I,J,0))/A12
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
     &                    = (W1*QC(I,J,M-1)-W2*QC(I,J,M-2)
     &                    -W6*QC(I,J,M-3))/A21
                     IF ( A22.NE.0 ) QC(I,J,M)
     &                    = (-W3*PC(I,J,M-1)+W4*PC(3-I,J,M-1)
     &                    +W5*PC(I,J,M-2)+W7*PC(I,J,M-3))/A22
                  END DO
               END DO
            END DO
         END IF
CMBE
C
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
C
C ======================================================================
C                                                     == EMPTY SPHERE ==
      ELSE
C                     assume constant pot: V=V(1)   ignore coupling: B=0
C
         DO N = 1,NABM
            ZZ = CDSQRT(E-V(1))*R(N)
            EFAC = (ZZ/R(N))*C/(E+CSQR)
C
            DO J = 1,NSOL
               I = 3 - J
               PR(J,J,N) = CJLZ(L,ZZ)*R(N)
               QR(J,J,N) = EFAC*SK(J)*CJLZ(LB(J),ZZ)*R(N)*C
               DP(J,J,N) = (DBLE(L+1)*CJLZ(L,ZZ)-ZZ*CJLZ(L+1,ZZ))
     &                     *DRDI(N)
               M = LB(J)
               DQ(J,J,N) = EFAC*SK(J)
     &                     *(DBLE(M+1)*CJLZ(M,ZZ)-ZZ*CJLZ(M+1,ZZ))
     &                     *DRDI(N)*C
C
               PR(I,J,N) = C0
               QR(I,J,N) = C0
               DP(I,J,N) = C0
               DQ(I,J,N) = C0
            END DO
         END DO
C
      END IF
C
C =============================================================== n ====
C     DO 400 J=1,NSOL
C
      NFY = 0
      DO J = 1,NSOL
         DO I = 1,NSOL
            FY(NFY+1) = PR(I,J,1)
            FY(NFY+2) = QR(I,J,1)
            DY(NFY+1) = DP(I,J,1)
            DY(NFY+2) = DQ(I,J,1)
            NFY = NFY + 2
         END DO
      END DO
      X = 1.0D0
      DIX = 1.0D0

      DO N = 2,NMESH
C
         NRADBS = N
C
         CALL DIRBSSTP(FY,DY,NFY,X,DIX,EPSBS,FY,B,V,R,DRDI,NMESH)
C
C
         NFY = 0
         DO J = 1,NSOL
            DO I = 1,NSOL
               PR(I,J,N) = FY(NFY+1)
               QR(I,J,N) = FY(NFY+2)
               NFY = NFY + 2
            END DO
         END DO
C
      END DO
C
C
      IF ( .NOT.GETIRRSOL ) RETURN
C
C
C =============================================================== n ====
C
C ######################################################################
C                     irregular solution
C ######################################################################
C
C         calculate the initial values of the wavefunction
C                     at the sphere boundary
C
      N = NMESH
C
      ZZ = PIS*R(N)
C
      DO J = 1,NSOL
         PI(J,J,N) = CJLZ(L,ZZ)*R(N)
         QI(J,J,N) = CFAC*SK(J)*CJLZ(LB(J),ZZ)*R(N)*C
         DP(J,J,N) = (DBLE(L+1)*CJLZ(L,ZZ)-ZZ*CJLZ(L+1,ZZ))
     &               *DRDI(N)
C
         M = LB(J)
         DQ(J,J,N) = CFAC*SK(J)
     &               *(DBLE(M+1)*CJLZ(M,ZZ)-ZZ*CJLZ(M+1,ZZ))
     &               *DRDI(N)*C
C
         I = 3 - J
         PI(I,J,N) = C0
         QI(I,J,N) = C0
         DP(I,J,N) = C0
         DQ(I,J,N) = C0
      END DO
C
C =============================================================== n ====
      HBS = -1.0D0
C
      NFY = 0
      DO J = 1,NSOL
         DO I = 1,NSOL
            FY(NFY+1) = PI(I,J,NMESH)
            FY(NFY+2) = QI(I,J,NMESH)
            DY(NFY+1) = DP(I,J,NMESH)
            DY(NFY+2) = DQ(I,J,NMESH)
            NFY = NFY + 2
         END DO
      END DO
      X = DBLE(NMESH)
C
      NRADBS = NMESH
      CALL DIRBSRAD(X,FY,DY,DRDI,B,V,R,NMESH)
C
      DO N = NMESH - 1,1, - 1
         NRADBS = N
C
         CALL DIRBSSTP(FY,DY,NFY,X,HBS,EPSBS,FY,B,V,R,DRDI,NMESH)
C
         NFY = 0
         DO J = 1,NSOL
            DO I = 1,NSOL
               PI(I,J,N) = FY(NFY+1)
               QI(I,J,N) = FY(NFY+2)
C              DP(I,J,N) = DY(NFY+1)
C              DQ(I,J,N) = DY(NFY+2)
               NFY = NFY + 2
            END DO
         END DO
C
      END DO
C
C =============================================================== n ====
      END 
