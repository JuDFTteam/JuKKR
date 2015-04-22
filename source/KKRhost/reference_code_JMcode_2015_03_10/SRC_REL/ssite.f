C*==ssite.f    processed by SPAG 6.05Rc at 17:31 on 29 Apr 2001
      SUBROUTINE SSITE(IWRREGWF,IWRIRRWF,NFILCBWF,CALCINT,GETIRRSOL,
     &                 SOCTL,CTL,ERYD,P,IHYPER,IPRINT,IKM1LIN,IKM2LIN,
     &                 NLQ,NKMQ,NLINQ,NT,NKM,IQAT,TSST,MSST,TSSTLIN,DZZ,
     &                 DZJ,SZZ,SZJ,OZZ,OZJ,BZZ,BZJ,QZZ,QZJ,TZZ,TZJ,VT,
     &                 BT,AT,Z,NUCLEUS,R,DRDI,R2DRDI,JWS,IMT,AMEOPC,
     &                 AMEOPO,LOPT,SOLVER,CGC,OZZS,OZJS,NLMAX,NQMAX,
     &                 LINMAX,NRMAX,NMMAX,NTMAX,NKMMAX,NKMPMAX,NLAMAX)
C   ********************************************************************
C   *                                                                  *
C   * ASSIGN QUANTUM NUMBERS AND CALL ROUTINE TO SOLVE                 *
C   * 8 COUPLED PARTIAL DIFFERENTIAL RADIAL DIRAC EQUATIONS. THE       *
C   * RESULTING WAVEFUNCTIONS ARE USED TO CALCULATE T-MATRICES IN      *
C   * THE KAPPA-MU REPRESENTATION                                      *
C   *                                                                  *
C   * + CALCULATION OF THE RADIAL INTEGRALS                            *
C   *   [ G1*G2 + F1*F2 ] R**2 DR                                      *
C   *                                                                  *
C   * FOR IHYPER <> 0 :                                                *
C   * CALCULATION OF THE HYPERFINE MATRIX ELEMENTS                     *
C   *                                                                  *
C   * RYD-UNITS USED THROUGHOUT                                        *
C   *                                                                  *
C   * NOTE: to save storage force  JG/JF  and  PR/QR  to share the     *
C   *       same storage by corresponding argument list in CALL ....   *
C   *                                                                  *
C   * 28/10/94  HE  tidy up,  P,Q used in <DIRAC> instead of g,f       *
C   * 05/10/96  HE  workspace for wavefunctions and matrices           *
C   *               is allocated dynamically !!!                       *
C   * 07/02/05  VP  few changes connected to the calculation of orbital*
C   *               polarisation                                       *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      COMPLEX*16 CI
      PARAMETER ( CI =(0.0D0, 1.0D0) )

      REAL*8 F1,E0,A0,CAUTOG
C
C conversion factor for hyperfine fields from A.U. to GAUSS
C                                      electron charge     in esu
C                                      Bohr-radius         in cm
C
      PARAMETER (F1=1.0D0,E0=1.6021892D-19*2.997930D+09,
     &           A0=0.52917706D-08,CAUTOG=E0/(A0*A0))
C
C Dummy arguments
C
      LOGICAL CALCINT,GETIRRSOL
      COMPLEX*16 ERYD,P
      INTEGER IHYPER,IPRINT,IWRIRRWF,IWRREGWF,LINMAX,NFILCBWF,NKM,
     &        NKMMAX,NLAMAX,NLMAX,NMMAX,NQMAX,NRMAX,NT,NTMAX,NUCLEUS
      INTEGER NKMPMAX
      CHARACTER*10 SOLVER
      REAL*8 AMEOPC(NKMMAX,NKMMAX,NLAMAX,3),
     &       AMEOPO(NKMMAX,NKMMAX,NLAMAX,3),AT(NRMAX,NLAMAX,3,NTMAX),
     &       BT(NRMAX,NTMAX),CTL(NTMAX,NLMAX),DRDI(NRMAX,NMMAX),
     &       R(NRMAX,NMMAX),R2DRDI(NRMAX,NMMAX),SOCTL(NTMAX,NLMAX),
     &       VT(NRMAX,NTMAX)
      REAL*8 CGC(NKMPMAX,2)
      COMPLEX*16 OZJS(LINMAX,NTMAX,2),OZZS(LINMAX,NTMAX,2)
      COMPLEX*16 BZJ(LINMAX,NTMAX),BZZ(LINMAX,NTMAX),DZJ(LINMAX,NTMAX),
     &           DZZ(LINMAX,NTMAX),MSST(NKMMAX,NKMMAX,NTMAX),
     &           OZJ(LINMAX,NTMAX),OZZ(LINMAX,NTMAX),QZJ(LINMAX,NTMAX),
     &           QZZ(LINMAX,NTMAX),SZJ(LINMAX,NTMAX),SZZ(LINMAX,NTMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX),TSSTLIN(LINMAX,NTMAX),
     &           TZJ(LINMAX,NTMAX),TZZ(LINMAX,NTMAX)
      INTEGER IKM1LIN(LINMAX),IKM2LIN(LINMAX),IMT(NTMAX),
     &        IQAT(NQMAX,NTMAX),JWS(NMMAX),LOPT(NTMAX),NKMQ(NQMAX),
     &        NLINQ(NQMAX),NLQ(NQMAX),Z(NTMAX),MUEM05
C
C Local variables
C
      COMPLEX*16 A(2,2),ARG,B1,B2,CINT(NRMAX),CRSQ,CSUM,DET,DXP(2,2),
     &           F11,F12,F21,F22,G11,G11P,G12,G12P,G21,G21P,G22,G22P,
     &           GAM(2,2),GAMINV(2,2),HL,HLB1,HLB2,JF(NRMAX,2,2),
     &           JG(NRMAX,2,2),JL,JLB1,JLB2,JLP,MAUX(NKMMAX,NKMMAX),
     &           MSST2(2,2),NL,NLB1,NLB2,NLP,NORM,PFAC,PI(2,2,NRMAX),
     &           PR(2,2,NRMAX),QI(2,2,NRMAX),QR(2,2,NRMAX),RMEHF(2,2),
     &           RMEHF1(2,2),RMEHF2(2,2),SIG(2,2),TSST2(2,2),XSST2(2,2),
     &           ZF(NRMAX,2,2),ZFJF(2,2),ZFZF(2,2),ZG(NRMAX,2,2),
     &           ZGJG(2,2),ZGZG(2,2)
      REAL*8 AP(2,2,NRMAX),AQ(2,2,NRMAX),C,CFF(2,2),CFG(2,2),CG1,CG2,
     &       CG4,CG5,CG8,CGF(2,2),CGG(2,2),CH(2,2),CSQR,CTF(2,2),
     &       CTG(2,2),DOVR(NRMAX),DROVRN(NRMAX),MJ,R1M(2,2),RKD(2,2),
     &       RNUC,SK1,SK2,TDIA1,TDIA2,TOFF
      REAL*8 COG(2,2,2),COF(2,2,2)
      COMPLEX*16 CDJLZDZ,CDNLZDZ,CJLZ,CNLZ
      DOUBLE PRECISION DBLE,DSQRT
      INTEGER I,I1,I2,I3,I5,IKM1,IKM2,IL,IM,IN,INFO,IPIV(NKMMAX),IQ,
     &        ISK1,ISK2,IT,J,JLIM,JTOP,K,K1,K2,KAP1,KAP2,KC,L,L1,LB1,
     &        LB2,LIN,N,NSOL,IMKM1,IMKM2,IS
      INTEGER IKAPMUE
      INTEGER ISIGN,NINT
      REAL*8 RNUCTAB
      LOGICAL WRONSKI
C
C*** End of declarations rewritten by SPAG
C
      DATA R1M/1.0D0,0.0D0,0.0D0,1.0D0/
      DATA RKD/1.0D0,0.0D0,0.0D0, - 1.0D0/
C     DATA RKD / 1.0D0, 0.0D0, 0.0D0, 1.0D0 /
C
      CALL CINIT(NTMAX*LINMAX,DZZ)
      CALL CINIT(NTMAX*LINMAX,DZJ)
      CALL CINIT(NTMAX*LINMAX,SZZ)
      CALL CINIT(NTMAX*LINMAX,SZJ)
      CALL CINIT(NTMAX*LINMAX,OZZ)
      CALL CINIT(NTMAX*LINMAX,OZJ)
      CALL CINIT(NTMAX*LINMAX,BZZ)
      CALL CINIT(NTMAX*LINMAX,BZJ)
      CALL CINIT(NTMAX*LINMAX,QZZ)
      CALL CINIT(NTMAX*LINMAX,QZJ)
      CALL CINIT(NTMAX*LINMAX,TZZ)
      CALL CINIT(NTMAX*LINMAX,TZJ)
C
      WRONSKI = .TRUE.
      WRONSKI = .FALSE.
C------------------------------------------------------------------------
C
      C = CTL(1,1)
      CSQR = C*C
C
C     calculate relativistic momentum
C
      P = SQRT(ERYD*(1.0D0+ERYD/CSQR))
C
      PFAC = P/(1.0D0+ERYD/CSQR)
C
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = 1,NT
C
         IF ( IPRINT.GT.0 ) WRITE (*,'(A,I3,A,10F10.4)') ' SOLVER ',IT,
     &                             SOLVER,(SOCTL(IT,IL),IL=1,NLMAX)
C
         IQ = IQAT(1,IT)
         IM = IMT(IT)
         JTOP = JWS(IM)
         IF ( NUCLEUS.NE.0 ) THEN
            RNUC = RNUCTAB(Z(IT))
            IN = 1
            DO WHILE ( R(IN,IM).LT.RNUC )
               IN = IN + 1
            END DO
C        RLIM=R(IN,IM)
            JLIM = IN
            IF ( MOD(JLIM,2).EQ.0 ) JLIM = JLIM - 1
            RNUC = R(JLIM,IM)
         END IF
         DO I = 1,JTOP
            DOVR(I) = DRDI(I,IM)/R(I,IM)
            IF ( NUCLEUS.NE.0 ) DROVRN(I) = (R(I,IM)/RNUC)**3*DRDI(I,IM)
         END DO
C
         LIN = 0
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
         DO L = 0,(NLQ(IQ)-1)
            IL = L + 1
            C = CTL(IT,IL)
            CSQR = C*C
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            KAP1 = -L - 1
            KAP2 = L
            IF ( L.EQ.0 ) KAP2 = KAP1
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
            ISK1 = ISIGN(1,KAP1)
            ISK2 = ISIGN(1,KAP2)
            SK1 = DBLE(ISK1)
            SK2 = DBLE(ISK2)
            L1 = L
            LB1 = L - ISK1
            LB2 = L - ISK2
C
            ARG = P*R(JTOP,IM)
            JL = CJLZ(L1,ARG)
            JLB1 = CJLZ(LB1,ARG)
            JLB2 = CJLZ(LB2,ARG)
            NL = CNLZ(L1,ARG)
            NLB1 = CNLZ(LB1,ARG)
            NLB2 = CNLZ(LB2,ARG)
            HL = JL + CI*NL
            HLB1 = JLB1 + CI*NLB1
            HLB2 = JLB2 + CI*NLB2
C
            IF ( SOLVER(1:7).EQ.'ABM-SOC' ) THEN
               JLP = CDJLZDZ(L1,ARG,1)*P
               NLP = CDNLZDZ(L1,ARG,1)*P
            END IF
C
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
            DO MJ = -(DBLE(L)+0.5D0), + (DBLE(L)+0.5D0),1.0D0
C
C------------------------------------------------------------------------
C NO COUPLING FOR:  ABS(MUE)= J   +  J=L+1/2 == KAP=-L-1
               IF ( ABS(MJ).GE.DBLE(L) ) THEN
                  NSOL = 1
               ELSE
                  NSOL = 2
               END IF
C------------------------------------------------------------------------
               MUEM05 = NINT(MJ-0.5D0)
               IKM1 = IKAPMUE(KAP1,MUEM05)
               IKM2 = IKAPMUE(KAP2,MUEM05)
               IMKM1 = IKAPMUE(-KAP1,MUEM05)
               IMKM2 = IKAPMUE(-KAP2,MUEM05)
C------------------------------------------------------------------------
               I5 = NRMAX*2*2
               CALL CINIT(I5,ZG)
               CALL CINIT(I5,ZF)
               CALL CINIT(I5,JG)
               CALL CINIT(I5,JF)
               CALL CINIT(I5,PR)
               CALL CINIT(I5,QR)
               CALL CINIT(I5,PI)
               CALL CINIT(I5,QI)
C
               IF ( SOLVER(1:2).EQ.'BS' ) THEN
                  CALL DIRBS(GETIRRSOL,CTL(IT,IL),IT,ERYD,L,MJ,KAP1,
     &                       KAP2,P,CG1,CG2,CG4,CG5,CG8,VT(1,IT),
     &                       BT(1,IT),Z(IT),NUCLEUS,R(1,IM),DRDI(1,IM),
     &                       DOVR,JTOP,PR,QR,PI,QI,ZG,ZF)
               ELSE IF ( SOLVER.EQ.'ABM-BI' ) THEN
                  CALL DIRABMBI(GETIRRSOL,CTL(IT,IL),IT,ERYD,L,MJ,KAP1,
     &                          KAP2,P,CG1,CG2,CG4,CG5,CG8,AMEOPC,
     &                          AMEOPO,VT(1,IT),BT(1,IT),AT,Z(IT),
     &                          NUCLEUS,R(1,IM),DRDI(1,IM),DOVR,JTOP,PR,
     &                          QR,PI,QI,ZG,ZF,AP,AQ,NTMAX,NLAMAX,
     &                          NKMMAX,NRMAX)
               ELSE IF ( SOLVER(1:6).EQ.'ABM-OP' ) THEN
                  CALL DIRABMOP(GETIRRSOL,CTL(IT,IL),IT,ERYD,L,MJ,KAP1,
     &                          KAP2,P,CG1,CG2,CG4,CG5,CG8,AMEOPO,
     &                          VT(1,IT),BT(1,IT),AT,Z(IT),NUCLEUS,
     &                          R(1,IM),DRDI(1,IM),DOVR,JTOP,PR,QR,PI,
     &                          QI,ZG,ZF,AP,AQ,LOPT(IT),NTMAX,NLAMAX,
     &                          NKMMAX,NRMAX)
               ELSE IF ( SOLVER.EQ.'ABM-SOC   ' ) THEN
                  CALL DIRABMSOC(GETIRRSOL,CTL(IT,IL),SOCTL(IT,IL),IT,
     &                           ERYD,L,MJ,KAP1,KAP2,P,CG1,CG2,CG4,CG5,
     &                           CG8,VT(1,IT),BT(1,IT),Z(IT),NUCLEUS,
     &                           R(1,IM),DRDI(1,IM),DOVR,JTOP,DXP,PR,QR,
     &                           PI,QI,ZG,ZF,NRMAX)
               ELSE IF ( SOLVER.EQ.'ABM-SOC-II' ) THEN
                  CALL DIRABMSOC2(GETIRRSOL,CTL(IT,IL),SOCTL(IT,IL),IT,
     &                            ERYD,L,MJ,KAP1,KAP2,P,CG1,CG2,CG4,CG5,
     &                            CG8,VT(1,IT),BT(1,IT),Z(IT),NUCLEUS,
     &                            R(1,IM),DRDI(1,IM),DOVR,JTOP,DXP,PR,
     &                            QR,PI,QI,ZG,ZF,NRMAX)
               ELSE
                  WRITE (6,*) 'No solver found for: ',SOLVER
                  STOP
               END IF
C
C  wavefunctions at the muffin-tin-radius
C
               N = JTOP
C
               G11 = PR(1,1,N)/R(N,IM)
               F11 = QR(1,1,N)/(R(N,IM)*C)
               G21 = PR(2,1,N)/R(N,IM)
               F21 = QR(2,1,N)/(R(N,IM)*C)
               G22 = PR(2,2,N)/R(N,IM)
               F22 = QR(2,2,N)/(R(N,IM)*C)
               G12 = PR(1,2,N)/R(N,IM)
               F12 = QR(1,2,N)/(R(N,IM)*C)
C
               IF ( SOLVER(1:7).EQ.'ABM-SOC' ) THEN
                  G11P = (DXP(1,1)/DRDI(N,IM)-G11)/R(N,IM)
                  G21P = (DXP(2,1)/DRDI(N,IM)-G21)/R(N,IM)
                  G12P = (DXP(1,2)/DRDI(N,IM)-G12)/R(N,IM)
                  G22P = (DXP(2,2)/DRDI(N,IM)-G22)/R(N,IM)
               END IF
C
C ------- the minor component for the soc-manipulated wf is meaningless
C
               IF ( SOLVER(1:7).EQ.'ABM-SOC' ) THEN
                  CALL CINIT(2*2*NRMAX,QR)
                  CALL CINIT(2*2*NRMAX,QI)
               END IF
C
C      COSD= NL * C * F11 - PFAC * SK1 * NLB1 * G11
C      SIND= JL * C * F11 - PFAC * SK1 * JLB1 * G11
C
C -------------------------------------------------------------------
C       T-SS  CONSTRUCTED USING EXPRESSIONS FROM H.E. + B.L.G. (1988)
C -------------------------------------------------------------------
C
               NL = (HL-JL)/CI
               NLB1 = (HLB1-JLB1)/CI
               NLB2 = (HLB2-JLB2)/CI
C
C
               IF ( SOLVER(1:7).EQ.'ABM-SOC' ) THEN
                  GAM(1,1) = JL*G11P - JLP*G11
                  GAM(2,1) = JL*G21P - JLP*G21
                  GAM(1,2) = JL*G12P - JLP*G12
                  GAM(2,2) = JL*G22P - JLP*G22
C
                  SIG(1,1) = NL*G11P - NLP*G11
                  SIG(2,1) = NL*G21P - NLP*G21
                  SIG(1,2) = NL*G12P - NLP*G12
                  SIG(2,2) = NL*G22P - NLP*G22
               ELSE
                  GAM(1,1) = JL*C*F11 - PFAC*SK1*JLB1*G11
                  GAM(2,1) = JL*C*F21 - PFAC*SK2*JLB2*G21
                  GAM(1,2) = JL*C*F12 - PFAC*SK1*JLB1*G12
                  GAM(2,2) = JL*C*F22 - PFAC*SK2*JLB2*G22
C
                  SIG(1,1) = NL*C*F11 - PFAC*SK1*NLB1*G11
                  SIG(2,1) = NL*C*F21 - PFAC*SK2*NLB2*G21
                  SIG(1,2) = NL*C*F12 - PFAC*SK1*NLB1*G12
                  SIG(2,2) = NL*C*F22 - PFAC*SK2*NLB2*G22
               END IF
C
               CALL ZCOPY(NSOL*NSOL,GAM,1,GAMINV,1)
               CALL ZGETRF(NSOL,NSOL,GAMINV,2,IPIV,INFO)
               CALL ZGETRI(NSOL,GAMINV,2,IPIV,MAUX,2*2,INFO)
C
               DO I2 = 1,NSOL
                  DO I1 = 1,NSOL
                     CSUM = 0.0D0
                     DO I3 = 1,NSOL
                        CSUM = CSUM + SIG(I1,I3)*GAMINV(I3,I2)
                     END DO
                     XSST2(I1,I2) = P*CSUM
                  END DO
               END DO
C
               DO I1 = 1,NSOL
                  DO I2 = 1,NSOL
                     MSST2(I1,I2) = -XSST2(I1,I2)
                  END DO
                  MSST2(I1,I1) = MSST2(I1,I1) + CI*P
               END DO
C
               CALL ZCOPY(NSOL*NSOL,MSST2,1,TSST2,1)
C
               CALL ZGETRF(NSOL,NSOL,TSST2,2,IPIV,INFO)
               CALL ZGETRI(NSOL,TSST2,2,IPIV,MAUX,2*2,INFO)
C
               IF ( IPRINT.GE.3 ) WRITE (6,99001) IT,L,MJ,
     &              ((TSST2(I1,I2),I2=1,NSOL),I1=1,NSOL)
C------------------------------------------------------------------------
C
C   COEFFICIENTS TO CALCULATE THE SPIN  MAGNETISATION
C
               CGG(1,1) = CG1
               CGG(1,2) = CG2
               CGG(2,1) = CG2
               CGG(2,2) = CG4
               CALL RINIT(4,CGF)
               CGF(1,1) = CG5
               CGF(2,2) = CG8
C
C
C   COEFFICIENTS TO CALCULATE THE SPIN  DIPOLAR MOMENT TZ
C
               TDIA1 = 2*MJ/DBLE((2*L1+1)*(2*LB1+1))
               TDIA2 = 2*MJ/DBLE((2*L1+1)*(2*LB2+1))
               TOFF = -SQRT((L1+0.5D0)**2-MJ**2)/DBLE(2*L1+1)
C
               CTG(1,1) = 0.5D0*(CG1-3.0D0*TDIA1)
               CTG(1,2) = 0.5D0*(CG2-3.0D0*TOFF)
               CTG(2,1) = 0.5D0*(CG2-3.0D0*TOFF)
               CTG(2,2) = 0.5D0*(CG4-3.0D0*TDIA2)
               CALL RINIT(4,CTF)
C
C
C   COEFFICIENTS TO CALCULATE THE ORBITAL MAGNETISATION
C
               CFG(1,1) = MJ*(KAP1+1.0D0)/(KAP1+0.5D0)
               CFG(2,2) = MJ*(KAP2+1.0D0)/(KAP2+0.5D0)
               CFG(1,2) = 0.5D0*DSQRT(1.0D0-(MJ/(KAP1+0.5D0))**2)
               CFG(2,1) = CFG(1,2)
               CALL RINIT(4,CFF)
               CFF(1,1) = MJ*(-KAP1+1.0D0)/(-KAP1+0.5D0)
               CFF(2,2) = MJ*(-KAP2+1.0D0)/(-KAP2+0.5D0)
C
C-----------------------------------------------------------------------
C   COEFFICIENTS TO CALCULATE THE ORBITAL POLARISATION
C
               CALL RINIT(4*2,COG)
               CALL RINIT(4*2,COF)
               DO IS = 1,2
                  COG(1,1,IS) = CGC(IKM1,IS)*CGC(IKM1,IS)
     &                          * DBLE(MUEM05-IS+2)
                  COF(1,1,IS) = CGC(IMKM1,IS)*CGC(IMKM1,IS)
     &                          * DBLE(MUEM05-IS+2)
               END DO
C
               IF ( NSOL.EQ.2 ) THEN
                  DO IS = 1,2
                     COG(2,2,IS) = CGC(IKM2,IS)*CGC(IKM2,IS)
     &                          * DBLE(MUEM05-IS+2)
                     COF(2,2,IS) = CGC(IMKM2,IS)*CGC(IMKM2,IS)
     &                          * DBLE(MUEM05-IS+2)
C
                     COG(1,2,IS) = CGC(IKM1,IS)*CGC(IKM2,IS)
     &                          * DBLE(MUEM05-IS+2)
                     COG(2,1,IS) = COG(1,2,IS)
                  END DO
               END IF
C-----------------------------------------------------------------------
C   ANGULAR HYPERFINE MATRIX ELEMENTS   SEE E.G.  E.M.ROSE
C        THE FACTOR  I  HAS BEEN OMITTED
C
               CH(1,1) = 4.0D0*KAP1*MJ/(4.0D0*KAP1*KAP1-1.0D0)
               CH(2,2) = 4.0D0*KAP2*MJ/(4.0D0*KAP2*KAP2-1.0D0)
               IF ( NSOL.EQ.2 ) THEN
                  CH(1,2) = DSQRT(0.25D0-(MJ/DBLE(KAP1-KAP2))**2)
                  CH(2,1) = CH(1,2)
               END IF
C-----------------------------------------------------------------------
CALCULATE RADIAL INTEGRALS  UP TO   OR RWS(JTOP=JWS)
C
C
               IF ( NSOL.EQ.1 ) THEN
C====================================================================
C NO COUPLING TO OTHER SCATTERING CHANNELS
C REGULAR PART    Z*Z    Z = (GRA,FRA)
C
                  NORM = (P*NL-JL*XSST2(1,1))/G11
C
                  DO I = 1,JTOP
                     ZG(I,1,1) = (PR(1,1,I)/R(I,IM))*NORM
                     ZF(I,1,1) = (QR(1,1,I)/R(I,IM)/C)*NORM
                     JG(I,1,1) = PI(1,1,I)/R(I,IM)
                     JF(I,1,1) = QI(1,1,I)/R(I,IM)/C
                  END DO
C
C
                  IF ( IWRREGWF.NE.0 )
     &                 WRITE (NFILCBWF,REC=IKM1+(IT-1)*NKM) IT,L,MJ,
     &                 NSOL,'REG',KAP1,IKM1,
     &                 (ZG(I,1,1),ZF(I,1,1),I=1,JTOP)
C
                  IF ( IWRIRRWF.NE.0 )
     &                 WRITE (NFILCBWF,REC=IKM1+(IT-1+NT)*NKM) IT,L,MJ,
     &                 NSOL,'IRR',KAP1,IKM1,
     &                 (JG(I,1,1),JF(I,1,1),I=1,JTOP)
C
C============================================== NO COUPLING = END ===
               ELSE
C====================================================================
C COUPLING OF TWO SCATTERING CHANNELS
C   Z(K1,K2):  INDEX 1: SPIN-ANGULAR CHARACTER
C              INDEX 2: BOUNDARY CONDITION
C
C
                  DET = G11*G22 - G12*G21
C
COEFFICIENTS TO GET:   Z(K1,K1)  Z(K2,K1)
                  B1 = P*NL - XSST2(1,1)*JL
                  B2 = -XSST2(2,1)*JL
                  A(1,1) = (G22*B1-G12*B2)/DET
                  A(2,1) = (G11*B2-G21*B1)/DET
C
COEFFICIENTS TO GET:   Z(K1,K2)  Z(K2,K2)
                  B1 = -XSST2(1,2)*JL
                  B2 = P*NL - XSST2(2,2)*JL
                  A(1,2) = (G22*B1-G12*B2)/DET
                  A(2,2) = (G11*B2-G21*B1)/DET
C
CALCULATE FUNCTIONS: Z(K1,K1), Z(K2,K1), Z(K1,K2), Z(K2,K2)
                  DO K = 1,NSOL
                     DO I = 1,JTOP
                        ZG(I,1,K) = (PR(1,1,I)*A(1,K)+PR(1,2,I)*A(2,K))
     &                              /R(I,IM)
                        ZF(I,1,K) = (QR(1,1,I)*A(1,K)+QR(1,2,I)*A(2,K))
     &                              /R(I,IM)/C
C
                        ZG(I,2,K) = (PR(2,1,I)*A(1,K)+PR(2,2,I)*A(2,K))
     &                              /R(I,IM)
                        ZF(I,2,K) = (QR(2,1,I)*A(1,K)+QR(2,2,I)*A(2,K))
     &                              /R(I,IM)/C
                     END DO
                  END DO
                  DO K = 1,NSOL
                     KC = 3 - K
                     DO I = 1,JTOP
                        JG(I,K,K) = PI(K,K,I)/R(I,IM)
                        JF(I,K,K) = QI(K,K,I)/R(I,IM)/C
                        JG(I,KC,K) = PI(KC,K,I)/R(I,IM)
                        JF(I,KC,K) = QI(KC,K,I)/R(I,IM)/C
                     END DO
                  END DO
C
C-----------------------------------------------------------------------
                  IF ( IWRREGWF.NE.0 ) THEN
C solution 1
                     WRITE (NFILCBWF,REC=IKM1+(IT-1)*NKM) IT,L,MJ,NSOL,
     &                      'REG',KAP1,IKM1,
     &                      (ZG(I,1,1),ZF(I,1,1),I=1,JTOP),KAP2,IKM2,
     &                      (ZG(I,2,1),ZF(I,2,1),I=1,JTOP)
C
C solution 2
                     WRITE (NFILCBWF,REC=IKM2+(IT-1)*NKM) IT,L,MJ,NSOL,
     &                      'REG',KAP2,IKM2,
     &                      (ZG(I,2,2),ZF(I,2,2),I=1,JTOP),KAP1,IKM1,
     &                      (ZG(I,1,2),ZF(I,1,2),I=1,JTOP)
                  END IF
C
                  IF ( IWRIRRWF.NE.0 ) THEN
C solution 1
                     WRITE (NFILCBWF,REC=IKM1+(IT-1+NT)*NKM) IT,L,MJ,
     &                      NSOL,'IRR',KAP1,IKM1,
     &                      (JG(I,1,1),JF(I,1,1),I=1,JTOP),KAP2,IKM2,
     &                      (JG(I,2,1),JF(I,2,1),I=1,JTOP)
C
C solution 2
                     WRITE (NFILCBWF,REC=IKM2+(IT-1+NT)*NKM) IT,L,MJ,
     &                      NSOL,'IRR',KAP2,IKM2,
     &                      (JG(I,2,2),JF(I,2,2),I=1,JTOP),KAP1,IKM1,
     &                      (JG(I,1,2),JF(I,1,2),I=1,JTOP)
C
                  END IF
C================================================= COUPLING = END ===
               END IF
C
C
C
CALCULATE SUM OF INTEGRALS TO BE MULTIPLIED TO   TAU(K1,K2)
               DO K1 = 1,NSOL
                  DO K2 = 1,NSOL
C
                     LIN = LIN + 1
                     TSSTLIN(LIN,IT) = TSST2(K1,K2)
                     IF ( CALCINT ) THEN
C REGULAR PART    Z*Z
C
                        CALL CINTABR(ZG(1,1,K1),ZG(1,1,K2),ZGZG,
     &                               ZF(1,1,K1),ZF(1,1,K2),ZFZF,
     &                               R2DRDI(1,IM),NSOL,NSOL,JTOP,NRMAX)
C
                        CALL SUMUPINT(DZZ(LIN,IT),F1,ZGZG,R1M,F1,ZFZF,
     &                                R1M,NSOL)
                        CALL SUMUPINT(SZZ(LIN,IT),F1,ZGZG,CGG,-F1,ZFZF,
     &                                CGF,NSOL)
                        CALL SUMUPINT(OZZ(LIN,IT),F1,ZGZG,CFG,-F1,ZFZF,
     &                                CFF,NSOL)
C ----------------------------------------------------------------------
                        DO IS = 1,2
                           CALL SUMUPINT(OZZS(LIN,IT,IS),F1,ZGZG,
     &                                   COG(1,1,IS),-F1,ZFZF,
     &                                   COF(1,1,IS),NSOL)
                        END DO
C ----------------------------------------------------------------------
                        CALL SUMUPINT(QZZ(LIN,IT),F1,ZGZG,RKD,F1,ZFZF,
     &                                RKD,NSOL)
                        CALL SUMUPINT(TZZ(LIN,IT),F1,ZGZG,CTG,-F1,ZFZF,
     &                                CTF,NSOL)
C
C       write(66,'(3I3,2e16.7)') it,nsol,lin,DZZ(LIN,IT)
C       write(66,'(4e16.7)') ((ZGZG(ii,jj),ii=1,nsol),jj=1,nsol)
C       write(66,'(4e16.7)') ((ZFZF(ii,jj),ii=1,nsol),jj=1,nsol)
C
C
C-----------------------------------------------------------------------
                        IF ( IHYPER.EQ.1 ) THEN
                           CALL CINTHFF(ZG(1,1,K1),ZF(1,1,K1),ZG(1,1,K2)
     &                                  ,ZF(1,1,K2),RMEHF,NSOL,NSOL,
     &                                  JTOP,CINT,R(1,IM),DRDI(1,IM),
     &                                  NRMAX)
C
                           IF ( NUCLEUS.NE.0 ) THEN
C calculates integrals inside nucleus but up to now only
C approximately because jlim is not the nuclear radius
C the same arguments are valid for the irregular parts below
                              CALL CINTHFF(ZG(1,1,K1),ZF(1,1,K1),
     &                           ZG(1,1,K2),ZF(1,1,K2),RMEHF1,NSOL,NSOL,
     &                           JLIM,CINT,R(1,IM),DRDI(1,IM),NRMAX)
                              CALL CINTHFF(ZG(1,1,K1),ZF(1,1,K1),
     &                           ZG(1,1,K2),ZF(1,1,K2),RMEHF2,NSOL,NSOL,
     &                           JLIM,CINT,R(1,IM),DROVRN,NRMAX)
                              DO I = 1,NSOL
                                 DO J = 1,NSOL
                                    RMEHF(I,J) = RMEHF(I,J)
     &                                 - RMEHF1(I,J) + RMEHF2(I,J)
                                 END DO
                              END DO
                           END IF
C                !end of nucleus.eq.0
                           CALL SUMUPINT(BZZ(LIN,IT),CAUTOG,RMEHF,CH,
     &                        0.0D0,RMEHF,CH,NSOL)
                        END IF
C-----------------------------------------------------------------------
C
C IRREGULAR PART    Z*J
C THE  ENDING  A (B)  STANDS FOR THE DOMINATING (DIMINATED)
C SET OF SPIN-ANGULAR-CHAR:  I.E.  J==J(A,A)  FOR R>RMT
C
                        IF ( K1.EQ.K2 ) THEN
C
                           CALL CINTABR(ZG(1,1,K1),JG(1,1,K1),ZGJG,
     &                                  ZF(1,1,K1),JF(1,1,K1),ZFJF,
     &                                  R2DRDI(1,IM),NSOL,NSOL,JTOP,
     &                                  NRMAX)
C
                           CALL SUMUPINT(DZJ(LIN,IT),F1,ZGJG,R1M,F1,
     &                        ZFJF,R1M,NSOL)
                           CALL SUMUPINT(SZJ(LIN,IT),F1,ZGJG,CGG,-F1,
     &                        ZFJF,CGF,NSOL)
                           CALL SUMUPINT(OZJ(LIN,IT),F1,ZGJG,CFG,-F1,
     &                        ZFJF,CFF,NSOL)
C ----------------------------------------------------------------------
                           DO IS = 1,2
                              CALL SUMUPINT(OZJS(LIN,IT,IS),F1,ZGJG,
     &                                      COG(1,1,IS),-F1,ZFJF,
     &                                      COF(1,1,IS),NSOL)
                           END DO
C ----------------------------------------------------------------------
                           CALL SUMUPINT(QZJ(LIN,IT),F1,ZGJG,RKD,F1,
     &                        ZFJF,RKD,NSOL)
                           CALL SUMUPINT(TZJ(LIN,IT),F1,ZGJG,CTG,-F1,
     &                        ZFJF,CTF,NSOL)
C
C-----------------------------------------------------------------------
                           IF ( IHYPER.EQ.1 ) THEN
                              CALL CINTHFF(ZG(1,1,K1),ZF(1,1,K1),
     &                           JG(1,1,K1),JF(1,1,K1),RMEHF,NSOL,NSOL,
     &                           JTOP,CINT,R(1,IM),DRDI(1,IM),NRMAX)
                              IF ( NUCLEUS.NE.0 ) THEN
C calculates integrals inside nucleus but up to now only
C approximately because jlim is not the nuclear radius
C the same arguments are valid for the irregular parts below
                                 CALL CINTHFF(ZG(1,1,K1),ZF(1,1,K1),
     &                              JG(1,1,K1),JF(1,1,K1),RMEHF1,NSOL,
     &                              NSOL,JLIM,CINT,R(1,IM),DRDI(1,IM),
     &                              NRMAX)
                                 CALL CINTHFF(ZG(1,1,K1),ZF(1,1,K1),
     &                              JG(1,1,K1),JF(1,1,K1),RMEHF2,NSOL,
     &                              NSOL,JLIM,CINT,R(1,IM),DROVRN,NRMAX)
                                 DO I = 1,NSOL
                                    DO J = 1,NSOL
                                       RMEHF(I,J) = RMEHF(I,J)
     &                                    - RMEHF1(I,J) + RMEHF2(I,J)
                                    END DO
                                 END DO
                              END IF
C                !end of nucleus.eq.0
C
                              CALL SUMUPINT(BZJ(LIN,IT),CAUTOG,RMEHF,CH,
     &                           0.0D0,RMEHF,CH,NSOL)
                           END IF
                        END IF
C-----------------------------------------------------------------------
                     END IF
C           ! OF IF (.CALCINT.)
                  END DO
               END DO
C
Check WRONSKI-relationship
C
               IF ( WRONSKI ) THEN
                  DO I = 1,JTOP,40
                     CRSQ = C*R(I,IM)**2
                     WRITE (6,99002) IT,L,NINT(2*MJ),I,R(I,IM),
     &                               1.0D0 - (ZF(I,1,1)*JG(I,1,1)
     &                               -ZG(I,1,1)*JF(I,1,1)+ZF(I,2,1)
     &                               *JG(I,2,1)-ZG(I,2,1)*JF(I,2,1))
     &                               *CRSQ
                     IF ( NSOL.EQ.2 ) THEN
                        WRITE (6,99003) 1.0D0 - (ZF(I,1,2)*JG(I,1,2)-ZG(
     &                                  I,1,2)*JF(I,1,2)+ZF(I,2,2)
     &                                  *JG(I,2,2)-ZG(I,2,2)*JF(I,2,2))
     &                                  *CRSQ,
     &                                  1.0D0 - (ZF(I,1,2)*JG(I,1,1)
     &                                  -ZG(I,1,2)*JF(I,1,1)+ZF(I,2,2)
     &                                  *JG(I,2,1)-ZG(I,2,2)*JF(I,2,1))
     &                                  *CRSQ,
     &                                  1.0D0 - (ZF(I,1,1)*JG(I,1,2)
     &                                  -ZG(I,1,1)*JF(I,1,2)+ZF(I,2,1)
     &                                  *JG(I,2,2)-ZG(I,2,1)*JF(I,2,2))
     &                                  *CRSQ
                     ELSE
                        WRITE (6,*)
                     END IF
                  END DO
               END IF
C
            END DO
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
         END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
         CALL CINIT(NKMMAX*NKMMAX,TSST(1,1,IT))
C
         DO LIN = 1,NLINQ(IQ)
            I1 = IKM1LIN(LIN)
            I2 = IKM2LIN(LIN)
            TSST(I1,I2,IT) = TSSTLIN(LIN,IT)
         END DO
C
         DO J = 1,NKMQ(IQ)
            CALL ZCOPY(NKMQ(IQ),TSST(1,J,IT),1,MSST(1,J,IT),1)
         END DO
C
         CALL ZGETRF(NKMQ(IQ),NKMQ(IQ),MSST(1,1,IT),NKMMAX,IPIV,INFO)
         CALL ZGETRI(NKMQ(IQ),MSST(1,1,IT),NKMMAX,IPIV,MAUX,
     &               NKMMAX*NKMMAX,INFO)
C
         IF ( IPRINT.GE.4 ) THEN
            DO LIN = 1,NLINQ(IQ)
               I1 = IKM1LIN(LIN)
               I2 = IKM2LIN(LIN)
               WRITE (6,99004) IT,LIN
               WRITE (6,99004) IT,I1,I2,' DZZ ',DZZ(LIN,IT),DZJ(LIN,IT)
               WRITE (6,99004) IT,I1,I2,' SZZ ',SZZ(LIN,IT),SZJ(LIN,IT)
               WRITE (6,99004) IT,I1,I2,' OZZ ',OZZ(LIN,IT),OZJ(LIN,IT)
            END DO
         END IF
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
99001 FORMAT ('  t-ss(TLM)',2I3,F5.1,2E14.5,2X,2E14.5,:,/,22X,2E14.5,2X,
     &        2E14.5)
99002 FORMAT (' IT=',I2,' L=',I2,' MJ=',I2,'/2',I7,F10.6,1(2X,2F9.6),$)
99003 FORMAT (3(2X,2F9.6))
99004 FORMAT (' IT=',I2,2I3,A,2X,2E14.5,2X,2E14.5)
      END
C*==readwfun.f    processed by SPAG 6.05Rc at 17:31 on 29 Apr 2001
      SUBROUTINE READWFUN(NFIL,IT,L,MJ,NSOL,SREG,SIRR,IKM1,KAP1,IKM2,
     &                    KAP2,NT,NKM,ZG,ZF,JG,JF,JTOP,NRMAX)
C   ********************************************************************
C   *                                                                  *
C   *  reread the wave functions written by  <SSITE>  or  <CORE>       *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER IKM1,IKM2,IT,JTOP,KAP1,KAP2,L,NFIL,NKM,NRMAX,NSOL,NT
      REAL*8 MJ
      CHARACTER*3 SIRR,SREG
      COMPLEX*16 JF(NRMAX,2,2),JG(NRMAX,2,2),ZF(NRMAX,2,2),ZG(NRMAX,2,2)
C
C Local variables
C
      INTEGER I,IFLAG,IKMIN(2),ITP,K,KAPIN(2),LP,NSOLP
      REAL*8 MJP
      CHARACTER*3 STRP
C
C*** End of declarations rewritten by SPAG
C
      IFLAG = 0
C-----------------------------------------------------------------------
C                                    REGULAR wave function -- solution 1
      IF ( SREG.EQ.'REG' .OR. SREG.EQ.'COR' ) THEN
         READ (NFIL,REC=IKM1+(IT-1)*NKM) ITP,LP,MJP,NSOLP,STRP,
     &         (KAPIN(K),IKMIN(K),(ZG(I,K,1),ZF(I,K,1),I=1,JTOP),K=1,
     &         NSOL)
         IF ( ITP.NE.IT .OR. LP.NE.L .OR. ABS(MJP-MJ).GT.0.001D0 .OR. 
     &        NSOLP.NE.NSOL .OR. STRP.NE.'REG' ) IFLAG = IFLAG + 1
         IF ( KAP1.NE.KAPIN(1) ) IFLAG = IFLAG + 1
         IF ( IKM1.NE.IKMIN(1) ) IFLAG = IFLAG + 1
         IF ( NSOL.GT.1 ) THEN
            IF ( KAP2.NE.KAPIN(2) ) IFLAG = IFLAG + 1
            IF ( IKM2.NE.IKMIN(2) ) IFLAG = IFLAG + 1
         END IF
      END IF
C
C-----------------------------------------------------------------------
C                                  IRREGULAR wave function -- solution 1
      IF ( SIRR.EQ.'IRR' ) THEN
         READ (NFIL,REC=IKM1+(IT-1+NT)*NKM) ITP,LP,MJP,NSOLP,STRP,
     &         (KAPIN(K),IKMIN(K),(JG(I,K,1),JF(I,K,1),I=1,JTOP),K=1,
     &         NSOL)
         IF ( ITP.NE.IT .OR. LP.NE.L .OR. ABS(MJP-MJ).GT.0.001D0 .OR. 
     &        NSOLP.NE.NSOL .OR. STRP.NE.'IRR' ) IFLAG = IFLAG + 1
         IF ( KAP1.NE.KAPIN(1) ) IFLAG = IFLAG + 1
         IF ( IKM1.NE.IKMIN(1) ) IFLAG = IFLAG + 1
         IF ( NSOL.GT.1 ) THEN
            IF ( KAP2.NE.KAPIN(2) ) IFLAG = IFLAG + 1
            IF ( IKM2.NE.IKMIN(2) ) IFLAG = IFLAG + 1
         END IF
      END IF
C
      IF ( NSOL.EQ.2 ) THEN
C-----------------------------------------------------------------------
C                                    REGULAR wave function -- solution 2
         IF ( SREG.EQ.'REG' .OR. SREG.EQ.'COR' ) THEN
            READ (NFIL,REC=IKM2+(IT-1)*NKM) ITP,LP,MJP,NSOLP,STRP,
     &            (KAPIN(K),IKMIN(K),(ZG(I,K,2),ZF(I,K,2),I=1,JTOP),K=2,
     &            1,-1)
            IF ( ITP.NE.IT .OR. LP.NE.L .OR. ABS(MJP-MJ).GT.0.001D0 .OR. 
     &           NSOLP.NE.NSOL .OR. STRP.NE.'REG' ) IFLAG = IFLAG + 1
         END IF
C
C-----------------------------------------------------------------------
C                                  IRREGULAR wave function -- solution 2
         IF ( SIRR.EQ.'IRR' ) THEN
            READ (NFIL,REC=IKM2+(IT-1+NT)*NKM) ITP,LP,MJP,NSOLP,STRP,
     &            (KAPIN(K),IKMIN(K),(JG(I,K,2),JF(I,K,2),I=1,JTOP),K=2,
     &            1,-1)
            IF ( ITP.NE.IT .OR. LP.NE.L .OR. ABS(MJP-MJ).GT.0.001D0 .OR. 
     &           NSOLP.NE.NSOL .OR. STRP.NE.'IRR' ) IFLAG = IFLAG + 1
         END IF
C
      END IF
C
C-----------------------------------------------------------------------
C
      IF ( IFLAG.GT.0 ) THEN
         WRITE (*,99001) IFLAG,IT,L,MJ,SREG,SIRR
         STOP ' in <READWFUN>'
      END IF
99001 FORMAT (//,1x,79('*'),/,10X,'error reading the wave functions',
     &        ' IFLAG = 1',/,10X,'for  IT=',I2,' L=',I2,' MJ=',F4.1,
     &        ' KEYS=',A,A)
      END
