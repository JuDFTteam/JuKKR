      SUBROUTINE DIRABMBI(GETIRRSOL,C,IT,E,L,MJ,KAP1,KAP2,PIS,CG1,CG2,
     &                    CG4,CG5,CG8,AMEBI1,AMEBI2,V,B,AT,Z,NUCLEUS,R,
     &                    DRDI,DOVR,NMESH,PR,QR,PI,QI,DP,DQ,AP,AQ,NTMAX,
     &                    NLAMAX,NKMMAX,NRMAX)
C   ********************************************************************
C   *                                                                  *
C   *   ROUTINE TO SOLVE THE SPIN-POLARISED RADIAL DIRAC EQUATIONS     *
C   *                                                                  *
C   *            including a vector potential A(L,M)                   *
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
C   *  26/01/95  HE                                                    *
C   ********************************************************************
      IMPLICIT COMPLEX*16(A-H,O-Z)
C
C
C PARAMETER definitions
C
      INTEGER NLMAXLOC,NTMAXLOC,NRMAXLOC,MPSMAX,NPEMAX,NABM
      PARAMETER (NLMAXLOC=4,NTMAXLOC=2,NRMAXLOC=750,MPSMAX=40,NPEMAX=4,
     &           NABM=4)
      COMPLEX*16 CZ
      PARAMETER (CZ=(0.0D0,0.0D0))
      REAL*8 TOL
      PARAMETER (TOL=1.0D-6)
      INTEGER ITMAX
      PARAMETER (ITMAX=50)
C
C Dummy arguments
C
      REAL*8 C,CG1,CG2,CG4,CG5,CG8,MJ
      COMPLEX*16 E,PIS
      LOGICAL GETIRRSOL
      INTEGER IT,KAP1,KAP2,L,NKMMAX,NLAMAX,NMESH,NRMAX,NTMAX,NUCLEUS,Z
      REAL*8 AMEBI1(NKMMAX,NKMMAX,NLAMAX,-1:+1),
     &       AMEBI2(NKMMAX,NKMMAX,NLAMAX,-1:+1),AP(2,2,NRMAXLOC),
     &       AQ(2,2,NRMAXLOC),AT(NRMAXLOC,NLAMAX,-1:+1,NTMAXLOC),
     &       B(NRMAX),DOVR(NRMAX),DRDI(NRMAX),R(NRMAX),V(NRMAX)
      COMPLEX*16 DP(2,2,NRMAX),DQ(2,2,NRMAX),PI(2,2,NRMAX),PR(2,2,NRMAX)
     &           ,QI(2,2,NRMAX),QR(2,2,NRMAX)
C
C Local variables
C
      COMPLEX*16 A11,A12,A21,A22,AA11,AA12,AA21,AA22,ALPHA,ARG,BB1,BB2,
     &           BETA,BPP,BQQ,CFAC,CGO,D14,DETD,DH,DIFFA,DIFFB,EFAC,
     &           EMVPP,EMVQQ,MP1(2,2),MP2(2,2),MP3(2,2),MP4(2,2),
     &           MQ1(2,2),MQ2(2,2),MQ3(2,2),MQ4(2,2),P1(2,2),P2(2,2),
     &           P3(2,2),P4(2,2),PC(2,2,-NPEMAX:MPSMAX),PNEW(2,2),
     &           POLD(2,2),Q1(2,2),Q2(2,2),Q3(2,2),Q4(2,2),
     &           QC(2,2,-NPEMAX:MPSMAX),QNEW(2,2),QOLD(2,2),W1,W2,W3,W4,
     &           W5,W6,W7,ZZ
      REAL*8 ACORR(0:NABM-1),ACORR0(0:NABM-1),APRED(NABM),APRED0(NABM),
     &       ASTEP,B14,BC(0:NPEMAX),BH,BHLP(NABM+4),CGD(2),CGMD(2),
     &       CM(NPEMAX,NPEMAX),CMI(NPEMAX,NPEMAX),CSQR,DHLP(NABM+4),
     &       GAM(2),GPM,HLP(NABM+4),HLP1,KAP(2),R14,RH,RHLP(NABM+4),
     &       RPWGPM,RR,SK(2),SK1,SK2,SO2,SO6,SRK,TZ,V14,VC(0:NPEMAX),VH,
     &       VHLP(NABM+4),WR,X14,XH
      COMPLEX*16 CJLZ
      DOUBLE PRECISION DABS,DBLE,DSQRT
      INTEGER I,IC,IKM(2),IKMI,IKMJ,ILA,IP,IRK,ISK1,ISK2,IV,J,JCORR,K,
     &        LB(2),LB1,LB2,M,MPS,N,NACORR,NDIV,NHLP,NM,NPE,NSOL,NTOP
      INTEGER IKAPMUE
      INTEGER INT,ISIGN,NINT
      REAL*8 YLAG
C
C     DATA APRED0 / 1901.0D0, -2774.0D0, 2616.0D0, -1274.0D0, 251.0D0 /
C     DATA ACORR0 /  251.0D0,  +646.0D0, -264.0D0,  +106.0D0, -19.0D0 /
C     DATA ASTEP  /  720.0D0 /
      DATA APRED0/55.0D0, - 59.0D0, + 37.0D0, - 9.0D0/
      DATA ACORR0/9.0D0, + 19.0D0, - 5.0D0, + 1.0D0/
      DATA ASTEP/24.0D0/
C
C#######################################################################

      stop ' < DIRABMBI > : Not implemented. Set SOLVER=BS in inputcard'

      END
