      module ShapeIntegrationHelpers_mod
      implicit none

      contains

C-----------------------------------------------------------------------
C>    THIS ROUTINE  ACCOMPLISHES THE  FI-INTEGRATION  OF REAL  SPHERICAL
C>    HARMONICS BY THE REPEATED SIMPSON'S METHOD , OR ANALYTICALLY ACCOR
C>    DING TO THE VALUE OF ITYPE. THE OBTAINED RESULTS HAVE TO BE MULTI-
C>    PLIED BY THE APPROPRIATE EXPANSION COEFFICIENTS.
C-----------------------------------------------------------------------
      SUBROUTINE PINTG(X1,X2,DLT,S,LMAX,ISI,ARG,FD,ITYPE)
      USE shape_constants_mod, ONLY: LMAXD1, NDIM
      implicit none

C
C     .. PARAMETER STATEMENTS ..
C
C     include 'inc.geometry'
c      INTEGER LMAXD,NDIM
c      PARAMETER (LMAXD=25,NDIM=1000)

C
C     .. SCALAR ARGUMENTS ..
C
      REAL*8 X1,X2,DLT,ARG,FD
      INTEGER   LMAX,ISI,ITYPE
C
C     .. ARRAY ARGUMENTS ..
C
      REAL*8 S(-LMAXD1:LMAXD1,0:LMAXD1)
C
C     .. LOCAL SCALARS ..
C
      INTEGER   I,M,N,K
      REAL*8    X,THETA,W
C
C     .. LOCAL ARRAYS ..
C
      REAL*8    XX(NDIM),WW(NDIM)
C
C     .. INTRINSIC FUNCTIONS ..
C
      INTRINSIC DFLOAT,ACOS,ATAN,COS,IABS

C-----------------------------------------------------------------------
      IF(LMAX.LE.LMAXD1) GO TO 1
      WRITE(6,200)LMAX,LMAXD1
  200 FORMAT(3X,'FROM PINTG: LMAX=',I4,' GREATER THAN DIMENSIONED',I4)
      STOP
    1 CONTINUE
      DO 6 I=0,LMAX
    6 S( 0,I)=0.D0
      DO 5 M=1,LMAX
      DO 5 I=0,LMAX-M
      S(-M,I)=0.D0
    5 S( M,I)=0.D0
      IF(ITYPE.NE.0)      GO TO 10
      THETA=ACOS(ARG)
      CALL RECUR0(LMAX,X1,THETA,-DFLOAT(ISI),S)
      CALL RECUR0(LMAX,X2,THETA, DFLOAT(ISI),S)
      RETURN
C                         E N D    I F
   10 CONTINUE
      N=(X2-X1)/DLT+3
      IF(N.GT.NDIM) STOP 'INCREASE NDIM'
      CALL GAULEG(X1,X2,XX,WW,N)
      DO 2 K=1,N
      X=XX(K)
      W=DFLOAT(ISI)*WW(K)
      THETA=ATAN(ARG/COS(X-FD))
      CALL RECUR(LMAX,X,THETA,W,S)
    2 CONTINUE
      RETURN
      END SUBROUTINE

C======================================================================

C-----------------------------------------------------------------------
C>    THIS ROUTINE IS USED TO PERFORM THE FI-INTEGRATION OF REAL SPHE-
C>    RICAL HARMONICS .THE THETA-INTEGRATION IS PERFORMED ANALYTICALLY
C>    USING RECURRENCE RELATIONS.
C-----------------------------------------------------------------------
      SUBROUTINE RECUR(LMAX,X,THETA,FAC,S)
      USE shape_constants_mod, ONLY: LMAXD1
      implicit none

C
C     .. PARAMETER STATEMENTS ..
C
C      include 'inc.geometry'
C
C     .. SCALAR ARGUMENTS ..
C
      INTEGER   LMAX
      REAL*8 X,THETA,FAC
C
C     .. ARRAY ARGUMENTS ..
C
      REAL*8 S(-LMAXD1:LMAXD1,0:LMAXD1)
C
C     .. LOCAL SCALARS ..
C
      INTEGER   M,I
      REAL*8 OL0,OL,EL0,EL,C1,C2,SS,CC
      REAL*8 C01(LMAXD1),C02(LMAXD1),SSA(LMAXD1+2),CCA(LMAXD1+2)
C
C     .. INTRINSIC FUNCTIONS ..
C
      INTRINSIC COS,SIN,DFLOAT
C-----------------------------------------------------------------------
      SS=SIN(THETA)
      CC=COS(THETA)
      DO 13 I=1,LMAX
      C01(I)=FAC*SIN(DFLOAT(I)*X)
   13 C02(I)=FAC*COS(DFLOAT(I)*X)
      DO 11 I=1,LMAX+2
      SSA(I)=SS**I
      CCA(I)=CC**I
   11 CONTINUE
      OL0=(THETA-SS*CC)/2.D0
      EL0=0D0
      DO 1 M=1,LMAX,2
      OL=OL0
      EL=EL0
      C1=C01(M)
      C2=C02(M)
      DO 2 I=0,LMAX-M,2
      S(-M,I)=S(-M,I)+C1*(OL-EL)
      S( M,I)=S( M,I)+C2*(OL-EL)
      EL= DFLOAT(I+1)*EL/DFLOAT(I+M+3)
    2 OL=(DFLOAT(I+1)*OL+(SSA(M+2))*(CCA(I+1)))/DFLOAT(I+M+3)
      EL0= DFLOAT(M+2)*EL0/DFLOAT(M+3)
      OL0=(DFLOAT(M+2)*OL0-(SSA(M+2))*CC)/DFLOAT(M+3)
    1 CONTINUE
      OL0=SSA(3)/3.D0
      EL0=0D0
      DO 3 M=1,LMAX,2
      OL=OL0
      EL=EL0
      C1=C01(M)
      C2=C02(M)
      DO 4 I=1,LMAX-M,2
      S(-M,I)=S(-M,I)+C1*(OL-EL)
      S( M,I)=S( M,I)+C2*(OL-EL)
      OL=(DFLOAT(I+1)*OL+SSA(M+2)*CCA(I+1))/DFLOAT(I+M+3)
    4 EL= DFLOAT(I+1)*EL/DFLOAT(I+M+3)
      OL0=(DFLOAT(M+2)*OL0-(SSA(M+2))*CCA(2))/DFLOAT(M+4)
      EL0= DFLOAT(M+2)*EL0/DFLOAT(M+4)
    3 CONTINUE
      OL0=-CC
      EL0=-1D0
      OL=OL0
      EL=EL0
      C2=FAC
      DO 9 I=0,LMAX,2
      S(0,I)=S(0,I)+C2*(OL-EL)
      OL=(DFLOAT(I+1)*OL+SSA(2)*CCA(I+1))/DFLOAT(I+3)
    9 EL= DFLOAT(I+1)*EL/DFLOAT(I+3)
      OL0=(2.D0*OL0-SSA(2)*CC)/3.D0
      EL0= 2.D0*EL0/3.D0
      DO 5 M=2,LMAX,2
      OL=OL0
      EL=EL0
      C1=C01(M)
      C2=C02(M)
      DO 6 I=0,LMAX-M,2
      S(-M,I)=S(-M,I)+C1*(OL-EL)
      S( M,I)=S( M,I)+C2*(OL-EL)
      OL=(DFLOAT(I+1)*OL+SSA(M+2)*CCA(I+1))/DFLOAT(I+M+3)
    6 EL= DFLOAT(I+1)*EL/DFLOAT(I+M+3)
      OL0=(DFLOAT(M+2)*OL0-SSA(M+2)*CC)/DFLOAT(M+3)
      EL0= DFLOAT(M+2)*EL0/DFLOAT(M+3)
    5 CONTINUE
      OL0=-CCA(2)/2D0
      EL0=-0.5D0
      OL=OL0
      EL=EL0
      C2=FAC
      DO 10 I=1,LMAX,2
      S(0,I)=S(0,I)+C2*(OL-EL)
      OL=(DFLOAT(I+1)*OL+SS*SS*CCA(I+1))/DFLOAT(I+3)
   10 EL= DFLOAT(I+1)*EL/DFLOAT(I+3)
      OL0=(2.D0*OL0-SSA(2)*CCA(2))/4.D0
      EL0= 2.D0*EL0/4.D0
      DO 7 M=2,LMAX,2
      OL=OL0
      EL=EL0
      C1=C01(M)
      C2=C02(M)
      DO 8 I=1,LMAX-M,2
      S(-M,I)=S(-M,I)+C1*(OL-EL)
      S( M,I)=S( M,I)+C2*(OL-EL)
      OL=(DFLOAT(I+1)*OL+SSA(M+2)*CCA(I+1))/DFLOAT(I+M+3)
    8 EL= DFLOAT(I+1)*EL/DFLOAT(I+M+3)
      OL0=(DFLOAT(M+2)*OL0-SSA(M+2)*CCA(2))/DFLOAT(M+4)
      EL0= DFLOAT(M+2)*EL0/DFLOAT(M+4)
    7 CONTINUE
      RETURN
      END SUBROUTINE

C======================================================================

C-----------------------------------------------------------------------
C>    THIS ROUTINE IS USED TO PERFORM  THE  FI - INTEGRATION  OF REAL SP
C>    RICAL HARMONICS ANALYTICALLY.  THE  THETA-INTEGRATION  IS   PERFOR
C>    ALSO ANALYTICALLY USING RECURRENCE RELATIONS.(THETA IS FI-INDEPEND
C-----------------------------------------------------------------------
      SUBROUTINE RECUR0(LMAX,X,THETA,FAC,S)
      USE shape_constants_mod, ONLY: LMAXD1
      implicit none
C
C     .. PARAMETER STATEMENTS ..
C
C      include 'inc.geometry'
c      INTEGER LMAXD
c      PARAMETER (LMAXD=25)
C
C     .. SCALAR ARGUMENTS ..
C
      INTEGER   LMAX
      REAL*8 X,THETA,FAC
C
C     .. ARRAY ARGUMENTS ..
C
      REAL*8 S(-LMAXD1:LMAXD1,0:LMAXD1)
C
C     .. LOCAL SCALARS ..
C
      INTEGER   M,I
      REAL*8 OL0,OL,EL0,EL,C1,C2,SS,CC
C
C     .. INTRINSIC FUNCTIONS ..
C
      INTRINSIC COS,SIN,DFLOAT
C-----------------------------------------------------------------------
      SS=SIN(THETA)
      CC=COS(THETA)
      OL0=(THETA-SS*CC)/2D0
      EL0=0D0
      DO 1 M=1,LMAX,2
      OL=OL0
      EL=EL0
      C1=-FAC*COS(DFLOAT(M)*X)/DFLOAT(M)
      C2= FAC*SIN(DFLOAT(M)*X)/DFLOAT(M)
      DO 2 I=0,LMAX-M,2
      S(-M,I)=S(-M,I)+C1*(OL-EL)
      S( M,I)=S( M,I)+C2*(OL-EL)
      EL= DFLOAT(I+1)*EL/DFLOAT(I+M+3)
    2 OL=(DFLOAT(I+1)*OL+(SS**(M+2))*(CC**(I+1)))/DFLOAT(I+M+3)
      EL0= DFLOAT(M+2)*EL0/DFLOAT(M+3)
      OL0=(DFLOAT(M+2)*OL0-(SS**(M+2))*CC)/DFLOAT(M+3)
    1 CONTINUE
      OL0=SS*SS*SS/3D0
      EL0=0D0
      DO 3 M=1,LMAX,2
      OL=OL0
      EL=EL0
      C1=-FAC*COS(DFLOAT(M)*X)/DFLOAT(M)
      C2= FAC*SIN(DFLOAT(M)*X)/DFLOAT(M)
      DO 4 I=1,LMAX-M,2
      S(-M,I)=S(-M,I)+C1*(OL-EL)
      S( M,I)=S( M,I)+C2*(OL-EL)
      OL=(DFLOAT(I+1)*OL+(SS**(M+2))*(CC**(I+1)))/DFLOAT(I+M+3)
    4 EL= DFLOAT(I+1)*EL/DFLOAT(I+M+3)
      OL0=(DFLOAT(M+2)*OL0-(SS**(M+2))*CC*CC)/DFLOAT(M+4)
      EL0= DFLOAT(M+2)*EL0/DFLOAT(M+4)
    3 CONTINUE
      OL0=-CC
      EL0=-1D0
      OL=OL0
      EL=EL0
      C2=FAC*X
      DO 9 I=0,LMAX,2
      S(0,I)=S(0,I)+C2*(OL-EL)
      OL=(DFLOAT(I+1)*OL+SS*SS*(CC**(I+1)))/DFLOAT(I+3)
    9 EL= DFLOAT(I+1)*EL/DFLOAT(I+3)
      OL0=(2.D0*OL0-SS*SS*CC)/3.D0
      EL0= 2.D0*EL0/3.D0
      DO 5 M=2,LMAX,2
      OL=OL0
      EL=EL0
      C1=-FAC*COS(DFLOAT(M)*X)/DFLOAT(M)
      C2 =FAC*SIN(DFLOAT(M)*X)/DFLOAT(M)
      DO 6 I=0,LMAX-M,2
      S(-M,I)=S(-M,I)+C1*(OL-EL)
      S( M,I)=S( M,I)+C2*(OL-EL)
      OL=(DFLOAT(I+1)*OL+(SS**(M+2))*(CC**(I+1)))/DFLOAT(I+M+3)
    6 EL= DFLOAT(I+1)*EL/DFLOAT(I+M+3)
      OL0=(DFLOAT(M+2)*OL0-(SS**(M+2))*CC)/DFLOAT(M+3)
      EL0= DFLOAT(M+2)*EL0/DFLOAT(M+3)
    5 CONTINUE
      OL0=-CC*CC/2D0
      EL0=-0.5D0
      OL=OL0
      EL=EL0
      C2=FAC*X
      DO 10 I=1,LMAX,2
      S(0,I)=S(0,I)+C2*(OL-EL)
      OL=(DFLOAT(I+1)*OL+SS*SS*(CC**(I+1)))/DFLOAT(I+3)
   10 EL= DFLOAT(I+1)*EL/DFLOAT(I+3)
      OL0=(2.D0*OL0-SS*SS*CC*CC)/4.D0
      EL0= 2.D0*EL0/4.D0
      DO 7 M=2,LMAX,2
      OL=OL0
      EL=EL0
      C1=-FAC*COS(DFLOAT(M)*X)/DFLOAT(M)
      C2= FAC*SIN(DFLOAT(M)*X)/DFLOAT(M)
      DO 8 I=1,LMAX-M,2
      S(-M,I)=S(-M,I)+C1*(OL-EL)
      S( M,I)=S( M,I)+C2*(OL-EL)
      OL=(DFLOAT(I+1)*OL+(SS**(M+2))*(CC**(I+1)))/DFLOAT(I+M+3)
    8 EL= DFLOAT(I+1)*EL/DFLOAT(I+M+3)
      OL0=(DFLOAT(M+2)*OL0-(SS**(M+2))*CC*CC)/DFLOAT(M+4)
      EL0= DFLOAT(M+2)*EL0/DFLOAT(M+4)
    7 CONTINUE
      RETURN
      END SUBROUTINE

C======================================================================

C     ----------------------------------------------------------------
C>    GINEN THE LOWER AND UPPER LIMITS OF INTEGRATION  X1 AND X2, AND
C>    GIVEN N, THIS SUBROUTINE RETURNS THE  ARRAYS X(1:N) AND  W(1:N)
C>    OF LENGTH N, CONTAINING THE ABSCISSAS AND WEIGHTS OF THE  GAUSS
C>    LEGENDRE N-POINT QUADRATURE FORMULA (NUMERICAL RECIPES,2ND ED.).
C     ----------------------------------------------------------------
      SUBROUTINE GAULEG(X1,X2,X,W,N)
      IMPLICIT NONE
C
C     .. SCALAR ARGUMENTS ..
C
      INTEGER   N
      REAL*8    X1,X2
C
C     .. ARRAY ARGUMENTS ..
C
      REAL*8    X(1),W(1)
C
C     .. LOCAL SCALARS ..
C
      INTEGER   I,J,M
      REAL*8    P1,P2,P3,PP,XL,XM,Z,Z1,PI314
C     ----------------------------------------------------------------
      PI314 = 4.D0*DATAN(1.D0)

      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
      DO 12 I=1,M
      Z=DCOS(PI314*(I-.25D0)/(N+.5D0))
    1 CONTINUE
      P1=1.D0
      P2=0.D0
      DO 11 J=1,N
      P3=P2
      P2=P1
      P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
   11 CONTINUE
      PP=N*(Z*P1-P2)/(Z*Z-1.D0)
      Z1=Z
      Z=Z1-P1/PP
      IF(ABS(Z-Z1).GT.3.D-14) GO TO 1
      X(I)=XM-XL*Z
      X(N+1-I)=XM+XL*Z
      W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
      W(N+1-I)=W(I)
   12 CONTINUE
      RETURN
      END SUBROUTINE

C-----------------------------------------------------------------------
C>    THIS ROUTINE CALCULATES THE COEFFICIENTS OF A POLYNOMIAL EXPANSION
C>    OF RENORMALIZED LEGENDRE FUNCTIONS IN POWERS OF COSINES.
C>    THE POSSIBILITY OF OVERFLOW (HIGH LMAX) IS AVOIDED BY USING FACTO-
C>    RIZED FORMS FOR THE NUMBERS.
C-----------------------------------------------------------------------
C     CHANGED: get constants from module shape_constants_mod instead
C              of inc.geometry
      SUBROUTINE CCOEF(LMAX,CL,COE)
      USE shape_constants_mod, ONLY: ICD, ICED, LMAXD1
      implicit none

C
C     .. PARAMETER STATEMENTS ..
C
C     include 'inc.geometry'
C     INTEGER ICD,ICED
      INTEGER IFMX,LMA2D
C     PARAMETER (ICD=1729,ICED=((LMAXD1+1)*(LMAXD1+2))/2)
      PARAMETER (IFMX=25,LMA2D=LMAXD1/2+1)
C
C     .. SCALAR ARGUMENTS ..
C
      INTEGER   LMAX
C
C     .. ARRAY ARGUMENTS ..
C
      REAL*8    CL(ICD),COE(ICED)
C
C     .. LOCAL SCALARS ..
C
      INTEGER   ICMAX,L,LI,ICE,IC,I1,L2P,M,K,K0,ISI,IRE,IR,IC1,IC2
      INTEGER   LA,LB,IEUPSQ,IEINT,IEMOD
      REAL*8    UP,DOWN,UPSQ
C
C     .. LOCAL ARRAYS ..
C
      INTEGER   IE(IFMX,LMA2D),IED(IFMX),IFI(IFMX)
      INTEGER   L1ST(IFMX),L2ST(IFMX),L1(IFMX),L2(IFMX),JM0(IFMX)
      INTEGER   IEA(IFMX),IEB(IFMX),IL2P(IFMX)
C
C     .. INTRINSIC FUNCTIONS ..
C
      INTRINSIC MOD,SQRT,DFLOAT
C
C     .. DATA STATEMENTS ..
C
      DATA IFI/2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,
     &         61,67,71,73,79,83,89,97/
C-----------------------------------------------------------------------
      ICMAX=0
      DO 9 L=0,LMAX
      LI=L/2+1
    9 ICMAX=ICMAX+(LMAX+1-L)*LI
      IF(LMAX.GT.LMAXD1.OR.ICMAX.GT.ICD) GO TO 100
C     if (TEST('SHAPE   ')) WRITE(6,203) ICMAX
      ICE=0
      IC=1
      L=0
      DO 11 I1=1,IFMX
      L1ST(I1)=0
   11 L2ST(I1)=0
    1 CONTINUE
      L2P=2*L+1
      CALL REDUCE(L2P,IFMX,IFI,IL2P)
      IL2P(1)=IL2P(1)+1
      M=L
      DO 12 I1=1,IFMX
      L1(I1)=L1ST(I1)
      L2(I1)=L2ST(I1)
   12 JM0(I1)=0
    2 CONTINUE
      ICE=ICE+1
      ISI=1
C  THIS IS CHANGED
C     ISI=1-2*MOD(M,2)
C
      K0=(L+M+1)/2
      K=L
      IRE=1
      IC1=IC
      DO 13 I1=1,IFMX
   13 IE(I1,IRE)=JM0(I1)
    3 CONTINUE
      IF((K-1). LT .K0)  GO TO 30
      IRE=IRE+1
      IC=IC+1
      LA=(2*K-L-M)*(2*K-L-M-1)
      LB=2*(2*K-1)*(L-K+1)
      CALL REDUCE(LA,IFMX,IFI,IEA)
      CALL REDUCE(LB,IFMX,IFI,IEB)
      DO 14 I1=1,IFMX
   14 IE(I1,IRE)=IE(I1,IRE-1)+IEA(I1)-IEB(I1)
      K=K-1
      GO TO 3
   30 CONTINUE
      IC2=IC
      DO 4 I1=1,IFMX
      IED(I1)=IE(I1,1)
      DO 5 IR=2,IRE
      IF(IE(I1,IR). LT .IED(I1))  IED(I1)=IE(I1,IR)
    5 CONTINUE
      DO 6 IR=1,IRE
    6 IE(I1,IR)=IE(I1,IR)-IED(I1)
    4 CONTINUE
      IR=0
      DO 7 IC=IC1,IC2
      IR=IR+1
      CL(IC)=1.D0
      DO 8 I1=1,IFMX
    8 CL(IC)=CL(IC)*IFI(I1)**IE(I1,IR)
      CL(IC)=ISI*CL(IC)
      ISI=-ISI
    7 CONTINUE
      IF(M.EQ.0) IL2P(1)=IL2P(1)-1
      UP  =1.D0
      UPSQ=1.D0
      DOWN=1.D0
      DO 17 I1=1,IFMX
      IEUPSQ=2*IED(I1)+IL2P(I1)+L1(I1)
      IEINT =IEUPSQ/2-L2(I1)
      IEMOD =MOD(IEUPSQ,2)
      UPSQ =UPSQ*IFI(I1) **IEMOD
      IF(IEINT.GE.0)                   T H E N
      UP   =UP  *IFI(I1) **IEINT
                                       E L S E
      DOWN =DOWN*IFI(I1) **(-IEINT)
                                       E N D    I F
   17 CONTINUE
      COE(ICE)=SQRT(UPSQ)* UP / DOWN
C     if (TEST('SHAPE   '))
C    &      WRITE(6,201) L,M,UP,UPSQ,DOWN,(CL(IC),IC=IC1,IC2)
      IF(M. EQ .0)  GO TO 20
      LA=L+M
      LB=L-M+1
      CALL REDUCE(LA,IFMX,IFI,IEA)
      CALL REDUCE(LB,IFMX,IFI,IEB)
      DO 15 I1=1,IFMX
      JM0(I1)=JM0(I1)+IEA(I1)-IEB(I1)
   15 L1 (I1)=L1 (I1)-IEA(I1)+IEB(I1)
      M=M-1
      GO TO 2
   20 CONTINUE
C     if (TEST('SHAPE   ')) WRITE(6,202)
      IF(L. EQ .LMAX)   GO TO 10
      LA=(2*L+1)*(2*L+2)
      LB=(L+1)*2
      CALL REDUCE(LA,IFMX,IFI,IEA)
      CALL REDUCE(LB,IFMX,IFI,IEB)
      DO 16 I1=1,IFMX
      L1ST(I1)=L1ST(I1)+IEA(I1)
   16 L2ST(I1)=L2ST(I1)+IEB(I1)
      L=L+1
      GO TO 1
   10 CONTINUE
      RETURN
  100 WRITE(6,204) LMAX,LMAXD1,ICMAX,ICD
      STOP
  201 FORMAT(2X,'L=',I2,' M=',I2,F10.3,' *SQRT(',F16.2,')/',F10.3/
     *       2X,'CL  :',6F14.2)
  202 FORMAT(80('*'))
  203 FORMAT(13X,'THERE ARE',I5,'  COEFFICIENTS'/)
  204 FORMAT(13X,'FROM CCOEF: INCONSISTENCY DATA-DIMENSION'/
     *       14X,'LMAX:',2I5/13X,'ICMAX:',2I5)
      END SUBROUTINE


C-----------------------------------------------------------------------
C>    THIS ROUTINE REDUCES A POSITIVE INTEGER   INPUT NUMBER 'NMBR'
C>    TO A PRODUCT  OF  FIRST  NUMBERS 'IFI' , AT POWERS  'IEXP'.
C-----------------------------------------------------------------------
      SUBROUTINE REDUCE(NMBR,IFMX,IFI,IEXP)
      implicit none

C
C     .. SCALAR ARGUMENTS ..
C
      INTEGER   NMBR,IFMX
C
C     .. ARRAY  ARGUMENTS ..
C
      INTEGER   IEXP(1),IFI(1)
C
C     .. LOCAL SCALARS ..
C
      INTEGER   I,NMB
C
C     .. INTRINSIC FUNCTIONS ..
C
      INTRINSIC MOD
C-----------------------------------------------------------------------
      IF(NMBR.LE.0) GO TO 4
      DO 5 I=1,IFMX
    5 IEXP(I)=0
      IF(NMBR.EQ.1) RETURN
      NMB=NMBR
      DO 1 I=1,IFMX
      IEXP(I)=0
    2 IF(MOD(NMB,IFI(I)).NE.0) GO TO 3
      NMB=NMB/IFI(I)
      IEXP(I)=IEXP(I)+1
      GO TO 2
    3 CONTINUE
      IF(NMB.EQ.1) RETURN
    1 CONTINUE
      WRITE(6,100) NMBR
      STOP
    4 WRITE(6,101) NMBR
      STOP
  100 FORMAT(3X,I15,'  CANNOT BE REDUCED IN THE BASIS OF FIRST NUMBERS G
     *IVEN'/20X,'INCREASE THE BASIS OF FIRST NUMBERS')
  101 FORMAT(3X,I15,'  NON POSITIVE NUMBER')
      END SUBROUTINE

C------------------------------------------------------------------
C>    THIS ROUTINE COMPUTES TRANSFORMATION MATRICES ASSOCIATED TO
C>    THE ROTATION THROUGH THE EULER ANGLES ALPHA,BETA,GAMMA  FOR
C>    REAL SPHERICAL HARMONICS UP TO QUANTUM NUMBER LMAX. THE RE-
C>    SULTS ARE STORED IN DMATL(ISUMD)
C------------------------------------------------------------------
      SUBROUTINE DREAL(LMAX,ALPHA,BETA,GAMMA,
     &                 DMATL, ISUMD,
     &                 LMAXD1)
C     new parameters: DMATL, ISUMD
C     LMAXD1 from inc.geometry
      implicit none
C
C     .. PARAMETER STATEMENTS ..
C
C     include 'inc.geometry'
c      INTEGER LMAXD,ISUMD
c     PARAMETER (LMAXD=25,ISUMD=23426)
c      PARAMETER (LMAXD=25,ISUMD=100000)
C
C     .. SCALAR ARGUMENTS ..
C
      INTEGER   LMAXD1
      INTEGER   ISUMD
      INTEGER   LMAX
      REAL*8    ALPHA,BETA,GAMMA
C
C     .. LOCAL SCALARS ..
C
      INTEGER   L,M,MP,I,IMAX,IP,IPMAX,ISU,ISUM
      REAL*8    SQR2,FAC,FAC1,FAC2,D,D1,D2
C
C     .. LOCAL ARRAYS ..
C
      REAL*8 DMN(LMAXD1+1,LMAXD1+1),DPL(LMAXD1+1,LMAXD1+1),DMATL(ISUMD)
C
C     .. EXTERNAL ROUTINES ..
C
C     REAL*8 DROT
C
C     .. INTRINSIC FUNCTIONS ..
C
      INTRINSIC COS,SIN,MOD,SQRT
C-----------------------------------------------------------------------
      SQR2=SQRT(2.D0)
      ISU=0
      DO 2 L=0,LMAX
      FAC2=1.D0
      M=0
    1 FAC1=1.D0
      MP=0
    3 CONTINUE
      FAC=FAC1*FAC2/2.D0
      D1=DROT(L,MP ,M,BETA)
      D2=DROT(L,MP,-M,BETA)
      IF ( MOD ( M , 2       ) .NE. 0 )  D2 = - D2
      DPL(MP+1,M+1)=(D1+D2)*FAC
      DMN(MP+1,M+1)=(D1-D2)*FAC
      IF ( MOD ( M + MP , 2 ) .NE. 0 )  GOTO  4
      DPL(M+1,MP+1)=DPL(MP+1,M+1)
      DMN(M+1,MP+1)=DMN(MP+1,M+1)
                                             GOTO  5
    4 DMN(M+1,MP+1)=-DMN(MP+1,M+1)
      DPL(M+1,MP+1)=-DPL(MP+1,M+1)
    5 CONTINUE
      FAC1=SQR2
      MP=1+MP
      IF (MP .LE. M) GOTO 3
      FAC2=SQR2
      M=1+M
      IF (M .LE. L) GOTO 1
      IMAX=1
      M=0
   12 CONTINUE
      DO 13 I=1,IMAX
      IPMAX=1
      MP=0
   11 CONTINUE
      DO 6 IP=1,IPMAX
      IF ( I   .EQ. 2 )                      GOTO   7
      IF ( IP .EQ. 2 )                      GOTO  10
      D= COS(MP*ALPHA)*COS(M*GAMMA)*DPL(MP+1,M+1)
     *  -SIN(MP*ALPHA)*SIN(M*GAMMA)*DMN(MP+1,M+1)
                                                 GOTO   9
    7 IF ( IP .EQ. 2 )                      GOTO   8
      D=-COS(MP*ALPHA)*SIN(M*GAMMA)*DPL(MP+1,M+1)
     *  -SIN(MP*ALPHA)*COS(M*GAMMA)*DMN(MP+1,M+1)
                                                 GOTO   9
    8 D=-SIN(MP*ALPHA)*SIN(M*GAMMA)*DPL(MP+1,M+1)
     *  +COS(MP*ALPHA)*COS(M*GAMMA)*DMN(MP+1,M+1)
                                                 GOTO   9
   10 D= SIN(MP*ALPHA)*COS(M*GAMMA)*DPL(MP+1,M+1)
     *  +COS(MP*ALPHA)*SIN(M*GAMMA)*DMN(MP+1,M+1)
    9 CONTINUE
C THIS IS CHANGED
      IF(MOD(M+MP,2) .NE. 0)  D=-D
C
      ISU=ISU+1
      DMATL(ISU)=D
    6 CONTINUE
      IPMAX=2
      MP=MP+1
      IF(MP.LE.L) GOTO 11
   13 CONTINUE
      IMAX=2
      M=M+1
      IF(M .LE.L) GOTO 12
    2 CONTINUE
      ISUM=ISU
C     WRITE(ITEMP) (DMATL(ISU),ISU=1,ISUM)
      RETURN
      END SUBROUTINE

C-----------------------------------------------------------------------
C>    CALCULATION OF D COEFFICIENT ACCORDING TO ROSE, ELEMENTARY THEORY
C>    ANGULAR MOMENTUM,J.WILEY & SONS ,1957 , EQ. (4.13).
C-----------------------------------------------------------------------
      FUNCTION DROT(L,MP,M,BETA)
      implicit none
C
C     .. SCALAR ARGUMENTS ..
C
      INTEGER   L,M,MP
      REAL*8    BETA,DROT
C
C     .. LOCAL SCALARS ..
C
      INTEGER   L0,M0,MP0,N1,N2,N3,N4,NN,I,KMIN,KMAX,LTRM,N,K
      REAL*8    SINB,TERM,FF,BETA2,COSB
C
C     .. LOCAL ARRAYS ..
C
      INTEGER        NF(4)
      EQUIVALENCE (N1,NF(1)),(N2,NF(2)),(N3,NF(3)),(N4,NF(4))
C
C     .. INTRINSIC FUNCTIONS ..
C
      INTRINSIC IABS,ABS,SQRT,COS,SIN,MOD,MIN0,MAX0
C-----------------------------------------------------------------------
      DATA L0,M0,MP0/-1,1,1/
         L0=-1
         M0=1
         MP0=1
   10 IF(L.NE.L0) GO TO 2
                           WRITE(6,300) L0,M0,MP0
  300                      FORMAT(3X,'L0,M0,MP0=',3I3)
      IF((IABS(M ).EQ.IABS(M0).AND.IABS(MP).EQ.IABS(MP0)).OR.
     1   (IABS(MP).EQ.IABS(M0).AND.IABS(M ).EQ.IABS(MP0))) GO TO 1
    2 FF=1.D0
      IF(IABS(M).LE.L.AND.IABS(MP).LE.L) GO TO 3
      WRITE(6,100) L,M,MP
  100 FORMAT('     L=',I5,'    M=',I5,'    MP =',I5)
      RETURN
    3 N1   =L+M
      N2   =L-M
      N3   =L+MP
      N4   =L-MP
      L0=L
      M0=M
      MP0=MP
      DO 4 N=1,4
      NN=NF(N)
      IF(NN.EQ.0) GO TO 4
      DO 5 I=1,NN
    5 FF=FF*I
    4 CONTINUE
      FF=SQRT(FF)
    1 BETA2=BETA/2.D0
      COSB= COS(BETA2)
      SINB=-SIN(BETA2)
      IF(ABS(COSB).LT.1.E-4) GO TO 9
      IF(ABS(SINB).LT.1.E-4) GO TO 11
      KMAX=MIN0(L-MP,L+M)
      KMIN=MAX0(M-MP,0)
      TERM=COSB**(2*L+M-MP-2*KMIN)*SINB**(MP-M+2*KMIN)*FF
      GO TO 12
    9 LTRM=L
      TERM=FF
      IF(SINB.LT.0.D0.AND.MOD(MP-M,2).NE.0) TERM=-TERM
      GO TO 14
   11 LTRM=0
      TERM=FF
      IF(COSB.LT.0.D0.AND.MOD(MP-M,2).NE.0) TERM=-TERM
   14 KMAX=M-MP
      IF(MOD(KMAX,2).NE.0) GO TO 13
      KMAX=LTRM+KMAX/2
      IF(KMAX.LT.MAX0(M-MP,0)) GO TO 13
      IF(KMAX.GT.MIN0(L-MP,L+M)) GO TO 13
      KMIN=KMAX
   12 IF(MOD(KMIN,2).NE.0) TERM=-TERM
      N1   =L-MP-KMIN
      N2   =L+M-KMIN
      N3   =KMIN+MP-M
      N4   =KMIN
      DO 6 N=1,4
      NN=NF(N)
      IF(NN.EQ.0) GO TO 6
      DO 7 I=1,NN
    7 TERM=TERM/I
    6 CONTINUE
      DROT=TERM
      IF(KMIN.EQ.KMAX) RETURN
      KMIN=KMIN+1
      COSB=COSB**2
      SINB=SINB**2
      N3=N3   +1
      DO 8 K=KMIN,KMAX
      TERM=-N1*N2*TERM*SINB/(COSB*K*N3)
      DROT=DROT+TERM
      N1=N1-1
      N2=N2-1
    8 N3=N3+1
      RETURN
   13 DROT=0.D0
      RETURN
      END FUNCTION

      end module ShapeIntegrationHelpers_mod
