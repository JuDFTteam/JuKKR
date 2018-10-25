      SUBROUTINE DIRBSRAD(XBS,Y,DYDX,DRDI,B,V,R,NMESH)
C
C   ********************************************************************
C   *                                                                  *
C   *   supply the derivatives for the coupled set of                  *
C   *   radial dirac equation in case of a spin-dependent potential    *
C   *                                                                  *
C   ********************************************************************
C

      IMPLICIT NONE

      INCLUDE 'sprkkr_rmesh.dim'

C
C PARAMETER definitions
C
      INTEGER NLAG
      PARAMETER (NLAG=3)
C
C COMMON variables
C
      REAL*8 CGD(2),CGMD(2),CGO,KAP(2)
      REAL*8 CSQR
      COMPLEX*16 EBS
      INTEGER NRADBS,NSOLBS
      COMMON /COMMBS/ EBS,CSQR,CGD,CGMD,CGO,KAP,NSOLBS,NRADBS
C
C Dummy arguments
C
      INTEGER NMESH
      REAL*8 XBS
      REAL*8 B(NRMAX),DRDI(NRMAX),R(NRMAX),V(NRMAX)
      COMPLEX*16 DYDX(NCFMAX),Y(NCFMAX)
C
C Local variables
C
      REAL*8 BBS,BPP,BQQ,DRDIBS,RBS,VBS
      COMPLEX*16 EMVPP,EMVQQ
      INTEGER I,J,K,M
C
      CALL DIRBSLAG(XBS,VBS,BBS,RBS,DRDIBS,V,B,R,DRDI,NRADBS,NLAG,NMESH)
C
      EMVQQ = (EBS-VBS+CSQR)*DRDIBS/CSQR
      EMVPP = -(EBS-VBS)*DRDIBS
      BQQ = BBS*DRDIBS/CSQR
      BPP = BBS*DRDIBS
C
C
      M = 0
      DO J = 1,NSOLBS
         DO I = 1,NSOLBS
            K = 3 - 4*(I-1)
            DYDX(M+1) = -KAP(I)*Y(M+1)/RBS*DRDIBS + (EMVQQ+BQQ*CGMD(I))
     &                  *Y(M+2)
            DYDX(M+2) = KAP(I)*Y(M+2)/RBS*DRDIBS + (EMVPP+BPP*CGD(I))
     &                  *Y(M+1) + BPP*CGO*Y(M+K)
            M = M + 2
         END DO
      END DO
      END 
