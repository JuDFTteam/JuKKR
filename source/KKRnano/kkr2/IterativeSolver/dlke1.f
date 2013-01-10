      SUBROUTINE DLKE1(ALAT,NACLS,RR,EZOA,
     &                 BZKP,IC,EIKRM,EIKRP,
C                      new input parameters after inc.p removal
     &                 nrd, naclsd)
      IMPLICIT NONE
c ----------------------------------------------------------------------
c
c     Fourier transformation of the cluster Greens function
c     Prepares the calculation (calculates Fourier factors) for dlke0
c ----------------------------------------------------------------------


C     >> Input parameters
C     ALAT .. lattice constant a
C     NACLS . array that contains number of atoms in each cluster
C     RR .... array of real space vectors
C     EZOA ..
C     BZKP ..
C     IC .... reference cluster index
C     nrd ... There are nrd+1 real space vectors in RR
C     naclsd. maximal number of atoms in cluster

C     << Output parameters
C     EIKRM . Fourier exponential factor with minus sign
C     EIKRP . Fourier exponential factor with plus sign

      INTEGER nrd
      INTEGER naclsd

      DOUBLE COMPLEX CI,CONE
      PARAMETER (CI= (0.0D0,1.0D0),CONE=(1.D0,0.D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALAT
      INTEGER IC
C     ..
C     .. Array Arguments ..
      INTEGER EZOA(*),NACLS(*)
      DOUBLE COMPLEX EIKRP(NACLSD),EIKRM(NACLSD)
      DOUBLE PRECISION BZKP(*),RR(3,0:NRD)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CONVPU,TPI
      INTEGER M
      DOUBLE COMPLEX TT
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX ARG(3)
C     ..
C     .. External Subroutines ..
      EXTERNAL CINIT,TEST,OPT,ZAXPY
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,EXP
C     ..
C     .. Save statement ..
      SAVE
C     ..
C
      TPI = 8.0D0*ATAN(1.0D0)         
      CONVPU = ALAT/TPI

      DO 90 M = 1,NACLS(IC)

c
c     
c     Here we do   --                  nn'
c                  \                   ii'          ii'
c                  /  exp(+ik(x  -x ))G   (E)  =   G   (k,E)
c                  --          n'  n   LL'          LL'
c                  n'
c  Be careful a minus sign must be included here. RR is not
c  symmetric around each atom. The minus comes from the fact that
c  the repulsive potential GF is calculated for 0n and not n0!                   
c  and that is why we need a minus sign extra!
c  

           ARG(1) = -CI*TPI*RR(1,EZOA(M))
           ARG(2) = -CI*TPI*RR(2,EZOA(M))
           ARG(3) = -CI*TPI*RR(3,EZOA(M))
c
        TT = BZKP(1)*ARG(1)+BZKP(2)*ARG(2)+BZKP(3)*ARG(3)
c
c  convert to p.u. and multiply with 1/2
        EIKRP(M) = EXP(TT) * CONVPU * 0.5D0
        EIKRM(M) = EXP(-TT) * CONVPU * 0.5D0
c
 90   CONTINUE                    

      RETURN
      END
