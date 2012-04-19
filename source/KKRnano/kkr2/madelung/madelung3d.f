      SUBROUTINE MADELUNG3D(LPOT,YRG,WG,ALAT,
     &                      RMAX,GMAX,BRAVAIS,RECBV,
     &                      LMXSPD,LASSLD,LPOTD,LMPOTD,
     &                      NMAXD,ISHLD,
     &                      LMPOT,CLEB,ICLEB,IEND,
     &                      NCLEBD,LOFLM,DFAC,
     &                      NGMAX,NRMAX,NSG,NSR,NSHLG,NSHLR,GN,RM)
C **********************************************************************
C *                                                                    *
C * This subroutine calculates the Madelung potential coefficients     *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER LPOT
      INTEGER LMXSPD,LASSLD,LPOTD,LMPOTD,NMAXD,ISHLD
      DOUBLE PRECISION ALAT,RMAX,GMAX
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION YRG(LASSLD,0:LASSLD,0:LASSLD),WG(LASSLD)
      DOUBLE PRECISION BRAVAIS(3,3),RECBV(3,3)
C     ..
C     .. Local Scalars ..
      INTEGER I,IEND,IPRINT,L,M,L1,L2,LMPOT
      INTEGER NCLEBD
      INTEGER NGMAX,NRMAX,NSHLG,NSHLR

      DOUBLE PRECISION PI,FPI

C     ..
C     .. Local Arrays ..
C     .. Attention: Dimension LMXSPD*LMPOTD appears sometimes as NCLEB1
      DOUBLE PRECISION CLEB(LMXSPD*LMPOTD)
      DOUBLE PRECISION GN(3,NMAXD),RM(3,NMAXD)
      DOUBLE PRECISION DFAC(0:LPOTD,0:LPOTD)
      INTEGER ICLEB(LMXSPD*LMPOTD,3)
      INTEGER NSG(ISHLD),NSR(ISHLD)
      INTEGER LOFLM(LMXSPD)
C     ..
C     .. External subroutines
      EXTERNAL LATTICE3D,MADELGAUNT
C ......................................................................
      IPRINT = 0
      NCLEBD = LMXSPD*LMPOTD
      PI = 4.0D0*ATAN(1.0D0)
      FPI = 4.0D0*PI
      LMPOT = (LPOT+1)**2
      I = 1
C
C --> determine the l-value for given lm
C
      DO L = 0,2*LPOT
         DO M = -L,L
            LOFLM(I) = L
            I = I + 1
         END DO
      END DO
C
C --> calculate:                             (2*(l+l')-1)!!
C                 dfac(l,l') = 4pi**2 *  ----------------------
C                                        (2*l+1)!! * (2*l'+1)!!
C
      DFAC(0,0) = FPI*FPI
      DO L1 = 1,LPOT
         DFAC(L1,0) = DFAC(L1-1,0)*DBLE(2*L1-1)/DBLE(2*L1+1)
         DFAC(0,L1) = DFAC(L1,0)
         DO L2 = 1,L1
            DFAC(L1,L2) = DFAC(L1,L2-1)*DBLE(2*(L1+L2)-1)/DBLE(2*L2+1)
            DFAC(L2,L1) = DFAC(L1,L2)
         END DO
      END DO
C
C --> calculate the gaunt coefficients
C
      CALL MADELGAUNT(LPOT,YRG,WG,CLEB,ICLEB,IEND,LASSLD,NCLEBD)
C
C ======================================================================
      CALL LATTICE3D(ALAT,BRAVAIS,RECBV,NGMAX,NRMAX,NSHLG,NSHLR,NSG,NSR,
     &               RMAX,GMAX,GN,RM,IPRINT,NMAXD,ISHLD)
C
C ======================================================================
C
      END
