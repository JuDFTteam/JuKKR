C*==madelcoef.f    processed by SPAG 6.05Rc at 12:05 on 17 May 2004
      SUBROUTINE MADELCOEF(LINTERFACE,LPOT,A,B,SMAT,CLEB,ICLEB,IEND,
     &                     LPOTD,LMPOTD,LMXSPD,NCLEBD)
      IMPLICIT NONE
C     ..
C     .. Scalar arguments
      INTEGER LPOT,IEND,LPOTD,LMPOTD,LMXSPD,NCLEBD
      LOGICAL LINTERFACE
C     ..
C     .. Array arguments
      DOUBLE PRECISION A(LMPOTD,LMPOTD),B(LMPOTD)
      DOUBLE PRECISION SMAT(LMXSPD),CLEB(NCLEBD)
      INTEGER ICLEB(NCLEBD,3)
C     ..
C     .. Local scalars
      DOUBLE PRECISION PI,FPI
      INTEGER I,L,L1,L2,LM1,LM2,LM3,LMPOT,LOFLM(LMXSPD),M
      INTEGER ICALL
C     ..
C     .. Local arrays
      DOUBLE PRECISION DFAC(0:LPOTD,0:LPOTD)
C     ..
C     .. Data statements
      DATA ICALL /0/
C     ..
C     .. Intrinsic functions
      INTRINSIC ABS,DBLE
C     ..
C     .. Save statements
      SAVE ICALL,LMPOT,PI,FPI
C     ..................................................................
      ICALL = ICALL + 1
C
      IF ( ICALL.EQ.1 ) THEN
         PI = 4.0D0*ATAN(1.0D0)
         FPI = 4.0D0*PI
         LMPOT = (LPOT+1)**2
      END IF
C
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
C --> initialize
C
      DO LM1 = 1,LMPOT
         DO LM2 = 1,LMPOT
            A(LM1,LM2) = 0.0D0
         END DO
      END DO
C
C --> calculate a(lm1,lm2)
C
      DO I = 1,IEND
         LM1 = ICLEB(I,1)
         LM2 = ICLEB(I,2)
         LM3 = ICLEB(I,3)
         L1 = LOFLM(LM1)
         L2 = LOFLM(LM2)
C
C --> this loop has to be calculated only for l1+l2=l3
C
         A(LM1,LM2) = A(LM1,LM2) + 2.0D0*DFAC(L1,L2)*SMAT(LM3)*CLEB(I)
      END DO
C
      IF ( LINTERFACE ) RETURN
C
C --> initialize
C
      DO LM1 = 1,LMPOT
         B(LM1) = 0.0D0
      END DO
C
C --> calculate b(lm1)
C
      DO LM1 = 1,LMPOT
         L1 = LOFLM(LM1)
         B(LM1) = B(LM1) - 2.0D0*FPI/DBLE(2*L1+1)*SMAT(LM1)
      END DO
      END
