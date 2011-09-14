C ************************************************************************
      SUBROUTINE GAUNT(LMAX,LPOT,W,YR,CLEB,LOFLM,ICLEB,IEND,JEND,
     &                 NCLEB)
C     NCLEB is a new argument after inc.p remove
C ************************************************************************
c
c   - fills the array cleb with the gaunt coeffients ,i.e.
c      the integral of y(l1,m1)*y(l2,m2)*y(l3,m3)
c      but only for lm2.le.lm1 and lm3>1
c   - calculate the pointer array jend  to project the indices
c      array cleb with the same lm3,l1,l2 values - because of
c      the special ordering of array cleb only the last index
c      has to be determined .
c     (the parameter n has to be chosen that l1+l2+l3 .lt. 2*n)
c     using gaussian quadrature as given by
c     m. abramowitz and i.a. stegun, handbook of mathematical functions,
c     nbs applied mathematics series 55 (1968), pages 887 and 916
c     m. weinert and e. wimmer
c     northwestern university march 1980
c
c     an index array -icleb- is used to save storage place .
c     fills the array loflm which is used to determine the
c     l-value of a given lm-value .
c     this subroutine has to be called only once !
c
c                               b.drittler   november 1987
c
c     modified gaunt coefficients are als calculated defined by
c     the integral of y(l1,m1)*y(l2,m2)*y(l3,m3)*i**(l2-l1+l3)
c-----------------------------------------------------------------------
c
c---> attention : ncleb is an empirical factor - it has to be optimized
c
C     .. Parameters ..
c      include 'inc.p'
C     ..

      IMPLICIT NONE

      DOUBLE COMPLEX CI
      PARAMETER (CI= (0.0D0,1.0D0))
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT(OUT) :: IEND
      INTEGER, INTENT(IN) :: LMAX
      INTEGER, INTENT(IN) :: LPOT

      INTEGER, INTENT(IN) :: NCLEB
C     ..
C     .. Array Arguments ..

C      DOUBLE PRECISION CLEB(NCLEB,2),W(*),YR(N,0:N,0:N)
C      INTEGER ICLEB(NCLEB,3),JEND(LMPOTD,0:LMAXD,0:LMAXD),LOFLM(*)
      DOUBLE PRECISION CLEB(NCLEB,2)
      DOUBLE PRECISION W(*)
      DOUBLE PRECISION YR(4*LMAX,0:4*LMAX,0:4*LMAX)
C     INTEGER ICLEB(NCLEB,3)
      INTEGER ICLEB(NCLEB,3)
C      INTEGER JEND(LMPOTD,0:LMAXD,0:LMAXD)
      INTEGER JEND((LPOT+1)**2,0:LMAX,0:LMAX)
      INTEGER LOFLM(*)

C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CLECG,FACTOR,FCI,S
      INTEGER I,J,L,L1,L1P,L2,L2P,L3,LM1,LM2,LM3,LM3P,LMPOT,M,M1,M1A,
     +        M1S,M2,M2A,M2S,M3,M3A,M3S
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,MOD,REAL,SIGN
C     ..
C     .. External Subroutines ..
      EXTERNAL RCSTOP
C     ..

      INTEGER N

      N=4*LMAX

      I = 1
      DO 20 L = 0,2*LMAX
        DO 10 M = -L,L
          LOFLM(I) = L
          I = I + 1
   10   CONTINUE
   20 CONTINUE
c
      IF (LPOT.EQ.0) THEN
        IEND = 1
        ICLEB(1,1) = (LMAX+1)**2
        ICLEB(1,3) = 1
      END IF
c
      IF (LPOT.NE.0) THEN
c
c---> set up of the gaunt coefficients with an index field
c
        I = 1
        DO 90 L3 = 1,LPOT
          DO 80 M3 = -L3,L3
c
            DO 70 L1 = 0,LMAX
              DO 60 L2 = 0,L1
c
                IF (MOD((L1+L2+L3),2).NE.1 .AND. (L1+L2-L3).GE.0 .AND.
     +              (L1-L2+L3).GE.0 .AND. (L2-L1+L3).GE.0) THEN
c
                  FCI = DBLE(CI** (L2-L1+L3))
                  DO 50 M1 = -L1,L1
                    DO 40 M2 = -L2,L2
c
c---> store only gaunt coeffients for lm2.le.lm1
c
                      LM1 = L1* (L1+1) + M1 + 1
                      LM2 = L2* (L2+1) + M2 + 1
                      IF (LM2.LE.LM1) THEN
c
                        M1S = SIGN(1,M1)
                        M2S = SIGN(1,M2)
                        M3S = SIGN(1,M3)
c
                        IF (M1S*M2S*M3S.GE.0) THEN
c
                          M1A = ABS(M1)
                          M2A = ABS(M2)
                          M3A = ABS(M3)
c
                          FACTOR = 0.0
c
                          IF (M1A+M2A.EQ.M3A) FACTOR = FACTOR +
     +                        REAL(3*M3S+SIGN(1,-M3))/8.0d0
                          IF (M1A-M2A.EQ.M3A) FACTOR = FACTOR +
     +                        REAL(M1S)/4.0d0
                          IF (M2A-M1A.EQ.M3A) FACTOR = FACTOR +
     +                        REAL(M2S)/4.0d0
c
                          IF (FACTOR.NE.0.0) THEN
c
                            IF (M1S*M2S.NE.1 .OR. M2S*M3S.NE.1 .OR.
     +                          M1S*M3S.NE.1) FACTOR = -FACTOR
c
                            S = 0.0
                            DO 30 J = 1,N
                              S = S + W(J)*YR(J,L1,M1A)*YR(J,L2,M2A)*
     +                            YR(J,L3,M3A)
   30                       CONTINUE
                            CLECG = S*FACTOR
                            IF (ABS(CLECG).GT.1.D-10) THEN
                              CLEB(I,1) = CLECG
                              CLEB(I,2) = FCI*CLECG
                              ICLEB(I,1) = LM1
                              ICLEB(I,2) = LM2
                              ICLEB(I,3) = L3* (L3+1) + M3 + 1
                              I = I + 1
                            END IF

                          END IF

                        END IF

                      END IF

   40               CONTINUE
   50             CONTINUE
                END IF

   60         CONTINUE
   70       CONTINUE
   80     CONTINUE
   90   CONTINUE
        IEND = I - 1
        IF (NCLEB.LT.IEND) THEN
          WRITE (6,FMT=9000) NCLEB,IEND
          CALL RCSTOP('33      ')

        ELSE
c
c---> set up of the pointer array jend,use explicitly
c     the ordering of the gaunt coeffients
c
          LMPOT = (LPOT+1)* (LPOT+1)
          DO 120 L1 = 0,LMAX
            DO 110 L2 = 0,L1
              DO 100 LM3 = 2,LMPOT
                JEND(LM3,L1,L2) = 0
  100         CONTINUE
  110       CONTINUE
  120     CONTINUE
c
          LM3 = ICLEB(1,3)
          L1 = LOFLM(ICLEB(1,1))
          L2 = LOFLM(ICLEB(1,2))
c
          DO 130 J = 2,IEND
            LM3P = ICLEB(J,3)
            L1P = LOFLM(ICLEB(J,1))
            L2P = LOFLM(ICLEB(J,2))
c
            IF (LM3.NE.LM3P .OR. L1.NE.L1P .OR. L2.NE.L2P) THEN
              JEND(LM3,L1,L2) = J - 1
              LM3 = LM3P
              L1 = L1P
              L2 = L2P
            END IF

  130     CONTINUE
          JEND(LM3,L1,L2) = IEND
c
c
        END IF

      END IF
c
c

 9000 FORMAT (13x,'error stop in gaunt : dimension of NCLEB = ',i10,
     +       ' too small ',/,
     +       13x,'change NCLEB to ',i6)
      END
