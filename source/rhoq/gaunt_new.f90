! ************************************************************************
      SUBROUTINE GAUNT_NEW(LMAX,LMAX2,LPOT,W,YR,CLEB,LOFLM,ICLEB,IEND,JEND,NCLEB,LMAXD,LMPOTD)
! ************************************************************************
!
!   - fills the array cleb with the gaunt coeffients ,i.e.
!      the integral of y(l1,m1)*y(l2,m2)*y(l3,m3)
!      but only for lm2.le.lm1 and lm3>1
!   - calculate the pointer array jend  to project the indices
!      array cleb with the same lm3,l1,l2 values - because of
!      the special ordering of array cleb only the last index
!      has to be determined .
!     (the parameter n has to be chosen that l1+l2+l3 .lt. 2*n)
!     using gaussian quadrature as given by
!     m. abramowitz and i.a. stegun, handbook of mathematical functions,
!     nbs applied mathematics series 55 (1968), pages 887 and 916
!     m. weinert and e. wimmer
!     northwestern university march 1980
!
!     an index array -icleb- is used to save storage place .
!     fills the array loflm which is used to determine the
!     l-value of a given lm-value .
!     this subroutine has to be called only once !
!
!                               b.drittler   november 1987
!
!     modified gaunt coefficients are als calculated defined by
!     the integral of y(l1,m1)*y(l2,m2)*y(l3,m3)*i**(l2-l1+l3)
!-----------------------------------------------------------------------
!---> attention : ncleb is an empirical factor - it has to be optimized
!
      use mod_rcstop, only: rcstop
      IMPLICIT NONE
!     ..
      DOUBLE COMPLEX CI
      PARAMETER (CI= (0.0D0,1.0D0))
!     ..
      INTEGER LMPOTD,LMAXD,NCLEB
!     .. Scalar Arguments ..
      INTEGER IEND,LMAX,LPOT, LMAX2
!     .. Array Arguments ..
      DOUBLE PRECISION CLEB(NCLEB),W(*),YR(4*LMAXD,0:4*LMAXD,0:4*LMAXD)
      INTEGER ICLEB(NCLEB,3),JEND(LMPOTD,0:LMAX,0:LMAX2),LOFLM((2*LPOT+1)**2)
!     .. Local Scalars ..
      DOUBLE PRECISION CLECG,FACTOR,FCI,S
      INTEGER I,J,L,L1,L1P,L2,L2P,L3,LM1,LM2,LM3,LM3P,LMPOT,M,M1,M1A,M1S,M2,M2A,M2S,M3,M3A,M3S
!     .. Intrinsi! Functions ..
      INTRINSIC ABS,DBLE,MOD,REAL,SIGN

      I = 1
      DO 20 L = 0,LPOT
        DO 10 M = -L,L
          LOFLM(I) = L
          I = I + 1
   10   CONTINUE
   20 CONTINUE

      ICLEB=0
      CLEB=0d0
      IF (LPOT.EQ.0) THEN
        IEND = 1
        ICLEB(1,1) = (LMAX+1)**2
        ICLEB(1,3) = 1
      END IF

      IF (LPOT.NE.0) THEN
!---> set up of the gaunt coefficients with an index field
        I = 1
        DO 90 L3 = 1,LPOT
          DO 80 M3 = -L3,L3
            DO 70 L1 = 0,LMAX
              DO 60 L2 = 0,LMAX2
                IF (MOD((L1+L2+L3),2).NE.1 .AND. (L1+L2-L3).GE.0 .AND. (L1-L2+L3).GE.0 .AND. (L2-L1+L3).GE.0) THEN
                  FCI = DBLE(CI** (L2-L1+L3))
                  DO 50 M1 = -L1,L1
                    DO 40 M2 = -L2,L2
!---> store only gaunt coeffients for lm2.le.lm1
                      LM1 = L1* (L1+1) + M1 + 1
                      LM2 = L2* (L2+1) + M2 + 1
                      IF (LM2.LE.LM1) THEN
                        M1S = SIGN(1,M1)
                        M2S = SIGN(1,M2)
                        M3S = SIGN(1,M3)
                        IF (M1S*M2S*M3S.GE.0) THEN
                          M1A = ABS(M1)
                          M2A = ABS(M2)
                          M3A = ABS(M3)
                          FACTOR = 0.0
                          IF (M1A+M2A.EQ.M3A) FACTOR = FACTOR + REAL(3*M3S+SIGN(1,-M3))/8.0d0
                          IF (M1A-M2A.EQ.M3A) FACTOR = FACTOR + REAL(M1S)/4.0d0
                          IF (M2A-M1A.EQ.M3A) FACTOR = FACTOR + REAL(M2S)/4.0d0
                          IF (FACTOR.NE.0.0) THEN
                            IF (M1S*M2S.NE.1 .OR. M2S*M3S.NE.1 .OR. M1S*M3S.NE.1) FACTOR = -FACTOR
                            S = 0.0
                            DO 30 J = 1,4*LMAXD
                              S = S + W(J)*YR(J,L1,M1A)*YR(J,L2,M2A)*YR(J,L3,M3A)
   30                       CONTINUE
                            CLECG = S*FACTOR
                            IF (ABS(CLECG).GT.1.D-10) THEN
                              CLEB(I) = CLECG
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

        ELSE ! (NCLEB.LT.IEND)

!---> set up of the pointer array jend,use explicitly
!     the ordering of the gaunt coeffients

          LMPOT = (LPOT+1)* (LPOT+1)
          DO 120 L1 = 0,LMAX
            DO 110 L2 = 0,LMAX2
              DO 100 LM3 = 2,LMPOT
                JEND(LM3,L1,L2) = 0
  100         CONTINUE
  110       CONTINUE
  120     CONTINUE

          LM3 = ICLEB(1,3)
          L1 = LOFLM(ICLEB(1,1))
          L2 = LOFLM(ICLEB(1,2))

          DO 130 J = 2,IEND
            LM3P = ICLEB(J,3)
            L1P = LOFLM(ICLEB(J,1))
            L2P = LOFLM(ICLEB(J,2))

            IF (LM3.NE.LM3P .OR. L1.NE.L1P .OR. L2.NE.L2P) THEN
              JEND(LM3,L1,L2) = J - 1
              LM3 = LM3P
              L1 = L1P
              L2 = L2P
            END IF

  130     CONTINUE
          JEND(LM3,L1,L2) = IEND

        END IF ! (NCLEB.LT.IEND)

      END IF ! (LPOT.NE.0)

 9000 FORMAT (13x,'error stop in gaunt : dimension of NCLEB = ',i10,' too small ',/,13x,'change NCLEB to ',i6)
      END
