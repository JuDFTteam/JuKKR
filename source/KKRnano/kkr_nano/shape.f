c ************************************************************************
      SUBROUTINE SHAPE(LPOT,NAEZ,GSH,ILM,IMAXSH,LMSP,NTCELL,W,YR)
c ************************************************************************
c   - prepares shape corrections
c     (the parameter n has to be chosen that l1+l2+l3 .le. 2*n)
c     using gaussian quadrature as given by
c     m. abramowitz and i.a. stegun, handbook of mathematical functions,
c     nbs applied mathematics series 55 (1968), pages 887 and 916
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.p'
      INTEGER N,LASSLD,LMPOTD,LMXSPD
      PARAMETER (N=4*LMAXD,LASSLD=N,LMPOTD= (LPOTD+1)**2)
      PARAMETER (LMXSPD= (2*LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER LPOT,NAEZ
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION GSH(*),W(N),YR(N,0:LASSLD,0:LASSLD)
      INTEGER ILM(NGSHD,3),IMAXSH(0:LMPOTD),LMSP(LMXSPD,NAEZD),NTCELL(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION FACTOR,GAUNT,S
      INTEGER I,IAT,ICELL,ISUM,J,L1,L2,L3,LM1,LM2,LM3,M1,M1A,M1S,M2,M2A,
     +        M2S,M3,M3A,M3S
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MOD,REAL,SIGN
C     ..
C     .. External Subroutines ..
      EXTERNAL RCSTOP
C     ..
c
c---> set up of the gaunt coefficients with an index field
c     so c(lm,lm',lm'') is mapped to c(i)
      I = 1
      DO 80 L1 = 0,LPOT
        DO 70 M1 = -L1,L1
          LM1 = L1*L1 + L1 + M1 + 1
          IMAXSH(LM1-1) = I - 1
          DO 60 L3 = 0,LPOT*2
            DO 50 M3 = -L3,L3
              LM3 = L3*L3 + L3 + M3 + 1
              ISUM = 0
              DO 10 IAT = 1,NAEZ
                ICELL = NTCELL(IAT)
                ISUM = ISUM + LMSP(LM3,ICELL)
   10         CONTINUE
              IF (ISUM.GT.0) THEN
                DO 40 L2 = 0,LPOT
                  IF (MOD((L1+L2+L3),2).NE.1 .AND. (L1+L2-L3).GE.0 .AND.
     +                (L1-L2+L3).GE.0 .AND. (L2-L1+L3).GE.0) THEN
                    DO 30 M2 = -L2,L2
                      LM2 = L2* (L2+1) + M2 + 1
c---> use the m-conditions for the gaunt coefficients not to be 0
                      M1S = SIGN(1,M1)
                      M2S = SIGN(1,M2)
                      M3S = SIGN(1,M3)
c
                      IF (M1S*M2S*M3S.GE.0) THEN
                        M1A = ABS(M1)
                        M2A = ABS(M2)
                        M3A = ABS(M3)
c
                        FACTOR = 0.0D0
c
                        IF (M1A+M2A.EQ.M3A) FACTOR = FACTOR +
     +                      REAL(3*M3S+SIGN(1,-M3))/8.0D0
                        IF (M1A-M2A.EQ.M3A) FACTOR = FACTOR +
     +                      REAL(M1S)/4.0D0
                        IF (M2A-M1A.EQ.M3A) FACTOR = FACTOR +
     +                      REAL(M2S)/4.0D0
c
                        IF (FACTOR.NE.0.0D0) THEN
c
                          IF (M1S*M2S.NE.1 .OR. M2S*M3S.NE.1 .OR.
     +                        M1S*M3S.NE.1) FACTOR = -FACTOR
                          S = 0.0D0
                          DO 20 J = 1,N
                            S = S + W(J)*YR(J,L1,M1A)*YR(J,L2,M2A)*
     +                          YR(J,L3,M3A)
   20                     CONTINUE
                          GAUNT = S*FACTOR
                          IF (ABS(GAUNT).GT.1.D-10) THEN
                            GSH(I) = GAUNT
                            ILM(I,1) = LM1
                            ILM(I,2) = LM2
                            ILM(I,3) = LM3
                            I = I + 1
                          END IF

                        END IF


                      END IF

   30               CONTINUE
                  END IF

   40           CONTINUE



              END IF

   50       CONTINUE

   60     CONTINUE

   70   CONTINUE

   80 CONTINUE
      IMAXSH(LM1) = I - 1

      WRITE (*,FMT=9000) IMAXSH(LM1),NGSHD
      IF (IMAXSH(LM1).GT.NGSHD) CALL RCSTOP('SHAPE   ')
c
 9000 FORMAT(' >>> SHAPE : IMAXSH(',I4,'),NGSHD :',2I6)
c
      END
