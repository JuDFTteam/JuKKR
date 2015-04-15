      SUBROUTINE SHAPE(LPOT,NATYP,GSH,ILM,IMAXSH,LMSP,NTCELL,W,YR,
     &                 LASSLD,LMPOTD,NATYPD,NGSHD)
C **********************************************************************
C *  Prepares shape corrections using gaussian quadrature as given by  *
C *  m. abramowitz and i.a. stegun, handbook of mathematical functions *
C *  nbs applied mathematics series 55 (1968), pages 887 and 916       *
C *                                                                    *
C *  the parameter LASSLD has to be chosen such that                   *
C *                        l1+l2+l3 .le. 2*LASSLD                      *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER LASSLD,LMPOTD,NATYPD,NGSHD
      INTEGER LPOT,NATYP
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION GSH(*),W(LASSLD),YR(LASSLD,0:LASSLD,0:LASSLD)
      INTEGER ILM(NGSHD,3),IMAXSH(0:LMPOTD),LMSP(NATYPD,*),NTCELL(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION FACTOR,GAUNT,S
      INTEGER I,IAT,ICELL,ISUM,J,L1,L2,L3,LM1,LM2,LM3,M1,M1A,M1S,M2,M2A,
     +        M2S,M3,M3A,M3S
      LOGICAL TRIANGLE
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,SIGN
C     ..
C     .. External Subroutines ..
      EXTERNAL RCSTOP,TRIANGLE
C     ..
C
C -> set up of the gaunt coefficients with an index field
C    so that  c(lm,lm',lm'') is mapped to c(i)
      I = 1
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
      DO L1 = 0,LPOT
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
         DO M1 = -L1,L1
C
            LM1 = L1*L1 + L1 + M1 + 1
            IMAXSH(LM1-1) = I - 1
C llllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll
            DO L3 = 0,LPOT*2
C mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
               DO M3 = -L3,L3
C
                  LM3 = L3*L3 + L3 + M3 + 1
                  ISUM = 0
C     
                  DO IAT = 1,NATYP
                     ICELL = NTCELL(IAT)
                     ISUM = ISUM + LMSP(ICELL,LM3)
                  END DO
C     
C ======================================================================
                  IF ( ISUM.GT.0 ) THEN
                     DO L2 = 0,LPOT
C ----------------------------------------------------------------------
                        IF ( TRIANGLE(L1,L2,L3) ) THEN
                           DO M2 = -L2,L2
C
                              LM2 = L2*L2 + L2 + M2 + 1
C
C -> use the m-conditions for the gaunt coefficients not to be 0
C
                              M1S = SIGN(1,M1)
                              M2S = SIGN(1,M2)
                              M3S = SIGN(1,M3)
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                              IF ( M1S*M2S*M3S.GE.0 ) THEN
                                 M1A = ABS(M1)
                                 M2A = ABS(M2)
                                 M3A = ABS(M3)
                                 FACTOR = 0.0D0
C
                                 IF (M1A+M2A.EQ.M3A) FACTOR = FACTOR +
     +                                   DBLE(3*M3S+SIGN(1,-M3))/8.0D0
C
                                 IF (M1A-M2A.EQ.M3A) FACTOR = FACTOR +
     +                                   DBLE(M1S)/4.0D0
C
                                 IF (M2A-M1A.EQ.M3A) FACTOR = FACTOR +
     +                                   DBLE(M2S)/4.0D0
C ......................................................................
                                 IF (FACTOR.NE.0.0D0) THEN
C
                                    IF ( M1S*M2S.NE.1 .OR. M2S*M3S.NE.1 
     &                               .OR.M1S*M3S.NE.1 ) FACTOR = -FACTOR
C
                                    S = 0.0D0
                                    DO J = 1,LASSLD
                                       S = S + W(J) * YR(J,L1,M1A)
     &                                     * YR(J,L2,M2A) * YR(J,L3,M3A)
                                    END DO
C
                                    GAUNT = S*FACTOR
                                    IF ( ABS(GAUNT).GT.1D-10 ) THEN
                                       GSH(I) = GAUNT
                                       ILM(I,1) = LM1
                                       ILM(I,2) = LM2
                                       ILM(I,3) = LM3
                                       I = I + 1
                                    END IF
                                 END IF
C ......................................................................
                              END IF
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                           END DO
                        END IF
C ----------------------------------------------------------------------
                     END DO
                  END IF
C ======================================================================
               END DO
C mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
            END DO
C llllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll
         END DO
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
      END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
      IMAXSH(LM1) = I - 1
      WRITE (*,FMT=9000) IMAXSH(LM1),NGSHD
      IF ( IMAXSH(LM1).GT.NGSHD ) CALL RCSTOP('SHAPE   ')
C
 9000 FORMAT(' >>> SHAPE : IMAXSH(',I4,'),NGSHD :',2I6)
C
      END
C
      FUNCTION TRIANGLE(L1,L2,L3)
      IMPLICIT NONE
      INTEGER L1,L2,L3
      LOGICAL TRIANGLE
      INTRINSIC MOD
C     ..
      TRIANGLE = (L1.GE.ABS(L3-L2)) .AND. (L1.LE.(L3+L2))
     &     .AND. (MOD((L1+L2+L3),2).EQ.0)
      END
