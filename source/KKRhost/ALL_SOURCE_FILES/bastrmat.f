C*==bastrmat.f    processed by SPAG 6.05Rc at 11:40 on 21 Jul 2001
      SUBROUTINE BASTRMAT(LMAX,CGC,RC,CREL,RREL,NKMMAX,NKMPMAX)
C   ********************************************************************
C   *                                                                  *
C   *    INITIALIZE TRANSFORMATION MATRIX THAT TAKES MATRICES FROM     *
C   *    RELATIVISTIC  TO  REAL SPERICAL HARM.  REPRESENTATION         *
C   *                                                                  *
C   *    this is a special version of <STRSMAT> passing the            *
C   *    full BASis TRansformation MATrices  RC, CREL and RREL         *
C   *                                                                  *
C   * 13/01/98  HE                                                     *
C   ********************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      COMPLEX*16 CI,C1,C0
      PARAMETER (CI=(0.0D0,1.0D0),C1=(1.0D0,0.0D0),C0=(0.0D0,0.0D0))
C
C Dummy arguments
C
      INTEGER LMAX,NKMMAX,NKMPMAX
      REAL*8 CGC(NKMPMAX,2)
      COMPLEX*16 CREL(NKMMAX,NKMMAX),RC(NKMMAX,NKMMAX),
     &           RREL(NKMMAX,NKMMAX)
C
C Local variables
C
      INTEGER I,IKM,J,JP05,K,L,LM,LNR,M,MUEM05,MUEP05,NK,NKM,NLM
      REAL*8 W
C
C*** End of declarations rewritten by SPAG
C
      NK = 2*(LMAX+1) + 1
      NLM = (LMAX+1)**2
      NKM = 2*NLM
C     ===================================================
C     INDEXING:
C     IKM  = L*2*(J+1/2) + J + MUE + 1
C     LM   = L*(L+1)     +     M   + 1
C     ===================================================
C
C ----------------------------------------------------------------------
C CREL  transforms from  COMPLEX (L,M,S)  to  (KAP,MUE) - representation
C                 |LAM> = sum[LC] |LC> * CREL(LC,LAM)
C ----------------------------------------------------------------------
      CALL CINIT(NKMMAX*NKMMAX,CREL)
C
      LM = 0
      DO LNR = 0,LMAX
         DO M = -LNR,LNR
            LM = LM + 1
C
            IKM = 0
            DO K = 1,NK
               L = K/2
               IF ( 2*L.EQ.K ) THEN
                  JP05 = L
               ELSE
                  JP05 = L + 1
               END IF
C
               DO MUEM05 = -JP05,(JP05-1)
                  MUEP05 = MUEM05 + 1
                  IKM = IKM + 1
C
                  IF ( L.EQ.LNR ) THEN
                     IF ( MUEP05.EQ.M ) CREL(LM,IKM) = CGC(IKM,1)
                     IF ( MUEM05.EQ.M ) CREL(LM+NLM,IKM) = CGC(IKM,2)
                  END IF
C
               END DO
            END DO
C
         END DO
      END DO
C
C ----------------------------------------------------------------------
C    RC  transforms from  REAL to  COMPLEX (L,M,S) - representation
C                 |LC> = sum[LR] |LR> * RC(LR,LC)
C ----------------------------------------------------------------------
      CALL CINIT(NKMMAX*NKMMAX,RC)
C
      W = 1.0D0/SQRT(2.0D0)
C
      DO L = 0,LMAX
         DO M = -L,L
            I = L*(L+1) + M + 1
            J = L*(L+1) - M + 1
C
            IF ( M.LT.0 ) THEN
               RC(I,I) = -CI*W
               RC(J,I) = W
               RC(I+NLM,I+NLM) = -CI*W
               RC(J+NLM,I+NLM) = W
            END IF
            IF ( M.EQ.0 ) THEN
               RC(I,I) = C1
               RC(I+NLM,I+NLM) = C1
            END IF
            IF ( M.GT.0 ) THEN
               RC(I,I) = W*(-1.0D0)**M
               RC(J,I) = CI*W*(-1.0D0)**M
               RC(I+NLM,I+NLM) = W*(-1.0D0)**M
               RC(J+NLM,I+NLM) = CI*W*(-1.0D0)**M
            END IF
         END DO
      END DO
C
C ----------------------------------------------------------------------
C RREL  transforms from   REAL (L,M,S)  to  (KAP,MUE) - representation
C                 |LAM> = sum[LR] |LR> * RREL(LR,LAM)
C ----------------------------------------------------------------------
C
      CALL ZGEMM('N','N',NKM,NKM,NKM,C1,RC,NKMMAX,CREL,NKMMAX,C0,RREL,
     &           NKMMAX)
C
      END

