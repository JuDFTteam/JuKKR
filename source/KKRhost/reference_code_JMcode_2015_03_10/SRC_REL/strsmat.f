      SUBROUTINE STRSMAT(LMAX,CGC,SRREL,NRREL,IRREL,NKMMAX,NKMPMAX)
C   ********************************************************************
C   *                                                                  *
C   *    INITIALIZE TRANSFORMATION MATRIX THAT TAKES MATRICES FROM     *
C   *    RELATIVISTIC  TO  REAL SPERICAL HARM.  REPRESENTATION         *
C   *                                                                  *
C   *    ONLY THE NON-0 ELEMENTS OF THE MATRIX ARE STORED              *
C   *                                                                  *
C   * 25/10/95  HE  proper convention of trans. matrix introduced      *
C   ********************************************************************
C
      IMPLICIT NONE
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
      INTEGER IRREL(2,2,NKMMAX),NRREL(2,NKMMAX)
      COMPLEX*16 SRREL(2,2,NKMMAX)
C
C Local variables
C
      COMPLEX*16 CREL(NKMMAX,NKMMAX),RC(NKMMAX,NKMMAX),
     &           RREL(NKMMAX,NKMMAX)
      INTEGER I,IKM,J,JP05,K,L,LAM,LM,LNR,LR,M,MUEM05,MUEP05,NK,NKM,NLM,
     &        NS1,NS2
      REAL*8 W
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
      CALL ZGEMM('N','N',NKM,NKM,NKM,C1,RC,NKMMAX,CREL,NKMMAX,C0,RREL,
     &              NKMMAX)
C
C     ---------------------------------------------------
C     store the elements of  RREL
C     ---------------------------------------------------
      DO LAM = 1,NKM
         NS1 = 0
         NS2 = 0
C
         DO LR = 1,2*NLM
            IF ( CDABS(RREL(LR,LAM)).GT.1D-6 ) THEN
               IF ( LR.LE.NLM ) THEN
                  NS1 = NS1 + 1
                  IF ( NS1.GT.2 ) STOP ' IN <STRSMAT>   NS1 > 2'
                  SRREL(NS1,1,LAM) = RREL(LR,LAM)
                  IRREL(NS1,1,LAM) = LR
               ELSE
                  NS2 = NS2 + 1
                  IF ( NS2.GT.2 ) STOP ' IN <STRSMAT>   NS2 > 2'
                  SRREL(NS2,2,LAM) = RREL(LR,LAM)
                  IRREL(NS2,2,LAM) = LR - NLM
               END IF
            END IF
         END DO
C
         NRREL(1,LAM) = NS1
         NRREL(2,LAM) = NS2
      END DO
C
      END
