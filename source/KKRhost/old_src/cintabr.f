      SUBROUTINE CINTABR(AG,BG,AGBG,AF,BF,AFBF,RPW,NKA,NKB,JTOP,NRMAX)
C   ********************************************************************
C   *                                                                  *
C   *  SIMPSON - INTERGRATION FOR COMPLEX INTEGRAND  FX FROM 1 TO JTOP *
C   *  AND EQUIDISTANT MESH    I                                       *
C   *   INT = [ F1 + 4*F2 + 2*F3 + .... + 4F(N-1) + FN ]/3     N:ODD   *
C   *                                                                  *
C   *            FX = AG*AG*RPW   and   FX = AF*AF*RPW                 *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C
C Dummy arguments
C
      INTEGER JTOP,NKA,NKB,NRMAX
      COMPLEX*16 AF(NRMAX,2),AFBF(2,2),AG(NRMAX,2),AGBG(2,2),BF(NRMAX,2)
     &           ,BG(NRMAX,2)
      REAL*8 RPW(NRMAX)
C
C Local variables
C
      REAL*8 F,SIMP
      INTEGER I,KA,KB
C
      DO KB = 1,NKB
         DO KA = 1,NKA
            AGBG(KA,KB) = AG(1,KA)*BG(1,KB)*RPW(1)
            AFBF(KA,KB) = AF(1,KA)*BF(1,KB)*RPW(1)
         END DO
      END DO
C
      IF ( MOD(JTOP,2).EQ.0 ) STOP '<CINTABR>  JTOP is even !!!'
C
      SIMP = -1.0D0
      DO I = 2,JTOP - 1
         SIMP = -SIMP
         F = (3.0D0+SIMP)*RPW(I)
         DO KB = 1,NKB
            DO KA = 1,NKA
               AGBG(KA,KB) = AGBG(KA,KB) + AG(I,KA)*BG(I,KB)*F
               AFBF(KA,KB) = AFBF(KA,KB) + AF(I,KA)*BF(I,KB)*F
            END DO
         END DO
      END DO
C
      DO KB = 1,NKB
         DO KA = 1,NKA
            AGBG(KA,KB) = (AGBG(KA,KB)+AG(JTOP,KA)*BG(JTOP,KB)*RPW(JTOP)
     &                    )/3.0D0
            AFBF(KA,KB) = (AFBF(KA,KB)+AF(JTOP,KA)*BF(JTOP,KB)*RPW(JTOP)
     &                    )/3.0D0
         END DO
      END DO
C
      END
