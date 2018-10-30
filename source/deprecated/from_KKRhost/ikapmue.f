      FUNCTION IKAPMUE(KAPPA,MUEM05)
C
C   ********************************************************************
C   *                                                                  *
C   *  INDEXING OF MATRIX-ELEMENTS:                                    *
C   *                                                                  *
C   *  I = 2*L*(J+1/2) + J + MUE + 1                                   *
C   *                                                                  *
C   ********************************************************************
C

      IMPLICIT NONE

C
C Dummy arguments
C
      INTEGER KAPPA,MUEM05
      INTEGER IKAPMUE
C
C Local variables
C
      INTEGER IABS
      INTEGER JP05,L
C
      JP05 = IABS(KAPPA)
C
      IF ( KAPPA.LT.0 ) THEN
         L = -KAPPA - 1
      ELSE
         L = KAPPA
      END IF
C
      IKAPMUE = 2*L*JP05 + JP05 + MUEM05 + 1
C
      END
