      SUBROUTINE CALCCGC(LTAB,KAPTAB,NMUETAB,CGC,NKMAX,NMUEMAX,NKMPMAX)
C   ********************************************************************
C   *                                                                  *
C   *   CLEBSCH-GORDON-COEFFICIENTS     CGC(IKM,IS)                    *
C   *                                                                  *
C   *   IKM NUMBERS  CGC  FOR INCREASING  K  AND  MUE                  *
C   *   IKM  = L*2*(J+1/2) + J + MUE + 1                               *
C   *   IS= 1/2  SPIN DOWN/UP                                          *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C
C Dummy arguments
C
      INTEGER NKMAX,NKMPMAX,NMUEMAX
      REAL*8 CGC(NKMPMAX,2)
      INTEGER KAPTAB(NMUEMAX),LTAB(NMUEMAX),NMUETAB(NMUEMAX)
C
C Local variables
C
      INTEGER IKM,K,KAPPA,M
      REAL*8 J,L,MUE,TWOLP1
C
      IKM = 0
      DO K = 1,(NKMAX+1)
         L = LTAB(K)
         KAPPA = KAPTAB(K)
         J = ABS(KAPPA) - 0.5D0
         MUE = -J - 1.0D0
         TWOLP1 = 2.0D0*L + 1.0D0
C
         IF ( KAPPA.LT.0 ) THEN
C
C     J = L + 1/2
            DO M = 1,NMUETAB(K)
C
               MUE = MUE + 1.0D0
               IKM = IKM + 1
               CGC(IKM,1) = DSQRT((L-MUE+0.5D0)/TWOLP1)
               CGC(IKM,2) = DSQRT((L+MUE+0.5D0)/TWOLP1)
            END DO
         ELSE
C     J = L - 1/2
            DO M = 1,NMUETAB(K)
C
               MUE = MUE + 1.0D0
               IKM = IKM + 1
               CGC(IKM,1) = DSQRT((L+MUE+0.5D0)/TWOLP1)
               CGC(IKM,2) = -DSQRT((L-MUE+0.5D0)/TWOLP1)
C
            END DO
         END IF
C
C
      END DO
C
      END
