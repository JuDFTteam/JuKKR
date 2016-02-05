c ************************************************************************
      SUBROUTINE GLLCOPYA(GLLKE0,M,NSPBLOCK,J,D1,D2)
c ************************************************************************
C     .. Parameters ..
      include 'inc.p'
      INTEGER LMAXSQ
      PARAMETER (LMAXSQ= (LMAXD+1)**2)
      INTEGER ALM
      PARAMETER (ALM=LMAXSQ*NAEZD)
c
c     .. arguments
      INTEGER NSPBLOCK
      DOUBLE COMPLEX GLLKE0(ALM,*),
     +               M(LMAXSQ*NSPBLOCK,*)
      INTEGER J,D1,D2
C
c     .. Local
      INTEGER D1LM,D2LM,JLM,LM1,LM2
c ------------------------------------------------------------------------
      D1LM = (D1-1)*LMAXSQ
      D2LM = (D2-1)*LMAXSQ
      JLM  = (J-1)*LMAXSQ

      DO 10 LM2 = 1,LMAXSQ
        DO 20 LM1 = 1,LMAXSQ
          M(D1LM+LM1,D2LM+LM2) = GLLKE0(JLM+LM1,LM2)
 20     END DO
 10   END DO
      RETURN
      END
