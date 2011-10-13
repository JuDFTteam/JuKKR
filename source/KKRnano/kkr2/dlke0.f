c **********************************************************************
      SUBROUTINE DLKE0(I,GLLH,EIKRM,EIKRP,IC,NACLS,
     +                  ATOM,NUMN0,INDN0,GINP)
c **********************************************************************
      IMPLICIT NONE
C     .. Parameters ..
      INCLUDE 'inc.p'
      INCLUDE 'inc.cls'
C
      INTEGER LMGF0D
      PARAMETER (LMGF0D= (LMAXD+1)**2)
      INTEGER ALMGF0
      PARAMETER (ALMGF0=LMGF0D*NAEZD)
      INTEGER NGTBD
      PARAMETER (NGTBD = NACLSD*LMGF0D)
C     ..
      DOUBLE COMPLEX GINP(LMGF0D,LMGF0D,NACLSD)
      DOUBLE COMPLEX GLLH(LMGF0D,NGTBD,*)
      DOUBLE COMPLEX EIKRM(NACLSD),EIKRP(NACLSD)
      INTEGER ATOM(NACLSD),NACLS(*),INDN0(NAEZD,*),NUMN0(*)
C     ..
      INTEGER AM,AN,I,IC,IN,J,LM1,LM2,M,N1,N2,IND1,IND2
C     ..
c ----------------------------------------------------------------------

      DO M = 1,NACLS(IC)
        DO N1 = 1,NUMN0(I)
        IND1 = INDN0(I,N1)
        IF(ATOM(M).EQ.IND1) THEN
           AM = (N1-1)*LMGF0D
           DO LM1 = 1,LMGF0D
           DO LM2 = 1,LMGF0D
             GLLH(LM1,AM+LM2,I) = GLLH(LM1,AM+LM2,I)
     +                         + EIKRM(M)*GINP(LM2,LM1,M)
           ENDDO
           ENDDO
        ENDIF
        ENDDO
        J = ATOM(M)
        DO N2 = 1,NUMN0(J)
        IND2 = INDN0(J,N2)
        IF(I.EQ.IND2) THEN
           AN = (N2-1)*LMGF0D
           DO LM1 = 1,LMGF0D
           DO LM2 = 1,LMGF0D
             GLLH(LM1,AN+LM2,J) = GLLH(LM1,AN+LM2,J)
     +                         + EIKRP(M)*GINP(LM1,LM2,M)
           ENDDO
           ENDDO
        ENDIF
        ENDDO
      ENDDO

      RETURN

      END
