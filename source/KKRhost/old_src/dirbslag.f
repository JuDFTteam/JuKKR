      SUBROUTINE DIRBSLAG(XI,Y1I,Y2I,Y3I,Y4I,Y1,Y2,Y3,Y4,IND1,N,IMAX)
C
C   ********************************************************************
C   *                                                                  *
C   *      lagrangian interpolation of Y(X) at position XI             *
C   *                                                                  *
C   *      XI      entry into x-array                                  *
C   *              for regular solution:   X(IND1-1) < XI <=X(IND1)    *
C   *              for irregular solution: X(IND1)   < XI <=X(IND1+1)  *
C   *      X/Y     X/Y-arrays                                          *
C   *      N       order of lagrangian interpolation                   *
C   *      IND     min-I for which  X(I) > XI                          *
C   *      IMAX    max index of X/Y-arrays                             *
C   *                                                                  *
C   ********************************************************************
C
 
      IMPLICIT NONE
      
C
C
C Dummy arguments
C
      INTEGER IMAX,IND1,N
      REAL*8 XI
      REAL*8 Y1I,Y2I,Y3I,Y4I
      REAL*8 Y1(IMAX),Y2(IMAX),Y3(IMAX),Y4(IMAX)
C
C Local variables
C
      REAL*8 D,P,XD
C      
      INTEGER I,IND,INL,INU,J
C
      IND = IND1
      IF ( ABS(XI-DBLE(IND)).LT.1.0D-12 ) THEN
         Y1I = Y1(IND)
         Y2I = Y2(IND)
         Y3I = Y3(IND)
         Y4I = Y4(IND)
         RETURN
      END IF
C ------------------------------------- shift IND for irregular solution
      IF ( XI.GT.DBLE(IND) ) IND = IND + 1
C
      INL = MAX(1,IND-(N+1)/2)
      INU = INL + N - 1
C
      IF ( INU.GT.IMAX ) THEN
         INL = IMAX - N + 1
         INU = IMAX
      END IF
C
      Y1I = 0.0D0
      Y2I = 0.0D0
      Y3I = 0.0D0
      Y4I = 0.0D0
      P = 1.0D0
      DO J = INL,INU
         P = P*(XI-DBLE(J))
         D = 1.0D0
         DO I = INL,INU
            IF ( I.NE.J ) THEN
               XD = DBLE(J)
            ELSE
               XD = XI
            END IF
            D = D*(XD-DBLE(I))
         END DO
C
         Y1I = Y1I + Y1(J)/D
         Y2I = Y2I + Y2(J)/D
         Y3I = Y3I + Y3(J)/D
         Y4I = Y4I + Y4(J)/D
C
      END DO
      Y1I = Y1I*P
      Y2I = Y2I*P
      Y3I = Y3I*P
      Y4I = Y4I*P
      END
