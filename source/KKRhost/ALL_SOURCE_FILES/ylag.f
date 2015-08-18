      FUNCTION YLAG(XI,X,Y,IND1,N1,IMAX)
C   ********************************************************************
C   *                                                                  *
C   * lagrangian interpolation                                         *
C   * xi is interpolated entry into x-array                            *
C   * n is the order of lagrangran interpolation                       *
C   * y is array from which ylag is obtained by interpolation          *
C   * ind is the min-i for x(i).gt.xi                                  *
C   * if ind=0,x-array will be searched                                *
C   * imax is max index of x-and y-arrays                              *
C   *                                                                  *
C   * 07/12/94  HE  arg. IEX removed                                   *
C   ********************************************************************
      IMPLICIT NONE
C
C
C Dummy arguments
C
      INTEGER IMAX,IND1,N1
      REAL*8 XI
      REAL*8 X(IMAX),Y(IMAX)
      REAL*8 YLAG
C
C Local variables
C
      REAL*8 D,P,S,XD
      INTEGER I,IND,INL,INU,J,N
      SAVE D,I,IND,INL,INU,J,N,P,S,XD
C
      IND = IND1
      N = N1
      IF ( N.GT.IMAX ) N = IMAX
      IF ( IND.GT.0 ) GOTO 200
      DO J = 1,IMAX
         IF ( ABS(XI-X(J)).LT.1.0D-12 ) GOTO 600
C      IF( ABS(XI-X(J)).LT.1.0e-06) GO TO 130
         IF ( XI.LT.X(J) ) GOTO 100
         IF ( XI.EQ.X(J) ) GOTO 600
      END DO
      GOTO 300
 100  CONTINUE
      IND = J
 200  CONTINUE
      IF ( IND.GT.1 ) THEN
      END IF
      INL = IND - (N+1)/2
      IF ( INL.LE.0 ) INL = 1
      INU = INL + N - 1
      IF ( INU.LE.IMAX ) GOTO 400
 300  CONTINUE
      INL = IMAX - N + 1
      INU = IMAX
 400  CONTINUE
      S = 0.0D0
      P = 1.0D0
      DO J = INL,INU
         P = P*(XI-X(J))
         D = 1.0D0
         DO I = INL,INU
            IF ( I.NE.J ) THEN
               XD = X(J)
            ELSE
               XD = XI
            END IF
            D = D*(XD-X(I))
         END DO
         S = S + Y(J)/D
      END DO
      YLAG = S*P
 500  CONTINUE
      RETURN
 600  CONTINUE
      YLAG = Y(J)
      GOTO 500
      END
C
C
