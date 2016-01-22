c***********************************************************************
      SUBROUTINE SORTVERTICES(N,SINFI,X,Y,Z)
c Sorts the array SINFI(N) in ascending order using straight insertion. 
c The arrays Z(N), Y(N), and Z(N) follow.
c On output, arrays SINFI, X, Y, and Z return sorted.
      implicit none
      INTEGER N,I,J
      REAL*8           SINFI(*),X(*),Y(*),Z(*),TMPS,TMPX,TMPY,TMPZ

      DO J = 2,N
         TMPS = SINFI(J)
         TMPX = X(J)
         TMPY = Y(J)
         TMPZ = Z(J)
         DO I = J-1,1,-1
            IF (SINFI(I).LE.TMPS) GOTO 10
            SINFI(I+1) = SINFI(I)
            X(I+1) = X(I)
            Y(I+1) = Y(I)
            Z(I+1) = Z(I)
         ENDDO
         I = 0
 10      SINFI(I+1) = TMPS
         X(I+1) = TMPX
         Y(I+1) = TMPY
         Z(I+1) = TMPZ
      ENDDO

      RETURN
      END
