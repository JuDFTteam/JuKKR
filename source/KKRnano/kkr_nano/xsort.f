c ************************************************************************
      SUBROUTINE XSORT (W,IND,MAX,POS)
c ************************************************************************
c     p.zahn, april 96
c     W   is the original array returned unchanged
c     IND is an array that holds the new positions 
c     max number of ellements to be sorted
c     pos the position where the first element is found
c ------------------------------------------------------------------------
      INTEGER MAX,POS
      DOUBLE PRECISION  W(*)
      INTEGER IND(*)

      INTEGER I,II,J,JJ,K
      DATA BOUND /1.0D-12/
c ------------------------------------------------------------------------
      DO 10 I = 1,MAX
        IND(I) = I
 10   END DO

      J = MAX
      J = 1
      DO 60 WHILE (J.LT.MAX/3)
        J = 3*J+1
 60   END DO

      DO 20 WHILE (J.GT.1)
        J = J/3
        JJ = 1
        DO 30 WHILE (JJ.EQ.1)
          JJ = 0
          DO 40 K=1,MAX-J
            DIFF = ABS( W(IND(K)) - W(IND(K+J)) )
            IF ( W(IND(K)) .GT. W(IND(K+J)) .AND.
     +           DIFF.GT.BOUND ) THEN
              II       = IND(K)
              IND(K)   = IND(K+J)
              IND(K+J) = II
              JJ = 1
            END IF
 40       END DO                    ! K=1,MAX-J
 30     END DO                      ! WHILE (JJ.EQ.1)
 20   END DO
      
      DO 50 I=1,MAX
        IF (IND(I) .EQ. 1) POS=I
 50   END DO

      RETURN
      END


