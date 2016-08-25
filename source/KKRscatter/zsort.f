c ************************************************************************
      SUBROUTINE ZSORT (W,IND,MAX,POS)
c ************************************************************************
c     p.zahn, april 96
c     modified quick sort algorithm,
c     sorts due to the real part of the double complex array W
c ------------------------------------------------------------------------
      INTEGER MAX,POS
      DOUBLE COMPLEX W(*)
      INTEGER IND(*)

      INTEGER I,II,J,JJ,K
      data bound /1.0d-12/
c ------------------------------------------------------------------------
      DO 10 I=1,MAX
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
            diff=abs(DREAL(W(IND(K)))-DREAL(W(IND(K+j))))
            IF ( DREAL(W(IND(K))) .GT. DREAL(W(IND(K+J))) ) THEN

c            IF ( (DREAL(W(IND(K))).GT.DREAL(W(IND(K+J))) .and.
c     +           diff.gt.bound ) .or.
c     +           ( diff.le.bound .and. ind(k).gt.ind(k+j) ) ) THEN

              II       = IND(K)
              IND(K)   = IND(K+J)
              IND(K+J) = II
              JJ = 1
            END IF
 40       END DO                    ! K=1,MAX-J
 30     END DO                      ! WHILE (JJ.EQ.1)
 20   END DO
      
      DO 50 I=1,MAX
        if (IND(I) .eq. 1) pos=i
 50   END DO

      RETURN
      END
