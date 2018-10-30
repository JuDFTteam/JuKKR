c ************************************************************************
      SUBROUTINE DSORT (W,IND,MAX,POS)
      implicit none
c#@# KKRtags: VORONOI
c ************************************************************************
c     p.zahn, april 96
c     W   is the original array returned unchanged
c     IND is an array that holds the new positions 
c     max number of ellements to be sorted
c     pos the position where the first element is found
c ------------------------------------------------------------------------
      INTEGER MAX,POS
      REAL*8            W(*)
      REAL*8            BOUND, DIFF
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
            DIFF = DABS( W(IND(K)) - W(IND(K+J)) )
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



c ************************************************************************
      SUBROUTINE DSORT_NCOMP(W,NCOMP,LTEST,IND,MAX,POS)
      implicit none
c#@# KKRtags: VORONOI undefined
c ************************************************************************
!     Phivos Mavropoulos 2014
!     Sorting according to multiple components
c     W   is the original array returned unchanged
c     IND is an array that holds the new positions 
c     max number of ellements to be sorted
c     pos the position where the first element is found
!     NCOMP is the number of components:
!     First sort according to the 1st component, then according to the 2nd etc.
!     Uses integer function CMPR
c ------------------------------------------------------------------------
      INTEGER MAX,POS,NCOMP
      REAL*8  W(NCOMP,*)
      LOGICAL LTEST
      INTEGER IND(*)

      INTEGER I,II,J,JJ,K,DIFF_INT,CMPR
c ------------------------------------------------------------------------
      DO 10 I = 1,MAX
        IND(I) = I
 10   END DO

      DO I = 1,MAX
         IND(I) = 1
         DO J = 1,MAX
            DIFF_INT = CMPR( W(1,I),W(1,J),NCOMP )
            IF (DIFF_INT.GT.0) THEN 
               IND(I) = IND(I) + 1
            ENDIF
            IF (IND(I).GT.MAX) STOP 'Error 1 in DSORT_NCOMP'
         ENDDO
      ENDDO

      ! Test
      IF (LTEST) THEN
         DO I = 1,MAX
            JJ = 0
            DO J = 1,MAX
               IF (IND(I).EQ.J) JJ = JJ + 1
            ENDDO
            IF (JJ.NE.1) STOP 'Error 2 in DSORT_NCOMP'
         ENDDO
      ENDIF

      
      DO 50 I=1,MAX
        IF (IND(I) .EQ. 1) POS=I
 50   END DO

      RETURN
      END
c ************************************************************************


      FUNCTION CMPR(A1,A2,NCOMP)
c#@# KKRtags: VORONOI undefined
      IMPLICIT NONE
      ! Input
      INTEGER NCOMP
      REAL*8 A1(NCOMP),A2(NCOMP)
      !Output
      INTEGER CMPR
      ! Local
      REAL*8 TOL
      PARAMETER(TOL=0.D0)
      REAL*8 DIFF
      INTEGER ICOMP

      CMPR = 0
      
      DO ICOMP = 1,NCOMP
         DIFF = A1(ICOMP)-A2(ICOMP)
         IF (DIFF.GT.TOL) THEN
            CMPR = 1
            RETURN
         ENDIF
         IF (-DIFF.GT.TOL) THEN
            CMPR = -1
            RETURN
         ENDIF
      ENDDO
      

      END FUNCTION CMPR
