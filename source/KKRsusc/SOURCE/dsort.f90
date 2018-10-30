MODULE MOD_DSORT
CONTAINS

! ************************************************************************
      SUBROUTINE DSORT (W,IND,IMAX,POS)
! ************************************************************************
!     p.zahn, april 96
!     W   is the original array returned unchanged
!     IND is an array that holds the new positions 
!     max number of ellements to be sorted
!     pos the position where the first element is found
! ------------------------------------------------------------------------
      INTEGER           :: IMAX,POS
      DOUBLE PRECISION  :: W(:)
      INTEGER           :: IND(:)

      INTEGER           :: I,II,J,JJ,K
      INTEGER,PARAMETER :: BOUND =1.0D-12
! ------------------------------------------------------------------------

    IF(imax>5) then

      DO I = 1,IMAX
        IND(I) = I
      END DO

      J = IMAX
      J = 1
      DO WHILE (J.LT.IMAX/3)
        J = 3*J+1
      END DO

      DO WHILE (J.GT.1)
        J = J/3
        JJ = 1
        DO WHILE (JJ.EQ.1)
          JJ = 0
          DO K=1,IMAX-J
            DIFF = ABS( W(IND(K)) - W(IND(K+J)) )
            IF ( W(IND(K)) .GT. W(IND(K+J)) .AND. &
                 DIFF.GT.BOUND ) THEN
              II       = IND(K)
              IND(K)   = IND(K+J)
              IND(K+J) = II
              JJ = 1
            END IF
          END DO                    ! K=1,IMAX-J
        END DO                      ! WHILE (JJ.EQ.1)
      END DO
      
      DO I=1,IMAX
        IF (IND(I) .EQ. 1) POS=I
      END DO

    ELSE

      DO I = 1,IMAX
        IND(I) = I
      END DO

      DO I=1,IMAX
        DO J=I,IMAX
          IF(W(IND(I))>W(IND(J))) THEN
            II=IND(J)
            IND(J)=IND(I)
            IND(I)=II
          END IF
        END DO !J
      END DO !I

    END IF
    END SUBROUTINE DSORT
END MODULE MOD_DSORT


