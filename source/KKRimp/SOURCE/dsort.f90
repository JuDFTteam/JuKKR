!-------------------------------------------------------------------------------
!> Summary: Sort double precision array returning sorted index array
!> Author: P. Zahn
!> Date: April 96
!> Sort double precision array returning sorted index array
!-------------------------------------------------------------------------------
MODULE MOD_DSORT
CONTAINS

      !-------------------------------------------------------------------------------
      !> Summary: Sort double precision array returning sorted index array
      !> Author: P. Zahn
      !> Date: April 96
      !> Category: KKRimp, numerical-tools
      !> Deprecated: False
      !> Sort double precision array returning sorted index array
      !-------------------------------------------------------------------------------
      SUBROUTINE DSORT (W,IND,IMAX,POS)
      IMPLICIT NONE
      INTEGER           :: IMAX         !! Number of elements to be sorted
      INTEGER           :: POS          !! Position where the first element is found
      DOUBLE PRECISION  :: W(:)         !! Original array returned unchanged
      INTEGER           :: IND(:)       !! Array that holds the new positions

      INTEGER           :: I,II,J,JJ,K
      INTEGER,PARAMETER :: BOUND =1.0D-12
      DOUBLE PRECISION  :: DIFF
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


