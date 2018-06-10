! ************************************************************************
SUBROUTINE dsort (w,ind,MAX,pos)
! ************************************************************************
!     p.zahn, april 96
!     W   is the original array returned unchanged
!     IND is an array that holds the new positions
!     max number of ellements to be sorted
!     pos the position where the first element is found
! ------------------------------------------------------------------------
IMPLICIT NONE
INTEGER MAX,POS
DOUBLE PRECISION  W(*)
INTEGER IND(*)

INTEGER I,II,J,JJ,K
DOUBLE PRECISION BOUND, DIFF
DATA BOUND /1.0D-12/
! ------------------------------------------------------------------------
DO  i = 1,MAX
  ind(i) = i
END DO

j = MAX
j = 1
DO  WHILE (j < MAX/3)
  j = 3*j+1
END DO

DO  WHILE (j > 1)
  j = j/3
  jj = 1
  DO  WHILE (jj == 1)
    jj = 0
    DO  k=1,MAX-j
      diff = ABS( w(ind(k)) - w(ind(k+j)) )
      IF ( w(ind(k)) > w(ind(k+j)) .AND. diff > bound ) THEN
        ii       = ind(k)
        ind(k)   = ind(k+j)
        ind(k+j) = ii
        jj = 1
      END IF
    END DO                    ! K=1,MAX-J
  END DO
END DO                      ! WHILE (JJ.EQ.1)

DO  i=1,MAX
  IF (ind(i) == 1) pos=i
END DO

RETURN
END SUBROUTINE dsort


