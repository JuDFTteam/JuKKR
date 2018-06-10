! **********************************************************************
DOUBLE PRECISION FUNCTION dclock()
! **********************************************************************
!     .. External Functions ..
REAL :: etime,tarry(2)
EXTERNAL etime

!     ..
dclock = DBLE(etime(tarry))
RETURN
END FUNCTION dclock
