! ************************************************************************
SUBROUTINE rcstop(c)
! ************************************************************************
!     .. Scalar Arguments ..

CHARACTER (LEN=8), INTENT(IN) :: c

!     ..
PRINT *,'ERROR: STOP AT POSITION ',c
STOP

END SUBROUTINE rcstop
