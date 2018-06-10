! ***********************************************************************
LOGICAL FUNCTION opt(string)
! ***********************************************************************

!     OPT = 'STRING  ' IS CONTAINED IN /OPTC/.

! ------------------------------------------------------------------------
use mod_wunfiles, only: t_params

IMPLICIT NONE
CHARACTER (LEN=8), INTENT(IN)            :: string
INTEGER :: i
CHARACTER (LEN=8) :: optc(32)

optc = t_params%optc

opt=.false.
DO i=1,32
  IF(string == optc(i)) opt=.true.
END DO
RETURN
END FUNCTION opt
