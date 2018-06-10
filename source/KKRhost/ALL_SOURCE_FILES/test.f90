! *********************************************************** 17.05.91 **
LOGICAL FUNCTION test(string)
! ***********************************************************************

!     TEST = 'STRING  ' IS CONTAINED IN /TESTC/.

! ------------------------------------------------------------------------
use mod_wunfiles, only: t_params

IMPLICIT NONE
CHARACTER (LEN=8), INTENT(IN)            :: string
INTEGER :: i
CHARACTER (LEN=8) :: testc(32)

testc = t_params%testc

test=.false.
DO i=1,32
  IF(string == testc(i)) test=.true.
END DO
RETURN
END FUNCTION test
