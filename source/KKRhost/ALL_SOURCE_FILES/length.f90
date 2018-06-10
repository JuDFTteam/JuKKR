! ************************************************************************
INTEGER FUNCTION length(s,MAX)
! ************************************************************************

CHARACTER (LEN=1), INTENT(IN OUT)        :: s(*)
INTEGER, INTENT(IN)                      :: MAX

INTEGER :: i
! ------------------------------------------------------------------------
i = MAX

DO  WHILE (s(i) == ' ')
  i = i - 1
END DO

length = i

RETURN
END FUNCTION length
