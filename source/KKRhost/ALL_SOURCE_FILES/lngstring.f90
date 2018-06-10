INTEGER FUNCTION lngstring(string,lstrmax)
!   ********************************************************************
!   *                                                                  *
!   *  find position of last non-blank character in STRING(1:LSTRMAX)  *
!   *                                                                  *
!   ********************************************************************
IMPLICIT NONE

! Dummy arguments
INTEGER LSTRMAX
CHARACTER*(*) STRING

! Local variables
CHARACTER C
INTEGER I,ICHAR

lngstring = 0
DO i = lstrmax,1, - 1
  c = string(i:i)
  IF ( c /= ' ' .AND. ICHAR(c) > 0 ) THEN
    lngstring = i
    RETURN
  END IF
END DO
END FUNCTION lngstring
