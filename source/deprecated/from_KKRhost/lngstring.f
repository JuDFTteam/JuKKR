      INTEGER FUNCTION LNGSTRING(STRING,LSTRMAX)
C   ********************************************************************
C   *                                                                  *
C   *  find position of last non-blank character in STRING(1:LSTRMAX)  *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C Dummy arguments
      INTEGER LSTRMAX
      CHARACTER*(*) STRING
C
C Local variables
      CHARACTER C
      INTEGER I,ICHAR
C
      LNGSTRING = 0
      DO I = LSTRMAX,1, - 1
         C = STRING(I:I)
         IF ( C.NE.' ' .AND. ICHAR(C).GT.0 ) THEN
            LNGSTRING = I
            RETURN
         END IF
      END DO
      END
