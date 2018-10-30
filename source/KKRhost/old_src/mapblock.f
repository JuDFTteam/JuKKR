      INTEGER FUNCTION MAPBLOCK(IE,IE1,NE,ITERSTEP,
     +                          NODEFIRST,NODELAST)
C **********************************************************************
C *                                                                    *
C *                                                                    *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C ..
C ..  Arguments ..
      INTEGER IE,IE1,NE,ITERSTEP
      INTEGER NODEFIRST,NODELAST
C ..
C ..  Locals ..
      INTEGER INC,IP,IPP,IPROC,JE,KE
      INTEGER IESORT(NE),IPROCE(NE)
C ......................................................................
      IPP = ITERSTEP            !         dummy use of argument iterstep
      DO JE = IE1,NE
         IESORT(JE) = JE
         IPROCE(JE) = 0
      END DO
C     
      IPP=0
      DO IP=NODEFIRST,NODELAST
         IPP=IPP+1
      END DO
C ----------------------------------------------------------------------
      IF ( IPP.GT.1 ) THEN
         IPROC = 0
         INC = 1
         DO JE = IE1,NE - 1
            KE = IESORT(JE)
            IPROC = IPROC + INC
C
            IF ( IPROC.EQ.IPP ) THEN
               IPROC = 0
               INC =  1
            ELSE IF ( IPROC.EQ.-1 ) THEN
               IPROC = 0
               INC = 1
            END IF
C
            IPROCE(KE) = IPROC
         END DO
         MAPBLOCK=IPROCE(IE)
C ----------------------------------------------------------------------
      ELSE
C ----------------------------------------------------------------------
         MAPBLOCK=0
      END IF
C ----------------------------------------------------------------------
C
      RETURN
      END
