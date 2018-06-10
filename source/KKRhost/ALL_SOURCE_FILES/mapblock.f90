INTEGER FUNCTION mapblock(ie,ie1,NE,iterstep,  &
        nodefirst,nodelast)
! **********************************************************************
! *                                                                    *
! *                                                                    *
! *                                                                    *
! **********************************************************************
IMPLICIT NONE

!Arguments ..
INTEGER IE,IE1,NE,ITERSTEP
INTEGER NODEFIRST,NODELAST

!Locals ..
INTEGER INC,IP,IPP,IPROC,JE,KE
INTEGER IESORT(NE),IPROCE(NE)
! ......................................................................
ipp = iterstep            !         dummy use of argument iterstep
DO je = ie1,NE
  iesort(je) = je
  iproce(je) = 0
END DO

ipp=0
DO ip=nodefirst,nodelast
  ipp=ipp+1
END DO
! ----------------------------------------------------------------------
IF ( ipp > 1 ) THEN
  iproc = 0
  inc = 1
  DO je = ie1,NE - 1
    ke = iesort(je)
    iproc = iproc + inc
    
    IF ( iproc == ipp ) THEN
      iproc = 0
      inc =  1
    ELSE IF ( iproc == -1 ) THEN
      iproc = 0
      inc = 1
    END IF
    
    iproce(ke) = iproc
  END DO
  mapblock=iproce(ie)
! ----------------------------------------------------------------------
ELSE
! ----------------------------------------------------------------------
  mapblock=0
END IF
! ----------------------------------------------------------------------

RETURN
END FUNCTION mapblock
