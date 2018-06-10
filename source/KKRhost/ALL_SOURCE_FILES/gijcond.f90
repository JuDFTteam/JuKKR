SUBROUTINE gijcond(ido,naez,rbasis,iqat,natomimp,rclsimp,atomimp,  &
    ijtabcalc,natomimpd)
! **********************************************************************
! *                                                                    *
! * In case of tasks requiring Gij blocks calculation, set variables:  *
! *                                                                    *
! * NATOMIMP, RCLSIMP(3,1..NATOMIMP), ATOMIMP(1..NATOMIMP)             *
! * IJTABCALC flag to which pair is needed: I,J --> (I-1)*NATOMIMP + J *
! * IDO takes on the value 1 or 0 if setting up process was OK or not  *
! *                                                                    *
! * CONDUCTANCE calculation case                                       *
! *             still to be implemente the correct read in             *
! **********************************************************************
IMPLICIT NONE

!Parameters
INTEGER NCPAIRD
PARAMETER (NCPAIRD=10)

!Arguments
INTEGER IDO,NAEZ,NATOMIMP,NATOMIMPD
INTEGER ATOMIMP(*),IJTABCALC(*),IQAT(*)
DOUBLE PRECISION RBASIS(3,*),RCLSIMP(3,*)

!Locals
INTEGER I,IAT,IATCONDL(NCPAIRD),IATCONDR(NCPAIRD),J,JAT,NCONDPAIR, &
        NN

ido = 0
DO i = 1,ncpaird
  iatcondl(i) = 0
  iatcondr(i) = 0
END DO

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
WRITE (1337,99001)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

!     ---------------------------------------------------- dummy
!     settings so far, need to be replaced by conductance input
!     and some output
ncondpair = 4
IF ( ncondpair > ncpaird ) THEN
  WRITE (6,99002) 'local','NCPAIRD',ncondpair
  STOP
END IF
iatcondl(1) = 1
iatcondr(1) = 2
iatcondl(2) = 1
iatcondr(2) = 2
iatcondl(3) = 2
iatcondr(3) = 1
iatcondl(4) = 2
iatcondr(4) = 1

!     ---------------------------------------------------- dummy
IF ( ncondpair == 0 ) RETURN
DO i = 1,ncondpair
  IF ( (iatcondl(i) <= 0) .OR. (iatcondl(i) > naez) ) RETURN
  IF ( (iatcondr(i) <= 0) .OR. (iatcondr(i) > naez) ) RETURN
END DO

natomimp = 2*ncondpair
IF ( natomimp > natomimpd ) THEN
  WRITE (6,99002) 'global','NATOMIMPD',natomimp
  STOP
END IF

DO i = 1,natomimp
  nn = (i-1)*natomimp
  DO j = 1,natomimp
    ijtabcalc(nn+j) = 0
  END DO
END DO

nn = 0
DO i = 1,ncondpair
  iat = iqat(iatcondl(i)) ! left lead
  nn = nn + 1
  DO j = 1,3
    rclsimp(j,nn) = rbasis(j,iat)
  END DO
  atomimp(nn) = iat
  iat = nn
  
  jat = iqat(iatcondr(i)) ! right lead
  nn = nn + 1
  DO j = 1,3
    rclsimp(j,nn) = rbasis(j,jat)
  END DO
  atomimp(nn) = jat
  jat = nn
  ijtabcalc((iat-1)*natomimp+jat) = 1
END DO
IF ( natomimp /= nn ) THEN
  WRITE (6,'(6X,A,/,6X,A,/)')  &
      'ERROR: Found some inconsistencies in IATCOND arrays'  &
      ,'       Please check your CONDUCTANCE input'
  STOP
END IF
ido = 1
99001 FORMAT (5X,'< GIJCOND > : Conductance/conductivity calculation',/)
99002 FORMAT (6X,'Dimension ERROR: please increase the ',a,' parameter',  &
    /,6X,a,' to a value >=',i5,/)
END SUBROUTINE gijcond
