FUNCTION cdjlzdz(l,z,mode)
!   ********************************************************************
!   *                                                                  *
!   *     d j(L,Z) / dz    analytically                                *
!   *                                                                  *
!   ********************************************************************
IMPLICIT NONE

! Dummy arguments
INTEGER L,MODE
COMPLEX*16 Z
COMPLEX*16 CDJLZDZ

! Local variables
COMPLEX*16 CJLZ

IF ( mode == 1 ) THEN
  
  IF ( l == 0 ) THEN
    
    cdjlzdz = l*cjlz(l,z)/z - cjlz(l+1,z)
  ELSE
    cdjlzdz = (l*cjlz(l-1,z)-(l+1)*cjlz(l+1,z))/DBLE(2*l+1)
    RETURN
  END IF
ELSE IF ( mode == 2 ) THEN
  
  IF ( l == 0 ) THEN
    cdjlzdz = l*cjlz(l,z)/z - cjlz(l+1,z)
  ELSE
    cdjlzdz = cjlz(l-1,z) - (l+1)*cjlz(l,z)/z
    RETURN
  END IF
ELSE
  cdjlzdz = l*cjlz(l,z)/z - cjlz(l+1,z)
END IF
END FUNCTION cdjlzdz
