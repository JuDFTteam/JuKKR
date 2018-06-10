FUNCTION cdnlzdz(l,z,mode)
!   ********************************************************************
!   *                                                                  *
!   *     d n(L,Z) / dz    analytically                                *
!   *                                                                  *
!   ********************************************************************
IMPLICIT NONE

! Dummy arguments
      INTEGER L,MODE
      COMPLEX*16 Z
      COMPLEX*16 CDNLZDZ

! Local variables
      COMPLEX*16 CNLZ

IF ( mode == 1 ) THEN
  
  IF ( l == 0 ) THEN
    
    cdnlzdz = l*cnlz(l,z)/z - cnlz(l+1,z)
  ELSE
    cdnlzdz = (l*cnlz(l-1,z)-(l+1)*cnlz(l+1,z))/DBLE(2*l+1)
    RETURN
  END IF
ELSE IF ( mode == 2 ) THEN
  
  IF ( l == 0 ) THEN
    cdnlzdz = l*cnlz(l,z)/z - cnlz(l+1,z)
  ELSE
    cdnlzdz = cnlz(l-1,z) - (l+1)*cnlz(l,z)/z
    RETURN
  END IF
ELSE
  cdnlzdz = l*cnlz(l,z)/z - cnlz(l+1,z)
END IF
END FUNCTION cdnlzdz
