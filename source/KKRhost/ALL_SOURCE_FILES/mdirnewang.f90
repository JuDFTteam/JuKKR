SUBROUTINE mdirnewang(it,nmvec,mvevi,mvphi,mvtet,mvgam,  &
        natypd,lmaxd,nmvecmax)
!   ********************************************************************
!   *                                                                  *
!   *  this routine has been build up from the last part of the        *
!   *  original Munich CALCMVEC routine.                               *
!   *  After correcting MVEVI with the Fermi energy value MVEVIEF      *
!   *  (outside this routine) it calculates the new angles of the      *
!   *  LOCAL FRAME quantisation axis with respect to the GLOBAL FRAME  *
!   *                                                                  *
!   ********************************************************************
IMPLICIT NONE

!Parameter definitions
INTEGER LMAXDLOC
PARAMETER (LMAXDLOC=8)

!Scalar Arguments
INTEGER IT,NMVEC,NATYPD,LMAXD,NMVECMAX

!Array Arguments
DOUBLE COMPLEX MVEVI(NATYPD,3,NMVECMAX)
DOUBLE PRECISION MVPHI(NATYPD,NMVECMAX),MVTET(NATYPD,NMVECMAX), &
                 MVGAM(NATYPD,NMVECMAX)

!Local Scalars
DOUBLE PRECISION MV,MVX,MVXY,MVY,MVZ,PI
INTEGER I,IMV,ICALL

!Local Arrays
DOUBLE PRECISION MVGLO(3,NMVECMAX)

!Intrinsic Functions
INTRINSIC ABS,ATAN

!Data Statements
DATA ICALL /0/

!Save Statements
SAVE ICALL,PI

icall = icall + 1
!=======================================================================
IF (icall == 1) THEN
  
  IF (lmaxd > lmaxdloc) THEN
    WRITE (6,*)
    WRITE (6,*) ' Please increase parameter LMAXDLOC to ',lmaxd
    WRITE (6,*) ' in the < MVECGLOBAL > routine.'
    STOP ' < TBKKR2 > '
  END IF
  
  pi = 4.d0 * ATAN(1.d0)
  
END IF
!=======================================================================

DO imv = 1,nmvec
  
  DO i = 1,3
    mvglo(i,imv) = DIMAG(mvevi(it,i,imv))
  END DO
  
  mvx = mvglo(1,imv)
  mvy = mvglo(2,imv)
  mvz = mvglo(3,imv)
  
  mv = SQRT(mvx**2+mvy**2+mvz**2)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  IF ( mv < 1D-8 ) THEN
    mvphi(it,imv) = 0D0
    mvtet(it,imv) = 0D0
    mvgam(it,imv) = 0D0
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ELSE
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    mvxy = SQRT(mvx**2+mvy**2)
! ----------------------------------------------------------------------
    IF ( ABS(mvxy) < 1D-8 ) THEN
      mvphi(it,imv) = 0D0
! ----------------------------------------------------------------------
    ELSE
! ----------------------------------------------------------------------
      IF ( mvy >= 0D0 ) THEN
        mvphi(it,imv) = ACOS(mvx/mvxy)
      ELSE IF ( mvx < 0D0 ) THEN
        mvphi(it,imv) = pi + ACOS(-mvx/mvxy)
      ELSE
        mvphi(it,imv) = 2*pi - ACOS(mvx/mvxy)
      END IF
      mvphi(it,imv) = mvphi(it,imv)*180D0/pi
      IF ( ABS(mvphi(it,imv)-360.0D0) < 1D-8 ) mvphi(it,imv) = 0D0
    END IF
! ----------------------------------------------------------------------
    IF (mvphi(it,imv) >= 345.d0) mvphi(it,imv) = 360.d0 - mvphi(it,imv)
    mvtet(it,imv) = ACOS(mvz/mv)*180D0/pi
    mvgam(it,imv) = 0D0
  END IF
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
END DO
!=======================================================================
END SUBROUTINE mdirnewang
