SUBROUTINE impcheck(atomimp,natomimp,naez,rclsimp,  &
        rbasis,bravais,ndim)
! **********************************************************************
! * Checking the coordinates and site-index assignments of an impurity *
! * cluster read in from a file                                        *
! * The size of the Bravais lattice to be generated is determined      *
! * dynamically using the radius of the input cluster as reference     *
! *                                                                    *
! * For an input site not belonging to the Bravais lattice the program *
! * stops.                                                             *
! * A wrong site-index (unit cell) assignment is corrected and the     *
! * execution continues.                                               *
! **********************************************************************

      IMPLICIT NONE
!..
!.. Scalar arguments
      INTEGER NAEZ,NATOMIMP,NDIM
!.. Array arguments
      INTEGER ATOMIMP(*)
      DOUBLE PRECISION BRAVAIS(3,3),RBASIS(3,*),RCLSIMP(3,*)
!..
!.. Local scalars
      INTEGER I,IATOK,IPOSOK,IQ,J,N1,N2,N3,NMAX,NMAXZ
      DOUBLE PRECISION DIFF,RMAXCLUS,RMAXGEN
      CHARACTER*6 STRAT,STRPOS
!..
!.. Local arrays 
      INTEGER AIN(NATOMIMP),NBR(3)
      DOUBLE PRECISION RCLSNEW(3,NATOMIMP),VEC1(3),VEC2(3)
      LOGICAL LATOM(NATOMIMP),LPOS(NATOMIMP),LABSCORD
      DOUBLE PRECISION DIFFMIN(NATOMIMP)
!..
!.. External subroutine
      EXTERNAL GETCLUSNXYZ
!..

! ----------------------------------------------------------------------
!     initialize diffmin array with high value
diffmin(:) = 1D+5


!     LABSCORD - cluster coordinates are absolute atomic positions

labscord = .false.
j = 0
DO WHILE ( j < 3.AND..NOT.labscord )
  j = j + 1
  IF ( ABS(rclsimp(j,1)) > 1D-8 ) labscord = .true.
END DO

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! --> set cluster coordinates in absolute atomic positions
!     RCLSNEW(3,*)

DO i = 1,natomimp
  CALL dcopy(3,rclsimp(1,i),1,rclsnew(1,i),1)
END DO
IF ( .NOT.labscord ) THEN
  iq = atomimp(1)
  DO i = 1,natomimp
    CALL daxpy(3,1D0,rbasis(1,iq),1,rclsnew(1,i),1)
  END DO
END IF

! --> determine the maximum radius of the input cluster
!     this will be then compared to the maximum generated radius
!     when testing the positions -- setting NMAX for generating the
!     lattice

rmaxclus = 0D0
DO i = 2,natomimp
  diff = 0D0
  DO j = 1,3
    diff = diff + (rclsnew(j,i)-rclsnew(j,1))**2
  END DO
  diff = SQRT(diff)
  rmaxclus = MAX(rmaxclus,diff)
END DO
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

DO i = 1,3
  nbr(i) = 0
END DO
rmaxclus = 3.0*rmaxclus
!       RMAXCLUS = 1.5*RMAXCLUS
CALL getclusnxyz(rmaxclus,bravais,ndim,diff,nbr)
nmax = MAX(nbr(1),nbr(2),nbr(3))
nmaxz = nmax
IF ( ndim == 2 ) nmaxz = 0
rmaxgen = 0D0
WRITE(1337,*) nmax,nmaxz

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
WRITE(1337,*) 'rbasis of impurity (in cartesian coordinate)'
DO i = 1,natomimp
  WRITE(1337,*) (rclsnew(j,i),j=1,3)
  ain(i) = atomimp(i)
  lpos(i) = .false.
  latom(i) = .true.
!=======================================================================
  DO n1 = -nmax,nmax
    DO n2 = -nmax,nmax
      DO n3 = -nmaxz,nmaxz
        
        DO j = 1,3
          vec1(j) = DBLE(n1)*bravais(j,1) + DBLE(n2)  &
              *bravais(j,2) + DBLE(n3)*bravais(j,3)
        END DO
        
!-----------------------------------------------------------------------
        DO iq = 1,naez
          diff = 0D0
          DO j = 1,3
            vec2(j) = vec1(j) + rbasis(j,iq)
            diff = diff + vec2(j)**2
          END DO
          rmaxgen = MAX(rmaxgen,SQRT(diff))
          
          diff = SQRT((rclsnew(1,i)-vec2(1))**2 +(rclsnew(2,i)-vec2(2))**2  &
              +(rclsnew(3,i)-vec2(3))**2)
          
          IF ( diff <= (1D-5) ) THEN
            atomimp(i) = iq
            IF ( ain(i) /= iq ) latom(i) = .false.
            lpos(i) = .true.
            GO TO 100
          END IF
          
          IF (diff <= diffmin(i)) diffmin(i)=diff
          
        END DO
!-----------------------------------------------------------------------
      END DO
    END DO
  END DO
!=======================================================================
100  END DO
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
WRITE (1337,99001) rmaxclus,rmaxgen
WRITE (1337,99002) 'Input data for impurity sites - consistency check'


! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
iatok = 0
iposok = 0
OPEN(58,FILE='rimp.dat')
WRITE(58,*) natomimp,' NATOMIMP'
DO i = 1,natomimp
  IF ( lpos(i) ) THEN
    WRITE (strpos,'(A6)') 'OK'
  ELSE
    WRITE (strpos,'(A6)') 'neq BL'
    iposok = iposok + 1
    WRITE(*,*) 'minimal difference for atom',i,'=',diffmin(i)
  END IF
  WRITE (strat,'(I3,A3)') atomimp(i),' <?'
  IF ( .NOT.latom(i) ) THEN
    iatok = iatok + 1
  ELSE IF ( lpos(i) ) THEN
    WRITE (strat,'(I3)') atomimp(i)
  END IF
  WRITE (1337,99003) i,(rclsimp(j,i),j=1,3),ain(i),strpos,strat
  WRITE (58,FMT='(I6,3E16.8,I6)') i,(rclsimp(j,i),j=1,3),ain(i)
END DO
CLOSE(58)
WRITE (1337,99004)

IF ( iposok /= 0 ) THEN
  WRITE (6,99005)
  STOP
END IF

DO i = 1,natomimp
  IF ( (atomimp(i) > naez) .OR. (atomimp(i) == 0) ) THEN
    WRITE (6,99008) i,atomimp(i)
    STOP
  END IF
END DO

IF ( iatok /= 0 ) THEN
  WRITE (1337,99006)
ELSE
  WRITE (1337,99007)
END IF
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

99001 FORMAT (12X,'input-cluster R  : ',f11.6,/,12X,  &
    'test-cluster  R  : ',f11.6,/)
99002 FORMAT (13X,63('-'),/,15X,a,/,13X,63('-'),/,13X,  &
    ' imp |               READ IN DATA         host  |',  &
    ' CHECKED DATA ',/,13X,  &
    'index|       x           y           z    site  |',  &
    '  pos.   site ',/,13X,63('-'))
99003 FORMAT (13X,i3,2X,'|',3(f12.6),1X,i4,1X,'|',a6,2X,a6)
99004 FORMAT (13X,63('-'))
99005 FORMAT (/,6X,'ERROR: At least one of your input sites does not',  &
    ' belong to the Bravais ',/,13X,'lattice (neq BL). ',  &
    'Please check your input file',/)
99006 FORMAT (13X,'WARNING: At least one inconsistent assignment of ',  &
    'site indices',/,13X, '         was found in your input. The program will',  &
    ' override the',/,13X, '         input data. Crosscheck?  ',/,13X,63('-'),/)
99007 FORMAT (13X,'Your cluster data is consistent',/,13X,63('-'),/)
99008 FORMAT (/,6X,'ERROR: Wrong assignment of impurity site ',i3,  &
    ' to the unit-cell site ',i3,/,13X, 'Please check your input file',/)
END SUBROUTINE impcheck
