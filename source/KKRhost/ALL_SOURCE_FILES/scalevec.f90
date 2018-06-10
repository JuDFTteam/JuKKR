SUBROUTINE scalevec(lcartesian,rbasis,abasis,bbasis,cbasis,  &
    nlbasis,nrbasis,nleft,nright, zperleft,zperight,tleft,tright,  &
    linterface,naez,nemb,bravais,kaoez,noq, naezd,natypd,nembd)

      IMPLICIT NONE
!..
!.. Arguments ..
      INTEGER NAEZD,NATYPD,NEMBD
      INTEGER NAEZ,NEMB,NLBASIS,NLEFT,NRBASIS,NRIGHT
      DOUBLE PRECISION ABASIS,BBASIS,CBASIS
      LOGICAL LINTERFACE
      INTEGER KAOEZ(NATYPD,*),NOQ(NAEZD)
      DOUBLE PRECISION BRAVAIS(3,3),RBASIS(3,*),TLEFT(3,*),TRIGHT(3,*), & 
             ZPERIGHT(3),ZPERLEFT(3)
!..
!.. Locals ..
      INTEGER I,I1,J
      LOGICAL LCARTESIAN
      DOUBLE PRECISION RBASIS1(3,NAEZD+NEMBD),TEMP(3),TX,TY,TZ

WRITE (1337,'(79(1H=))')
WRITE (1337,'(23X,A)') 'SCALEVEC: scale site coordinates'
WRITE (1337,'(23X,A)') '          bring all to CARTESIAN system'
WRITE (1337,'(79(1H=))')
WRITE (1337,*)

! -->   normalization of basis vectors
!       multiplication instead of division 04/2004

DO i = 1,naez + nemb
  rbasis1(1,i) = rbasis(1,i)*abasis
  rbasis1(2,i) = rbasis(2,i)*bbasis
  rbasis1(3,i) = rbasis(3,i)*cbasis
END DO

IF (linterface) THEN
  DO i = 1,nlbasis
    tleft(1,i) = tleft(1,i)*abasis
    tleft(2,i) = tleft(2,i)*bbasis
    tleft(3,i) = tleft(3,i)*cbasis
  END DO
  zperleft(1) = zperleft(1)*abasis
  zperleft(2) = zperleft(2)*bbasis
  zperleft(3) = zperleft(3)*cbasis
  
  DO i = 1,nrbasis
    tright(1,i) = tright(1,i)*abasis
    tright(2,i) = tright(2,i)*bbasis
    tright(3,i) = tright(3,i)*cbasis
  END DO
  zperight(1) = zperight(1)*abasis
  zperight(2) = zperight(2)*bbasis
  zperight(3) = zperight(3)*cbasis
END IF

IF ( abasis /= 1D0 .OR. bbasis /= 1D0 .OR. cbasis /= 1D0 ) THEN
  WRITE (1337,'(5X,A,2(/,34X,F12.8,A))') 'Scaling site coordinates with:',  &
      abasis,'  x',bbasis,'  y'
  IF ( .NOT.linterface ) WRITE (1337,'(34X,F12.8,A)') cbasis,'  z'
  WRITE (1337,'(5X,44(1H-))')
ELSE
  WRITE (1337,'(5X,A)') 'Site coordinates will not be scaled'
END IF

! ---> normalization of atomic positions in the unit cell

!      if lcartesian is true cartesian coordinates are used
!      else the basis atoms are in units of the lattice vectors

IF ( lcartesian ) THEN
  WRITE (1337,'(A)') ' CARTESIAN coordinates'
ELSE
  WRITE (1337,'(A)') ' LATTICE VECTOR coordinates will be',  &
      ' changed to CARTESIAN coordinates'
END IF

!**********************************************************************
! Change to cartesian coordinates
IF ( linterface ) THEN
!======================================================================
  IF ( .NOT.lcartesian ) THEN
!----------------------------------------------------------------------
    WRITE (1337,*)
    WRITE(1337,'(12X,49(1H-))')
    WRITE(1337,'(13X,A)') 'Input positions transformed to CARTESIAN system'
    WRITE(1337,'(12X,49(1H-),/,13X,A,/,12X,49(1H-))')  &
        'IQ        x             y             z        IT'
    DO i = 1,naez + nemb
      DO j = 1,2
        rbasis(j,i) = (rbasis1(1,i)*bravais(j,1) +rbasis1(2,i)*bravais(j,2))
      END DO
      rbasis(3,i) = rbasis1(3,i)
      
      IF ( i <= naez ) THEN
        WRITE (1337,99005) i,(rbasis(j,i),j=1,3), (kaoez(j,i),j=1,noq(i))
      ELSE
        WRITE (1337,99005) i,(rbasis(j,i),j=1,3),kaoez(1,i)
      END IF
      IF ( i == naez ) WRITE(1337,'(12X,49(1H.))')
    END DO
    WRITE(1337,'(12X,49(1H-),/)')
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
! -->  Do the same for the boundary vectors
    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      left side
    
    DO i = 1,nlbasis
      DO i1 = 1,2
        temp(i1) = tleft(i1,i)
      END DO
      DO j = 1,2
        tleft(j,i) = (temp(1)*bravais(j,1)+ temp(2)*bravais(j,2))
      END DO
    END DO
    
    DO i1 = 1,2
      temp(i1) = zperleft(i1)
    END DO
    DO j = 1,2
      zperleft(j) = (temp(1)*bravais(j,1)+temp(2)*bravais(j,2))
    END DO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      right side
    
    DO i = 1,nrbasis
      DO i1 = 1,2
        temp(i1) = tright(i1,i)
      END DO
      DO j = 1,2
        tright(j,i) = (temp(1)*bravais(j,1) +temp(2)*bravais(j,2))
      END DO
    END DO
    
    DO i1 = 1,2
      temp(i1) = zperight(i1)
    END DO
    DO j = 1,2
      zperight(j) = (temp(1)*bravais(j,1)+temp(2)*bravais(j,2))
    END DO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!----------------------------------------------------------------------
  ELSE
!----------------------------------------------------------------------
    WRITE(1337,'(42X,A)') '---> No transformation required'
    DO i = 1,3
      DO j = 1,naez + nemb
        rbasis(i,j) = rbasis1(i,j)
      END DO
    END DO
!----------------------------------------------------------------------
  END IF                 ! IF (.NOT.LCARTESIAN)
!======================================================================
  
  WRITE (1337,99002)
  DO i = nleft,1, - 1
    DO i1 = nlbasis,1, - 1
      tx = tleft(1,i1) + (i-1)*zperleft(1)
      ty = tleft(2,i1) + (i-1)*zperleft(2)
      tz = tleft(3,i1) + (i-1)*zperleft(3)
      WRITE (1337,99001) (i-1)*nlbasis + i1,tx,ty,tz, kaoez(1,naez+i1)
    END DO
  END DO
  
  WRITE (1337,99003)
  DO i = 1,naez
    WRITE (1337,99001) i,(rbasis(i1,i),i1=1,3), (kaoez(i1,i),i1=1,noq(i))
  END DO
  
  WRITE (1337,99004)
  DO i = 1,nright
    DO i1 = 1,nrbasis
      tx = tright(1,i1) + (i-1)*zperight(1)
      ty = tright(2,i1) + (i-1)*zperight(2)
      tz = tright(3,i1) + (i-1)*zperight(3)
      WRITE (1337,99001) (i-1)*nrbasis + i1,tx,ty,tz, kaoez(1,naez+nlbasis+i1)
    END DO
  END DO
  WRITE(1337,'(14X,45(1H-),/)')
!======================================================================
ELSE IF ( .NOT.lcartesian ) THEN ! Rescale lattice
!----------------------------------------------------------------------
  WRITE (1337,*)
  WRITE(1337,'(12X,49(1H-))')
  WRITE(1337,'(13X,A)') 'Input positions transformed to CARTESIAN system'
  WRITE(1337,'(12X,49(1H-),/,13X,A,/,12X,49(1H-))')  &
      'IQ        x             y             z        IT'
  DO i = 1,naez + nemb
    DO j = 1,3
      rbasis(j,i) = (rbasis1(1,i)*bravais(j,1) +rbasis1(2,i)*bravais(j,2)  &
          +rbasis1(3,i)*bravais(j,3))
    END DO
    
    IF ( i <= naez ) THEN
      WRITE (1337,99005) i,(rbasis(j,i),j=1,3), (kaoez(j,i),j=1,noq(i))
    ELSE
      WRITE (1337,99005) i,(rbasis(j,i),j=1,3),kaoez(1,i)
    END IF
    IF ( i == naez .AND. nemb > 0 ) WRITE(1337,'(12X,49(1H.))')
  END DO
  WRITE(1337,'(12X,49(1H-),/)')
!----------------------------------------------------------------------
ELSE
!----------------------------------------------------------------------
  WRITE(1337,'(42X,A,/)') '---> No transformation required'
  WRITE(1337,99006)
!     changed by v.Bellini 21/10/99
  DO j = 1,naez + nemb
    DO i = 1,3
      rbasis(i,j) = rbasis1(i,j)
    END DO
    
    IF ( j <= naez ) THEN
      WRITE (1337,99001) j,(rbasis(i,j),i=1,3), (kaoez(i,j),i=1,noq(j))
    ELSE
      WRITE (1337,99001) j,(rbasis(i,j),i=1,3),kaoez(1,j)
    END IF
    IF ( i == naez .AND. nemb > 0 ) WRITE(1337,'(12X,51(1H.))')
  END DO
!     end of the change
  WRITE(1337,'(12X,51(1H-),/)')
!----------------------------------------------------------------------
!======================================================================
END IF                    !  IF (.NOT.LINTERFACE )
!**********************************************************************

! FROM NOW ON after < SCALEVEC > RBASIS are the basis vectors
! in units of au/alat in (xyz) reference

!**********************************************************************
99001 FORMAT (13X,i5,3F12.6,10I3)
99002 FORMAT (14X,45(1H-),/,  &
    15X,'     Positions of ALL generated sites ',/,  &
    15X,'   in CARTESIAN coordinates (ALAT units)',/, 14X,45(1H-),/,  &
    15X,'IQ       x           y           z       IT',/,  &
    15X,'**************** Left  Host ***************')
99003 FORMAT (15X,'****************   S L A B  ***************')
99004 FORMAT (15X,'**************** Right Host ***************')
99005 FORMAT (12X,i3,3F14.8,10I3)
99006 FORMAT (12X,51(1H-),/,  &
    16X,'    Positions of (ALL) generated sites',/,  &
    16X,'   in CARTESIAN coordinates (ALAT units)',/, 12X,51(1H-),/,  &
    15X,'IQ       x           y           z       IT',/, 12X,51(1H-))
END SUBROUTINE scalevec
