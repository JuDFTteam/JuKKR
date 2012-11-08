      SUBROUTINE SCALEVEC(RBASIS,ABASIS,BBASIS,CBASIS,
     &                    NAEZ,BRAVAIS, LCARTESIAN)
      IMPLICIT NONE
C      INCLUDE 'inc.p'
C
C PARAMETER definitions
C
C Dummy arguments
C
      DOUBLE PRECISION ABASIS,BBASIS,CBASIS
      INTEGER NAEZ
      DOUBLE PRECISION BRAVAIS(3,3),RBASIS(3,*)
      LOGICAL LCARTESIAN
C
C Local variables
C
      INTEGER I,J
      DOUBLE PRECISION RBASIS1(3,NAEZ)
C
      WRITE(6,'(79(1H=))')
      WRITE (6,'(23X,A)') 'SCALEVEC: scale site coordinates'
      WRITE (6,'(23X,A)') '          bring all to CARTESIAN system'
      WRITE (6,'(79(1H=))')
      WRITE (6,*)
C
C -->   normalization of basis vectors
C
      DO I = 1,NAEZ
         RBASIS1(1,I) = RBASIS(1,I)*ABASIS
         RBASIS1(2,I) = RBASIS(2,I)*BBASIS
         RBASIS1(3,I) = RBASIS(3,I)*CBASIS
      END DO
C
      IF ( ABASIS.NE.1D0 .OR. BBASIS.NE.1D0 .OR. CBASIS.NE.1D0 ) THEN
         WRITE (6,'(5X,A,2(/,34X,F12.8,A))') 
     &        'Scaling site coordinates with:',
     &        ABASIS,'  x',BBASIS,'  y'
         WRITE (6,'(34X,F12.8,A)') CBASIS,'  z'
         WRITE (6,'(5X,44(1H-))')
      ELSE
         WRITE (6,'(5X,A)') 
     &        'Site coordinates will not be scaled'
      END IF
C
C ---> normalization of atomic positions in the unit cell
C
C      if lcartesian is true cartesian coordinates are used
C      else the basis atoms are in units of the lattice vectors
C
C
      WRITE (6,'(5X,A,$)')
     &     'Position of atoms in the unit cell READ IN as:'
      IF ( LCARTESIAN ) THEN
         WRITE (6,'(A)') ' CARTESIAN coordinates'
      ELSE
         WRITE (6,'(A)') ' LATTICE VECTORS units'
      END IF
C
C**********************************************************************
      IF ( .NOT.LCARTESIAN ) THEN ! Rescale lattice
C----------------------------------------------------------------------
         WRITE (6,*)
         WRITE(6,'(12X,49(1H-))') 
         WRITE(6,'(13X,A)') 
     &        'Input positions transformed to CARTESIAN system'
         WRITE(6,'(12X,49(1H-),/,13X,A,/,12X,49(1H-))') 
     &        'IQ        x             y             z        IT'
         DO I = 1,NAEZ
            DO J = 1,3
               RBASIS(J,I) = (RBASIS1(1,I)*BRAVAIS(J,1)
     &                       +RBASIS1(2,I)*BRAVAIS(J,2)
     &                       +RBASIS1(3,I)*BRAVAIS(J,3))
            END DO
C
               WRITE (6,99005) I,(RBASIS(J,I),J=1,3),I
         END DO
         WRITE(6,'(12X,49(1H-),/)')
C----------------------------------------------------------------------
      ELSE
C----------------------------------------------------------------------
         WRITE(6,'(42X,A,/)') 
     &           '---> No transformation required'
         WRITE(6,99006) 
C     changed by v.Bellini 21/10/99
         DO J = 1,NAEZ
            DO I = 1,3
               RBASIS(I,J) = RBASIS1(I,J)
            END DO
C
               WRITE (6,99001) J,(RBASIS(I,J),I=1,3),J
         END DO
C     end of the change
         WRITE(6,'(12X,51(1H-),/)')
C----------------------------------------------------------------------
C======================================================================
      END IF
C**********************************************************************
C
C FROM NOW ON after < SCALEVEC > RBASIS are the basis vectors
C in units of au/alat in (xyz) reference
C
C**********************************************************************
99001 FORMAT (13X,I5,3F12.6,10I3)
99005 FORMAT (12X,I3,3F14.8,10I3)
99006 FORMAT (12X,51(1H-),/,
     &        16X,'    Positions of (ALL) generated sites',/,
     &        16X,'   in CARTESIAN coordinates (ALAT units)',/,
     &        12X,51(1H-),/,
     &        15X,'IQ       x           y           z       IT',/,
     &        12X,51(1H-))
      END
