C*==scalevec.f    processed by SPAG 6.05Rc at 11:25 on 26 Apr 2004
      SUBROUTINE SCALEVEC(LCARTESIAN,RBASIS,ABASIS,BBASIS,CBASIS,
     &                    NLBASIS,NRBASIS,NLEFT,NRIGHT,
     &                    ZPERLEFT,ZPERIGHT,TLEFT,TRIGHT,
     &                    LINTERFACE,NAEZ,NEMB,BRAVAIS,KAOEZ,NOQ,
     &                    NAEZD,NATYPD,NEMBD)
      IMPLICIT NONE
C     ..
C     .. Arguments ..
      INTEGER NAEZD,NATYPD,NEMBD
      INTEGER NAEZ,NEMB,NLBASIS,NLEFT,NRBASIS,NRIGHT
      DOUBLE PRECISION ABASIS,BBASIS,CBASIS
      LOGICAL LINTERFACE
      INTEGER KAOEZ(NATYPD,*),NOQ(NAEZD)
      DOUBLE PRECISION BRAVAIS(3,3),RBASIS(3,*),TLEFT(3,*),TRIGHT(3,*),
     &       ZPERIGHT(3),ZPERLEFT(3)
C     ..
C     .. Locals ..
      INTEGER I,I1,J
      LOGICAL LCARTESIAN
      DOUBLE PRECISION RBASIS1(3,NAEZD+NEMBD),TEMP(3),TX,TY,TZ
C
      WRITE (1337,'(79(1H=))')
      WRITE (1337,'(23X,A)') 'SCALEVEC: scale site coordinates'
      WRITE (1337,'(23X,A)') '          bring all to CARTESIAN system'
      WRITE (1337,'(79(1H=))')
      WRITE (1337,*)
C
C -->   normalization of basis vectors
C       multiplication instead of division 04/2004
C
      DO I = 1,NAEZ + NEMB
         RBASIS1(1,I) = RBASIS(1,I)*ABASIS
         RBASIS1(2,I) = RBASIS(2,I)*BBASIS
         RBASIS1(3,I) = RBASIS(3,I)*CBASIS
      END DO
C
      IF (LINTERFACE) THEN
         DO I = 1,NLBASIS
            TLEFT(1,I) = TLEFT(1,I)*ABASIS
            TLEFT(2,I) = TLEFT(2,I)*BBASIS
            TLEFT(3,I) = TLEFT(3,I)*CBASIS
         END DO
         ZPERLEFT(1) = ZPERLEFT(1)*ABASIS
         ZPERLEFT(2) = ZPERLEFT(2)*BBASIS
         ZPERLEFT(3) = ZPERLEFT(3)*CBASIS
C
         DO I = 1,NRBASIS
            TRIGHT(1,I) = TRIGHT(1,I)*ABASIS
            TRIGHT(2,I) = TRIGHT(2,I)*BBASIS
            TRIGHT(3,I) = TRIGHT(3,I)*CBASIS
         END DO
         ZPERIGHT(1) = ZPERIGHT(1)*ABASIS
         ZPERIGHT(2) = ZPERIGHT(2)*BBASIS
         ZPERIGHT(3) = ZPERIGHT(3)*CBASIS
      END IF
C
      IF ( ABASIS.NE.1D0 .OR. BBASIS.NE.1D0 .OR. CBASIS.NE.1D0 ) THEN
         WRITE (1337,'(5X,A,2(/,34X,F12.8,A))') 
     &        'Scaling site coordinates with:',
     &        ABASIS,'  x',BBASIS,'  y'
         IF ( .NOT.LINTERFACE ) WRITE (1337,'(34X,F12.8,A)') 
     &                               CBASIS,'  z'
         WRITE (1337,'(5X,44(1H-))')
      ELSE
         WRITE (1337,'(5X,A)') 
     &        'Site coordinates will not be scaled'
      END IF
C
C ---> normalization of atomic positions in the unit cell
C
C      if lcartesian is true cartesian coordinates are used
C      else the basis atoms are in units of the lattice vectors
C
      IF ( LCARTESIAN ) THEN
         WRITE (1337,'(A)') ' CARTESIAN coordinates'
      ELSE
         WRITE (1337,'(A)') ' LATTICE VECTOR coordinates will be', 
     &                   ' changed to CARTESIAN coordinates'
      END IF
C
C**********************************************************************
      ! Change to cartesian coordinates
      IF ( LINTERFACE ) THEN
C======================================================================
         IF ( .NOT.LCARTESIAN ) THEN
C----------------------------------------------------------------------
            WRITE (1337,*)
            WRITE(1337,'(12X,49(1H-))') 
            WRITE(1337,'(13X,A)') 
     &           'Input positions transformed to CARTESIAN system'
            WRITE(1337,'(12X,49(1H-),/,13X,A,/,12X,49(1H-))') 
     &           'IQ        x             y             z        IT'
            DO I = 1,NAEZ + NEMB
               DO J = 1,2
                  RBASIS(J,I) = (RBASIS1(1,I)*BRAVAIS(J,1)
     &                          +RBASIS1(2,I)*BRAVAIS(J,2))
               END DO
               RBASIS(3,I) = RBASIS1(3,I)
C
               IF ( I.LE.NAEZ ) THEN
                  WRITE (1337,99005) I,(RBASIS(J,I),J=1,3),
     &                            (KAOEZ(J,I),J=1,NOQ(I))
               ELSE
                  WRITE (1337,99005) I,(RBASIS(J,I),J=1,3),KAOEZ(1,I)
               END IF
               IF ( I.EQ.NAEZ ) WRITE(1337,'(12X,49(1H.))')
            END DO
            WRITE(1337,'(12X,49(1H-),/)')
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C -->  Do the same for the boundary vectors
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C      left side
C
            DO I = 1,NLBASIS
               DO I1 = 1,2
                  TEMP(I1) = TLEFT(I1,I)
               END DO
               DO J = 1,2
                  TLEFT(J,I) = (TEMP(1)*BRAVAIS(J,1)+
     &                          TEMP(2)*BRAVAIS(J,2))
               END DO
            END DO
C
            DO I1 = 1,2
               TEMP(I1) = ZPERLEFT(I1)
            END DO
            DO J = 1,2
               ZPERLEFT(J) = (TEMP(1)*BRAVAIS(J,1)+TEMP(2)*BRAVAIS(J,2))
            END DO
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C      right side
C
            DO I = 1,NRBASIS
               DO I1 = 1,2
                  TEMP(I1) = TRIGHT(I1,I)
               END DO
               DO J = 1,2
                  TRIGHT(J,I) = (TEMP(1)*BRAVAIS(J,1)
     &                          +TEMP(2)*BRAVAIS(J,2))
               END DO
            END DO
C
            DO I1 = 1,2
               TEMP(I1) = ZPERIGHT(I1)
            END DO
            DO J = 1,2
               ZPERIGHT(J) = (TEMP(1)*BRAVAIS(J,1)+TEMP(2)*BRAVAIS(J,2))
            END DO
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C----------------------------------------------------------------------
         ELSE
C----------------------------------------------------------------------
            WRITE(1337,'(42X,A)') 
     &           '---> No transformation required'
            DO I = 1,3
               DO J = 1,NAEZ + NEMB
                  RBASIS(I,J) = RBASIS1(I,J)
               END DO
            END DO
C----------------------------------------------------------------------
         END IF                 ! IF (.NOT.LCARTESIAN)
C======================================================================
C
         WRITE (1337,99002)
         DO I = NLEFT,1, - 1
            DO I1 = NLBASIS,1, - 1
               TX = TLEFT(1,I1) + (I-1)*ZPERLEFT(1)
               TY = TLEFT(2,I1) + (I-1)*ZPERLEFT(2)
               TZ = TLEFT(3,I1) + (I-1)*ZPERLEFT(3)
               WRITE (1337,99001) (I-1)*NLBASIS + I1,TX,TY,TZ,
     &                         KAOEZ(1,NAEZ+I1)
            END DO
         END DO
C
         WRITE (1337,99003)
         DO I = 1,NAEZ
            WRITE (1337,99001) I,(RBASIS(I1,I),I1=1,3),
     &                      (KAOEZ(I1,I),I1=1,NOQ(I))
         END DO
C
         WRITE (1337,99004)
         DO I = 1,NRIGHT
            DO I1 = 1,NRBASIS
               TX = TRIGHT(1,I1) + (I-1)*ZPERIGHT(1)
               TY = TRIGHT(2,I1) + (I-1)*ZPERIGHT(2)
               TZ = TRIGHT(3,I1) + (I-1)*ZPERIGHT(3)
               WRITE (1337,99001) (I-1)*NRBASIS + I1,TX,TY,TZ,
     &                         KAOEZ(1,NAEZ+NLBASIS+I1)
            END DO
         END DO
         WRITE(1337,'(14X,45(1H-),/)')
C======================================================================
      ELSE IF ( .NOT.LCARTESIAN ) THEN ! Rescale lattice
C----------------------------------------------------------------------
         WRITE (1337,*)
         WRITE(1337,'(12X,49(1H-))') 
         WRITE(1337,'(13X,A)') 
     &        'Input positions transformed to CARTESIAN system'
         WRITE(1337,'(12X,49(1H-),/,13X,A,/,12X,49(1H-))') 
     &        'IQ        x             y             z        IT'
         DO I = 1,NAEZ + NEMB
            DO J = 1,3
               RBASIS(J,I) = (RBASIS1(1,I)*BRAVAIS(J,1)
     &                       +RBASIS1(2,I)*BRAVAIS(J,2)
     &                       +RBASIS1(3,I)*BRAVAIS(J,3))
            END DO
C
            IF ( I.LE.NAEZ ) THEN
               WRITE (1337,99005) I,(RBASIS(J,I),J=1,3),
     &                         (KAOEZ(J,I),J=1,NOQ(I))
            ELSE
               WRITE (1337,99005) I,(RBASIS(J,I),J=1,3),KAOEZ(1,I)
            END IF
            IF ( I.EQ.NAEZ .AND. NEMB.GT.0 ) WRITE(1337,'(12X,49(1H.))')
         END DO
         WRITE(1337,'(12X,49(1H-),/)')
C----------------------------------------------------------------------
      ELSE
C----------------------------------------------------------------------
         WRITE(1337,'(42X,A,/)') 
     &           '---> No transformation required'
         WRITE(1337,99006) 
C     changed by v.Bellini 21/10/99
         DO J = 1,NAEZ + NEMB
            DO I = 1,3
               RBASIS(I,J) = RBASIS1(I,J)
            END DO
C
            IF ( J.LE.NAEZ ) THEN
               WRITE (1337,99001) J,(RBASIS(I,J),I=1,3),
     &                         (KAOEZ(I,J),I=1,NOQ(J))
            ELSE
               WRITE (1337,99001) J,(RBASIS(I,J),I=1,3),KAOEZ(1,J)
            END IF
            IF ( I.EQ.NAEZ .AND. NEMB.GT.0 ) WRITE(1337,'(12X,51(1H.))')
         END DO
C     end of the change
         WRITE(1337,'(12X,51(1H-),/)')
C----------------------------------------------------------------------
C======================================================================
      END IF                    !  IF (.NOT.LINTERFACE )
C**********************************************************************
C
C FROM NOW ON after < SCALEVEC > RBASIS are the basis vectors
C in units of au/alat in (xyz) reference
C
C**********************************************************************
99001 FORMAT (13X,I5,3F12.6,10I3)
99002 FORMAT (14X,45(1H-),/,
     &        15X,'     Positions of ALL generated sites ',/,
     &        15X,'   in CARTESIAN coordinates (ALAT units)',/,
     &        14X,45(1H-),/,
     &        15X,'IQ       x           y           z       IT',/,
     &        15X,'**************** Left  Host ***************')
99003 FORMAT (15X,'****************   S L A B  ***************')
99004 FORMAT (15X,'**************** Right Host ***************')
99005 FORMAT (12X,I3,3F14.8,10I3)
99006 FORMAT (12X,51(1H-),/,
     &        16X,'    Positions of (ALL) generated sites',/,
     &        16X,'   in CARTESIAN coordinates (ALAT units)',/,
     &        12X,51(1H-),/,
     &        15X,'IQ       x           y           z       IT',/,
     &        12X,51(1H-))
      END
