      SUBROUTINE IMPCHECK(ATOMIMP,NATOMIMP,NAEZ,RCLSIMP,
     &                    RBASIS,BRAVAIS,NDIM)
C **********************************************************************
C * Checking the coordinates and site-index assignments of an impurity *
C * cluster read in from a file                                        *
C * The size of the Bravais lattice to be generated is determined      *
C * dynamically using the radius of the input cluster as reference     *
C *                                                                    *
C * For an input site not belonging to the Bravais lattice the program *
C * stops.                                                             *
C * A wrong site-index (unit cell) assignment is corrected and the     *
C * execution continues.                                               *
C **********************************************************************
      IMPLICIT NONE
C     ..
C     .. Scalar arguments
      INTEGER NAEZ,NATOMIMP,NDIM
C     .. Array arguments
      INTEGER ATOMIMP(*)
      DOUBLE PRECISION BRAVAIS(3,3),RBASIS(3,*),RCLSIMP(3,*)
C     ..
C     .. Local scalars
      INTEGER I,IATOK,IPOSOK,IQ,J,N1,N2,N3,NMAX,NMAXZ
      DOUBLE PRECISION DIFF,RMAXCLUS,RMAXGEN
      CHARACTER*6 STRAT,STRPOS
C     ..
C     .. Local arrays 
      INTEGER AIN(NATOMIMP),NBR(3)
      DOUBLE PRECISION RCLSNEW(3,NATOMIMP),VEC1(3),VEC2(3)
      LOGICAL LATOM(NATOMIMP),LPOS(NATOMIMP),LABSCORD
      DOUBLE PRECISION DIFFMIN(NATOMIMP)
C     ..
C     .. External subroutine
      EXTERNAL GETCLUSNXYZ
C     ..
C
C ----------------------------------------------------------------------
!     initialize diffmin array with high value
      DIFFMIN(:) = 1D+5

C
C     LABSCORD - cluster coordinates are absolute atomic positions 
C
      LABSCORD = .FALSE.
      J = 0
      DO WHILE ( J.LT.3.AND..NOT.LABSCORD )
         J = J + 1
         IF ( ABS(RCLSIMP(J,1)).GT.1D-8 ) LABSCORD = .TRUE.
      END DO
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C --> set cluster coordinates in absolute atomic positions
C     RCLSNEW(3,*)
C
      DO I = 1,NATOMIMP
         CALL DCOPY(3,RCLSIMP(1,I),1,RCLSNEW(1,I),1)
      END DO
      IF ( .NOT.LABSCORD ) THEN
         IQ = ATOMIMP(1)
         DO I = 1,NATOMIMP
            CALL DAXPY(3,1D0,RBASIS(1,IQ),1,RCLSNEW(1,I),1)
         END DO
      END IF
C
C --> determine the maximum radius of the input cluster
C     this will be then compared to the maximum generated radius
C     when testing the positions -- setting NMAX for generating the
C     lattice
C
      RMAXCLUS = 0D0
      DO I = 2,NATOMIMP
         DIFF = 0D0
         DO J = 1,3
            DIFF = DIFF + (RCLSNEW(J,I)-RCLSNEW(J,1))**2
         END DO
         DIFF = SQRT(DIFF)
         RMAXCLUS = MAX(RMAXCLUS,DIFF)
      END DO
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      DO I = 1,3
         NBR(I) = 0
      END DO
      RMAXCLUS = 3.0*RMAXCLUS
!       RMAXCLUS = 1.5*RMAXCLUS
      CALL GETCLUSNXYZ(RMAXCLUS,BRAVAIS,NDIM,DIFF,NBR)
      NMAX = MAX(NBR(1),NBR(2),NBR(3))
      NMAXZ = NMAX 
      IF ( NDIM.EQ.2 ) NMAXZ = 0
      RMAXGEN = 0D0
      WRITE(1337,*) NMAX,NMAXZ
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      WRITE(1337,*) 'rbasis of impurity (in cartesian coordinate)'
      DO I = 1,NATOMIMP
         WRITE(1337,*) (RCLSNEW(J,I),J=1,3)
         AIN(I) = ATOMIMP(I)
         LPOS(I) = .FALSE.
         LATOM(I) = .TRUE.
C=======================================================================
         DO N1 = -NMAX,NMAX
            DO N2 = -NMAX,NMAX
               DO N3 = -NMAXZ,NMAXZ
C
                  DO J = 1,3
                     VEC1(J) = DBLE(N1)*BRAVAIS(J,1) + DBLE(N2)
     &                         *BRAVAIS(J,2) + DBLE(N3)*BRAVAIS(J,3)
                  END DO
C
C-----------------------------------------------------------------------
                  DO IQ = 1,NAEZ
                     DIFF = 0D0
                     DO J = 1,3
                        VEC2(J) = VEC1(J) + RBASIS(J,IQ)
                        DIFF = DIFF + VEC2(J)**2
                     END DO
                     RMAXGEN = MAX(RMAXGEN,SQRT(DIFF))
C
                     DIFF = SQRT((RCLSNEW(1,I)-VEC2(1))**2
     &                          +(RCLSNEW(2,I)-VEC2(2))**2
     &                          +(RCLSNEW(3,I)-VEC2(3))**2)
C
                     IF ( DIFF.LE.(1D-5) ) THEN
                        ATOMIMP(I) = IQ
                        IF ( AIN(I).NE.IQ ) LATOM(I) = .FALSE.
                        LPOS(I) = .TRUE.
                        GOTO 100
                     END IF

                     IF (DIFF.LE.DIFFMIN(I)) DIFFMIN(I)=DIFF

                  END DO
C-----------------------------------------------------------------------
               END DO
            END DO
         END DO
C=======================================================================
 100  END DO
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      WRITE (1337,99001) RMAXCLUS,RMAXGEN
      WRITE (1337,99002) 
     &               'Input data for impurity sites - consistency check'
C
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      IATOK = 0
      IPOSOK = 0
      OPEN(58,FILE='rimp.dat')
      WRITE(58,*) NATOMIMP,' NATOMIMP'
      DO I = 1,NATOMIMP
         IF ( LPOS(I) ) THEN
            WRITE (STRPOS,'(A6)') 'OK'
         ELSE
            WRITE (STRPOS,'(A6)') 'neq BL'
            IPOSOK = IPOSOK + 1
            WRITE(*,*) 'minimal difference for atom',I,'=',DIFFMIN(I)
         END IF
         WRITE (STRAT,'(I3,A3)') ATOMIMP(I),' <?'
         IF ( .NOT.LATOM(I) ) THEN
            IATOK = IATOK + 1
         ELSE IF ( LPOS(I) ) THEN
            WRITE (STRAT,'(I3)') ATOMIMP(I)
         END IF
         WRITE (1337,99003) I,(RCLSIMP(J,I),J=1,3),AIN(I),STRPOS,STRAT
         WRITE (58,FMT='(I6,3E16.8,I6)') I,(RCLSIMP(J,I),J=1,3),AIN(I)
      END DO
      CLOSE(58)
      WRITE (1337,99004)
C
      IF ( IPOSOK.NE.0 ) THEN
         WRITE (6,99005)
         STOP
      END IF
C
      DO I = 1,NATOMIMP
         IF ( (ATOMIMP(I).GT.NAEZ) .OR. (ATOMIMP(I).EQ.0) ) THEN
            WRITE (6,99008) I,ATOMIMP(I)
            STOP
         END IF
      END DO
C
      IF ( IATOK.NE.0 ) THEN
         WRITE (1337,99006)
      ELSE
         WRITE (1337,99007)
      END IF
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
99001 FORMAT (12X,'input-cluster R  : ',F11.6,/,12X,
     &        'test-cluster  R  : ',F11.6,/)
99002 FORMAT (13X,63('-'),/,15X,A,/,13X,63('-'),/,13X,
     &        ' imp |               READ IN DATA         host  |',
     &        ' CHECKED DATA ',/,13X,
     &        'index|       x           y           z    site  |',
     &        '  pos.   site ',/,13X,63('-'))
99003 FORMAT (13X,I3,2X,'|',3(F12.6),1X,I4,1X,'|',A6,2X,A6)
99004 FORMAT (13X,63('-'))
99005 FORMAT (/,6X,'ERROR: At least one of your input sites does not',
     &        ' belong to the Bravais ',/,13X,'lattice (neq BL). ',
     &        'Please check your input file',/)
99006 FORMAT (13X,'WARNING: At least one inconsistent assignment of ',
     &        'site indices',/,13X,
     &        '         was found in your input. The program will',
     &        ' override the',/,13X,
     &        '         input data. Crosscheck?  ',/,13X,63('-'),/)
99007 FORMAT (13X,'Your cluster data is consistent',/,13X,63('-'),/)
99008 FORMAT (/,6X,'ERROR: Wrong assignment of impurity site ',I3,
     &        ' to the unit-cell site ',I3,/,13X,
     &        'Please check your input file',/)
      END
