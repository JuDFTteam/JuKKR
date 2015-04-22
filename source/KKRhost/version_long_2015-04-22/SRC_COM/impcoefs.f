C*==impcoefs.f    processed by SPAG 6.05Rc at 15:35 on 18 Oct 2004
      SUBROUTINE IMPCOEFS(NATOMIMP,NAEZ,ATOMIMP,RCLSIMP,NSHELL,NSH1,
     &                    NSH2,RATOM,NSYMAT,ISYMINDEX,ROTNAME,HOSTIMP,
     &                    NATYPD,LMAXD,NSHELD,NSIZE)
C **********************************************************************
C *                                                                    *
C * Writes out the auxiliary file impurity.coefs which is needed for   *
C * impurity calculations                                              *
C * Sets up the array HOSTIMP -- also needed for impurity case         *
C *                                adopted from N. Papanikolaou        *
C **********************************************************************
C
      IMPLICIT NONE
C     .. 
C     .. Scalar arguments
      INTEGER LMAXD,NAEZ,NATOMIMP,NATYPD,NSHELD,NSIZE,NSYMAT
C     ..
C     .. Array arguments
      INTEGER ATOMIMP(*),HOSTIMP(0:NATYPD),ISYMINDEX(*),NSH1(*),NSH2(*),
     &        NSHELL(0:NSHELD)
      DOUBLE PRECISION RATOM(3,NSHELD),RCLSIMP(3,*)
      CHARACTER*10 ROTNAME(*)
C     ..
C     .. Local scalars
      INTEGER AI,I,II,J,NB,NDIM,NHOST,NREP,NS
      DOUBLE PRECISION R1
C     ..
C     .. Local arrays
      INTEGER IMPHOST(NAEZ),NSHOUT(NATOMIMP)
      LOGICAL EXIST(NAEZ)
C     ..
C     
C -->  shells around atom icc are prepared for storing the 
C      cluster-gf in subroutine kkrmat in GMATLL(LMMAXD,LMMAXD,*) 
C     
      DO I = 1,NAEZ
         EXIST(I) = .FALSE.
      END DO
      DO I = 1,NATOMIMP
         EXIST(ATOMIMP(I)) = .TRUE.
      END DO
C
      NHOST = 0
      DO I = 1,NAEZ
         IMPHOST(I) = 0
         IF ( EXIST(I) ) THEN
            NHOST = NHOST + 1
            IMPHOST(I) = NHOST
            HOSTIMP(NHOST) = I
         END IF
      END DO
      HOSTIMP(0) = NHOST
      IF ( NHOST.NE.NAEZ ) WRITE (6,99001)
C
      NREP = 1
      NDIM = 1
C
      DO I = 1,NATOMIMP
         NSHOUT(I) = 1
      END DO
C
!       IF ( OPT('KKRFLEX ') ) THEN
!         OPEN (58,FILE='kkrflex_impurity.coefs',FORM='FORMATTED')
!         DO I = 1,NATOMIMP
!           WRITE (58,99004) (RCLSIMP(J,I),J=1,3)
!         END DO
!         CLOSE (58)
!       END DO
      OPEN (58,FILE='impurity.coefs',FORM='FORMATTED')
      WRITE (58,99002) NREP,NATOMIMP,LMAXD,NATOMIMP,
     &                 (NSHOUT(I),I=1,NATOMIMP)
      WRITE (58,99003)
C-----------------------------------------------------------------------
      DO I = 1,NATOMIMP
C
         R1 = SQRT(RCLSIMP(1,I)**2+RCLSIMP(2,I)**2+RCLSIMP(3,I)**2)
C
         IF ( NAEZ.EQ.NHOST ) THEN
            AI = ATOMIMP(I)
         ELSE
            AI = IMPHOST(ATOMIMP(I))
         END IF
C
         WRITE (58,99004) (RCLSIMP(J,I),J=1,3),AI,I,I,R1,ATOMIMP(I)
      END DO
C-----------------------------------------------------------------------
C
      NB = 0
      DO NS = 1,NSHELL(0)
         NB = NB + NSHELL(NS)
      END DO
C
      WRITE (58,99011) NSIZE,NB
      WRITE (58,99002) NDIM
      WRITE (58,99005) NHOST
      WRITE (58,99006) (HOSTIMP(I),I=1,NHOST)
      WRITE (58,99007) NSYMAT
      WRITE (58,99008) (ROTNAME(ISYMINDEX(I)),I=1,NSYMAT)
      WRITE (58,99009)
      WRITE (58,99002) NSHELL(0)
      WRITE (58,99010) (NS,NSH1(NS),NSH2(NS),(RATOM(II,NS),II=1,3),
     &                 NSHELL(NS),
     &                 SQRT(RATOM(1,NS)**2+RATOM(2,NS)**2+RATOM(3,NS)
     &                 **2),NS=1,NSHELL(0))
      CLOSE (58)
C ======================================================================
99001 FORMAT (8X,'WARNING: Some host atoms are missing in the ',
     &        'impurity cluster',/,8X,
     &        '         Indexing will be changed. Check ',
     &        'impurity.coefs file?',/)
99002 FORMAT (11I5)
99003 FORMAT ('     Position of Impurity            Host Imp Shell',
     &        '   Dist     Host id in Bulk')
99004 FORMAT (3F12.8,I4,I4,I5,F10.6,I5)
99005 FORMAT ('Host order, no of host sites: ',I5)
99006 FORMAT (12I4)
99007 FORMAT (I5,'    Symmetries for the Bulk')
99008 FORMAT (5A10)
99009 FORMAT ('Shells in the reduced format')
99010 FORMAT (3I5,3F12.8,I8,F10.5)
99011 FORMAT (11I20)
      END
