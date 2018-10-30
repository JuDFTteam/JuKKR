C*==decipothead.f    processed by SPAG 6.05Rc at 14:06 on  3 Dec 2004
      SUBROUTINE DECIPOTHEAD(IHOST,FILEHOST,ILHOST,NATHOST,VACFLAG,ALAT,
     &                       BRAVSYS,NQ,NT,BRAVAIS,EFERMI,INSH,KRELH,
     &                       NSPINH,INS,KREL,NSPIN,KMROT)
C **********************************************************************
C *                                                                    *
C * reads in the header of the host potential-file 'decimate.pot'      *
C * checking for consistencies with the actual 2D system               *
C *                                        v.popescu - munich, Dec 04  *
C *                                                                    *
C **********************************************************************
      use mod_version_info
      IMPLICIT NONE
C     ..
C     .. Arguments
      DOUBLE PRECISION ALAT,EFERMI
      CHARACTER*40 FILEHOST
      INTEGER IHOST,ILHOST,INS,INSH,KMROT,KREL,KRELH,NATHOST,
     &        NQ,NSPIN,NSPINH,NT
      DOUBLE PRECISION BRAVAIS(3,3),BRAVSYS(3,3)
      LOGICAL VACFLAG(2)
C     ..
C     .. Locals
      DOUBLE PRECISION ALATH
      DOUBLE PRECISION DABS
      INTEGER I,IH,IOS,KMROTH
C ----------------------------------------------------------------------
      IF (.NOT.(FILEHOST(1:7).EQ.'vacuum')) THEN
         OPEN (36+IHOST,FILE=FILEHOST,STATUS='OLD',IOSTAT=IOS)
         call version_check_header(36+IHOST)
C ......................................................................
         IF ( IOS.GT.0 ) THEN
            WRITE (6,'(/,5X,2A)') 'ERROR: Can not open host file ',
     &                            FILEHOST(1:ILHOST)
            WRITE (6,'(12X,A,A)') 
     &                           'vacuum    entry should be used to',
     &                           ' set a vacuum-host on this side'
            STOP
         END IF
C
         DO IH = 1,2
            READ (36+IHOST,*)
         END DO
         READ (36+IHOST,99002) KRELH,INSH,NSPINH,KMROTH
         READ (36+IHOST,'(2(7X,I3),7X,F12.8)') NQ,NT,ALATH
C
C --> non-spherical and/or CPA cases not implemented yet
C
         IF ( INSH.NE.0 ) STOP ' INS<>0 not implemented '
         IF ( NQ.NE.NT ) STOP ' CPA-host not implemented'
C ......................................................................
         IF ( (NQ.NE.NATHOST) .OR. (KMROTH.NE.KMROT) .OR. 
     &        (DABS(ALATH-ALAT).GT.1D-6) ) THEN
            WRITE (6,99006) FILEHOST(1:ILHOST)
            WRITE (6,99003) '  NAEZ KMROT ALAT '
            WRITE (6,99004) 'syst: ',NATHOST,KMROT,ALAT
            WRITE (6,99004) 'host: ',NQ,KMROTH,ALATH
            WRITE (6,*)
            STOP
         END IF
         IF ( (INSH.NE.INS) .OR. (NSPINH.NE.NSPIN) .OR. (KRELH.NE.KREL)
     &        ) THEN
            WRITE (6,99006) FILEHOST(1:ILHOST)
            WRITE (6,99003) '  KREL   INS NSPIN'
            WRITE (6,99005) 'syst: ',KREL,INS,NSPIN
            WRITE (6,99005) 'host: ',KRELH,INSH,NSPINH
            WRITE (6,*)
            STOP
         END IF
C ......................................................................
         READ (36+IHOST,'(10X,F10.6)') EFERMI
         READ (36+IHOST,*)
         DO IH = 1,3
            READ (36+IHOST,99001) (BRAVAIS(I,IH),I=1,3)
         END DO
         ALATH = 0D0
         DO IH = 1,2
            DO I = 1,2
               ALATH = ALATH + BRAVAIS(I,IH) - BRAVSYS(I,IH)
            END DO
         END DO
         IF ( DABS(ALATH).GT.1D-6 ) THEN
            WRITE (6,99006) FILEHOST(1:ILHOST)
            WRITE (6,99007)
            DO IH = 1,2
               WRITE (6,99008) IH,(BRAVAIS(I,IH),I=1,2),
     &                         (BRAVSYS(I,IH),I=1,2)
            END DO
            WRITE (6,*)
            STOP
         END IF
C ......................................................................
         READ (36+IHOST,*)
      ELSE
         VACFLAG(IHOST) = .TRUE.
      END IF
C ----------------------------------------------------------------------
99001 FORMAT (3F8.4)
99002 FORMAT (8(7X,I3))
99003 FORMAT (14X,A,/,8X,28('-'))
99004 FORMAT (8X,A6,2I6,F10.6)
99005 FORMAT (8X,A6,3I6)
99006 FORMAT (6X,'ERROR: host-potential ( ',A,' )',/,13X,
     &        'not compatible with your input/sytem')
99007 FORMAT (14X,'2D lattice',4X,'host',14X,'system',/,14X,38('-'))
99008 FORMAT (10X,'a_',I1,1X,2F9.5,2X,2F9.5)
      END
