SUBROUTINE decipothead(ihost,filehost,ilhost,nathost,vacflag,alat,  &
    bravsys,nq,nt,bravais,efermi,insh,krelh, nspinh,ins,krel,nspin,kmrot)
! **********************************************************************
! *                                                                    *
! * reads in the header of the host potential-file 'decimate.pot'      *
! * checking for consistencies with the actual 2D system               *
! *                                        v.popescu - munich, Dec 04  *
! *                                                                    *
! **********************************************************************
      use mod_version_info
      IMPLICIT NONE
!..
!.. Arguments
      DOUBLE PRECISION ALAT,EFERMI
      CHARACTER*40 FILEHOST
      INTEGER IHOST,ILHOST,INS,INSH,KMROT,KREL,KRELH,NATHOST, &
     &        NQ,NSPIN,NSPINH,NT
      DOUBLE PRECISION BRAVAIS(3,3),BRAVSYS(3,3)
      LOGICAL VACFLAG(2)
!..
!.. Locals
      DOUBLE PRECISION ALATH
      DOUBLE PRECISION DABS
      INTEGER I,IH,IOS,KMROTH
! ----------------------------------------------------------------------
IF (.NOT.(filehost(1:7) == 'vacuum')) THEN
  OPEN (36+ihost,FILE=filehost,STATUS='OLD',IOSTAT=ios)
  CALL version_check_header(36+ihost)
! ......................................................................
  IF ( ios > 0 ) THEN
    WRITE (6,'(/,5X,2A)') 'ERROR: Can not open host file ', filehost(1:ilhost)
    WRITE (6,'(12X,A,A)') 'vacuum    entry should be used to',  &
        ' set a vacuum-host on this side'
    STOP
  END IF
  
  DO ih = 1,2
    READ (36+ihost,*)
  END DO
  READ (36+ihost,99002) krelh,insh,nspinh,kmroth
  READ (36+ihost,'(2(7X,I3),7X,F12.8)') nq,nt,alath
  
! --> non-spherical and/or CPA cases not implemented yet
  
  IF ( insh /= 0 ) STOP ' INS<>0 not implemented '
  IF ( nq /= nt ) STOP ' CPA-host not implemented'
! ......................................................................
  IF ( (nq /= nathost) .OR. (kmroth /= kmrot) .OR.  &
        (DABS(alath-alat) > 1D-6) ) THEN
    WRITE (6,99006) filehost(1:ilhost)
    WRITE (6,99003) '  NAEZ KMROT ALAT '
    WRITE (6,99004) 'syst: ',nathost,kmrot,alat
    WRITE (6,99004) 'host: ',nq,kmroth,alath
    WRITE (6,*)
    STOP
  END IF
  IF ( (insh /= ins) .OR. (nspinh /= nspin) .OR. (krelh /= krel) ) THEN
    WRITE (6,99006) filehost(1:ilhost)
    WRITE (6,99003) '  KREL   INS NSPIN'
    WRITE (6,99005) 'syst: ',krel,ins,nspin
    WRITE (6,99005) 'host: ',krelh,insh,nspinh
    WRITE (6,*)
    STOP
  END IF
! ......................................................................
  READ (36+ihost,'(10X,F10.6)') efermi
  READ (36+ihost,*)
  DO ih = 1,3
    READ (36+ihost,99001) (bravais(i,ih),i=1,3)
  END DO
  alath = 0D0
  DO ih = 1,2
    DO i = 1,2
      alath = alath + bravais(i,ih) - bravsys(i,ih)
    END DO
  END DO
  IF ( DABS(alath) > 1D-6 ) THEN
    WRITE (6,99006) filehost(1:ilhost)
    WRITE (6,99007)
    DO ih = 1,2
      WRITE (6,99008) ih,(bravais(i,ih),i=1,2), (bravsys(i,ih),i=1,2)
    END DO
    WRITE (6,*)
    STOP
  END IF
! ......................................................................
  READ (36+ihost,*)
ELSE
  vacflag(ihost) = .true.
END IF
! ----------------------------------------------------------------------
99001 FORMAT (3F8.4)
99002 FORMAT (8(7X,i3))
99003 FORMAT (14X,a,/,8X,28('-'))
99004 FORMAT (8X,a6,2I6,f10.6)
99005 FORMAT (8X,a6,3I6)
99006 FORMAT (6X,'ERROR: host-potential ( ',a,' )',/,13X,  &
    'not compatible with your input/sytem')
99007 FORMAT (14X,'2D lattice',4X,'host',14X,'system',/,14X,38('-'))
99008 FORMAT (10X,'a_',i1,1X,2F9.5,2X,2F9.5)
END SUBROUTINE decipothead
