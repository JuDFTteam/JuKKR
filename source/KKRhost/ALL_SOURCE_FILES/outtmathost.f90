SUBROUTINE outtmathost(alat,ins,krel,kmrot,nspin,naez,lmmax,  &
        bravais,rbasis,qmtet,qmphi,e2in,tk,  &
        npol,npnt1,npnt2,npnt3)
! **********************************************************************
! *                                                                    *
! *  Writes out the header of the t-matrices decimation file           *
! *                                                                    *
! **********************************************************************
      use mod_version_info
      IMPLICIT NONE
!..
!.. Arguments ..
      INTEGER INS,KREL,KMROT,NSPIN,NAEZ,LMMAX,NPOL,NPNT1,NPNT2,NPNT3
      DOUBLE PRECISION ALAT,E2IN,TK
      DOUBLE PRECISION BRAVAIS(3,3),RBASIS(3,*),QMTET(*),QMPHI(*)
!..
!.. Locals ..
      INTEGER I,IH
! ----------------------------------------------------------------------
WRITE(1337,'(5X,A,/)') '< DECIOPT > : writing header of decimation file'

OPEN (37,FILE='decifile',STATUS='unknown')
CALL version_print_header(37)
WRITE (37,FMT=*) 'INVERSE T-MATRIX AND CMOMS'
WRITE (37,FMT=99001)
WRITE (37,FMT=99002) alat,nspin,naez,lmmax,ins,krel,kmrot
WRITE (37,FMT=99003) bravais
IF ( krel == 0 ) THEN
  WRITE (37,FMT=99004)
  DO ih = 1,naez
    WRITE (37,FMT=99006) (rbasis(i,ih),i=1,3)
  END DO
ELSE
  WRITE (37,FMT=99005)
  DO ih = 1,naez
    WRITE (37,FMT=99007) (rbasis(i,ih),i=1,3), qmtet(ih),qmphi(ih)
  END DO
END IF
WRITE (37,FMT=99008) e2in,tk
WRITE (37,FMT=99009) npnt1,npnt2,npnt3,npol
CLOSE (37)
! ----------------------------------------------------------------------
99001 FORMAT (' Vectors in lattice constant units',/,  &
    '                                 ')
99002 FORMAT ('ALAT=',f9.6,' NSPIN=',i2,'  NAEZ=',i3,' LMMAX=',i3,  &
    ' INS=',i1,' KREL=',i1,' KMROT=',i1)
99003 FORMAT ('BRAVAIS ',/,3F8.4,/,3F8.4,/,3F8.4)
99004 FORMAT ('RBASIS')
99005 FORMAT ('RBASIS',20X,'MAGNETISATION ANGLES THETA/PHI')
99006 FORMAT (3F8.4)
99007 FORMAT (3F8.4,2F9.4)
99008 FORMAT ('EF=',f10.6,' TEMP=',f10.4,' Kelvin')
99009 FORMAT ('N1=',i3,' N2=',i3,' N3=',i3,' NPOL=',i3)
END SUBROUTINE outtmathost
