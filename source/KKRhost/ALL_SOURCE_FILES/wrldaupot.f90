SUBROUTINE wrldaupot(itrunldau,lopt,ueff,jeff,  &
        erefldau,natyp,wldau,uldau,phildau,  &
        irmd,natypd,nspind,mmaxd,irws)
! **********************************************************************
! *                                                                    *
! * Writes out LDA+U arrays into formatted file 'ldaupot'              *
! *                                                                    *
! **********************************************************************
use mod_version_info
IMPLICIT NONE
!..
INTEGER IRMD,MMAXD,NATYPD,NSPIND,IRWS(NATYPD)
!..
!.. Arguments ..
INTEGER ITRUNLDAU,NATYP
INTEGER LOPT(NATYPD)
DOUBLE PRECISION UEFF(NATYPD),JEFF(NATYPD),EREFLDAU(NATYPD)
DOUBLE PRECISION WLDAU(MMAXD,MMAXD,NSPIND,NATYPD)
DOUBLE PRECISION ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD) 
DOUBLE COMPLEX PHILDAU(IRMD,NATYPD)
!..
!..  Locals 
INTEGER IR,M1,M2,M3,M4,IT,IS
! ======================================================================

OPEN (67,FILE='ldaupot_new',FORM='FORMATTED')
CALL version_print_header(67)
WRITE(1337,99001)
WRITE(67,99002) itrunldau,'    ITRUNLDAU'
WRITE(67,99002) natyp,'    NATYP'
WRITE(67,99003) natyp
WRITE(67,99004) (lopt(it),it=1,natyp)
WRITE(67,99005)
DO it = 1,natyp
  IF ( lopt(it)+1 /= 0 ) WRITE(67,99006) it,ueff(it),jeff(it),erefldau(it)
END DO
DO it = 1,natyp
  IF ( lopt(it)+1 /= 0 ) THEN
    WRITE(67,99002) it,'    WLDAU'
    DO is = 1,nspind
      DO m1 = 1,mmaxd
        WRITE (67,99007) (wldau(m1,m2,is,it),m2=1,mmaxd)
      END DO
    END DO
    WRITE(67,99002) it,'    ULDAU'
    WRITE (67,99007) ((((uldau(m1,m2,m3,m4,it),m4=1,mmaxd),  &
        m3=1,mmaxd),m2=1,mmaxd),m1=1,mmaxd)
    WRITE(67,99002) it,'    PHILDAU'
    WRITE (67,99007) (phildau(ir,it),ir=1,irws(it))
  END IF
END DO
CLOSE (67)

99001 FORMAT(/,5X,'< WRLDAUPOT > : ',  &
    'Writing out LDA+U potential (file ldaupot_new)',/)
99002 FORMAT(i6,a)
99003 FORMAT('LOPT 1..',i3)
99004 FORMAT(16I3)
99005 FORMAT('IAT',6X,'UEFF',12X,'JEFF',12X,'EREF')
99006 FORMAT(i3,3(1X,e15.8))
99007 FORMAT(5E16.8)
END SUBROUTINE wrldaupot
