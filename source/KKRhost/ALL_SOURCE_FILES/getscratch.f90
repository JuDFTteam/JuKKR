SUBROUTINE scratchdir(tmpdir,itmpdir,iltmp)
! **********************************************************************
! *                                                                    *
! *  This routine looks for a system variable SCRATCH which is used    *
! *  then as a prefix for tmat, gmat and gref files. Like this these   *
! *  files can be stored localy on temporary file system               *
! *                                                                    *
! **********************************************************************

      IMPLICIT NONE
      CHARACTER *(*) TMPDIR
      INTEGER ITMPDIR,ILTMP,LNGSTRING
!     EXTERNAL GETENV

!     CALL GETENV('SCRATCH',TMPDIR)  commented out 15.09.2006 by fivos for iff820 cluster (with gfortran)
IF (tmpdir == ' ') THEN
  itmpdir=0
  iltmp=0
ELSE
  itmpdir=1
  iltmp=lngstring(tmpdir,80)
END IF
RETURN
END SUBROUTINE scratchdir
! **********************************************************************


SUBROUTINE opendafile(iunit,basename,lbasename,lrec,tmpdir, itmpdir,iltmp)
! **********************************************************************
! *                                                                    *
! * This routine is ment to open DA file BASENAME with prefix TMPDIR   *
! * (/TMPDIR/BASENAME). If TMPDIR was not set (ITMPDIR=0), local       *
! * running directory will be used.                                    *
! *                                                                    *
! **********************************************************************
      IMPLICIT NONE
      INTEGER IUNIT,LBASENAME,LREC,ITMPDIR,ILTMP
      character(len=ILTMP) :: TMPDIR
      character(len=LBASENAME) :: BASENAME


IF (itmpdir == 1) THEN
  OPEN (iunit,ACCESS='direct',RECL=lrec,  &
      FILE=tmpdir(1:iltmp)//basename(1:lbasename), FORM='unformatted')
ELSE
  OPEN (iunit,ACCESS='direct',RECL=lrec,  &
      FILE=basename(1:lbasename),FORM='unformatted')
END IF
RETURN
END SUBROUTINE opendafile
