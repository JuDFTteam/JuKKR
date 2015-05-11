      SUBROUTINE SCRATCHDIR(TMPDIR,ITMPDIR,ILTMP)
C **********************************************************************
C *                                                                    *
C *  This routine looks for a system variable SCRATCH which is used    *
C *  then as a prefix for tmat, gmat and gref files. Like this these   *
C *  files can be stored localy on temporary file system               *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
      CHARACTER *(*) TMPDIR
      INTEGER ITMPDIR,ILTMP,LNGSTRING
c     EXTERNAL GETENV
C    
c     CALL GETENV('SCRATCH',TMPDIR)  commented out 15.09.2006 by fivos for iff820 cluster (with gfortran)
      IF (TMPDIR.EQ.' ') THEN
          ITMPDIR=0
          ILTMP=0
      ELSE
          ITMPDIR=1
          ILTMP=LNGSTRING(TMPDIR,80)
      ENDIF
      RETURN
      END
C **********************************************************************
      SUBROUTINE OPENDAFILE(IUNIT,BASENAME,LBASENAME,LREC,TMPDIR,
     &                      ITMPDIR,ILTMP)
C **********************************************************************
C *                                                                    *
C * This routine is ment to open DA file BASENAME with prefix TMPDIR   *
C * (/TMPDIR/BASENAME). If TMPDIR was not set (ITMPDIR=0), local       *
C * running directory will be used.                                    *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
      INTEGER IUNIT,LBASENAME,LREC,ITMPDIR,ILTMP
      CHARACTER*80 BASENAME(80),TMPDIR


      IF (ITMPDIR.EQ.1) THEN
          OPEN (IUNIT,ACCESS='direct',RECL=LREC,
     &         FILE=TMPDIR(1:ILTMP)//BASENAME(1:LBASENAME),
     &         FORM='unformatted')
      ELSE
          OPEN (IUNIT,ACCESS='direct',RECL=LREC,
     &          FILE=BASENAME(1:LBASENAME),FORM='unformatted')
      END IF
      RETURN
      END
