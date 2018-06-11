    Subroutine scratchdir(tmpdir, itmpdir, iltmp)
! **********************************************************************
! *                                                                    *
! *  This routine looks for a system variable SCRATCH which is used    *
! *  then as a prefix for tmat, gmat and gref files. Like this these   *
! *  files can be stored localy on temporary file system               *
! *                                                                    *
! **********************************************************************

      Use mod_datatypes, Only: dp
      Implicit None
      Character (Len=*) :: tmpdir
      Integer :: itmpdir, iltmp, lngstring
!     EXTERNAL GETENV

!     CALL GETENV('SCRATCH',TMPDIR)  commented out 15.09.2006 by fivos for iff820 cluster (with gfortran)
      If (tmpdir==' ') Then
        itmpdir = 0
        iltmp = 0
      Else
        itmpdir = 1
        iltmp = lngstring(tmpdir, 80)
      End If
      Return
    End Subroutine
! **********************************************************************


    Subroutine opendafile(iunit, basename, lbasename, lrec, tmpdir, itmpdir, &
      iltmp)
! **********************************************************************
! *                                                                    *
! * This routine is ment to open DA file BASENAME with prefix TMPDIR   *
! * (/TMPDIR/BASENAME). If TMPDIR was not set (ITMPDIR=0), local       *
! * running directory will be used.                                    *
! *                                                                    *
! **********************************************************************
      Use mod_datatypes, Only: dp
      Implicit None
      Integer :: iunit, lbasename, lrec, itmpdir, iltmp
      Character (Len=iltmp) :: tmpdir
      Character (Len=lbasename) :: basename


      If (itmpdir==1) Then
        Open (iunit, Access='direct', Recl=lrec, File=tmpdir(1:iltmp)// &
          basename(1:lbasename), Form='unformatted')
      Else
        Open (iunit, Access='direct', Recl=lrec, File=basename(1:lbasename), &
          Form='unformatted')
      End If
      Return
    End Subroutine
