module mod_getscratch

  private
  public :: scratchdir, opendafile

contains

  !-------------------------------------------------------------------------------
  !> Summary: Looks for environment variable SRATCH to store tmat, gmat, gref there
  !> Author: 
  !> Category: KKRhost, input-output
  !> Deprecated: True ! This needs to be set to True for deprecated subroutines
  !>
  !> This routine looks for a system variable SCRATCH which is used 
  !> then as a prefix for tmat, gmat and gref files. Like this these
  !> files can be stored localy on temporary file system
  !>
  !> @warning Call to GETENV commented out, probably not working as intended! @endwarning
  !-------------------------------------------------------------------------------
  subroutine scratchdir(tmpdir, itmpdir, iltmp)

    use :: mod_lngstring, only: lngstring
    implicit none
    character (len=*) :: tmpdir
    integer :: itmpdir, iltmp
    ! EXTERNAL GETENV

    ! CALL GETENV('SCRATCH',TMPDIR)  commented out 15.09.2006 by fivos for
    ! iff820 cluster (with gfortran)
    if (tmpdir==' ') then
      itmpdir = 0
      iltmp = 0
    else
      itmpdir = 1
      iltmp = lngstring(tmpdir, 80)
    end if
    return
  end subroutine scratchdir



  !-------------------------------------------------------------------------------
  !> Summary: Open direct-access (unformatted) file 
  !> Author: 
  !> Category: KKRhost, input-output
  !> Deprecated: True ! This needs to be set to True for deprecated subroutines
  !>
  !> This routine is ment to open DA file BASENAME with prefix TMPDIR
  !> (/TMPDIR/BASENAME). If TMPDIR was not set (ITMPDIR=0), local
  !> running directory will be used.
  !>
  !> @note 
  !> Can be replaced by a straight forward `open` statement since 
  !> the`scratchdir` routine is not used anymore
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine opendafile(iunit, basename, lbasename, lrec, tmpdir, itmpdir, iltmp)

    implicit none
    integer :: iunit, lbasename, lrec, itmpdir, iltmp
    character (len=iltmp) :: tmpdir
    character (len=lbasename) :: basename


    if (itmpdir==1) then
      open (iunit, access='direct', recl=lrec, file=tmpdir(1:iltmp)//basename(1:lbasename), form='unformatted')
    else
      open (iunit, access='direct', recl=lrec, file=basename(1:lbasename), form='unformatted')
    end if
    return
  end subroutine opendafile

end module mod_getscratch
