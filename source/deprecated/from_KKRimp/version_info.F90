! choose one of the following by uncommenting for use in different codes
!#define host
#define imp
!#define scatter
!#define pkkr
!#define voro


!-------------------------------------------------------------------------------
!> Summary: Module containing routines for version control and serial number generation
!> Author: Philipp Rüßmann
!> Category: KKRimp, version-control
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!>
!-------------------------------------------------------------------------------
module mod_version_info

  implicit none
  
  private
  public serialnr, construct_serialnr, version_check_header, version_print_header
  
#ifdef host
  character(len=5), parameter :: codename='kkrjm'
#endif
#ifdef imp
  character(len=6), parameter :: codename='kkrimp'
#endif
#ifdef scatter
  character(len=10), parameter :: codename='kkrscatter'
#endif
#ifdef pkkr
  character(len=4), parameter :: codename='pkkr'
#endif
#ifdef voro
  character(len=4), parameter :: codename='voronoi'
#endif
  
  character(len=:), allocatable :: serialnr
  

contains

!-------------------------------------------------------------------------------
!> Summary: Take information from version file and create serial number with time stamp
!> Author: Philipp Rüßmann
!> Category: KKRimp, version-control
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!>
!-------------------------------------------------------------------------------
  subroutine construct_serialnr()
    use mod_version, only: version
    implicit none
    integer,dimension(8) :: values
    character(len=500)     :: tmpname
    integer :: slength, ierr
    
    ! check date and time when program starts
    call date_and_time(VALUES=values)
    
    ! write codename, version info and time stamp to serialnr
    write(tmpname, '(6A,I4.4,5I2.2)') trim(codename), '_', trim(version(1)), '_', trim(version(2)), '_', values(1), values(2), values(3), values(5), values(6), values(7)
    slength = len_trim(tmpname)
  
    allocate( character(len=slength) :: serialnr, stat=ierr)
    if(ierr/=0) stop '[construct_serialnr] Error allocating serialnr'
    serialnr = trim( tmpname )
  
  end subroutine construct_serialnr


!-------------------------------------------------------------------------------
!> Summary: Print version header to file
!> Author: Philipp Rüßmann
!> Category: KKRimp, version-control
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!> This is called after an open statement of a file that is written
!> prints header line
!-------------------------------------------------------------------------------
  subroutine version_print_header(unit)
    implicit none
    integer, intent(in) :: unit
    ! write header:             code     version     compver   timestamp
    !               "# serial: kkrjm_v2.0-38-g6593f48_debug_20160907113604"
    write(unit, '(2A)') '# serial: ',serialnr

  end subroutine version_print_header


!-------------------------------------------------------------------------------
!> Summary: Checks for version header in file and rewinds to beginning if not found 
!> Author: Philipp Rüßmann
!> Category: KKRimp, version-control
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!>
!> This is called after an open statement of a file that is read
!> checks if a header with serial-number is in the first line
!> if not rewinds the file back to start
!-------------------------------------------------------------------------------
  subroutine version_check_header(unit)
    implicit none
    integer, intent(in) :: unit
    character(len=10) :: first_characters
    
    read(unit, '(A)') first_characters
    if(first_characters/='# serial: ') then
      rewind(unit)
    end if

  end subroutine version_check_header


end module mod_version_info
