!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Wrapper module for the generation of version and serial headers 
!> Author: People who wrote it
!> Wrapper module for the generation of version and serial headers. It helps to track
!> what data was generated with which version of the code.
!------------------------------------------------------------------------------------
!> @note This module can be exported to the other KKR codes, to use them one just
!> needs to comment out the correct pre-compilation flag defined as the begining of 
!> the file.
!> @endnote
!------------------------------------------------------------------------------------
module mod_version_info

  implicit none

  private
  public :: serialnr, construct_serialnr, version_check_header, version_print_header

  character (len=5), parameter :: codename = 'JuKKR'
  character (len=:), allocatable :: serialnr

contains

  !-------------------------------------------------------------------------------
  !> Summary: Take information from version file and create serial number with time stamp
  !> Author: 
  !> Category: input-output, sanity-check, version-control ,KKRhost 
  !> Deprecated: False 
  !> Take information from version file and create serial number with time stamp
  !-------------------------------------------------------------------------------
  subroutine construct_serialnr()
    ! take information from version file and create serial number with time stamp
    use :: mod_version, only: version1, version2
    implicit none

    integer, dimension (8) :: values
    character (len=500) :: tmpname
    integer :: slength, ierr

    ! check date and time when program starts
    call date_and_time(values=values)

    ! write codename, version info and time stamp to serialnr
    write (tmpname, '(6A,I4.4,5I2.2)') trim(codename), '_', trim(version1), '_', trim(version2), '_', values(1), values(2), values(3), values(5), values(6), values(7)
    ! write(tmpname, '(6A,I4.4,5I2.2)') trim(codename), '_', trim(version(1)), '_', trim(version(2)), '_', values(1), values(2), values(3), values(5), values(6), values(7)
    slength = len_trim(tmpname)

    allocate (character(len=slength) :: serialnr, stat=ierr)
    if (ierr/=0) stop '[construct_serialnr] Error allocating serialnr'
    serialnr = trim(tmpname)

  end subroutine construct_serialnr


  !-------------------------------------------------------------------------------
  !> Summary: This is called after an open statement of a file that is written, prints header line
  !> Author:
  !> Category: input-output, sanity-check, version-control, KKRhost 
  !> Deprecated: False 
  !> This is called after an open statement of a file that is written, 
  !> prints header line
  !-------------------------------------------------------------------------------
  subroutine version_print_header(unit, addition, disable_print)

    implicit none
    integer, intent (in) :: unit
    logical, optional, intent (in) :: disable_print
    character (len=*), optional, intent (in) :: addition
    logical :: print_version

    print_version = .true.

    if (present(disable_print)) then
      if (disable_print) print_version = .false.
    end if!


    if (print_version) then
      if (.not. present(addition)) then
        ! write header:             code     version     compver   timestamp
        ! "# serial: kkrjm_v2.0-38-g6593f48_debug_20160907113604"
        write (unit, '(2A)') '# serial: ', serialnr
      else
        write (unit, '(2A)') '# serial: ', serialnr // addition
      end if
    end if

  end subroutine version_print_header

  !-------------------------------------------------------------------------------
  !> Summary: Checks if a header with serial-number is in the first line
  !> Author:
  !> Category: input-output, sanity-check, version-control, KKRhost 
  !> Deprecated: False 
  !> This is called after an open statement of a file that is read
  !> checks if a header with serial-number is in the first line
  !> if not rewinds the file back to start
  !-------------------------------------------------------------------------------
  subroutine version_check_header(unit)

    implicit none

    integer, intent (in) :: unit
    character (len=10) :: first_characters

    read (unit, '(A)') first_characters
    if (first_characters/='# serial: ') then
      rewind (unit)
    end if

  end subroutine version_check_header

  ! subroutine version_potname(unit)

  ! implicit none

  ! integer, intent(in) :: unit



  ! end subroutine version_potname



  ! subroutine version_shapename(unit)

  ! implicit none



  ! end subroutine version_shapename

end module mod_version_info
