!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Wrapper routine to handle the measurements of the different processes 
!> Author: 
!> Wrapper routine to handle the measurements of the different processes 
!------------------------------------------------------------------------------------
module mod_timing

  ! taken from impurity code of David Bauer
  use :: mod_datatypes, only: dp

  implicit none

  private
  public :: timings_1a, timings_1b, load_imbalance, timing_init, timing_start, timing_stop, timing_pause, print_time_and_date
  ! for adaptive load imbalance tackling
  real (kind=dp), dimension(:,:), allocatable, save :: timings_1a
  real (kind=dp), dimension(:), allocatable, save :: timings_1b
  integer, dimension(:), allocatable, save :: load_imbalance


  integer, parameter :: nkeys = 20
  integer, parameter :: nkeylen = 40
  character (len=nkeylen), dimension(nkeys) :: timingkeys = ''
  integer, dimension(nkeys) :: start_time = 0
  real (kind=dp), dimension(nkeys) :: interm_time = 0.0d0
  integer, dimension(nkeys) :: ispaused = 0
  integer :: init = 0

contains

  !-------------------------------------------------------------------------------
  !> Summary: Initialize the printing of the timing information for each rank 
  !> Author: 
  !> Category: profiling, input-output, KKRhost 
  !> Deprecated: False
  !> Initialize the printing of the timing information for each rank 
  !-------------------------------------------------------------------------------
  subroutine timing_init(my_rank, disable_serial_number)
    use :: mod_types, only: t_inc
    use :: mod_version_info, only: version_print_header
    implicit none
    integer, intent(in) :: my_rank
    character (len=3) :: ctemp
    logical, optional, intent(in) :: disable_serial_number
    logical :: print_serial

    print_serial = .false.
    if (present(disable_serial_number)) then
      if (disable_serial_number) then
        print_serial = .true.
      else
        print_serial = .false.
      end if
    end if

    if (init/=0) stop '[mod_timing] timing already initilized'
    write (ctemp, '(I03.3)') my_rank
    if (t_inc%i_time>0) then
      open (unit=43234059, file='out_timing.'//trim(ctemp)//'.txt')
      call version_print_header(43234059, disable_print=print_serial)
    end if
    init = 1
  end subroutine timing_init

  !-------------------------------------------------------------------------------
  !> Summary: Start the measurement of a given process identified by `mykey2`  
  !> Author: 
  !> Category: profiling, KKRhost 
  !> Deprecated: False
  !> Start the measurement of a given process identified by `mykey2` 
  !-------------------------------------------------------------------------------
  subroutine timing_start(mykey2)
    use :: mod_types, only: t_inc
    implicit none

    character (len=*), intent(in) :: mykey2
    ! .. Local variables
    integer :: ikey
    character (len=nkeylen) :: mykey

    if (t_inc%i_time>0) then

      if (init/=1) stop '[timing_start] please run init first'
  
      mykey = mykey2
  
      ikey = timing_findkey(mykey, 'noerror')
  
      if (ikey==-1) then
        ikey = timing_setkey(mykey)
        interm_time(ikey) = 0.0d0
      else
        if (ispaused(ikey)==1) then
          if (ispaused(ikey)==0) stop '[timing_start] key already present unpaused'
        end if
      end if
  
      call system_clock(count=start_time(ikey)) ! Start timing

    end if

  end subroutine timing_start

  !-------------------------------------------------------------------------------
  !> Summary: Pause the timing of the process described by `mykey2` 
  !> Author: 
  !> Category: profiling, KKRhost 
  !> Deprecated: False
  !> Pause the timing of the process described by `mykey2` 
  !-------------------------------------------------------------------------------
  subroutine timing_pause(mykey2)
    use :: mod_types, only: t_inc
    implicit none

    character (len=*), intent(in) :: mykey2
    ! .. Local variables
    integer :: stop_time
    integer :: ikey
    character (len=nkeylen) :: mykey
    real (kind=dp) :: timing
    integer :: clock_rate

    if (t_inc%i_time>0) then

      mykey = mykey2
      ikey = timing_findkey(mykey)
     
      call system_clock(count_rate=clock_rate) ! Find the rate
      call system_clock(count=stop_time) ! Stop timing
     
      timing = (stop_time-start_time(ikey))/real(clock_rate)
      interm_time(ikey) = interm_time(ikey) + timing
      call system_clock(count=start_time(ikey))
      ispaused(ikey) = 1

    end if

  end subroutine timing_pause

  !-------------------------------------------------------------------------------
  !> Summary: Stop the timing of the process described by `mykey2` 
  !> Author: 
  !> Category: profiling, KKRhost 
  !> Deprecated: False
  !> Stop the timing of the process described by `mykey2` 
  !-------------------------------------------------------------------------------
  subroutine timing_stop(mykey2, save_out)
    use :: mod_types, only: t_inc
    implicit none
    character (len=*), intent (in) :: mykey2
    real (kind=dp), intent (out), optional :: save_out

    integer :: stop_time
    integer :: ikey
    character (len=nkeylen) :: mykey
    real (kind=dp) :: timing
    integer :: clock_rate

    if (t_inc%i_time>0) then

      mykey = mykey2
      ikey = timing_findkey(mykey)
  
      call system_clock(count_rate=clock_rate) ! Find the rate
      call system_clock(count=stop_time) ! Stop timing
  
      if (ispaused(ikey)==1) then
        timing = interm_time(ikey)
      else
        timing = (stop_time-start_time(ikey))/real(clock_rate) + interm_time(ikey)
      end if
  
      call timing_delkey(mykey)
  
      ispaused(ikey) = 0
  
      if (present(save_out)) then
        save_out = timing
        return
      else
        if (t_inc%i_time>0) write (43234059, *) mykey, '  ', timing
      end if

    end if

  end subroutine timing_stop

  !-------------------------------------------------------------------------------
  !> Summary: Find if the current process in in my predefined timing keywords 
  !> Author: 
  !> Category: profiling, KKRhost 
  !> Deprecated: False
  !> Find if the current process in in my predefined timing keywords 
  !-------------------------------------------------------------------------------
  integer function timing_findkey(char1, char2)
    implicit none
    character (len=nkeylen) , intent(in):: char1
    character (len=*), optional, intent(in) :: char2
    integer :: ival

    do ival = 1, nkeys
      if (char1==timingkeys(ival)) then
        timing_findkey = ival
        return
      end if
    end do
    if (.not. present(char2)) then
      write (*, *) timingkeys
      write (*, *) '[timing_findkey] timing key ', trim(char1), ' not found'
      stop
    else
      timing_findkey = -1
    end if
  end function timing_findkey

  !-------------------------------------------------------------------------------
  !> Summary: Set the keyword for the timing processes. 
  !> Author: 
  !> Category: profiling, KKRhost 
  !> Deprecated: False
  !> Set the keyword for the timing processes. Checks if the keyword has been defined
  !> more than once, and if there is enough space in the option array for the present
  !> keyword.
  !-------------------------------------------------------------------------------
  integer function timing_setkey(char1)
    implicit none
    character (len=nkeylen), intent(in) :: char1
    integer :: ival, ifound

    ifound = 0
    do ival = 1, nkeys
      if (timingkeys(ival)==char1) then
        write (*, *) '[timing_setkey] key=', char1
        stop '[timing_setkey] error used twice the same timing key'

      end if
      if (trim(timingkeys(ival))=='' .and. ifound==0) then
        timingkeys(ival) = char1
        timing_setkey = ival
        ifound = 1
      end if
    end do
    if (ifound==0) stop '[timing_setkey] not enough free spots to set timing values'
  end function timing_setkey

  !-------------------------------------------------------------------------------
  !> Summary: Sets the key corresponding to the present keyword to blank, and resets its timing. 
  !> Author: 
  !> Category: profiling, KKRhost 
  !> Deprecated: False
  !> Sets the key corresponding to the present keyword to blank, and resets
  !> its timing. 
  !-------------------------------------------------------------------------------
  subroutine timing_delkey(char1)
    implicit none
    character (len=nkeylen) :: char1
    integer :: ival

    do ival = 1, nkeys
      if (char1==timingkeys(ival)) then
        timingkeys(ival) = ''
        start_time(ival) = 0
        return
      end if
    end do
    stop '[timing_delkey] timing key not found'
  end subroutine timing_delkey

  !-------------------------------------------------------------------------------
  !> Summary: Helper routine that print 'message' with a time stamp to screen
  !> Author: 
  !> Category: profiling, input-output, KKRhost 
  !> Deprecated: False
  !> Helper routine that print 'message' with a time stamp to screen
  !-------------------------------------------------------------------------------
  subroutine print_time_and_date(message)

    implicit none

    character (len=*), intent (in) :: message
    integer, dimension (8) :: values

    call date_and_time(values=values)
    print '(a,a,i4,5(a,i2))', message, ' on ', values(1), '/', values(2), '/', values(3), ' at ', values(5), ':', values(6), ':', values(7)
  end subroutine print_time_and_date


end module mod_timing
