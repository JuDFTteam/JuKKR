module mod_timing
  integer,parameter      :: nkeys=20
  integer,parameter      :: nkeylen=40
  character(len=nkeylen) :: timingkeys(nkeys)=''
  integer                :: start_time(nkeys)=0
  double precision       :: interm_time(nkeys)=0.0D0
  integer                :: ispaused(nkeys)=0
  integer                :: writetiming=1
  integer                :: init=0
contains 

subroutine timing_init(my_rank)
  use mod_version_info
  implicit none
  integer  :: my_rank
  character(len=3) :: ctemp
  if (init/=0) stop '[mod_timing] timing already initilized'
  write(ctemp,'(I03.3)') my_rank
  open(unit=43234059 , file='out_timing.'//trim(ctemp)//'.txt')
  call version_print_header(43234059)
!   open(unit=432340591, file='out_timing.'//trim(ctemp)//'.txt.test')
  init=1
end subroutine timing_init


subroutine timing_start(mykey2)
  implicit none
  integer ikey
  character(len=*)       :: mykey2
  character(len=nkeylen) :: mykey

if (init/=1) stop '[timing_start] please run init first'

  mykey=mykey2

  ikey=timing_findkey(mykey,'noerror')

  if (ikey==-1) then
    ikey=timing_setkey(mykey)
    interm_time(ikey)=0.0D0
  else
    if  ( ispaused(ikey)==1) then
      if (ispaused(ikey)==0) stop '[timing_start] key already present unpaused'
    end if
  end if

   CALL SYSTEM_CLOCK(COUNT=start_time(ikey)) ! Start timing
   
end subroutine timing_start

subroutine timing_pause(mykey2)
  implicit none
  integer    :: stop_time
  integer ikey
  character(len=*)       :: mykey2
  character(len=nkeylen) :: mykey
  real*8                   :: timing
  integer                   :: clock_rate

  mykey=mykey2
  ikey=timing_findkey(mykey)


  CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) ! Find the rate
  CALL SYSTEM_CLOCK(COUNT=stop_time) ! Stop timing

  timing = (stop_time-start_time(ikey))/real(clock_rate)
  interm_time(ikey)=interm_time(ikey)+timing
  CALL SYSTEM_CLOCK(COUNT=start_time(ikey))
  ispaused(ikey)=1

end subroutine timing_pause

subroutine timing_stop(mykey2)
  implicit none
  integer    :: stop_time
  integer ikey
  character(len=*)       :: mykey2
  character(len=nkeylen) :: mykey
  real*8                   :: timing
  integer                   :: clock_rate 

  mykey=mykey2
  ikey=timing_findkey(mykey)


  CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) ! Find the rate
  CALL SYSTEM_CLOCK(COUNT=stop_time) ! Stop timing


  timing = (stop_time-start_time(ikey))/real(clock_rate)+interm_time(ikey)

  call timing_delkey(mykey)

  write(43234059,*)  mykey,'  ',timing
  ispaused(ikey)=0


end subroutine timing_stop

integer function timing_findkey(char1,char2)
  implicit none
  character(len=nkeylen) :: char1
  character(len=*),optional :: char2
  integer                :: ival
  do ival = 1, nkeys
  if (char1==timingkeys(ival)) then
    timing_findkey=ival
    return
  end if
  end do
  if (.not. present(char2)) then 
    write(*,*) timingkeys
    write(*,*) '[timing_findkey] timing key ',trim(char1),' not found'
    stop 
  else
    timing_findkey=-1
  end if
end function timing_findkey

integer function timing_setkey(char1)
  implicit none
  character(len=nkeylen) :: char1
  integer                :: ival,ifound
  ifound=0
  do ival = 1, nkeys
    if (timingkeys(ival)==char1) then
      write(*,*) '[timing_setkey] key=',char1
      stop '[timing_setkey] error used twice the same timing key'

    end if
    if (trim(timingkeys(ival))=='' .and. ifound==0) then
      timingkeys(ival)=char1
      timing_setkey=ival
      ifound=1
    end if
  end do
  if (ifound==0) stop '[timing_setkey] not enough free spots to set timing values'
end function timing_setkey

subroutine timing_delkey(char1)
  implicit none
  character(len=nkeylen) :: char1
  integer                :: ival

  do ival = 1, nkeys
  if (char1==timingkeys(ival)) then
    timingkeys(ival)=''
    start_time(ival)=0.0D0
    return
  end if
  end do
  stop '[timing_delkey] timing key not found'
end subroutine timing_delkey

end module mod_timing
