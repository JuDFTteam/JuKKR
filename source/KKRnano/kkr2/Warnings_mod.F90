#define iounit_t integer
#define status_t integer

module Warnings_mod
implicit none
  private ! default visibility

  public :: launch_warning
  public :: get_number_of_warnings
  public :: show_warning_lines
  public :: test

!   public :: ERROR
!   private :: WARNING
! character(len=*), parameter, public :: ERROR = 'ERROR! ' !! the constant string "ERROR! " for easy grepping

! #ifdef COUNT_WARNINGS
!   ! WARNING(i) is a character(len=10) function, that returns the
!   ! string ' WARNING! ' for i==0 and the count of warnings for i==1
!   integer, private, save :: nWarnings = 0 !! counter for warnings
!   ! Important: the warning counter is only correct, if the output unit is on (o/=0).
! #else
!   character(len=9), parameter :: WARNING(0:0) = ['WARNING! '] !! the constant string WARNING(0)
! #endif
  character(len=*), parameter :: WARNING = 'WARNING! ' !< the constant string WARNING!


  type, private :: warning_line
    character(len=128) :: text
    character(len=32)  :: file
    integer            :: line
    integer            :: number
    double precision   :: time
    integer            :: level
  endtype

  integer, parameter, private :: Max_archived = 8
  integer, private, save      :: num_archived = 0
  integer, private, save      :: num_launched = 0
  
  type(warning_line), private :: archive(0:Max_archived-1)
  
  
  contains

  subroutine launch_warning(text, unit, file, line, level)
    character(len=*), intent(in)           :: text ! text
    iounit_t,         intent(in), optional :: unit ! write warning to this unit if present and > 0
    character(len=*), intent(in), optional :: file ! source file (for debugging)
    integer,          intent(in), optional :: line ! source line (for debugging)
    integer,          intent(in), optional :: level ! severity level
    
    iounit_t          :: unt
    character(len=32) :: fil
    integer           :: lin, idx, lev
    double precision  :: now
    double precision, external :: MPI_Wtime
    
    fil = '?'; if (present(file)) fil = adjustl(file)
    lin = 0; if (present(line)) lin = line
    unt = 0; if (present(unit)) unt = unit
    lev = 0; if (present(level)) lev = level
    if(unt>0) write(unit=unt, fmt='(9(5A,I0))') WARNING, trim(text), ' --- launched at ',trim(fil),':',lin
    
    now = MPI_Wtime() ! now
      
    idx = modulo(num_archived, Max_archived)
    num_archived = num_archived + 1 ! count up
    archive(idx) = warning_line(text, fil, lin, num_archived, now, lev) ! store using default constructor

  endsubroutine ! launch_warning
  
  status_t function show_warning_lines(unit) result(ios)
    iounit_t, intent(in) :: unit ! write warning to this unit if > 0
    
    integer :: idx, num
    double precision :: now
    character(len=9) :: wrn_or_err
    double precision, external :: MPI_Wtime
    
    ios = 0; if (unit < 1) return ! silent

    if (num_archived < 1) then
      write(unit=unit, fmt='(/,A,/)', iostat=ios) 'No Warnings launched!'
      return
    endif
    
    now = MPI_Wtime() ! wall clock time now
    
    write(unit=unit, fmt='(/,A,I0,A,/)') 'Show last ',min(num_archived, Max_archived),' warnings:'
    do num = max(0, num_archived - Max_archived), num_archived - 1
      idx = modulo(num, Max_archived)
      wrn_or_err = 'Warning! ' ; if (archive(idx)%level > 0) wrn_or_err = 'Error!   ' ! for soft errors
      write(unit=unit, fmt='(I8,6A,I0)', iostat=ios) floor(now - archive(idx)%time),' seconds ago: ', &
        wrn_or_err, trim(archive(idx)%text), ' --- at ',trim(archive(idx)%file),':',archive(idx)%line
    enddo ! num
    write(unit=unit, fmt='(A)', iostat=ios) '' ! empty line
    
  endfunction ! show_warning_lines

  !! keep track of the total number of warning
  integer function get_number_of_warnings() result(nw)
    nw = num_launched
  endfunction ! get_number_of_warnings

  status_t function test()
    write(*,*,iostat=test) __FILE__,': start module test:'
    write(*,*,iostat=test) __FILE__,': so far ',get_number_of_warnings(),' warnings have been launched.'
    test = show_warning_lines(unit=6)
        
    call launch_warning('w01')
    call launch_warning('w02', file=__FILE__, line=__LINE__, unit=6)
    call launch_warning('w03')
    call launch_warning('w04')
    call launch_warning('w06', file=__FILE__, line=__LINE__, unit=6)
    call launch_warning('w07')
    test = show_warning_lines(unit=6)
    write(*,*,iostat=test) __FILE__,': by now ',get_number_of_warnings(),' warnings have been launched.'
  endfunction ! test

endmodule ! warnings