!--------------------------------------------------------------------------------
! Copyright (c) 2018 Forschungszentrum Juelich GmbH, Juelich, Germany
! This file is part of KKRnano and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module Errors_mod
!-------------------------------------------------------------------------------
!> Summary: Unified treatment of deadly errors
!> Author: Paul F Baumeister, Marcel Bornemann
!> Category: KKRnano, input-output
!-------------------------------------------------------------------------------
  use Warnings_mod, only: launch_warning
implicit none
#define iounit_t integer
  private ! default visibility

  public :: die
  public :: treat_errors_as_warnings ! this function sets soft_errors
!   public :: test

#define ALLOW_SOFT_ERRORS

#ifdef  ALLOW_SOFT_ERRORS
  logical, public, protected :: soft_errors = .false.
#else
  logical, public, parameter :: soft_errors = .false.
#endif

  contains

  logical function treat_errors_as_warnings() result(soft)
#ifdef  ALLOW_SOFT_ERRORS
    soft_errors = .true.
#endif
    soft = soft_errors
  endfunction ! treat_errors_as_warnings
  
  subroutine die(text, unit, file, line)
    include 'mpif.h'
    
    character(len=*), intent(in)           :: text ! text
    iounit_t,         intent(in), optional :: unit ! write warning to this unit if present and > 0
    character(len=*), intent(in), optional :: file ! source file (for debugging)
    integer,          intent(in), optional :: line ! source line (for debugging)
    
    iounit_t            :: unt
    character(len=32)   :: fil
    integer             :: lin, ierr
    iounit_t, parameter :: err=0
    
    lin = 0  ; if (present(line)) lin = line
    unt = 0  ; if (present(unit)) unt = unit
    fil = '?'; if (present(file)) fil = adjustl(file)
    
    if (unt > 0) write(unt,'(/,2A,I0,9A)') file,':',lin,' Error: ',trim(text)

    if (soft_errors) then
    
      call launch_warning(text, unt, fil, lin, level=1)
      
    else  ! soft_errors
    
      write( * ,'(/,2A,I0,9A)') file,':',lin,' Error: ',trim(text)
      write(err,'(/,2A,I0,9A)') file,':',lin,' Error: ',trim(text)
      
      call MPI_Abort(MPI_COMM_WORLD, lin, ierr) ! pass sourceline as errorcode
      stop
      
    endif ! soft_errors
    
  endsubroutine ! die

!   status_t function test()
!     write(*,*,iostat=test) __FILE__,' no module test implemented!'
!   endfunction ! test

endmodule ! errors
