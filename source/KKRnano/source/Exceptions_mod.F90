module Exceptions_mod
!-------------------------------------------------------------------------------
!> Summary: Unified treatment of exceptions
!> Author: Paul F Baumeister
!> Category: KKRnano
!-------------------------------------------------------------------------------
  ! this modules purpose is to make the following functionality public
  use Errors_mod, only: die
  use Warnings_mod, only: launch_warning
  use StringHelpers_mod, only: operator(-), operator(+)
  implicit none
  public !
  
endmodule Exceptions_mod
