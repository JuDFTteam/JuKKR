!-------------------------------------------------------------------------------
!> Summary: 
!> Author: 
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!>
!> 
!-------------------------------------------------------------------------------
module mod_datatypes
  use, intrinsic :: iso_fortran_env

  implicit none
  ! everything is private unless otherwise stated
  private
  public :: sp, dp, si, di

  integer, parameter :: sp = real32
  integer, parameter :: dp = real64
  integer, parameter :: si = int32
  integer, parameter :: di = int64


end module mod_datatypes
