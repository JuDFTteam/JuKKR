!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------
!> Summary: Defines single and double precision kinds
!> Author: P. Ruessmann
!> Date: 2018
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!> 
!> Taken from example in Fortran Modernization workshop of held by Wadud Miah (NAG)
!-------------------------------------------------------------------------------
module mod_datatypes
  use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int64

  implicit none
  ! everything is private unless otherwise stated
  private
  public :: sp, dp, si, di

  integer, parameter :: sp = real32
  integer, parameter :: dp = real64
  integer, parameter :: si = int32
  integer, parameter :: di = int64


end module mod_datatypes
