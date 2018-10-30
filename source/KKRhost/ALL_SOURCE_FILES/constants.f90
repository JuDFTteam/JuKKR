!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Physical and mathematical constants
!> Author: Jonathan Chico 
!> Physical and mathematical constants
!------------------------------------------------------------------------------------
module mod_constants
  use :: mod_datatypes, only: dp

  implicit none

  private :: dp

  real (kind=dp) :: kb_ev = 8.61734e-5_dp                       !! Boltzmann constant in eV
  real (kind=dp), parameter :: kb = 0.6333659e-5_dp             !! Boltzmann constant in Ry
  real (kind=dp), parameter :: pi = 3.141592653589793e0_dp      !! value of \( \pi \)
  real (kind=dp), parameter :: ryd = 13.6058e0_dp               !! Rydbergs in eV
  real (kind=dp), parameter :: cvlight = 274.0720442e0_dp       !! Speed of light divided by the fine structure constant
  complex (kind=dp), parameter :: ci = (0.0e0_dp, 1.0e0_dp)     !! complex number 'i'
  complex (kind=dp), parameter :: cone = (1.0e0_dp, 0.0e0_dp)   !! complex number '1'
  complex (kind=dp), parameter :: czero = (0.0e0_dp, 0.0e0_dp)  !! complex number '0', usually used for initializations
  integer, parameter :: nsymaxd = 48 !! maximal number of lattice symmetries

end module mod_constants
