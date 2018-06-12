! -------------------------------------------------------------------------------
! MODULE: Constants
! > @brief Physical and mathematical constants
! > @author Jonathan Chico
! > @date 09.01.2018
! -------------------------------------------------------------------------------
module constants
  use :: mod_datatypes, only: dp

  implicit none
  private :: dp
  ! .. Scalar parameters
  real (kind=dp) :: kb_ev = 8.61734e-5_dp ! < Boltzmann constant in eV
  integer, parameter :: nsymaxd = 48
  integer, parameter :: maxmshd = 30
  real (kind=dp), parameter :: kb = 0.6333659e-5_dp ! < Boltzmann constant in
                                                    ! Ry
  real (kind=dp), parameter :: pi = 3.141592653589793e0_dp
  real (kind=dp), parameter :: ryd = 13.6058e0_dp ! < Rydbergs in eV
  real (kind=dp), parameter :: cvlight = 274.0720442e0_dp ! < Speed of light
                                                          ! divided by the
                                                          ! fine structure
                                                          ! constant
  complex (kind=dp), parameter :: ci = (0.0e0_dp, 1.0e0_dp) ! < Unitary
                                                            ! imaginary
                                                            ! complex number
  complex (kind=dp), parameter :: cone = (1.0e0_dp, 0.0e0_dp) ! < Unitary real
                                                              ! complex number
  complex (kind=dp), parameter :: czero = (0.0e0_dp, 0.0e0_dp) ! < Complex
                                                               ! zero
                                                               ! initialization
end module constants
