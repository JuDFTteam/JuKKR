!-------------------------------------------------------------------------------
! MODULE: Constants
!> @brief Physical and mathematical constants
!> @author Jonathan Chico
!> @date 09.01.2018
!-------------------------------------------------------------------------------
module constants

  implicit none
!.. Scalar parameters
  double precision :: kb_ev = 8.61734d-5 !< Boltzmann constant in eV
  integer, parameter :: nsymaxd = 48
  integer, parameter :: maxmshd = 30
  double precision, parameter :: kb = 0.6333659d-5 !< Boltzmann constant in Ry
  double precision, parameter :: pi = 3.141592653589793d0
  double precision, parameter :: ryd = 13.6058d0 !< Rydbergs in eV
  double precision, parameter :: cvlight = 274.0720442d0 !< Speed of light divided by the fine structure constant
  double complex, parameter :: ci = (0.0d0, 1.0d0) !< Unitary imaginary complex number
  double complex, parameter :: cone = (1.0d0, 0.0d0) !< Unitary real complex number
  double complex, parameter :: czero = (0.0d0, 0.0d0) !< Complex zero initialization
end module
