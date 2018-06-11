!-------------------------------------------------------------------------------
! MODULE: Constants
!> @brief Physical and mathematical constants
!> @author Jonathan Chico
!> @date 09.01.2018
!-------------------------------------------------------------------------------
    Module constants
      Use mod_datatypes, Only: dp

      Implicit None
      Private :: dp
!.. Scalar parameters
      Real (Kind=dp) :: kb_ev = 8.61734E-5_dp !< Boltzmann constant in eV
      Integer, Parameter :: nsymaxd = 48
      Integer, Parameter :: maxmshd = 30
      Real (Kind=dp), Parameter :: kb = 0.6333659E-5_dp !< Boltzmann constant in Ry
      Real (Kind=dp), Parameter :: pi = 3.141592653589793E0_dp
      Real (Kind=dp), Parameter :: ryd = 13.6058E0_dp !< Rydbergs in eV
      Real (Kind=dp), Parameter :: cvlight = 274.0720442E0_dp !< Speed of light divided by the fine structure constant
      Complex (Kind=dp), Parameter :: ci = (0.0E0_dp, 1.0E0_dp) !< Unitary imaginary complex number
      Complex (Kind=dp), Parameter :: cone = (1.0E0_dp, 0.0E0_dp) !< Unitary real complex number
      Complex (Kind=dp), Parameter :: czero = (0.0E0_dp, 0.0E0_dp) !< Complex zero initialization
    End Module
