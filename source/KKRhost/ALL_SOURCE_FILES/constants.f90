!-------------------------------------------------------------------------------
! MODULE: Constants
!> @brief Physical and mathematical constants
!> @author Jonathan Chico
!> @date 09.01.2018
!-------------------------------------------------------------------------------
module Constants

   implicit none
   !.. Scalar parameters
   double precision :: KB_eV=8.61734d-5                     !< Boltzmann constant in eV
   integer, parameter :: NSYMAXD=48
   integer, parameter :: MAXMSHD=30
   double precision, parameter :: KB=0.6333659D-5           !< Boltzmann constant in Ry
   double precision, parameter :: PI=3.141592653589793d0
   double precision, parameter :: RYD=13.6058D0             !< Rydbergs in eV
   double precision, parameter :: CVLIGHT = 274.0720442D0   !< Speed of light divided by the fine structure constant
   double complex, parameter :: CI=(0.0D0,1.0D0)            !< Unitary imaginary complex number
   double complex, parameter :: CONE=(1.0D0,0.0D0)          !< Unitary real complex number
   double complex, parameter :: CZERO=(0.0D0,0.0D0)         !< Complex zero initialization
end module Constants
