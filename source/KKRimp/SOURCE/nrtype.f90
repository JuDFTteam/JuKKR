!------------------------------------------------------------------------------------
!> Summary: Contains mathematical constants, single/double real/complex precision kinds
!> Author: 
!> 
!------------------------------------------------------------------------------------
!> @note Notes on the code
!> @endnote
!> @todo things that must be checked
!> @endtodo
!> @warning Important precautions
!> @endwarning
!> @bug If nasty things are found
!> @endbug
!------------------------------------------------------------------------------------
MODULE nrtype
  
  INTEGER, PARAMETER :: WLENGTH = 1      ! For I/O in direct access files; =1 for ifort, =4 for gfort
  INTEGER, PARAMETER :: SP = KIND(1.0)
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
  INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
  REAL(DP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_dp
  REAL(DP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_dp
  REAL(DP), PARAMETER :: PI=3.141592653589793238462643383279502884197_dp
  REAL(DP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_dp
  REAL(DP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_dp
!                 DOUBLE PRECISION,parameter :: CVLIGHT = 274.0720442D0

END MODULE nrtype