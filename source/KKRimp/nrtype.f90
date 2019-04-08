!------------------------------------------------------------------------------------
!> Summary: Contains mathematical constants, single/double real/complex precision kinds
!> Author: 
!> 
!------------------------------------------------------------------------------------
MODULE nrtype
  
  use iso_fortran_env, only: real32, real64
  INTEGER, PARAMETER :: WLENGTH = 1      ! For I/O in direct access files; =1 for ifort, =4 for gfort
  INTEGER, PARAMETER :: SP = real32
  INTEGER, PARAMETER :: DP = real64
  INTEGER, PARAMETER :: SPC = real32
  INTEGER, PARAMETER :: DPC = real64
  REAL(DP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_dp
  REAL(DP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_dp
  REAL(DP), PARAMETER :: PI=3.141592653589793238462643383279502884197_dp
  REAL(DP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_dp
  REAL(DP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_dp
!                 DOUBLE PRECISION,parameter :: CVLIGHT = 274.0720442D0

END MODULE nrtype
