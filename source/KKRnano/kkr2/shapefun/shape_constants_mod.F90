!------------------------------------------------------------------------------
!> Hardcoded constants for shape function calculation
module shape_constants_mod
  implicit none
  public ! all module vars are constant

  !> Parameter to control diagnostic output, set to 0 for no output
  !! 1 for output as in old program, 2 for additional output about tetrahedra

#ifdef DEBUGSHAPEFUNCTIONS
  integer, parameter :: VERBOSITY = 2
#else
  integer, parameter :: VERBOSITY = 0
#endif

  !> Parameter to enable/disable geometry checks (routine POLCHK)

  logical, parameter :: CHECK_GEOMETRY = .true.

  double precision, parameter :: PI = 3.1415926535897932d0

  !integer, parameter :: NVERTD = 250 !< maximal number of cell vertices
  !integer, parameter :: NFACED = 200 !< maximal number of cell faces
  !integer, parameter :: NVTOTD = NFACED*NVERTD !< maximal number of all tetrahedra vertices
  !integer, parameter :: NVRTD  = 500 !< this constant is used in polchk

  ! TODO: remove
  !> highest quantum number for shape functions (from inc.geometry) why 25? - limits KKR calculation to l_max=6
  integer, parameter :: LMAXD1 = 25

  integer, parameter :: ICD = (((2*LMAXD1+15)*LMAXD1+34)*LMAXD1)/24+1 !< number of compress contraction coefficients
  integer, parameter :: ICED = ((LMAXD1+1)*(LMAXD1+2))/2 !< sum_l=0..lmaxd1 (l+1) number of compressed scaling coefficients
  integer, parameter :: ISUMD = (LMAXD1*(3+4*(LMAXD1+1)*(LMAXD1+2)))/3+1 !< number of compressed rotation matrix elements = sum_l=0..lmaxd1 (2*l+1)^2

  !> maximal number of Gauss-Legendre integration points, used in PINTG
  integer, parameter :: NDIM = 1000
  
end module


