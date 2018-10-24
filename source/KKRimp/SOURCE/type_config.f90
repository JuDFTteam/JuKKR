!------------------------------------------------------------------------------------
!> Summary: Config type
!> Author: 
!> Category: KKRimp, input-output
!> Deprecated: False 
!> Contains all information from the config.cfg file
!> Each option of the config file should be contained in this type
!------------------------------------------------------------------------------------
module type_config

! -------------------------
! test and run flags
! -------------------------
integer, parameter                 ::           dim_flags = 20 	!! dimension of testflag and runflag array
character(len=20),dimension(dim_flags)   ::     testflag  = ''  !! testflag array
character(len=20),dimension(dim_flags)   ::     runflag   = ''  !! runflag array

TYPE                              ::  CONFIG_TYPE

! -------------------------
! selfconsistency
! -------------------------
  integer                      ::  icst     = 4   !! number of born iterations
  integer                      ::  ins      = 1   !! non-spherical calculation =1 full pot, =0 asa
  integer                      ::  kvrel    = 1   !! option for the scalar relativistic approximation (=1 sets nsra to 2)
  integer                      ::  nsra     = 2   !! option for the scalar relativistic approximation (=2 for sra)
  integer                      ::  nspin    = 2   !! number of spins
  integer                      ::  kte      = 1   !! 
!   integer                      ::  kxc      = 2
  character(len=20)            ::  modeexcorr      = 'LDA' 	!! exchange correlation mode
  integer                      ::  kshape   = 1			!! = ins
  integer                      ::  kspinorbit = 0		!! spin-orbit coupling
  integer                      ::  ncoll = 0			!! non-collinear calculation
  ! changed default value to save the first 20 wavefunctions. This needs up to
  ! 1GB of additional memory, which should usually be available. By using the
  ! keyword 'WAVEFUNC_RECALC_THRESHHOLD' in the inputcard this can be modified. 
  !integer                      ::  wavefunc_recalc_threshhold=0
  integer                      ::  wavefunc_recalc_threshhold=20  !! Number of stored wavefunctions
! -------------------------
!
! -------------------------
  integer                      ::  npan_log      = 40 		!! number of panels in the log region
  integer                      ::  npan_eq       = 40		!! number of panels in the equidistant region
  integer                      ::  ncheb         = 16		!! probably number of chebyshev nodes
  double precision             ::  npan_logfac   = 2.0D0	!! factor for the generation of the log mesh
real(kind=8)                   ::  rmin   = -1.0D0      	!! first point of the new radial mesh
real(kind=8)                   ::  rlogpan      = 1.0D0         !! radius of the log panel


 INTEGER                       ::  SCFSTEPS = 1			!! number of iterations
! -------------------------
! switch for the calculation of different properties
! -------------------------

  integer                      ::  calcforce= 0			!! calculate the force
  integer                      ::  calcorbitalmoment= 0		!! calculate the orbital moments
  integer                      ::  calcJijmat = 0		!! calculate the magnetic exchange interactions


  integer                      ::  hfield_apply_niter=0		!! number of iterations for a magnetic field (should only be used to force a spin splitting)
  real(kind=8)                 ::  hfield=0.0D0			!! magnitude of the magnetic field

  integer                      ::  hfield_apply_niter2=0	!! number of iterations for a asymmetric magnetic field 
  real(kind=8)                 ::  hfield2(2)=0.0D0		!! asymmetric magnetic field for the two spin channels


! -------------------------
! mixing
! -------------------------
integer                        ::  imix     = 2          	!! mixing scheme, imix=2 (straight), 3 (broy1), 4 (broy2), 5 (anderson)
integer                        ::  NSIMPLEMIXFIRST = 0   	!! number of initial simple mixing steps
integer                        ::  IMIXSPIN   = 0		!! spin mixing (0 straight, 1 broyden) 
real(kind=8)                   ::  SPINMIXFAC = 1.0D0		!! spin mixing factor
integer                        ::  spinmixbound=99999

real(kind=8)                   ::  mixfac   = 0.1      		!! mixing factor 
real(kind=8)                   ::  fcm      = 2.0         	!! 
real(kind=8)                   ::  qbound   = 1d-8     		!! 
integer                        ::  itdbry   = 40          	!! number of iterations which are used for the broyden mixing
! -------------------------
! lattice relaxation
! -------------------------
integer                        ::  lattice_relax   = 0 		!! 



 END TYPE CONFIG_TYPE

end module type_config
