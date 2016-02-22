module type_config

! -------------------------
! test and run flags
! -------------------------
integer, parameter                 ::           dim_flags = 20 ! dimension of testflag array
character(len=20),dimension(dim_flags)   ::     testflag  = ''    ! testflag array
character(len=20),dimension(dim_flags)   ::     runflag   = ''     ! runflag array

TYPE                              ::  CONFIG_TYPE

! -------------------------
! selfconsistency
! -------------------------
  integer                      ::  icst     = 4   ! number of born iterations
  integer                      ::  ins      = 1    ! =1 full pot, =0 asa
  integer                      ::  kvrel    = 1    ! =1 full pot, =0 asa
  integer                      ::  nsra     = 2    ! =1 full pot, =0 asa
  integer                      ::  nspin    = 2    ! =1 full pot, =0 asa
  integer                      ::  kte      = 1
!   integer                      ::  kxc      = 2
  character(len=20)            ::  modeexcorr      = 'LDA'
  integer                      ::  kshape   = 1
  integer                      ::  kspinorbit = 0
  integer                      ::  ncoll = 0
  ! changed default value to save the first 20 wavefunctions. This needs up to
  ! 1GB of additional memory, which should usually be available. By using the
  ! keyword 'WAVEFUNC_RECALC_THRESHHOLD' in the inputcard this can be modified. 
  !integer                      ::  wavefunc_recalc_threshhold=0
  integer                      ::  wavefunc_recalc_threshhold=20
! -------------------------
!
! -------------------------
  integer                      ::  npan_log      = 40
  integer                      ::  npan_eq       = 40
  integer                      ::  ncheb         = 16
  double precision             ::  npan_logfac   = 2.0D0
real(kind=8)                   ::  rmin   = -1.0D0      ! cluster radius
real(kind=8)                   ::  rlogpan      = 1.0D0         ! cluster radius


 INTEGER                       ::  SCFSTEPS = 1
! -------------------------
! switch for the calculation of different properties
! -------------------------

  integer                      ::  calcforce= 0
  integer                      ::  calcorbitalmoment= 0
  integer                      ::  calcJijmat = 0


  integer                      ::  hfield_apply_niter=0
  real(kind=8)                 ::  hfield=0.0D0

  integer                      ::  hfield_apply_niter2=0
  real(kind=8)                 ::  hfield2(2)=0.0D0


! -------------------------
! mixing
! -------------------------
integer                        ::  imix     = 2          ! mix=2 (straight), 3 (broy1), 4 (broy2), 5 (anderson)
integer                        ::  NSIMPLEMIXFIRST = 0   ! number of initial simple mixing steps
integer                        ::  IMIXSPIN   = 0
real(kind=8)                   ::  SPINMIXFAC = 1.0D0
integer                        ::  spinmixbound=99999

real(kind=8)                   ::  mixfac   = 0.1      ! cluster radius
real(kind=8)                   ::  fcm      = 2.0         ! cluster radius
real(kind=8)                   ::  qbound   = 1d-8     ! cluster radius
integer                        ::  itdbry   = 40          ! mix=2 (straight), 3 (broy1), 4 (broy2), 5 (anderson)
! -------------------------
! lattice relaxation
! -------------------------
integer                        ::  lattice_relax   = 0 



 END TYPE CONFIG_TYPE

end module type_config
