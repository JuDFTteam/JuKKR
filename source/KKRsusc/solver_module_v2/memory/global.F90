  module global
! Global parameters and storage


  implicit none


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!                          Kind definitions
!        different compilers might handle kind values differently
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! short int
  integer, private :: idummy
  integer, public, parameter :: i4b = kind(idummy)
! single prec
  real,    private :: rdummy
  integer, public, parameter :: r4b = kind(rdummy)
  complex, private :: cdummy
  integer, public, parameter :: c4b = kind(cdummy)
! double prec
  double precision, private  :: ddummy
  integer, public, parameter :: r8b = kind(ddummy)
  double complex,   private  :: zdummy
  integer, public, parameter :: c8b = kind(zdummy)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!                    Unit numbers and input file
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! --> unit number for main I/O
  integer(kind=i4b), parameter   :: iomain = 99
! --> Gmat and tmat files
  integer(kind=i4b), parameter   :: iomain1 = 1333 
  integer(kind=i4b), parameter   :: iomain2 = 1666
! --> same, for debugging
  integer(kind=i4b), parameter :: iodb = 666
! --> unit number for opening, reading/writing and closing a file in a single subroutine
  integer(kind=i4b), parameter :: iofile  = 300
  integer(kind=i4b), parameter :: iofile2 = 301
  integer(kind=i4b), parameter :: iofile3 = 302
  integer(kind=i4b), parameter :: iofile4 = 303
! maximum number of characters in one line
  integer(kind=i4b), parameter :: nchars = 1000
! maximum number of lines in input file
  integer(kind=i4b) :: nlines
! storage for contents of input file
  character(len=nchars), allocatable :: inputfile(:)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!              Calculation options set by the input file
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!                 General parameters of calculation
!-----------------------------------------------------------------------
! --> input warnings
  logical           :: lwarn = .false.
! --> which KKR program: 1 for old and 2 for new
  integer(kind=i4b) :: ikkr = 1
! --> DOS calculation:
  logical           :: ldos = .false.
! --> SOC correction:
  logical           :: lsoc = .false.
! --> SOC from host:
  logical           :: lsoc_new = .false.
! --> LDA+U correction:
  logical           :: lldau = .false.
! --> external magnetic field
  logical           :: lbfield = .false.
! --> susceptibility calculation:
  logical           :: lsusc = .false.
! --> spin rotations:
  logical           :: lrot = .false.
! --> energy fit to the GF:
  logical           :: lfit = .false. 
! --> disk I/O (pots, wfns, outsusc.dat)
  logical           :: lhdio = .false.
! --> Lichtenstein formula
  logical           :: ljij = .false.
! --> Model self-energy
  logical           :: lsemodel = .false.
! --> Fermi energy
  real(kind=r8b)    :: efscf = 1.d99
! --> Restart mode 
  logical           :: lrestart = .false.
! --> Serial mode
  logical           :: lserial = .true.

!-----------------------------------------------------------------------
!                         Basis construction
!-----------------------------------------------------------------------
! --> basis type (see subroutine new_basis)
  integer(kind=i4b) :: ibasis = 2
! --> orthogonalization method: 1 for Gram-Schmidt, 2 for overlap matrix
  integer(kind=i4b) :: ibasismethod = 2
! --> tolerance to chop off basis functions
  real(kind=r8b)    :: basistol = 1.d-3
! --> whether to recompute the basis at each iteration or not
  logical           :: lbasisfreeze = .false.
! --> whether to project input t-matrix in basis
  logical           :: ltmatproj = .false.
!-----------------------------------------------------------------------
!                         Density of states
!-----------------------------------------------------------------------
! --> DOS options: 1 for projection energies, 2 for straight line
  integer(kind=i4b) :: idos = 1
! --> number of energy points: only effective if fitting is used
  integer(kind=i4b) :: nedos = 0
! --> start, end: only effective if fitting is used
  complex(kind=c8b) :: dose0 = 0.d0, dose1 = 0.d0
! --> analytical continuation to real axis
  logical           :: ldosacon = .false.
! --> how many energies to use for rational function extrapolation
  integer(kind=i4b) :: nedosacon = 0
! --> energy zero
  real(kind=r8b)    :: dosezero = 0.d0
! --> conversion factor for energy
  real(kind=r8b)    :: dosefac = 1.d0
! --> density matrix from (lm,l'm') blocks of GF
  logical           :: ldosdmat = .false.
!-----------------------------------------------------------------------
!                              LDA+U       
!-----------------------------------------------------------------------
! --> read in starting density matrices:
  logical           :: lrhomat = .false.
! --> straight mixing for density matrices:
  real(kind=r8b)    :: ldaumix = 0.d0
! --> after which iteration to start mixing:
  integer(kind=i4b) :: ldauitc = 0
!-----------------------------------------------------------------------
!                          Magnetic field  
!-----------------------------------------------------------------------
! --> handle constraining fields:
  logical           :: lbconst = .false.
!-----------------------------------------------------------------------
!                     Dynamical susceptibility     
!-----------------------------------------------------------------------
! --> dynamical susceptibility:
  logical           :: ldynsusc = .false.
! --> number of frequency points
  integer(kind=i4b) :: nomega = 0
! --> start, end and step of real frequency mesh
  complex(kind=c8b) :: omegamin = 0.d0, omegamax = 0.d0
  real(kind=r8b)    :: domega = 0.d0
! --> Hartree kernel:
  logical           :: lkha = .false.
! --> xc kernel:
  logical           :: lkxc = .false.
! --> cartesian or spin representation
  logical           :: lcartesian = .false.
! --> compute analytic and nonanalytic integrals
  logical           :: lanalytic  = .false., lnonanalytic = .false.
! --> KS or full susc
  logical           :: lenhanced  = .false.
! --> iterations to correct the denominator
  integer(kind=i4b) :: itermax = 0
! --> mixing for lambda min
  real(kind=r8b)    :: lambdamix = 0.d0
! --> spin sum rule
  logical           :: lsumrule = .false.
! --> spin sum rule for noncollinear case
  logical           :: lnewsumrule = .false.
! --> spin-current spin correlation function (added by Sascha)
  logical           :: lcurrcorr = .false.
! --> interpolation of spin-current spin correlation function (added by Sascha)
  logical           :: lcurrcorrint = .false.
! --> divergence of spin-current spin correlation function (added by Sascha)
  logical           :: lcurrcorrdiv = .false.
! --> number of x,y,z points on interpolation grid 
  integer(kind=i4b) :: n_int
! --> scalar relativisct mass correction
  logical           :: lscalarcorr = .false.
!-----------------------------------------------------------------------
!                            Spin rotations     
!-----------------------------------------------------------------------
! --> rotation type: 0 for global rotation, 1 for local rotations
  integer(kind=i4b) :: ispinrot = 0
! --> all spins rotated to this direction
  real(kind=r8b)    :: urot(3) = (/0.d0,0.d0,1.d0/)
! --> mixing for local spin directions from torques
  real(kind=r8b)    :: dirmix2 = 0.d0
! --> mixing for local spin directions from output orientations
  real(kind=r8b)    :: dirmix3 = 0.d0
! --> input GF in collinear ASA form?
  logical           :: lgrefsph = .false.
!-----------------------------------------------------------------------
!                       Energy fitting of the GF     
!-----------------------------------------------------------------------
! --> fit options: 1 for rational fit
  integer(kind=i4b) :: ifit = 0
! --> max degree of numerator and denominator in rational function
  integer(kind=i4b) :: numd = 0, dend = 0
! --> shift of energy argument
  complex(kind=c8b) :: eshift = (0.d0,0.d0)
! --> constant shift of Re GF; fudge factor
  logical           :: lregf = 0
  real(kind=r8b)    :: fudge = 0.d0
!-----------------------------------------------------------------------
!                          Jij tensor
!-----------------------------------------------------------------------
! --> T for Ebert's formula, F for Lichtenstein formula
  logical           :: ljijtensor = .false.
! --> correction from EF due to constant Ne
  logical           :: ljijef = .false.
! --> add onsite correction
  logical           :: ljionsite = .false.
!-----------------------------------------------------------------------
!                          Groundstate switches     
!-----------------------------------------------------------------------
! --> diagonalize density matrix
  logical           :: lrhodiag = .false.
! --> use onsite GF in groundstate
  logical           :: lgsonsite = .true.
! --> use structural GF in groundstate
  logical           :: lgsstruct = .true.
! --> update charge density
  logical           :: lnewdensity = .true.
! --> energy mesh for integration from file
  logical           :: lscfmesh = .false.
! --> replace structural GF with free space GF
  logical           :: lfreegf = .false.
! --> exclude Fermi energy from band energy definition
  logical           :: lebandnoef = .false.
! --> read atomic positions from file
  logical           :: lpositions = .false.
! --> kill spin density in selected atoms
  logical           :: lnobxc = .false.
! --> calculate groundstate current
  logical           :: lcurrent = .false.
  logical           :: lcurrentsoc = .false.
  logical           :: lcurrentzeeman = .false.
  logical           :: lcurrentint = .false.
  logical           :: lcurrentoutput = .false.
!-----------------------------------------------------------------------
!                       Atom-dependent parameters     
!-----------------------------------------------------------------------
! --> number of energy points for susc output
  integer(kind=i4b) :: nesusc = 0
! --> number of atoms for projection
  integer(kind=i4b) :: nasusc = 0
! --> number of atoms for susceptibility
  integer(kind=i4b) :: nasusc2 = 0
! --> max number of radial basis functions per l, s and atom, for GF; eventually reduced
  integer(kind=i4b) :: nbmax = 0
! --> max angular momentum
  integer(kind=i4b) :: nlmax = -1
! --> small components of SRA: 0 for no, 1 for yes -- NOT CODED
  integer(kind=i4b) :: isra  = -1
! --> max number of points in radial meshes
  integer(kind=i4b) :: nrmax = 0
! --> max number of spin components
  integer(kind=i4b) :: nsmax = 0
! --> kill bxc option: 0 for no, 1 for yes
  integer(kind=i4b), allocatable :: inobxc(:)
! --> atom-dependent density matrix switch: see ldosdmat above
  integer(kind=i4b), allocatable :: iadmat(:)
! --> orbital and spin-dependent density matrix switch: see ldosdmat above
  integer(kind=i4b), allocatable :: ildmat(:,:,:)
! --> atom-dependent SOC switch: see lsoc above
  integer(kind=i4b), allocatable :: isoc(:)
! --> atom-dependent SOC scaling
  real(kind=r8b),    allocatable :: socscaling(:)
! --> atom-dependent SOC control flag
  logical,           allocatable :: soc_applied(:)
! --> atom-dependent LDA+U switch: see lldau above
  integer(kind=i4b), allocatable :: ildau(:)
! --> atom-dependent and l-dependent effective U parameters
  real(kind=r8b),    allocatable :: ueff(:,:)
! --> atom-dependent and l-dependent effective J parameters
  real(kind=r8b),    allocatable :: jeff(:,:)
! --> atom-dependent Jij switch: see ljij above
  integer(kind=i4b), allocatable :: ijij(:)
! --> atom-dependent external field: 0 no field, 1 spin Zeeman, 2 orb Zeeman, 3 full Zeeman
  integer(kind=i4b), allocatable :: ibfield(:)
! --> atom-dependent strength and direction of external field
  real(kind=r8b),    allocatable :: blen(:), bdir(:,:)
! --> atom-dependent strength and direction of constraining field
  real(kind=r8b),    allocatable :: bconlen(:), bcondir(:,:)
! --> susceptibility option: static if < 0, dynamic if > 0
!     1 transverse, 2 longitudinal, 3 full
  integer(kind=i4b), allocatable :: isusc(:)
! --> Hartree kernel option: 1 yes, 0 no
  integer(kind=i4b), allocatable :: ikha(:)
! --> xc kernel option: 1 transverse, 2 longitudinal, 3 both, 0 none
  integer(kind=i4b), allocatable :: ikxc(:)
! --> spin-dependent radial basis: 1 for same ref energies, 2 for independent ref energies
  integer(kind=i4b), allocatable :: issusc(:)
! --> which wfns to use for projection (s,p,d,f) for each atom and spin channel
!     n = 0 means don't project
!     n > 0 means project using n regular solutions computed at ewsusc
  integer(kind=i4b), allocatable :: iwsusc(:,:,:)   ! l only
! --> energy for wfn calculation for projection
  complex(kind=c8b), allocatable :: ewsusc(:,:,:,:)
! --> this checks if the parameters were initialized -- set by init_param
  logical           :: noparameters = .true.
! --> Number of panels > 1
  integer(kind=i4b)              :: npanmax  
  integer(kind=i4b), allocatable :: npanat(:)
  integer(kind=i4b), allocatable :: ircutat(:,:)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!              Variables and storage for the basis module
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! --> non-zero Gaunt numbers
  real(kind=r8b),    parameter :: ylmtol = 1.d-10
! --> checks if the Gram-Schmidt process has found a zero
  real(kind=r8b),    parameter :: gstol = 1.d-8
! --> radial mesh and weights for radial integration with reference wfns
! --> the weight must ensure that the initial wnfs are integrable
  integer(kind=i4b), allocatable :: nrpts(:), nrpts0(:), nrpts1(:)
  real(kind=r8b),    allocatable :: rmesh(:,:), rsmesh(:,:,:), drmesh(:,:), drproj(:,:)
  logical,           allocatable :: normesh(:)
! --> energy mesh for projection of the GF
  complex(kind=c8b), allocatable :: ensusc(:)
! --> 2*nbmax, eventually reduced
  integer(kind=i4b) :: nbmax2 = 0
! --> 2*nsmax
  integer(kind=i4b) :: nsmax2 = 0
! --> maximum angular momentum in Lebedev Laikov mesh
  integer(kind=i4b) :: nlmaxll = 0, lmmaxll = 0
! --> number of angular mesh points
  integer(kind=i4b) :: nll = 0
! --> 2*lmax, 4*lmax
  integer(kind=i4b) :: nlmax2 = -1, nlmax4 = -1
! --> (lmax+1)^2, etc.
  integer(kind=i4b) :: lmmax = 0, lmmax2 = 0, lmmax4 = 0
! --> total number of angular momentum components for susceptibility
  integer(kind=i4b) :: nlmax0 = -1, lmmax0 = 0
! --> angular mesh, spherical harmonics, gaunt numbers
  real(kind=r8b),    allocatable :: ull(:,:), wll(:), cthetall(:), phill(:), ylm(:,:), rgaunt(:,:,:)
! --> angular momentum matrices
  complex(kind=c8b), allocatable :: lorb(:,:,:)
! --> matrix elements for L.S
  complex(kind=c8b), allocatable :: ldots(:,:,:)
! --> angular momentum and spin
  integer(kind=i4b) :: nlms  = 0
! --> angular momentum, spin and basis functions
  integer(kind=i4b) :: nlmsb = 0
! --> dimensions for storage of structural GF
  integer(kind=i4b) :: nalms = 0
! --> storage for reference regular wfns
  real(kind=r8b),    allocatable :: phiref(:,:,:,:,:), overlap(:,:,:)
  logical,           allocatable :: nowfns(:,:,:), nobasis(:,:,:)
! --> basis functions for the density
  integer(kind=i4b), allocatable :: iwsusc2(:,:,:)  ! lm
! --> storage for basis functions for density
  real(kind=r8b),    allocatable :: suscbasis(:,:,:,:,:), suscnorm(:)
! --> gradient basis function
  complex(kind=c8b), allocatable :: gradnorm(:,:,:)
  complex(kind=c8b), allocatable :: gradbasis_lm(:,:,:,:)
  complex(kind=c8b), allocatable :: gfnorm(:,:)
! --> Gaunt-like coefficients for density: overlaps * Gaunt numbers
  real(kind=r8b),    allocatable :: dengaunt(:,:,:)
! --> pointers from spin labels to storage and back
  integer(kind=i4b) :: is2i(2,2), i2is(2,4)
! --> pointers from lm to storage and back (gaunt numbers)
  integer(kind=i4b), allocatable :: lm2i(:,:), i2lm(:,:)
! --> pointers from lms to storage and back
  integer(kind=i4b), allocatable :: lms2i(:,:), i2lms(:,:), lms2i_new(:,:), i2lms_new(:,:)
! --> pointers from lmsb to storage and back
  integer(kind=i4b), allocatable :: lmsb2i(:,:,:,:), i2lmsb(:,:,:)
! --> sizes of blocks of projected GF
  integer(kind=i4b), allocatable :: nlmsba(:)
! --> pointers from alms to storage and back
  integer(kind=i4b), allocatable :: alms2i(:,:), i2alms(:,:)
! --> pointers for almsb (product basis for susceptibility)
  integer(kind=i4b) :: ngfsum, ngfmax
  integer(kind=i4b), allocatable :: nalmsbgf(:), i2almsbgf(:,:), almsbgf2i(:,:,:)
! --> pointers for almsb (product basis with gradient for correlation function)
  integer(kind=i4b) :: ngradsum, ngradmax
  integer(kind=i4b), allocatable :: nalmsbgrad(:), i2almsbgrad(:,:), almsbgrad2i(:,:,:)
! --> pointers for almsb (density basis for susceptibility)
  integer(kind=i4b) :: ndensum, ndenmax
  integer(kind=i4b), allocatable :: nalmsbden(:), i2almsbden(:,:)
! --> Pauli matrices
  complex(kind=c8b) :: pauli(2,2,4)
! --> transformation matrices from cartesian to spin labels and back
! potential from cartesian to spin
  complex(kind=c8b) :: pc2s(2,2,4)
! potential from spin to cartesian
  complex(kind=c8b) :: ps2c(4,2,2)
! density from cartesian to spin
  complex(kind=c8b) :: dc2s(2,2,4)
! density from spin to cartesian
  complex(kind=c8b) :: ds2c(4,2,2)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!           Variables and storage for the projection module
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! --> checks if the structural GF element is to be set to zero
  real(kind=r8b),    parameter :: gfilter = 1.d-6
! --> energy mesh
  complex(kind=c8b), allocatable :: esusc(:), eksusc(:), desusc(:), wsusc(:)
! --> coefficients of the projected scattering solutions
  complex(kind=c8b), allocatable :: pzc(:,:,:,:), fzc(:,:,:,:)
! --> the coefficients of the products of regular and irregular wfns
! --> these are computed in the impurity program, in ASA mode
  complex(kind=c8b), allocatable :: pqc(:,:,:,:,:), psc(:,:,:,:,:), fqc(:,:,:,:,:), fsc(:,:,:,:,:)
  logical,           allocatable :: noregcoeffs(:,:), noirrcoeffs(:,:)
  logical,           allocatable :: noregcoeffs_soc(:,:,:), noirrcoeffs_soc(:,:,:)
! --> storage for the t-matrices
  complex(kind=c8b), allocatable :: tmatrix(:,:,:,:)
! --> storage for the structural GF
  complex(kind=c8b), allocatable :: gstruct(:,:,:)
! --> storage for fitted GF
  complex(kind=c8b), allocatable :: gffit(:,:,:,:,:)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!         Variables and storage for the groundstate module      
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! --> non-zero values in the observables
  real(kind=r8b),    parameter :: atol = 1.d-8
! This drops small elements of the SOC wfn
  real(kind=r8b),    parameter :: soctol = 1.d-8
! This drops small elements of the LDA+U wfn
  real(kind=r8b),    parameter :: ldautol = 1.d-8
! This drops small elements of the self-energy wfn
  real(kind=r8b),    parameter :: sigmatol = 1.d-8
! --> atomic number
  real(kind=r8b),    allocatable :: zat(:)
! --> atomic positions
  real(kind=r8b),    allocatable :: ri(:,:)
! --> assumed direction of magnetization
  real(kind=r8b),    allocatable :: magdir(:,:)
! --> computed direction of magnetization
  real(kind=r8b),    allocatable :: newdir(:,:)
! --> how to update local spin directions
  integer(kind=i4b), allocatable :: iarot(:)
! --> spin projections
  complex(kind=c8b), allocatable :: spinproj(:,:,:)
! --> scalar and magnetic potentials (nuclear pot not included)
  real(kind=r8b),    allocatable :: vr(:,:), br(:,:)
! --> potential in the lmsb product and density basis
  complex(kind=c8b), allocatable :: vlmsbgf(:,:,:), vlmsbden(:)
! --> core charge and magnetisation densities
  real(kind=r8b),    allocatable :: nrc(:,:), mrc(:,:)
! --> l and spin decomposed valence band energy
  real(kind=r8b),    allocatable :: ebandv(:,:,:)
! --> l and spin decomposed LDA+U energy
  real(kind=r8b),    allocatable :: eldau(:,:,:)
! --> l and spin chemical potential shift
  real(kind=r8b),    allocatable :: vshift(:,:,:)
! --> l and spin decomposed torque anisotropy
  real(kind=r8b),    allocatable :: etorque(:,:,:)
! --> matrix elements of [Sx,Vsoc]  [Sy,Vsoc]  [Sz,Vsoc]
  complex(kind=c8b), allocatable :: torque(:,:,:,:,:)
! --> spherical part of valence charge and magnetisation densities
  real(kind=r8b),    allocatable :: nrv(:,:,:)
! --> all non-zero components of groundstate charge and magnetization densities
  real(kind=r8b),    allocatable :: old_rho2ns(:,:,:,:), gs_qlm(:,:), gs_mlm(:,:), new_rho2ns(:,:,:,:)
! --> charge and magnetization densities (not projected)
  real(kind=r8b),    allocatable :: rho_lm(:,:,:,:)
! --> density matrix integrated over WS sphere
  complex(kind=c8b), allocatable :: rhomat(:,:,:)
! --> xc kernel
  real(kind=r8b),    allocatable :: kxclm(:,:,:,:,:)
! --> total magnetization
  complex(kind=c8b), allocatable :: mtotsusc(:,:,:,:,:)
! --> magnetization from xc potential
  complex(kind=c8b), allocatable :: mxcsusc(:,:,:,:,:)
! --> magnetization from SOC+Bext potential
  complex(kind=c8b), allocatable :: msocsusc(:,:,:,:,:)
! --> number of non-vanishing multipoles
  integer(kind=i4b), allocatable :: nlmpot(:)
! --> pointers from lmpot to storage and back
  integer(kind=i4b), allocatable :: lmpot2i(:,:,:), i2lmpot(:,:)
! --> to compute SCF-like quantities when using the fitted GF
  integer(kind=i4b) :: nescf
  complex(kind=c8b), allocatable :: escf(:), ekscf(:), descf(:)
! --> coefficients of the projected scattering solutions of the Lippmann-Schwinger equation
  complex(kind=c8b), allocatable :: pzl(:,:,:,:), pzr(:,:,:,:), fzl(:,:,:,:), fzr(:,:,:,:)
! --> storage for single-site Green functions
  complex(kind=c8b), allocatable :: gfpq(:,:,:,:), gfps(:,:,:,:), gffq(:,:,:,:), gffs(:,:,:,:)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!         Variables and storage for the susceptibility module      
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This checks if the KS susc element is to be set to zero
  real(kind=r8b),    parameter   :: susctol = 1.d-7
! --> use onsite GF in susc
  logical                        :: lsusconsite = .false.
! --> use structural GF in susc
  logical                        :: lsuscstruct = .false.
! --> storage for kronecker form of KS susc
  complex(kind=c8b), allocatable :: kssusc(:,:)
  complex(kind=c8b), allocatable :: kssusc0(:,:)
  complex(kind=c8b), allocatable :: kssusc1(:,:)
  complex(kind=c8b), allocatable :: kssusc2(:,:)
! --> storage for spin-current spin correlation function
  complex(kind=c8b), allocatable :: kscurrcorr(:,:)
  complex(kind=c8b), allocatable :: kscurrcorr0(:,:)
  complex(kind=c8b), allocatable :: kscurrcorr1(:,:)
  complex(kind=c8b), allocatable :: kscurrcorr2(:,:)
! --> storage for kronecker form of xc kernel
  complex(kind=c8b), allocatable :: kernel(:,:)
! --> storage for kronecker form of susc denominator
  complex(kind=c8b), allocatable :: denominator(:,:)
! Gradient of scalar relativistic mass (added by Sascha)
  complex(kind=c8b), allocatable :: grad_mass(:,:,:,:,:)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!         Variables and storage for the input/output module      
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Numbering of atoms:
!   # atoms = # host atoms + # atoms in cluster
!   natypd  = natref       + natomd
! --> which atoms for projection, as listed in the impurity program
  integer(kind=i4b), allocatable :: iasusc(:)
! --> which atoms for susceptibility, as listed in inpsusc.dat
  integer(kind=i4b), allocatable :: iasusc2(:)
! --> useful arrays to read and write susc data
!     ngroup <= nasusc says how many groups of similar atoms there are
!     igroup(1:ngroup) says how many atoms there are in each group
  integer(kind=i4b) :: ngroup
  integer(kind=i4b), allocatable :: igroup(:)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! All done!
  end module global
