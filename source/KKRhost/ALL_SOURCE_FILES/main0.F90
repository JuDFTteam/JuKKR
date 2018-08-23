!-------------------------------------------------------------------------------
! MODULE: MOD_MAIN0
!> @brief Wrapper module for the reading and setup of the JM-KKR program
!> @author Philipp Rüssmann, Bernd Zimmermann, Phivos Mavropoulos, R. Zeller,
!> and many others ...
!> @todo JC: NATOMIMP and NATOMIMPD seem to be the same variable, however, right
!> now find no way to eliminate one of them.
!> @todo JC: There seem to be several repeated variables doing the same, e.g. INS,
!> KNOSPH, KWS and KSHAPE, all seem to dictate whether one has ASA or FP.
!> Maybe it would be good to consolidate and eliminate any unnecessary variables.
!> @todo JC: Several variables such as IRMD and IRNSD are actually determined in
!> the startb1 subroutine, maybe change the allocations such that they are done
!> there instead
!> @note
!> - Jonathan Chico Jan. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
#ifdef CPP_HYBRID
#define CPP_OMPSTUFF
#endif
#ifdef CPP_OMP
#define CPP_OMPSTUFF
#endif

module mod_main0

   use mod_wunfiles
   use mod_types, only: t_imp
#ifdef CPP_TIMING
   use mod_timing
#endif
   use Constants
   use memoryhandling
   use global_variables
   use mod_create_newmesh
   use mod_rhoqtools, only: rhoq_save_rmesh
   use rinput
   Use mod_datatypes, Only: dp

   use mod_addviratoms14
   use mod_bzkint0
   use mod_calcrotmat
   use mod_changerep
   use mod_cinit
   use mod_clsgen_tb
   use mod_convol
   use mod_deciopt
   use mod_drvbastrans
   use mod_epathtb
   use mod_gaunt2
   use mod_gaunt
   use mod_generalpot
   use mod_getbr3
   use mod_gfmask
   use mod_lattix99
   use mod_madelung2d
   use mod_madelung3d
   use mod_opt
   use mod_outpothost
   use mod_outtmathost
   use mod_readimppot
   use mod_relpotcvt
   use mod_rinit
   use mod_scalevec
   use mod_setgijtab
   use mod_shape_corr
   use mod_startb1
   use mod_startldau
   use mod_testdim
   use mod_write_tbkkr_files
   use mod_writehoststructure

   implicit none

   integer :: KTE       !< Calculation of the total energy On/Off (1/0)
   integer :: KWS       !< 0 (MT), 1(ASA)
   integer :: KXC       !< Type of xc-potential 0=vBH 1=MJW 2=VWN 3=PW91
   integer :: IGF       !< Do not print or print (0/1) the KKRFLEX_* files
   integer :: ICC       !< Enables the calculation of off-diagonal elements of the GF.(0=SCF/DOS; 1=cluster; -1=custom)
   integer :: INS       !< 0 (MT), 1(ASA), 2(Full Potential)
   integer :: IRM       !< Maximum number of radial points
   integer :: IPE       !< Not real used, IPFE should be 0
   integer :: IPF       !< Not real used, IPFE should be 0
   integer :: IPFE      !< Not real used, IPFE should be 0
   integer :: KCOR
   integer :: KEFG
   integer :: KHYP
   integer :: KPRE
   integer :: nprinc
   integer :: NSRA
   integer :: LPOT      !< Maximum l component in potential expansion
   integer :: IMIX      !< Type of mixing scheme used (0=straight, 4=Broyden 2nd, 5=Anderson)
   integer :: IEND      !< Number of nonzero gaunt coefficients
   integer :: ICST      !< Number of Born approximation
   integer :: NAEZ      !< Number of atoms in unit cell
   integer :: NEMB      !< Number of 'embedding' positions
   integer :: LMAX      !< Maximum l component in wave function expansion
   integer :: NCLS      !< Number of reference clusters
   integer :: NREF      !< Number of diff. ref. potentials
   integer :: NPOL      !< Number of Matsubara Poles (EMESHT)
   integer :: NPNT1     !< number of E points (EMESHT) for the contour integration
   integer :: NPNT2     !< number of E points (EMESHT) for the contour integration
   integer :: NPNT3     !< number of E points (EMESHT) for the contour integration
   integer :: LMMAX     !< (LMAX+1)^2
   integer :: NVIRT
   integer :: LMPOT     !< (LPOT+1)**2
   integer :: KVMAD
   integer :: ITSCF
   integer :: NCHEB     !< Number of Chebychev pannels for the new solver
   integer :: NINEQ     !< Number of ineq. positions in unit cell
   integer :: NATYP     !< Number of kinds of atoms in unit cell
   integer :: IFILE     !< Unit specifier for potential card
   integer :: KVREL     !< 0,1 : non / scalar relat. calculation
   integer :: NSPIN     !< Counter for spin directions
   integer :: NLEFT     !< Number of repeated basis for left host to get converged electrostatic potentials
   integer :: NRIGHT    !< Number of repeated basis for right host to get converged electrostatic potentials
   integer :: INVMOD    !< Inversion scheme
   integer :: KHFELD    !< 0,1: no / yes external magnetic field
   integer :: ITDBRY    !< Number of SCF steps to remember for the Broyden mixing
   integer :: INSREF    !< INS for reference pot. (usual 0)
   integer :: KSHAPE    !< Exact treatment of WS cell
   integer :: IELAST
   integer :: ISHIFT
   integer :: KFROZN
   integer :: NSYMAT
   integer :: NQCALC
   integer :: KFORCE    !< Calculation of the forces
   integer :: N1SEMI    !< Number of energy points for the semicore contour
   integer :: N2SEMI    !< Number of energy points for the semicore contour
   integer :: N3SEMI    !< Number of energy points for the semicore contour
   integer :: NLAYER    !< Number of principal layer
   integer :: NLBASIS   !< Number of basis layers of left host (repeated units)
   integer :: NRBASIS   !< Number of basis layers of right host (repeated units)
   integer :: INTERVX   !< Number of intervals in x-direction for k-net in IB of the BZ
   integer :: INTERVY   !< Number of intervals in y-direction for k-net in IB of the BZ
   integer :: INTERVZ   !< Number of intervals in z-direction for k-net in IB of the BZ
   integer :: MAXMESH
   integer :: NPAN_EQ   !< Number of intervals from [R_LOG] to muffin-tin radius Used in conjunction with runopt NEWSOSOL
   integer :: NPAN_LOG  !< Number of intervals from nucleus to [R_LOG] Used in conjunction with runopt NEWSOSOL
   integer :: NPOLSEMI  !< Number of poles for the semicore contour
   integer :: SCFSTEPS  !< number of scf iterations
   integer :: NATOMIMP  !< Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
   integer :: IESEMICORE
   integer :: IDOSEMICORE
   real (kind=dp) :: TK        !< Temperature
   real (kind=dp) :: FCM       !< Factor for increased linear mixing of magnetic part of potential compared to non-magnetic part.
   real (kind=dp) :: E2IN
   real (kind=dp) :: EMIN      !< Lower value (in Ryd) for the energy contour
   real (kind=dp) :: EMAX      !< Maximum value (in Ryd) for the DOS calculation Controls also [NPT2] in some cases
   real (kind=dp) :: ALAT      !< Lattice constant in a.u.
   real (kind=dp) :: RMAX      !< Ewald summation cutoff parameter for real space summation
   real (kind=dp) :: GMAX      !< Ewald summation cutoff parameter for reciprocal space summation
   real (kind=dp) :: R_LOG     !< Radius up to which log-rule is used for interval width. Used in conjunction with runopt NEWSOSOL
   real (kind=dp) :: RCUTZ     !< Parameter for the screening cluster along the z-direction
   real (kind=dp) :: RCUTXY    !< Parameter for the screening cluster along the x-y plane
   real (kind=dp) :: QBOUND    !< Convergence parameter for the potential
   real (kind=dp) :: VCONST    !< Potential shift in the first iteration
   real (kind=dp) :: HFIELD    !< External magnetic field, for initial potential shift in spin polarised case
   real (kind=dp) :: MIXING    !< Magnitude of the mixing parameter
   real (kind=dp) :: ABASIS    !< Scaling factors for rbasis
   real (kind=dp) :: BBASIS    !< Scaling factors for rbasis
   real (kind=dp) :: CBASIS    !< Scaling factors for rbasis
   real (kind=dp) :: EFERMI    !< Fermi energy
   real (kind=dp) :: ESHIFT
   real (kind=dp) :: TKSEMI    !< Temperature of semi-core contour
   real (kind=dp) :: TOLRDIF   !< For distance between scattering-centers smaller than [<TOLRDIF>], free GF is set to zero. Units are Bohr radii.
   real (kind=dp) :: ALATNEW
   real (kind=dp) :: VOLUME0
   real (kind=dp) :: EMUSEMI   !< Top of semicore contour in Ryd.
   real (kind=dp) :: EBOTSEMI  !< Bottom of semicore contour in Ryd
   real (kind=dp) :: FSEMICORE !< Initial normalization factor for semicore states (approx. 1.)
   real (kind=dp) :: LAMBDA_XC !< Scale magnetic moment (0 < Lambda_XC < 1, 0=zero moment, 1= full moment)
   character(len=10) :: SOLVER   !< Type of solver
   character(len=40) :: I12      !< File identifiers
   character(len=40) :: I13      !< Potential file name
   character(len=40) :: I19      !< Shape function file name
   character(len=40) :: I25      !< Scoef file name
   character(len=40) :: I40      !< File identifiers
   logical :: LRHOSYM
   logical :: LINIPOL    !< True: Initial spin polarization; false: no initial spin polarization
   logical :: LCARTESIAN !< True: Basis in cartesian coords; false: in internal coords

   !..
   !.. Local Arrays ..
   integer, dimension(NSYMAXD) :: ISYMINDEX
   integer, dimension(:), allocatable :: CLS    !< Cluster around atomic sites
   integer, dimension(:), allocatable :: IRC    !< R point for potential cutting
   integer, dimension(:), allocatable :: IMT    !< R point at MT radius
   integer, dimension(:), allocatable :: NFU    !< number of shape function components in cell 'icell'
   integer, dimension(:), allocatable :: NSH1   !< Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
   integer, dimension(:), allocatable :: NSH2   !< Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
   integer, dimension(:), allocatable :: LMXC
   integer, dimension(:), allocatable :: IPAN   !< Number of panels in non-MT-region
   integer, dimension(:), allocatable :: IRNS   !< Position of atoms in the unit cell in units of bravais vectors
   integer, dimension(:), allocatable :: IRWS   !< R point at WS radius
   integer, dimension(:), allocatable :: KMESH
   integer, dimension(:), allocatable :: IRMIN  !< Max R for spherical treatment
   integer, dimension(:), allocatable :: LOFLM  !< l of lm=(l,m) (GAUNT)
   integer, dimension(:), allocatable :: NACLS  !< Number of atoms in cluster
   integer, dimension(:), allocatable :: NCORE  !< Number of core states
   integer, dimension(:), allocatable :: IMAXSH
   integer, dimension(:), allocatable :: NSHELL !< Index of atoms/pairs per shell (ij-pairs); nshell(0) = number of shells
   integer, dimension(:), allocatable :: INIPOL !< Initial spin polarisation
   integer, dimension(:), allocatable :: IXIPOL !< Constraint of spin pol.
   integer, dimension(:), allocatable :: REFPOT !< Ref. pot. card  at position
   integer, dimension(:), allocatable :: NTCELL !< Index for WS cell
   integer, dimension(:), allocatable :: IQCALC
   integer, dimension(:), allocatable :: IOFGIJ !< Linear pointers, similar to NSH1/NSH2 but giving the actual index of sites I,J = 1,NATOMIMP in the cluster
   integer, dimension(:), allocatable :: JOFGIJ !< Linear pointers, similar to NSH1/NSH2 but giving the actual index of sites I,J = 1,NATOMIMP in the cluster
   integer, dimension(:), allocatable :: ATOMIMP
   integer, dimension(:), allocatable :: IJTABSH   !< Linear pointer, assigns pair (i,j) to a shell in the array GS(*,*,*,NSHELD)
   integer, dimension(:), allocatable :: IJTABSYM  !< Linear pointer, assigns pair (i,j) to the rotation bringing GS into Gij
   integer, dimension(:), allocatable :: NPAN_TOT
   integer, dimension(:), allocatable :: IJTABCALC !< Linear pointer, specifying whether the block (i,j) has to be calculated needs set up for ICC=-1, not used for ICC=1
   integer, dimension(:), allocatable :: NPAN_EQ_AT
   integer, dimension(:), allocatable :: NPAN_LOG_AT
   integer, dimension(:), allocatable :: IJTABCALC_I
   integer, dimension(:,:), allocatable :: ISH
   integer, dimension(:,:), allocatable :: JSH
   integer, dimension(:,:), allocatable :: ILM_MAP
   integer, dimension(:,:), allocatable :: KFG
   integer, dimension(:,:), allocatable :: ATOM    !< Atom at site in cluster
   integer, dimension(:,:), allocatable :: EZOA    !< EZ of atom at site in cluster
   integer, dimension(:,:), allocatable :: LMSP    !< 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
   integer, dimension(:,:), allocatable :: LCORE   !< Angular momentum of core states
   integer, dimension(:,:), allocatable :: ICLEB   !< Pointer array
   integer, dimension(:,:), allocatable :: IRCUT   !< R points of panel borders
   integer, dimension(:,:), allocatable :: LLMSP   !< lm=(l,m) of 'nfund'th nonvanishing component of non-spherical pot.
   integer, dimension(:,:), allocatable :: LMSP1
   integer, dimension(:,:), allocatable :: KAOEZ   !< Kind of atom at site in elem. cell
   integer, dimension(:,:), allocatable :: IFUNM
   integer, dimension(:,:), allocatable :: IFUNM1
   integer, dimension(:,:), allocatable :: ITITLE
   integer, dimension(:,:), allocatable :: ICHECK
   integer, dimension(:,:), allocatable :: IPAN_INTERVALL
   integer, dimension(:,:,:), allocatable :: JEND !< Pointer array for icleb()
   real (kind=dp), dimension(2) :: VBC       !< Potential constants
   real (kind=dp), dimension(3) :: ZPERIGHT  !< Vector to define how to repeat the basis of the right host
   real (kind=dp), dimension(3) :: ZPERLEFT  !< Vector to define how to repeat the basis of the left host
   real (kind=dp), dimension(3,3) :: RECBV   !< Reciprocal basis vectors
   real (kind=dp), dimension(3,3) :: BRAVAIS !< Bravais lattice vectors
   real (kind=dp), dimension(64,3,3) :: RSYMAT
   real (kind=dp), dimension(:), allocatable :: A   !< Constants for exponential R mesh
   real (kind=dp), dimension(:), allocatable :: B   !< Constants for exponential R mesh
   real (kind=dp), dimension(:), allocatable :: WG  !< Integr. weights for Legendre polynomials
   real (kind=dp), dimension(:), allocatable :: GSH
   real (kind=dp), dimension(:), allocatable :: ZAT !< Nuclear charge
   real (kind=dp), dimension(:), allocatable :: RMT !< Muffin-tin radius of true system
   real (kind=dp), dimension(:), allocatable :: RWS !< Wigner Seitz radius
   real (kind=dp), dimension(:), allocatable :: VREF
   real (kind=dp), dimension(:), allocatable :: MTFAC  !< Scaling factor for radius MT
   real (kind=dp), dimension(:), allocatable :: RMTNEW !< Adapted muffin-tin radius
   real (kind=dp), dimension(:), allocatable :: RMTREF !< Muffin-tin radius of reference system
   real (kind=dp), dimension(:), allocatable :: RMTREFAT
   real (kind=dp), dimension(:), allocatable :: FPRADIUS !< R point at which full-potential treatment starts
   real (kind=dp), dimension(:), allocatable :: SOCSCALE !< Spin-orbit scaling
   real (kind=dp), dimension(:,:), allocatable :: RMESH    !< Radial mesh ( in units a Bohr)
   real (kind=dp), dimension(:,:), allocatable :: S
   real (kind=dp), dimension(:,:), allocatable :: RR   !< Set of real space vectors (in a.u.)
   real (kind=dp), dimension(:,:), allocatable :: DRDI !< Derivative dr/di
   real (kind=dp), dimension(:,:), allocatable :: DROR
   real (kind=dp), dimension(:,:), allocatable :: CLEB !< GAUNT coefficients (GAUNT)
   real (kind=dp), dimension(:,:), allocatable :: VISP !< Spherical part of the potential
   real (kind=dp), dimension(:,:), allocatable :: CSCL !< Speed of light scaling
   real (kind=dp), dimension(:,:), allocatable :: RNEW
   real (kind=dp), dimension(:,:), allocatable :: RATOM
   real (kind=dp), dimension(:,:), allocatable :: ECORE   !< Core energies
   real (kind=dp), dimension(:,:), allocatable :: TLEFT   !< Vectors of the basis for the left host
   real (kind=dp), dimension(:,:), allocatable :: TRIGHT  !< Vectors of the basis for the right host
   real (kind=dp), dimension(:,:), allocatable :: SOCSCL
   real (kind=dp), dimension(:,:), allocatable :: RBASIS  !< Position of atoms in the unit cell in units of bravais vectors
   real (kind=dp), dimension(:,:), allocatable :: RCLSIMP
   real (kind=dp), dimension(:,:), allocatable :: CMOMHOST !< Charge moments of each atom of the (left/right) host
   real (kind=dp), dimension(:,:), allocatable :: RPAN_INTERVALL
   real (kind=dp), dimension(:,:,:), allocatable :: RS
   real (kind=dp), dimension(:,:,:), allocatable :: YRG  !< Spherical harmonics (GAUNT2)
   real (kind=dp), dimension(:,:,:), allocatable :: VINS !< Non-spherical part of the potential
   real (kind=dp), dimension(:,:,:), allocatable :: RCLS !< Real space position of atom in cluster
   real (kind=dp), dimension(:,:,:), allocatable :: RROT
   complex (kind=dp), dimension(:), allocatable :: EZ
   complex (kind=dp), dimension(:), allocatable :: DEZ
   complex (kind=dp), dimension(:), allocatable :: WEZ
   complex (kind=dp), dimension(:,:,:), allocatable :: DSYMLL
   complex (kind=dp), dimension(:,:,:), allocatable :: DSYMLL1
   complex (kind=dp), dimension(:,:,:,:,:), allocatable :: LEFTTINVLL
   complex (kind=dp), dimension(:,:,:,:,:), allocatable :: RIGHTTINVLL
   character(len=124), dimension(6) :: TXC
   logical, dimension(2) :: VACFLAG
   !
   !-------------------------------------------------------------------------
   !     Magnetisation angles -- description see RINPUT13
   !-------------------------------------------------------------------------
   integer :: KMROT !< 0: no rotation of the magnetisation; 1: individual rotation of the magnetisation for every site
   real (kind=dp), dimension(:), allocatable :: QMTET !< \f$ \theta\f$ angle of the agnetization with respect to the z-axis
   real (kind=dp), dimension(:), allocatable :: QMPHI !< \f$ \phi\f$ angle of the agnetization with respect to the z-axis
   !-------------------------------------------------------------------------
   !     CPA variables
   !-------------------------------------------------------------------------
   integer :: NCPA     !< NCPA = 0/1 CPA flag
   integer :: ITCPAMAX !< Max. number of CPA iterations
   integer, dimension(:), allocatable :: NOQ  !< Number of diff. atom types located
   integer, dimension(:), allocatable :: IQAT !< The site on which an atom is located on a given site
   integer, dimension(:), allocatable :: ICPA !< ICPA = 0/1 site-dependent CPA flag
   !
   !-------------------------------------------------------------------------
   !> @note ITERMDIR running option introduced Apr 2003 -- Munich
   !>              (H. Ebert + V. Popescu) allows a self-consistent
   !>              determination of the magnetic configuration in REL mode
   !-------------------------------------------------------------------------
   real (kind=dp), dimension(:), allocatable :: QMGAM
   real (kind=dp), dimension(:,:), allocatable :: QMGAMTAB
   real (kind=dp), dimension(:,:), allocatable :: QMPHITAB
   real (kind=dp), dimension(:,:), allocatable :: QMTETTAB
   !-------------------------------------------------------------------------
   !> @note changes for impurity 20/02/2004 -- v.popescu according to
   !>                                          n.papanikolaou VINS()
   !-------------------------------------------------------------------------
   integer, dimension(:), allocatable :: HOSTIMP
   real (kind=dp) :: CPATOL !< Convergency tolerance for CPA-cycle
   real (kind=dp), dimension(:), allocatable :: CONC !< Concentration of a given atom
   !-------------------------------------------------------------------------------
   complex (kind=dp), dimension(:,:), allocatable :: RC !< NREL REAL spher. harm. > CMPLX. spher. harm. NREL CMPLX. spher. harm. > REAL spher. harm.
   complex (kind=dp), dimension(:,:), allocatable :: CREL !< Non-relat. CMPLX. spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. CMPLX. spher. harm.
   complex (kind=dp), dimension(:,:), allocatable :: RREL !< Non-relat. REAL spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. REAL spher. harm.
   complex (kind=dp), dimension(:,:,:), allocatable :: SRREL
   complex (kind=dp), dimension(:,:,:), allocatable :: DROTQ !< Rotation matrices to change between LOCAL/GLOBAL frame of reference for magnetisation <> Oz or noncollinearity
   integer, dimension(:), allocatable :: ZREL    !< atomic number (cast integer)
   integer, dimension(:), allocatable :: JWSREL  !< index of the WS radius
   integer, dimension(:), allocatable :: IRSHIFT !< shift of the REL radial mesh with respect no NREL
   integer, dimension(:,:), allocatable :: NRREL
   integer, dimension(:,:,:), allocatable :: IRREL
   real (kind=dp), dimension(0:100) :: FACT
   real (kind=dp), dimension(:,:), allocatable :: VTREL !< potential (spherical part)
   real (kind=dp), dimension(:,:), allocatable :: BTREL !< magnetic field
   real (kind=dp), dimension(:,:), allocatable :: RMREL !< radial mesh
   real (kind=dp), dimension(:,:), allocatable :: DRDIREL !< derivative of radial mesh
   real (kind=dp), dimension(:,:), allocatable :: R2DRDIREL !< \f$ r^2 \frac{\partial}{\partial \mathbf{r}}\frac{\partial}{\partial i}\f$ (r**2 * drdi)
   real (kind=dp), dimension(:,:,:), allocatable :: THESME
   logical :: PARA
   logical, dimension(NSYMAXD) :: SYMUNITARY !< unitary/antiunitary symmetry flag
   !
   !-------------------------------------------------------------------------
   ! LDA+U LDA+U LDA+U
   !-------------------------------------------------------------------------
   !> @note ph. mavropoulos according to Munich SPR-KKR
   !>       h. ebert
   !>      input:
   !>            UEFF, JEFF : input U,J parameters for each atom
   !>            EREFLDAU(1..NATYP) : the energies of ggthe projector's wave
   !>                                  functions (REAL)
   !>            LOPT(1..NATYP): angular momentum QNUM for the atoms on
   !>                             which LDA+U should be applied (-1 to
   !>                             switch it OFF)
   !>      iteration index ITRUNLDAU
   !>      integer flag perform LDA+U IDOLDAU
   !>      integer flag LDA+U arrays available KREADLDAU
   !>      NTLDAU - number of atoms on which LDA+U is applied (<=NATYP)
   !>      arrays: ULDAU - calculated Coulomb matrix elements (EREFLDAU)
   !>              WLDAU - potential matrix
   !>              ITLDAU - integer pointer connecting the NTLDAU atoms to
   !>                       their corresponding index in the unit cell
   !
   !-------------------------------------------------------------------------
   integer :: NTLDAU    !< number of atoms on which LDA+U is applied
   integer :: IDOLDAU   !< flag to perform LDA+U
   integer :: ITRUNLDAU !< Iteration index for LDA+U
   integer :: KREADLDAU !< LDA+U arrays available
   integer, dimension(:), allocatable :: LOPT !< angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
   integer, dimension(:), allocatable :: ITLDAU !< integer pointer connecting the NTLDAU atoms to heir corresponding index in the unit cell
   real (kind=dp), dimension(:), allocatable :: UEFF !< input U parameter for each atom
   real (kind=dp), dimension(:), allocatable :: JEFF !< input J parameter for each atom
   real (kind=dp), dimension(:), allocatable :: EREFLDAU !< the energies of the projector's wave functions (REAL)
   !..
   !.. distinguish between spin-dependent and spin-independent
   !.. quantities
   real (kind=dp), dimension(:,:,:,:), allocatable :: WLDAU !< potential matrix
   real (kind=dp), dimension(:,:,:,:,:), allocatable :: ULDAU !< calculated Coulomb matrix elements (EREFLDAU)
   complex (kind=dp), dimension(:,:), allocatable :: PHILDAU
   !-------------------------------------------------------------------------
   ! LDA+U LDA+U LDA+U
   !-------------------------------------------------------------------------
   ! Lloyds formula
   integer :: LLY !< LLY <> 0 : apply Lloyds formula
   complex (kind=dp) :: DELTAE  !< Energy difference for numerical derivative

   ! SUSC (BEGIN: modifications by Manuel and Benedikt)             ! susc
   ! LOGICAL THAT CHECKS WHETHER ENERGY MESH FILE EXISTS            ! susc
   logical        :: emeshfile                                      ! susc
   ! SUSC (END:   modifications by Manuel and Benedikt)             ! susc

   ! ruess: IVSHIFT test option
   integer :: IVSHIFT

   !allocations:
   real (kind=dp), dimension(:,:,:), allocatable :: THETAS !< shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
   real (kind=dp), dimension(:,:,:), allocatable :: THETASNEW

   public :: main0, bshift_ns

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: main0
   !> @brief Main wrapper to handle input reading, allocation of arrays, and
   !> preparation of all the necessary data structures for a calculation.
   !> @author Philipp Rüssmann, Bernd Zimmermann, Phivos Mavropoulos, R. Zeller,
   !> and many others ...
   !> @note Jonathan Chico: Re-wrote the module from Fixed format Fortran to
   !> fortran90 Free Form. Also performed modifications to get rid of the inc.p
   !> file. 22.12.2017
   !----------------------------------------------------------------------------
   subroutine main0()

      use mod_types
#ifdef CPP_OMPSTUFF
      use omp_lib        ! necessary for omp functions
#endif
#ifdef CPP_MPI
      use mpi
#endif
      use mod_mympi, only: nranks
      use mod_version
      use mod_version_info
      use mod_md5sums

      implicit none
      !.. Local Scalars ..
      integer :: I
      integer :: J
      integer :: I1
      integer :: IE
      integer :: LM
      integer :: NS
      integer :: ISVATOM,NVATOM
      integer :: i_stat
      integer :: IREC
      integer :: LRECABMAD
      real (kind=dp) :: ZATTEMP
      integer :: ierr
      real (kind=dp), allocatable :: tmp_rr(:,:)
      ! for OPERATOR option
      logical :: lexist, operator_imp

#ifdef CPP_OMPSTUFF
      !     .. OMP ..
      integer nth, ith          ! total number of threads, thread number
#endif
      !     ..
      !     .. External Functions ..
      logical :: OPT,TEST
      external :: OPT,TEST
      !     ..
      !     .. External Subroutines ..
      external :: BZKINT0,CINIT,CLSGEN_TB,DECIOPT,EPATHTB,GAUNT,GAUNT2
      external :: GFMASK,LATTIX99,RINIT,SCALEVEC
      external :: STARTB1,STARTLDAU,TESTDIM,SHAPE_CORR
      external :: readimppot
      !     ..
      !     .. Intrinsic Functions ..
      intrinsic :: ATAN,DABS,DBLE,LOG,MAX,SQRT,product,shape
      !     ..
      !     ..
      !-------------------------------------------------------------------------
      ! Write version info:
      !-------------------------------------------------------------------------
      call print_versionserial (6, version1, version2, version3, version4, serialnr)
      call print_versionserial (1337, version1, version2, version3, version4, serialnr)

#ifdef CPP_OMPSTUFF
!$omp parallel shared(nth) private(ith)
      nth = omp_get_num_threads()
      ith = omp_get_thread_num()
      if(ith==0) then
         write(*,'(1X,A,I5//79("*")/)') 'Number of OpenMP threads used:',nth
         write(1337,'(1X,A,I5//79("*")/)') 'Number of OpenMP threads used:',nth
      endif
!$omp end parallel
#endif

#ifdef CPP_MPI
      write(*,'(1X,A,I5//79("*")/)') 'Number of MPI ranks used:',nranks
      write(1337,'(1X,A,I5//79("*")/)') 'Number of MPI ranks used:',nranks
#endif
      !-------------------------------------------------------------------------
      ! End write version info
      !-------------------------------------------------------------------------
      !

      ! allocate and initialize default values
      call init_all_wrapper()

      !
      !-------------------------------------------------------------------------
      ! Reading of the inputcard, and allocation of several arrays
      !> @note JC: have added reading calls for the parameters that used to be in
      !> the inc.p and can now be modified via the inputcard directly
      !-------------------------------------------------------------------------
      call RINPUT13(KTE,IGF,KXC,LLY,ICC,INS,KWS,IPE,IPF,IPFE,ICST,IMIX, &
         LPOT,NAEZ,NEMB,NREF,NCLS,NPOL,LMAX,KCOR,KEFG,KHYP,KPRE,  &
         KVMAD,LMMAX,LMPOT,NCHEB,NLEFT,IFILE,KVREL,NSPIN,NATYP,NINEQ,NPNT1,NPNT2,   &
         NPNT3,KFROZN,ISHIFT,      &
         N1SEMI,N2SEMI,N3SEMI,SCFSTEPS,INSREF,KSHAPE,ITDBRY,NRIGHT,KFORCE,  &
         IVSHIFT,KHFELD,NLBASIS,NRBASIS,INTERVX,INTERVY,INTERVZ,NPAN_EQ,    &
         NPAN_LOG,NPOLSEMI,TK,FCM,EMIN,EMAX,RMAX,GMAX,ALAT,R_LOG,RCUTZ,RCUTXY,      &
         ESHIFT,QBOUND,HFIELD,MIXING,ABASIS,BBASIS,CBASIS,VCONST,TKSEMI,TOLRDIF,    &
         EMUSEMI,EBOTSEMI,FSEMICORE,LAMBDA_XC,DELTAE,LRHOSYM,LINIPOL,LCARTESIAN,    &
         IMT,CLS,LMXC,IRNS,IRWS,NTCELL,REFPOT,INIPOL,IXIPOL,HOSTIMP,KFG, &
         VBC,ZPERLEFT,ZPERIGHT,BRAVAIS,RMT,ZAT,RWS,MTFAC,RMTREF,RMTNEW,RMTREFAT,    &
         FPRADIUS,TLEFT,TRIGHT,RBASIS,SOCSCALE,CSCL,SOCSCL,SOLVER,I12,I13,I19,I25,  &
         I40,TXC,DROTQ,NCPA,ITCPAMAX,CPATOL,NOQ,IQAT,ICPA,KAOEZ,CONC,KMROT,QMTET,   &
         QMPHI,KREADLDAU,LOPT,UEFF,JEFF,EREFLDAU)

      !Some consistency checks
      if ( (KREL.lt.0) .OR. (KREL.gt.1) ) &
         stop ' set KREL=0/1 (non/fully) relativistic mode in the inputcard'
      if ( (KREL.eq.1) .and. (NSPIND.eq.2) )&
         stop ' set NSPIND = 1 for KREL = 1 in the inputcard'

      ! Set the calculation of several parameters
      NREFD = naez!+nemb ! can be changed later on when it is determined in clsgen_tb

      IRM = IRMD

      NTOTD=IPAND+30
      NCELLD=NAEZ
      NRMAXD=NTOTD*(NCHEB+1)
      NCHEBD = NCHEB

      NATYPD = NATYP
      NAEZD = NAEZ

      LMAXD        = LMAX
      ALM          = NAEZD*LMMAXD
      ALMGF0       = NAEZD*LMGF0D
      NDIM_SLABINV = NPRINCD*LMMAXD
      NEMBD        = NEMBD1-1
      NEMBD2       = NEMBD+NAEZ
      IRMIND       = IRMD-IRNSD
      NOFGIJ       = NATOMIMPD*NATOMIMPD+1
      LPOTD = LPOT
      LMPOTD       = (LPOT+1)**2
      LMMAXSO      = lmmaxd ! lmmaxd already doubled in size! (KREL+1)*LMMAXD

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Allocation calls
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !> @note Jonathan Chico: The main idea here is to allocate all the needed arrays so that the inc.p
      !> file becomes irrelevant. In principle the philosophy would be to modularize
      !> the code such that each module has its own global variables and allocation routine
      !> e.g. a module called CPA_control could have defined all the needed CPA variables
      !> as well as the allocation calls, this module would be used in the needed routines
      !> and the arrays would only be allocated if a CPA calculation is actually performed
      !> in the current way ALL arrays are allocated which could cause an unnecessary memory
      !> consumption
      !
      ! Call to allocate the arrays associated with the potential
      call allocate_potential(1,IRMD,NATYPD,NPOTD,IPAND,NFUND,LMXSPD,    &
         LMPOTD,IRMIND,NSPOTD,NFU,IRC,NCORE,IRMIN,LMSP,LMSP1,IRCUT,LCORE,LLMSP,   &
         ITITLE,FPRADIUS,VISP,ECORE,VINS)
      ! Call to allocate the arrays associated with the LDA+U potential
      call allocate_ldau_potential(1,IRMD,NATYPD,MMAXD,NSPIND,ITLDAU,WLDAU,ULDAU,  &
         PHILDAU)
      ! Call to allocate the arrays associated with the energy
      call allocate_energies(1,IEMXD,EZ,DEZ,WEZ)
      ! Call to allocate the arrays associated with the relativistic corrections
      call allocate_relativistic(1,KREL,IRMD,NAEZD,NATYPD,ZREL,JWSREL,IRSHIFT,      &
         VTREL,BTREL,RMREL,DRDIREL,R2DRDIREL,QMGAM,QMGAMTAB,QMPHITAB,QMTETTAB)
      ! Call to allocate the arrays associated with the relativistic transformations
      call allocate_rel_transformations(1,LMMAXD,NRREL,IRREL,RC,CREL,RREL,SRREL)
      ! Call to allocate the arrays associated with the clusters
      call allocate_clusters(1,NAEZD,LMAXD,NCLEB,NCLSD,NEMBD1,NSHELD,NACLSD,LMPOTD, &
         NATOMIMPD,NSH1,NSH2,NACLS,NSHELL,ATOMIMP,ATOM,EZOA,ICLEB,JEND,RATOM,    &
         RCLSIMP,CMOMHOST,RCLS)
      ! Call to allocate the arrays associated with the expansion of the Green function
      call allocate_expansion(1,LM2D,IRID,NFUND,NTOTD,NCLEB,LASSLD,NCELLD,NCHEBD, &
         LOFLM,WG,CLEB,YRG,THETAS,THETASNEW)
      ! Call to allocate the arrays associated with the integration mesh
      call allocate_mesh(1,IRMD,NATYPD,A,B,RMESH,DRDI)
      ! Call to allocate the arrays associated with the pannels for the new solver
      call allocate_pannels(1,NATYPD,NTOTD,IPAN,NPAN_TOT,NPAN_EQ_AT,NPAN_LOG_AT,  &
         IPAN_INTERVALL,RPAN_INTERVALL)
      ! Call to allocate misc arrays
      call allocate_misc(1,NRD,IRMD,IRID,LMAXD,NAEZD,NATYPD,NFUND,NREFD,IEMXD,NTOTD,   &
         NSHELD,LMMAXD,NEMBD1,NCHEBD,NCELLD,LMXSPD,NSPINDD,NSYMAXD,NPRINCD,IFUNM, &
         IFUNM1,ICHECK,VREF,S,RR,DROR,RNEW,RS,RROT,THESME,DSYMLL,DSYMLL1,        &
         LEFTTINVLL,RIGHTTINVLL)
      ! Call to allocate the arrays associated with the Green function
      call allocate_green(1,NAEZD,IEMXD,NGSHD,NSHELD,LMPOTD,NOFGIJ,ISH,JSH,KMESH,  &
         IMAXSH,IQCALC,IOFGIJ,JOFGIJ,IJTABSH,IJTABSYM,IJTABCALC,IJTABCALC_I,ILM_MAP, &
         GSH)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End of allocation calls
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Deal with the lattice
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call LATTIX99(LINTERFACE,ALAT,NATYP,NAEZ,CONC,RWS,BRAVAIS,RECBV,VOLUME0,RR,&
         NRD,NATYP)

      ! rr has changed, fix allocation of array to new nrd size
      allocate(tmp_rr(3,0:nrd), stat=ierr)
      if ( ierr/=0 ) stop 'error allocating tmp_rr in main0'
      tmp_rr(:,:) = rr(1:3, 0:nrd)
      deallocate(rr, stat=ierr)
      if ( ierr/=0 ) stop 'error deallocating rr in main0'
      allocate(rr(3,0:nrd), stat=ierr)
      if ( ierr/=0 ) stop 'error reallocating rr in main0'
      rr(:,:) = tmp_rr(:,:)
      deallocate(tmp_rr, stat=ierr)
      if ( ierr/=0 ) stop 'error allocating tmp_rr in main0'

      call SCALEVEC(LCARTESIAN,RBASIS,ABASIS,BBASIS,CBASIS,NLBASIS,NRBASIS,NLEFT,&
         NRIGHT,ZPERLEFT,ZPERIGHT,TLEFT,TRIGHT,LINTERFACE,NAEZ,NEMB,BRAVAIS,     &
         KAOEZ,NOQ,NAEZ,NATYP,NEMB)
      ! After SCALEVEC all basis positions are in cartesian coords.

      NVIRT = 0
      if ( OPT('VIRATOMS') ) then
         write(1337,*) 'Calling ADDVIRATOMS'
         call ADDVIRATOMS14(LINTERFACE,NVIRT,NAEZ,NAEZ,NATYP,NEMB,NEMB,RBASIS,   &
         .TRUE.,BRAVAIS,NCLS,NINEQ,REFPOT,KAOEZ,NOQ,NREF,RMTREFAT,I25)
      endif

      call clsgen_tb(naez, nemb, nvirt, rr, rbasis, kaoez, zat, cls, ncls, &
         nacls, atom, ezoa, nlbasis, nrbasis, nleft, nright, zperleft, zperight, &
         tleft, tright, rmtref, rmtrefat, vref, refpot, nref, rcls, rcutz, rcutxy, &
         alat, natyp, nclsd, nrd, naclsd, nrefd, nembd, linterface, nprinc)

      write(*,*) 'after clsgen', nref, nrefd

      ! overwrite nprincd if chosen too small
      if ( nprincd<nprinc ) then
        nlayer = naez/nprinc
        if ( nlayer*nprinc /= naez ) nprinc = naez
        write(*,*) 'Automatically overwriting nprincd with ', nprinc
        write(1337,*) 'Automatically overwriting nprincd with ', nprinc
        nprincd = nprinc
        ! update parameter that depend on nprincd
        NDIM_SLABINV = NPRINCD*LMMAXD
        ! change allocations of arrays that have nprincd
        deallocate (icheck, stat=i_stat)
        call memocc(i_stat, product(shape(icheck))*kind(icheck), 'ICHECK', &
          'main0')
        allocate (icheck(naez/nprincd,naez/nprincd), stat=i_stat)
        call memocc(i_stat, product(shape(icheck))*kind(icheck), 'ICHECK', &
          'main0')
        icheck = 0
      end if


      ! Now the clusters, reference potentials and muffin-tin radii have been set.
      !-------------------------------------------------------------------------
      ! Consistency check
      !-------------------------------------------------------------------------
      if ( (KREL.EQ.1) .AND. (INS.NE.0) ) then
         write(6,*)' FULL-POTENTIAL RELATIVISTIC mode not implemented '
         stop ' set INS = 0 in the input'
      end if
      !
      if ( KVREL.LE.1 ) then
         if ( KREL.EQ.1 ) stop ' KVREL <= 1 in input, but relativistic program used'
      else
         if ( KREL.EQ.0 ) stop ' KVREL > 1 in input, but non-relativistic program used'
      end if
      !-------------------------------------------------------------------------
      !
      E2IN = EMAX
      NSRA = 1
      if (KVREL.GE.1) NSRA = 2
      if (KVREL.GE.2) NSRA = 3
      !
      call TESTDIM(NSPIN,NAEZ,NEMB,NATYP,INS,INSREF,NREF,IRNS,  &
         NLAYER,KREL,NSPIND,NPRINCD,KNOSPH,IRNSD,KORBIT)
      !
      if ( INS.GT.0 )    open (19,FILE=I19,STATUS='old',FORM='formatted')
      if ( IFILE.EQ.13 ) open (IFILE,FILE=I13,STATUS='old',FORM='formatted')
      if ( ICC.GT.0 )    open (25,FILE=I25,STATUS='unknown',FORM='formatted')
      !
      call STARTB1(IFILE,1337,1337,IPE,KREL,KWS,LMAX,1,NATYP,ALATNEW,RMTNEW, &
         RMT,ITITLE,IMT,IRC,VCONST,INS,IRNS,FPRADIUS,NSPIN,VINS,IRMIN,   &
         KSHAPE,NTCELL,IRCUT,IPAN,THETAS,IFUNM,NFU,LLMSP,LMSP,E2IN,VBC,       &
         DROR,RS,S,VISP,RWS,ECORE,LCORE,NCORE,DRDI,RMESH,ZAT,A,B,IRWS,1,LMPOT,    &
         IRMIND,IRMD,LMXSPD,IPAND,IRID,IRNSD,NATYP,NCELLD,NFUND,NSPOTD,IVSHIFT, &
         NPOTD)


      ! find md5sums for potential and shapefunction
      call get_md5sums(INS, I13, I19)
      write(1337,'(A,A)') 'Doing calculation with potential: ',md5sum_potential
      if(INS>0) then
         write(1337,'(A,A)') 'Doing calculation with shapefun: ',md5sum_shapefun
      end if

      IF(TEST('rhoqtest')) THEN
         call rhoq_save_rmesh(natyp,irm,ipand,irmin,irws,ipan,rmesh,ntcell,ircut,r_log,npan_log,npan_eq)

         ! check consistency
         if (OPT('GREENIMP') .or. OPT('WRTGREEN')) then
            write(*,*) 'warning! rhoqtest cannot be used together with '
            write(*,*) "'GREENIMP' or 'WRTGREEN' options"
            stop
         end if

         !! enforce MPIatom here
         !if (TEST('MPIenerg')) stop 'rhoqtest assumes MPIatom'
         !CALL ADDTEST('MPIatom')
         
         if (nranks>1) then
            write(*,*) 'at the moment rhoqtest does not work with MPI.'
            write(*,*) 'compile hybrid version and use OMP level only.'
            stop
         end if 

      ENDIF ! TEST('rhoqtest')

      if ( TEST('Vspher  ') ) then
         write(1337,*) 'TEST OPTION Vspher,',&
         'keeping only spherical component of potential.'
         VINS(IRMIND:IRMD,2:LMPOT,1:NSPOTD) = 0.D0
      endif

      if (OPT('zeropot ').OR.TEST('zeropot ')) then
         write(1337,*) 'Using OPT zeropot, setting potential to zero.'
         write(1337,*) 'Using OPT zeropot, setting nuclear charge to zero.'
         VINS(IRMIND:IRMD,1:LMPOT,1:NSPOTD) = 0.D0
         VISP(1:IRMD,1:NPOTD) = 0.D0
         ZAT(1:NATYP) = 0.D0
      endif
      !
      do I1 = 1,NATYP
         do LM = 1,LMXSPD
            IFUNM1(LM,I1) = IFUNM(I1,LM)
            LMSP1(LM,I1) = LMSP(I1,LM)
         end do
      end do
      !-------------------------------------------------------------------------
      ! update Fermi energy, adjust energy window according to running options
      !-------------------------------------------------------------------------
      if ( NPOL.EQ.0 ) EFERMI = E2IN
      if ( OPT('GF-EF   ') .OR. OPT('DOS-EF  ') ) then
         EMIN = E2IN
         if ( OPT('GF-EF   ') ) then
            write (1337,FMT=9070)
         else
            write (1337,FMT=9080)
         end if
      end if

      if ( DABS(E2IN-EMAX).GT.1D-10 .AND. NPOL.NE.0 ) EMAX = E2IN
      !-------------------------------------------------------------------------
      if (OPT('GENPOT  ')) then
         rewind(3)
         call GENERALPOT(3,1,NATYP,NSPIN,ZAT,ALAT,RMT,RMTNEW,RWS,RMESH,DRDI,VISP, &
            IRWS,A,B,INS,IRNS,LPOT,VINS,QBOUND,IRC,KSHAPE,E2IN,VBC,ECORE,     &
            LCORE,NCORE,LMPOT,IRMD,IRMIND)
         close(3)
      end if
      !-------------------------------------------------------------------------
      ! --> Apply external magnetic field
      !-------------------------------------------------------------------------
      !     from startb1 moved here
      if (KHFELD.eq.1) then
         !---> maybe apply a magnetic field
         call bshift_ns(IRMD,IRID,IPAND,LMPOT,NPOTD,NATYP,NSPIN,NGSHD,NFUND,NCELLD,  &
            IRMIND,LMXSPD,KSHAPE,IRC,IRMIN,INIPOL,NTCELL,IMAXSH,ILM_MAP,LMSP,IFUNM,     &
            IRCUT,HFIELD,GSH,RMESH,THESME,THETAS,VISP,VINS)
      end if
      if ( TEST('vpotout ') ) then !ruess
         open(unit=54633163,file='test_vpotout_bshift')
         do i1=1,natyp*nspin
            write(54633163,*) '# visp of atom ',i1
            write(54633163,'(50000E14.7)') visp(:,i1)
         end do !iatom
         do i1=1,natyp*nspin
            write(54633163,*) '# vins of atom ',i1
            write(54633163,'(50000E14.7)') vins(:,:,i1)
         end do !iatom
         close(54633163)
      end if
      !-------------------------------------------------------------------------
      ! Deal with the potential in the RELATIVISTIC CASE
      !-------------------------------------------------------------------------
      PARA = .TRUE.
      if (KREL+KORBIT.eq.1) then
         !----------------------------------------------------------------------
         if (NSPIN.eq.1) then
            !-------------------------------------------------------------------
            ! for paramagnetic (NSPIN=1) input potential fill up also the
            ! V(DOWN), ECORE(DOWN), LCORE(DOWN), NCORE(DOWN) and ITITLE(DOWN)
            ! arrays (needed)
            !-------------------------------------------------------------------
            do I = NATYP,1,-1
               J = 2*I - 1
               call DCOPY(IRMD,VISP(1,I),1,VISP(1,J),1)
               call DCOPY(IRMD,VISP(1,J),1,VISP(1,J+1),1)
               !
               call DCOPY(20,ECORE(1,I),1,ECORE(1,J),1)
               call DCOPY(20,ECORE(1,J),1,ECORE(1,J+1),1)
               !
               NCORE(J) = NCORE(I)
               NCORE(J+1) = NCORE(J)
               !
               do I1=1,20
                  LCORE(I1,J) = LCORE(I1,I)
                  LCORE(I1,J+1) = LCORE(I1,J)
                  ITITLE(I1,J ) = ITITLE(I1,I)
                  ITITLE(I1,J+1) = ITITLE(I1,J)
               end do
            end do
            !-------------------------------------------------------------------
         else !NSPIN.eq.1
            !-------------------------------------------------------------------
            ! --> check whether, although NSPIN=2 at input, the system is
            !     paramagnetic (useful for symmetry cosiderations)
            !-------------------------------------------------------------------
            do I = 1,2*NATYP-1,2
               do J=1,IRMD
                  if (ABS(VISP(J,I)-VISP(J,I+1)).GT.1D-5) PARA = .FALSE.
               end do
            end do
            if (PARA) then
               do I=1,2*NATYP-1,2
                  call DCOPY(IRMD,VISP(1,I),1,VISP(1,I+1),1)
               end do
            end if

         end if !NSPIN.eq.1
         !----------------------------------------------------------------------
         ! finally, convert input potential to the internal relativistic
         ! form VTREL,BTREL. Set up auxiliary arrays (needed in the REL
         ! routines) ZREL, JWSREL, RMREL, DRDIREL, R2DRDIREL, IRSHIFT
         !----------------------------------------------------------------------
         if(KREL.eq.1) then
            ! call this only if relativisitic solver is used
            call RELPOTCVT(1,VISP,ZAT,RMESH,DRDI,IRCUT,VTREL,BTREL,ZREL,RMREL,&
               JWSREL,DRDIREL,R2DRDIREL,IRSHIFT,IPAND,IRMD,NPOTD,NATYP)
         endif
      end if !KREL+KORBIT.EQ.1
      !-------------------------------------------------------------------------
      ! set up energy contour
      !-------------------------------------------------------------------------
      IDOSEMICORE = 0
      if ( OPT('SEMICORE') ) IDOSEMICORE = 1
      !
      call EPATHTB(EZ,DEZ,E2IN,IELAST,IESEMICORE,IDOSEMICORE,EMIN,EMAX,TK,NPOL,  &
         NPNT1,NPNT2,NPNT3,EBOTSEMI,EMUSEMI,TKSEMI,NPOLSEMI,N1SEMI,N2SEMI,       &
         N3SEMI,IEMXD)

      ! SUSC (BEGIN: modifications by Manuel and Benedikt)             ! susc
      !                                                                ! susc
      if ( opt('EMESH   ') ) then                                      ! susc
      ! write out the energy mesh and the corresponding                ! susc
      ! weights to a file called 'emesh.scf'                           ! susc
        write(*,'("main0: Runflag emesh is set.")')                    ! susc
        write(*,'("       File emesh.scf will be written!")')          ! susc
        write(*,*) 'writing emesh.scf file...'                         ! susc
        open(file='emesh.scf',unit=12111984,status='replace')          ! susc
          write(12111984,'(5x,i0)') ielast                             ! susc
          do ie = 1,ielast                                             ! susc
             write(12111984,'(4es16.8)') ez(ie),dez(ie)                ! susc
          end do                                                       ! susc
        close(12111984)                                                ! susc
        write(*,'("       Finished writing file emesh.scf.")')         ! susc
      end if                                                           ! susc
      !                                                                ! susc
      !                                                                ! susc
      if ( opt('KKRSUSC ') ) then                                      ! susc
      ! read in 'emesh.dat' with new energy mesh-points                ! susc
        inquire(file='emesh.dat',exist=emeshfile)                      ! susc
        write(*,'("main0: Runflag KKRSUSC is set.")')                  ! susc
        if (emeshfile) then                                            ! susc
          write(*,'("main0: File emesh.dat exists and will ")', advance='no')       ! susc
          write(*,                            '("be read in.")')       ! susc
          write(*,'("       Energy contour from inputcard ")', advance='no')        ! susc
          write(*,                 '("will be overwritten!")')         ! susc
          open(file='emesh.dat',unit=50)                               ! susc
            read(50,*) ielast                                          ! susc
            if (ielast > iemxd) stop 'ielast > iemxd!'                 ! susc
            do ie=1,ielast                                             ! susc
              read(50,'(4es16.8)') ez(ie), dez(ie)                     ! susc
              write(*,'(i8,4es16.8)') ie, ez(ie), dez(ie)              ! susc
            end do                                                     ! susc
          close(50)                                                    ! susc
          write(*,'("       Finished reading in file emesh.dat.")')    ! susc
        else                                                           ! susc
          stop'main0: Runflag KKRSUSC but cannot find file emesh.dat!' ! susc
        end if                                                         ! susc
      end if                                                           ! susc
      !                                                                ! susc
      !                                                                ! susc
      ! still missing: check here whether scfsteps is > 1              ! susc
      !   if scfsteps>1 --> option a) stop program here                ! susc
      !                     option b) set it to 1 and continue         ! susc
      if ( opt('KKRSUSC ') .and. scfsteps>1 ) then                     ! susc
        write(*,'("main0: Runflag KKRSUSC is set ")')                  ! susc
        write(*,            '("but scfsteps = ",i0)') scfsteps         ! susc
        write(*,'("       Here we enforce scfsteps = 1")')             ! susc
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! susc
        scfsteps = 1                                                   ! susc
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! susc
      end if                                                           ! susc
      !                                                                ! susc
      ! SUSC (END:   modifications by Manuel and Benedikt)             ! susc

      do IE = 1,IELAST
         WEZ(IE) = -2.D0/PI*DEZ(IE)
         if ( IE.LE.IESEMICORE ) WEZ(IE) = WEZ(IE)*FSEMICORE
      end do
      !-------------------------------------------------------------------------
      ! update energy contour for Fermi-surface generation                       ! fswrt=fermi-surface write
      !-------------------------------------------------------------------------
      if (OPT('FERMIOUT'))then                                                   ! fswrt
         if(AIMAG(EZ(1))>0d0) stop 'E has imaginary part'                        ! fswrt
         IELAST=3                                                                ! fswrt
         EZ(2) = EZ(1) + CMPLX(1.0D-03,0.0D0, kind=dp)                                    ! fswrt
         EZ(3) = EZ(1) - CMPLX(1.0D-03,0.0D0, kind=dp)                                    ! fswrt
      end if                                                                     ! fswrt
      !-------------------------------------------------------------------------
      ! update the value of NSPIN to be consistent with REL mode
      !-------------------------------------------------------------------------
      IF (KREL.EQ.1) NSPIN = 1
      !
      CALL GAUNT2(WG,YRG,4*LMAX)
      CALL GAUNT(LMAX,LPOT,WG,YRG,CLEB,LOFLM,ICLEB,IEND,JEND,NCLEB,LMAX,&
         LMGF0D,LMPOT)
      !
      !-------------------------------------------------------------------------
      ! set up of GAUNT coefficients C(l,m;l',m';l'',m'') for all
      ! nonvanishing (l'',m'')-components of the shape functions THETAS
      !-------------------------------------------------------------------------
      if (KSHAPE.NE.0) then
         CALL SHAPE_CORR(LPOT,NATYP,GSH,ILM_MAP,IMAXSH,LMSP,NTCELL,WG,YRG,LASSLD,&
            LMPOT,NATYP,NGSHD)
      endif
      !-------------------------------------------------------------------------
      ! calculate Madelung constants (needed only for SCF calculations)
      !-------------------------------------------------------------------------
      !fivos      IF ( SCFSTEPS.GT.1 .OR. ICC.GT.0 ) THEN
      if (NPOL.NE.0) then  ! No madelung calculation in case of DOS.
         !OPEN(99,FILE='madelinfo.txt')

         !> @note Use option 'ewald2d' if the madelung summation is to be carried out in
         !> single-slab mode, otherwise it is carried out in repeated (periodic)
         !> slab mode.
         !> Reason: the 2d-mode gives wrong results sometimes [e.g. in diamond
         !> structure (110)].
         if (LINTERFACE.AND.( OPT('ewald2d ').OR.OPT('DECIMATE'))) then ! ewald2d
            write(*,*) 'Calling MADELUNG2D'
            !-------------------------------------------------------------------
            ! 2D case
            !-------------------------------------------------------------------
            call MADELUNG2D(LPOT,YRG,WG,NAEZ,ALAT,VOLUME0,BRAVAIS,RECBV,RBASIS,  &
               RMAX,GMAX,NLBASIS,NLEFT,ZPERLEFT,TLEFT,NRBASIS,NRIGHT,ZPERIGHT,   &
               TRIGHT,LMXSPD,LASSLD,LPOT,LMPOT,NMAXD,ISHLD,NEMBD1,WLENGTH)
            write(*,*) 'Exited MADELUNG2D'
         else
            !-------------------------------------------------------------------
            ! 3D case
            !-------------------------------------------------------------------
            if ( LINTERFACE ) then
               call GETBR3(NEMBD1,NLBASIS,ALAT,TLEFT,NRBASIS,TRIGHT,BRAVAIS,&
                  RECBV,VOLUME0)
            end if

            write(*,*) 'Calling MADELUNG3D'
            call MADELUNG3D(LPOT,YRG,WG,NAEZ,ALAT,VOLUME0,BRAVAIS,RECBV,&
               RBASIS,RMAX,GMAX,NAEZ,LMXSPD,LASSLD,LPOT,LMPOT,NMAXD,    &
               ISHLD,NEMB,WLENGTH)
            write(*,*) 'Exited MADELUNG3D'
         end if

         !CLOSE(99)
      else !NPOL==0
         ! write dummy files

         !real (kind=dp) AVMAD(LMPOT,LMPOT),BVMAD(LMPOT)
         LRECABMAD = WLENGTH*2*LMPOT*LMPOT + WLENGTH*2*LMPOT
         open (69,ACCESS='direct',RECL=LRECABMAD,FILE='abvmad.unformatted', &
            FORM='unformatted')
         do I = 1,NAEZ
            do J = 1,NAEZ
               IREC = J + NAEZ*(I-1)
               write (69,REC=IREC) 0.0d0, 0.0d0!AVMAD,BVMAD
            end do
         end do
         close(69)

      endif !npol==0
      !-------------------------------------------------------------------------
      !fivos      END IF
      !-------------------------------------------------------------------------
      !
      !-------------------------------------------------------------------------
      ! Set up I,J pairs for ICC = -1
      !-------------------------------------------------------------------------
      if ( ICC.LT.0 ) then
         call SETGIJTAB(LINTERFACE,ICC,NAEZ,IQAT,RBASIS,BRAVAIS,NATOMIMP,  &
            ATOMIMP,RCLSIMP,NOFGIJ,IJTABCALC,IOFGIJ,JOFGIJ,NQCALC,IQCALC,  &
            NATOMIMPD,IJTABCALC_I)
      endif
      !
      !-------------------------------------------------------------------------
      !
      DSYMLL=(0d0,0d0)
      DSYMLL1=(0d0,0d0)

      call BZKINT0(NSHELL,NAEZ,NATYP,NOQ,RBASIS,KAOEZ,ICC,BRAVAIS,RECBV,ATOMIMP, &
         RSYMAT,ISYMINDEX,NSYMAT,I25,NATOMIMP,NSH1,NSH2,RCLSIMP,RATOM,IJTABSYM,  &
         IJTABSH,IJTABCALC,IOFGIJ,JOFGIJ,NOFGIJ,ISH,JSH,RROT,DSYMLL1,PARA,QMTET, &
         QMPHI,SYMUNITARY,HOSTIMP,INTERVX,INTERVY,INTERVZ,IELAST,EZ,KMESH,       &
         MAXMESH,MAXMSHD,NSYMAXD,KREL+KORBIT,LMAX,LMMAXD,KPOIBZ,NAEZ,NATYP,      &
         NATOMIMPD,NSHELD,NEMB)
      !
      !-------------------------------------------------------------------------
      !
      if ( OPT('KKRFLEX ') ) then

         call WRITEHOSTSTRUCTURE(BRAVAIS,NATYP,RBASIS,NAEZ,NEMB)

         open (58,FILE='kkrflex_atominfo',FORM='FORMATTED')
         call version_print_header(58,&
            '; '//md5sum_potential//'; '//md5sum_shapefun)
         NVATOM=0
         do I = 1,NATOMIMP
            if (KAOEZ(1,ATOMIMP(I))==-1) NVATOM=NVATOM+1
         end do
         write (58,'(500A)') '#NATOM   NTOTATOM'
         write (58,*) NATOMIMP,NATOMIMP-NVATOM
         write (58,'(500A)') '#Impurity positions x,y,z|Core Charge|Virtual Atom?|Remove Atom?|LMAX'
         do I = 1,NATOMIMP
            if (KAOEZ(1,ATOMIMP(I))==-1) then
               ZATTEMP=0.D0
               ISVATOM=1
               NVATOM=NVATOM+1
            else
               ISVATOM=0
               ZATTEMP=ZAT(KAOEZ(1,ATOMIMP(I)))
            end if
            write (58,'(3F14.7,F6.2,3I5)') (RCLSIMP(J,I),J=1,3), &
               ZATTEMP,ISVATOM,0,LMAX
         end do
         close (58)
      end if
      !
      !-------------------------------------------------------------------------
      ! fivos: write out nshell and nsh1,nsh2 into standard output and in file shells.dat
      !-------------------------------------------------------------------------
      if (ICC.NE.0 .and. .not.OPT('KKRFLEX ')) then
         open(58,FILE='shells.dat')
         write(1337,*) 'Writing out shells (also in shells.dat):'                ! fivos
         write(1337,*) 'itype,jtype,iat,jat,r(iat),r(jat)'                       ! fivos
         write(1337,*) NSHELL(0), 'NSHELL(0)'                                    ! fivos
         write(58,*) NSHELL(0), 'NSHELL(0)'                                      ! fivos
         do i1 = 1,NSHELL(0)                                                     ! fivos
            write(1337,*) i1,NSHELL(i1),'No. of shell, No. of atoms in shell'    ! fivos
            write(58,*) i1,NSHELL(i1),'No. of shell, No. of atoms in shell'      ! fivos
            do lm = 1,NSHELL(i1)                                                 ! fivos
               write(1337,*) 'ish(i1,lm)',ish(i1,lm)
               if(ISH(i1,lm)>0 .and. JSH(i1,lm)>0) then                          ! fix bernd
                  write(1337,8614) NSH1(i1),NSH2(i1),ISH(i1,lm),JSH(i1,lm),&     ! fivos
                     (RCLSIMP(i,ISH(i1,lm)),i=1,3),(RCLSIMP(i,JSH(i1,lm)),i=1,3) ! fivos
                  write(58,8614) NSH1(i1),NSH2(i1),ISH(i1,lm),JSH(i1,lm),  &     ! fivos
                     (RCLSIMP(i,ISH(i1,lm)),i=1,3),(RCLSIMP(i,JSH(i1,lm)),i=1,3) ! fivos
               else                                                              ! fix bernd
                  write(1337,8615) NSH1(i1),NSH2(i1),ISH(i1,lm),JSH(i1,lm)       ! fix bernd
                  write(58,8615) NSH1(i1),NSH2(i1),ISH(i1,lm),JSH(i1,lm)         ! fix bernd
               end if                                                            ! fix bernd
               8614          format(4i5,6f16.6)                                  ! fivos
               8615          format(4i5)                                         ! fix bernd
            enddo                                                                ! fivos
         enddo                                                                   ! fivos
         write(1337,*) '###################'
         close(58)
      endif
      !-------------------------------------------------------------------------
      ! end fivos
      !-------------------------------------------------------------------------

      call GFMASK(LINTERFACE,ICHECK,ICC,INVMOD,NSH1,NSH2,NAEZ,NSHELL(0),NAEZ,&
         NPRINCD)
      !-------------------------------------------------------------------------
      ! set up transformation matrices between REL/NREL representations
      !-------------------------------------------------------------------------
      if((KREL+KORBIT).EQ.1) then
         call DRVBASTRANS(RC,CREL,RREL,SRREL,NRREL,IRREL,LMAX+1,LMMAXD,&
            2*(LMAX+1),LMMAXD+2*(LMAX+1),MMAXD,2*(LMAX+1)*MMAXD)
      endif
      if (OPT('NEWSOSOL')) then
         do NS=1,NSYMAT
            call CHANGEREP(DSYMLL1(1,1,NS),'REL>RLM',DSYMLL(1,1,NS),LMMAXD,&
               LMMAXD,RC,CREL,RREL,'DSYMLL',0)
         enddo
         !       DSYMLL(:,:,:)=DSYMLL1(:,:,:)
      else
         DSYMLL(:,:,:)=DSYMLL1(:,:,:)
      endif
      !
      !-------------------------------------------------------------------------
      !  for the case that the magnetisation is rotated with respect to
      !  the (001)-direction (KMROT<>0) calculate the rotation matrices
      !  to switch between the CRYSTAL and LOCAL frames of reference
      !-------------------------------------------------------------------------
      call CINIT(LMMAXD*LMMAXD*NAEZ,DROTQ)

      if (KMROT.NE.0) then
         FACT(0) = 1.0D0
         DO I = 1,100
            FACT(I) = FACT(I-1)*DBLE(I)
         END DO
         !
         DO I1=1,NAEZ
            CALL CALCROTMAT(MMAXD,(KREL+KORBIT)*3,QMPHI(I1),QMTET(I1),&
            0.0D0,DROTQ(1,1,I1),FACT,LMMAXD)
         END DO
      END IF
      !-------------------------------------------------------------------------
      ! Treat decimation I/O cases
      !-------------------------------------------------------------------------
      if ( OPT('deci-pot') ) then
         CALL OUTPOTHOST(ALAT,INS,KREL+KORBIT,KMROT,NSPIN,NAEZ,NATYP,E2IN, &
            BRAVAIS,RBASIS,QMTET,QMPHI,NOQ,KAOEZ,IQAT,ZAT,CONC,IPAN,IRCUT, &
            SOLVER,SOCSCL,CSCL,IRWS,RMTNEW,RWS,RMESH,DRDI,VISP,IRSHIFT,RMREL,  &
            DRDIREL,VTREL,BTREL,LMAX,NATYP,NAEZ,IPAND,IRMD)
      endif
      if ( OPT('deci-out') ) then
         CALL OUTTMATHOST(ALAT,INS,KREL+KORBIT,KMROT,NSPIN,NAEZ,LMMAX,     &
            BRAVAIS,RBASIS,QMTET,QMPHI,E2IN,TK,NPOL,NPNT1,NPNT2,NPNT3)
      endif
      if ( OPT('DECIMATE') ) then
         CALL DECIOPT(ALAT,INS,KREL+KORBIT,KVREL,KMROT,NSPIN,NAEZ,LMMAX,   &
            BRAVAIS,TK,NPOL,NPNT1,NPNT2,NPNT3,EZ,IELAST,KAOEZ,LEFTTINVLL,  &
            RIGHTTINVLL,VACFLAG,NLBASIS,NRBASIS,CMOMHOST,VREF,RMTREF,NREF, &
            REFPOT(NAEZ),LMAX,LMGF0D,LMMAXD,LM2D,NEMBD1,IEMXD,NSPINDD,     &
            LMPOT,NATYP,IRMD,IPAND)
      endif
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      ! ITERMDIR  -- initialise
      !-------------------------------------------------------------------------
      IF (OPT('ITERMDIR')) THEN
         WRITE (1337,*)
         WRITE (1337,*) 'Angle mixing scheme will be applied '
         WRITE (1337,*)
         DO I1 = 1,NAEZ
            QMPHITAB(I1,1) = QMPHI(I1)
            QMTETTAB(I1,1) = QMTET(I1)
            QMGAMTAB(I1,1) = QMGAM(I1)
            DO I = 2,3
               QMPHITAB(I1,I) = 0D0
               QMTETTAB(I1,I) = 0D0
               QMGAMTAB(I1,I) = 0D0
            END DO
         END DO
      END IF
      !
      !-------------------------------------------------------------------------
      ! LDA+U -- initialise
      !-------------------------------------------------------------------------
      !
      if (OPT('LDA+U   ')) then
         CALL STARTLDAU(ITRUNLDAU,IDOLDAU,KREADLDAU,LOPT,UEFF,JEFF,EREFLDAU,  &
            NATYP,NSPIN,WLDAU,ULDAU,PHILDAU,IRWS,NTLDAU,ITLDAU,IRMD,NATYP,     &
            NSPIND,MMAXD)
      end if
      !-------------------------------------------------------------------------
      ! LDA+U
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !          write out information for the other program parts           =
      !-------------------------------------------------------------------------
      !
      ! new solver for full-potential, spin-orbit, initialise
      if (OPT('NEWSOSOL')) THEN
         call CREATE_NEWMESH(NATYP,IRMD,IPAND,IRID,NTOTD,NFUND,&
            NCHEB,NTOTD*(NCHEB+1),NSPIN,RMESH,IRMIN,IPAN,IRCUT,R_LOG,NPAN_LOG,    &
            NPAN_EQ,NPAN_LOG_AT,NPAN_EQ_AT,NPAN_TOT,RNEW,RPAN_INTERVALL,      &
            IPAN_INTERVALL,NCELLD,NTCELL,THETAS,THETASNEW)
      end if

      call WUNFILES(NPOL,NPNT1,NPNT2,NPNT3,IELAST,TK,EMIN,EMAX,EZ,WEZ,EFERMI, &
         NPOLSEMI,N1SEMI,N2SEMI,N3SEMI,IESEMICORE,TKSEMI,EBOTSEMI,EMUSEMI,FSEMICORE,&
         VINS,VISP,VBC,VTREL,BTREL,RMREL,DRDIREL,R2DRDIREL,ZREL,JWSREL,IRSHIFT,     &
         ITSCF,SCFSTEPS,CMOMHOST,ECORE,LCORE,NCORE,QMTET,QMPHI,QMPHITAB,QMTETTAB,   &
         QMGAMTAB,DROTQ,NSRA,INS,NATYPD,NAEZD,NINEQ,NREFD,NSPIN,NCLS,ICST,IPAN,IRCUT,  &
         ALAT,ZAT,RMESH,DRDI,REFPOT,RMTREF,VREF,IEND,JEND,CLEB,ICLEB,ATOM,CLS,RCLS,     &
         NACLS,LOFLM,SOLVER,SOCSCL,CSCL,ICC,IGF,NLBASIS,NRBASIS,NCPA,ICPA,ITCPAMAX, &
         CPATOL,RBASIS,RR,EZOA,NSHELL,NSH1,NSH2,IJTABCALC,IJTABCALC_I,ISH,JSH,      &
         IJTABSYM,IJTABSH,NOFGIJ,NQCALC,IQCALC,KMROT,KAOEZ,IQAT,NOQ,CONC,KMESH,     &
         MAXMESH,NSYMAT,SYMUNITARY,RROT,DSYMLL,INVMOD,ICHECK,NATOMIMP,RATOM,ATOMIMP,&
         RC,CREL,RREL,SRREL,NRREL,IRREL,LEFTTINVLL,RIGHTTINVLL,VACFLAG,A,B,IFUNM,   &
         IFUNM1,INTERVX,INTERVY,INTERVZ,ITITLE,LMSP1,NTCELL,THETAS,LPOTD,LMPOTD,      &
         NRIGHT,NLEFT,LINTERFACE,IMIX,MIXING,QBOUND,FCM,ITDBRY,IRNS,KPRE,KSHAPE,KTE,&
         KVMAD,KXC,LAMBDA_XC,TXC,ISHIFT,IXIPOL,LRHOSYM,KFORCE,LMSP,LLMSP,RMT,RMTNEW,&
         RWS,IMT,IRC,IRMIN,IRWS,NFU,HOSTIMP,GSH,ILM_MAP,IMAXSH,IDOLDAU,ITRUNLDAU,NTLDAU,&
         LOPT,ITLDAU,UEFF,JEFF,EREFLDAU,ULDAU,WLDAU,PHILDAU,IEMXD,IRMIND,IRMD,NSPOTD,&
         NPOTD,NEMBD1,LMMAXD,IPAND,NEMBD2,LMAX,NCLEB,NACLSD,NCLSD,LM2D,LMAX+1,MMAXD,&
         NRD,NSHELD,NSYMAXD,NAEZ/NPRINCD,NATOMIMPD,NSPIND,IRID,NFUND,NCELLD,LMXSPD,NGSHD, &
         KREL,NTOTD,NCHEB,NPAN_LOG,NPAN_EQ,NPAN_LOG_AT,NPAN_EQ_AT,R_LOG,NPAN_TOT,   &
         RNEW,RPAN_INTERVALL,IPAN_INTERVALL,NSPINDD,THETASNEW,SOCSCALE,TOLRDIF,LLY, &
         DELTAE,RCLSIMP)

      if (OPT('FERMIOUT'))then                                                   ! fswrt
         call WRITE_TBKKR_FILES(LMAX,NEMB,NCLS,NATYP,NAEZ,IELAST,INS,ALAT,    &  ! fswrt
            BRAVAIS,RECBV,RBASIS,CLS,NACLS,RCLS,EZOA,ATOM,RR,NSPIN,NRD,KORBIT, &  ! fswrt
            NCLSD,NACLSD)                                                        ! fswrt
      end if                                                                     ! fswrt
      !
      if (OPT('OPERATOR')) then
        ! check if impurity files are present (otherwise no imp.
        ! wavefunctions can be calculated)
        operator_imp = .true.
        inquire(file='potential_imp', exist=lexist)
        if (.not.lexist) operator_imp = .false.
        inquire(file='shapefun_imp', exist=lexist)
        if (.not.lexist) operator_imp = .false.
        inquire(file='scoef', exist=lexist)
        if (.not.lexist) operator_imp = .false.
      else
        operator_imp = .false.
      endif
      !
      if (OPT('GREENIMP') .or. operator_imp) then                      ! GREENIMP
         ! fill array dimensions and allocate arrays in t_imp          ! GREENIMP
         call init_params_t_imp(t_imp,IPAND,NATYP,IRMD,IRID,NFUND,     &! GREENIMP
                                NSPIN,IRMIND,LMPOT)                    ! GREENIMP
         call init_t_imp(t_inc,t_imp)                                  ! GREENIMP
                                                                       ! GREENIMP
         ! next read impurity potential and shapefunction              ! GREENIMP
         call readimppot(NATOMIMP,INS,1337,0,0,2,NSPIN,LPOT,          &! GREENIMP
                         t_imp%IPANIMP,t_imp%THETASIMP,t_imp%IRCUTIMP,&! GREENIMP
                         t_imp%IRWSIMP,KHFELD,HFIELD,t_imp%VINSIMP,   &! GREENIMP
                         t_imp%VISPIMP,t_imp%IRMINIMP,                &! GREENIMP
                         t_imp%RIMP,t_imp%ZIMP, IRMD, IRNSD,IRID,NFUND,IPAND) ! GREENIMP
      end if                                                           ! GREENIMP
      !
      !
      if (ISHIFT.EQ.2) then                                                      ! fxf
         open (67,FILE='vmtzero',FORM='formatted')                               ! fxf
         write (67,9090) VBC(1)                                                  ! fxf
         close(67)                                                               ! fxf
         9090    format(D20.12)                                                  ! fxf
      end if                                                                     ! fxf
      !
      ! Check for inputcard consistency in case of qdos option
      if (OPT('qdos    ')) then
         write(1337,*)
         write(1337,*) '     < QDOS > : consistency check '
         if ((NPOL.NE.0).AND.(NPNT1.EQ.0).AND.(NPNT3.EQ.0)) then
            stop 'For qdos calculation change enery contour to dos path'
         endif
         if (TK.GT.50.d0) write(*,*) 'WARNING:  high energy smearing due to high value of TEMPR for energy contour integration could not be of advantage. Consider changeing ''TEMPR'' to lower value'
         if (TK.GT.50.d0) write(1337,*) 'WARNING:  high energy smearing due to high value of TEMPR for energy contour integration could not be of advantage. Consider changeing ''TEMPR'' to lower value'
         write(1337,*) '       QDOS: consistecy check complete'
      endif
      !
      !-------------------------------------------------------------------------
      !
      write (1337,'(79("="),/,31X,"< KKR0 finished >",/,79("="),/)')
      9070 format (5X,'INFO:  Output of cluster Green function at E Fermi')
      9080 format (5X,'INFO:  Determination of DOS at E Fermi')

   end subroutine main0

   !----------------------------------------------------------------------------
   ! SUBROUTINE: BSHIFT_NS
   !> @brief Adds a constant (=VSHIFT) to the potentials of atoms
   !----------------------------------------------------------------------------
   subroutine bshift_ns(IRM,IRID,IPAND,LMPOT,NPOTD,NATYP,NSPIN,NGSHD,NFUND,NCELLD,  &
      IRMIND,LMXSPD,KSHAPE,IRC,IRMIN,INIPOL,NTCELL,IMAXSH,ILM_MAP,LMSP,IFUNM,IRCUT,     &
      HFIELD,GSH,RMESH,THESME,THETAS,VISP,VINS)
 
      use mod_DataTypes

      implicit none

      ! Adds a constant (=VSHIFT) to the potentials of atoms
      !
      ! Parameters:
      ! Input
      integer, intent(in) :: IRM
      integer, intent(in) :: IRID
      integer, intent(in) :: IPAND
      integer, intent(in) :: LMPOT
      integer, intent(in) :: NPOTD
      integer, intent(in) :: NATYP !< Number of kinds of atoms in unit cell
      integer, intent(in) :: NSPIN !< Counter for spin directions
      integer, intent(in) :: NGSHD
      integer, intent(in) :: NFUND
      integer, intent(in) :: NCELLD
      integer, intent(in) :: IRMIND
      integer, intent(in) :: LMXSPD
      integer, intent(in) :: KSHAPE !< exact treatment of WS cell
      integer, dimension(NATYP), intent(in) :: IRC !< r point for potential cutting
      integer, dimension(NATYP), intent(in) :: IRMIN !< max r for spherical treatment
      integer, dimension(NATYP), intent(in) :: INIPOL !< initial spin polarisation
      integer, dimension(NATYP), intent(in) :: NTCELL !< index for WS cell
      integer, dimension(0:LMPOT), intent(in) :: IMAXSH
      integer, dimension(NGSHD,3), intent(in) :: ILM_MAP
      integer, dimension(NATYP,LMXSPD), intent(in) :: LMSP !< 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
      integer, dimension(NATYP,LMXSPD), intent(in) :: IFUNM
      integer, dimension(0:IPAND,NATYP), intent(in) :: IRCUT !< r points of panel borders
      real (kind=dp), intent(in) :: HFIELD !< External magnetic field, for initial potential shift in spin polarised case
      real (kind=dp), dimension(NGSHD), intent(in) :: GSH
      real (kind=dp), dimension(IRM,NATYP), intent(in) :: RMESH
      real (kind=dp), dimension(IRID,NFUND,NCELLD), intent(in) :: THESME
      real (kind=dp), dimension(IRID,NFUND,NCELLD), intent(in) :: THETAS !< shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion

      ! Input/Output:
      real (kind=dp), dimension(IRM,NPOTD), intent(inout) :: VISP !< Spherical part of the potential
      real (kind=dp), dimension(IRMIND:IRM,LMPOT,NSPOTD), intent(inout) :: VINS !< Non-spherical part of the potential

      ! Inside
      integer :: ISPIN,IH,IPOT,IR,LM,IMT1,IRC1,IRMIN1
      real (kind=dp)  :: RFPI, VSHIFT
      real (kind=dp), dimension(IRM) :: PSHIFTR
      real (kind=dp), dimension(IRM,LMPOT) :: PSHIFTLMR

      RFPI = SQRT(16.0D0*ATAN(1.0D0))

      do IH = 1,NATYP

         IMT1 = IRCUT(1,IH)
         IRC1 = IRC(IH)
         IRMIN1 = IRMIN(IH)

         do ISPIN = 1,NSPIN
            ! shift potential spin dependent
            VSHIFT = -DBLE(2*ISPIN-3)*HFIELD*INIPOL(IH)

            write (1337,*) 'SHIFTING OF THE POTENTIALS OF ATOM',IH,  &
               'spin',ispin,' BY', VSHIFT, 'RY.'
            IPOT = NSPIN * (IH-1) + ISPIN

            call RINIT(IRM*LMPOT,PSHIFTLMR)
            call RINIT(IRM,PSHIFTR)
            do IR = 1,IRC1
               PSHIFTLMR(IR,1) = VSHIFT
            enddo

            if (KSHAPE.EQ.0) then ! ASA
               do IR = 1,IRC1
                  VISP(IR,IPOT) = VISP(IR,IPOT) + PSHIFTLMR(IR,1)
               end do
            else                ! Full-potential
               !
               call CONVOL(IMT1,IRC1,NTCELL(IH),IMAXSH(LMPOT),ILM_MAP,IFUNM,LMPOT,&
                  GSH,THETAS,THESME,0.d0,RFPI,RMESH(1,IH),PSHIFTLMR,PSHIFTR,  &
                  LMSP)
               !
               do IR = 1,IRC1
                  VISP(IR,IPOT) = VISP(IR,IPOT) + PSHIFTLMR(IR,1)
               enddo
               !
               do LM = 2,LMPOT
                  do IR = IRMIN1,IRC1
                     VINS(IR,LM,IPOT)=VINS(IR,LM,IPOT)+PSHIFTLMR(IR,LM)*RFPI
                  enddo
               enddo
            end if              ! (kshape.eq.0)
         end do
      end do

   end subroutine bshift_ns

   !-------------------------------------------------------------------------
   ! SUBROUTINE: init_misc_variables
   !> @brief Set default values for misc variables for the calculation
   !> @details The idea behind this kind of routine is to separate the initialization
   !> of variables such that one can use this to modularize the code
   !> @author Jonathan Chico
   !> @date 22.12.2017
   !-------------------------------------------------------------------------
   subroutine init_misc_variables()

      implicit none

      IPE         = 0
      IPF         = 0
      IPFE        = 0
      KCOR        = 0
      KEFG        = 0
      KHYP        = 0
      KPRE        = 0
      NSRA        = 1
      IEND        = 1
      NVIRT       = 0
      KVMAD       = 0
      ITSCF       = 0
      NSHELD      = 301 ! Number of blocks of the GF matrix that need to be calculated (NATYPD + off-diagonals in case of impurity)
      INVMOD      = 2   ! Corner band matrix inversion scheme
      IELAST      = 0
      ISHIFT      = 0
      KFROZN      = 0
      NSYMAT      = 0
      NQCALC      = 0
      KFORCE      = 0   ! Calculation of the forces
      PARA        = .True.
      LRHOSYM     = .False.

   end subroutine init_misc_variables

   !-------------------------------------------------------------------------
   ! subroutine: init_relativistic_variables
   !> @brief set default values for the relativistic variables
   !> @details The idea behind this kind of routine is to separate the initialization
   !> of variables such that one can use this to modularize the code
   !> @author Jonathan Chico
   !> @date 26.12.2017
   !-------------------------------------------------------------------------
   subroutine init_relativistic_variables()

      implicit none

      KREL        = 0      ! Switch for non- (or scalar-) relativistic/relativistic (Dirac) program (0/1). Attention: several other parameters depend explicitly on KREL, they are set automatically Used for Dirac solver in ASA
      KVREL       = 1      ! Scalar-relativistic calculation
      KORBIT      = 0      ! No spin-orbit coupling
      LNC         = .True. ! Coupled equations in two spins (switches true if KREL=1 or KORBIT=1 or KNOCO=1)

   end subroutine init_relativistic_variables

   !-------------------------------------------------------------------------
   ! subroutine: init_cluster_variables
   !> @brief set default values for the cluster variables
   !> @details The idea behind this kind of routine is to separate the initialization
   !> of variables such that one can use this to modularize the code
   !> @author Jonathan Chico
   !> @date 26.12.2017
   !-------------------------------------------------------------------------
   subroutine init_cluster_variables()

      implicit none

      NCLSD       = 2                                          ! NAEZD + NEMBD maximum number of different TB-clusters
      NACLSD      = 500                                        ! Maximum number of atoms in a TB-cluster
      NOFGIJ      = 2                                          ! NATOMIMPD*NATOMIMPD+1 probably the same variable than NOFGIJD
      NATOMIMP    = 0                                          ! Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
      NATOMIMPD   = 150                                        ! Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
      I25         = 'scoef                                   ' ! Default name of scoef file

   end subroutine init_cluster_variables

   !-------------------------------------------------------------------------
   ! subroutine: init_io_variables
   !> @brief set default values for the I/O variables
   !> @details The idea behind this kind of routine is to separate the initialization
   !> of variables such that one can use this to modularize the code
   !> @author Jonathan Chico
   !> @date 26.12.2017
   !-------------------------------------------------------------------------
   subroutine init_io_variables()

      implicit none

      IGF         = 0   ! Not printing the Green functions
      ICC         = 0   ! Not printing the Green functions
      WLENGTH     = 1   ! Word length for direct access files, compiler dependent ifort/others (1/4)
#ifdef __GFORTRAN__
     wlength = 4
#endif
      I12         = '                                        '
      I40         = '                                        '

   end subroutine init_io_variables

   !-------------------------------------------------------------------------
   ! subroutine: init_cell_variables
   !> @brief set default values for the unit cell variables
   !> @details The idea behind this kind of routine is to separate the initialization
   !> of variables such that one can use this to modularize the code
   !> @author Jonathan Chico
   !> @date 26.12.2017
   !-------------------------------------------------------------------------
   subroutine init_cell_variables()

      implicit none

      NAEZ        = 1         ! Number of atoms in the unit cell
      NEMB        = 1         ! Number of embedded atoms
      NCLS        = 1         ! Number of clusters
      NATYP       = 1         ! Number of kinds of atoms in the unit cell
      NINEQ       = 1
      ALAT        = 1.0D0     ! Lattice constant in a.u.
      ABASIS      = 1.0D0     ! Scaling factors for rbasis
      BBASIS      = 1.0D0     ! Scaling factors for rbasis
      CBASIS      = 1.0D0     ! Scaling factors for rbasis
      ALATNEW     = 1.0D0
      VOLUME0     = 1.0D0
      LCARTESIAN  = .False.   ! True: Basis in cartesian coords; false: in internal coords
      LINTERFACE  = .False.   ! If True a matching with semi-inifinite surfaces must be performed

   end subroutine init_cell_variables

   !-------------------------------------------------------------------------
   ! subroutine: init_slab_variables
   !> @brief set default values for the slab calculation variables
   !> @details The idea behind this kind of routine is to separate the initialization
   !> of variables such that one can use this to modularize the code
   !> @author Jonathan Chico
   !> @date 26.12.2017
   !-------------------------------------------------------------------------
   subroutine init_slab_variables()

      implicit none

      NLEFT       = 1
      NRIGHT      = 1
      NLAYER      = 1         ! Number of principal layer
      NLBASIS     = 0         ! Number of basis layers of left host (repeated units)
      NRBASIS     = 0         ! Number of basis layers of right host (repeated units)
      NPRINCD     = 1         ! Number of principle layers, set to a number >= NRPINC in output of main0
      NLAYERD     = 1         !(NAEZD/NPRINCD)

   end subroutine init_slab_variables

   !-------------------------------------------------------------------------
   ! subroutine: init_energy_variables
   !> @brief set default values for the energy variables
   !> @details The idea behind this kind of routine is to separate the initialization
   !> of variables such that one can use this to modularize the code
   !> @author Jonathan Chico
   !> @date 26.12.2017
   !-------------------------------------------------------------------------
   subroutine init_energy_variables()

      implicit none

      LLY         = 0         ! No Lloyds formula
      KTE         = 1         ! Calculate the total energy
      NPOL        = 7         ! Number of poles
      IEMXD       = 101       ! Dimension for energy-dependent arrays
      NPNT1       = 7         ! Number of points in the energy contour
      NPNT2       = 7         ! Number of points in the energy contour
      NPNT3       = 7         ! Number of points in the energy contour
      N1SEMI      = 0         ! Number of energy points for the semicore contour
      N2SEMI      = 0         ! Number of energy points for the semicore contour
      N3SEMI      = 0         ! Number of energy points for the semicore contour
      IVSHIFT     = 0
      NPOLSEMI    = 0
      IESEMICORE  = 0
      IDOSEMICORE = 0
      TK          = 800.D0    ! Temperature
      FCM         = 20.0D0
      E2IN        = 0.0D0
      EMIN        = -0.30D0   ! Energies needed in EMESHT
      EMAX        = 0.70D0    ! Energies needed in EMESHT
      EFERMI      = 0.D0      ! Setting the Fermi energy to zero
      ESHIFT      = 0.D0
      TKSEMI      = 800.D0    ! Temperature for the semi-core contour
      EMUSEMI     = 0.0D0
      EBOTSEMI    = -0.5D0
      FSEMICORE   = 0.0D0
      DELTAE      = (1.D-5,0.D0)

   end subroutine init_energy_variables

   !-------------------------------------------------------------------------
   ! subroutine: init_convergence_variables
   !> @brief set default values for the convergence and solver variables
   !> @details The idea behind this kind of routine is to separate the initialization
   !> of variables such that one can use this to modularize the code
   !> @author Jonathan Chico
   !> @date 26.12.2017
   !-------------------------------------------------------------------------
   subroutine init_convergence_variables()

      implicit none

      IMIX        = 0            ! Straight mixing
      ISHLD       = 200000       ! Paremeters for the Ewald summations
      NMAXD       = 2000000      ! Paremeters for the Ewald summations
      NCHEB       = 10           ! Number of Chebychev pannels for the new solver
      ITDBRY      = 10
      NTREFD      = 0            ! Parameter in broyden subroutine MUST BE 0 for the host program
      NTPERD      = 1            ! Parameter in broyden subroutines
      NPAN_EQ     = 30
      NPAN_LOG    = 30
      SCFSTEPS    = 100          ! Number of SCF steps
      RMAX        = 7.0D0        ! Ewald summation cutoff parameter for real space summation
      GMAX        = 100.D0       ! Ewald summation cutoff parameter for reciprocal space summation
      R_LOG       = 0.1D0
      RCUTZ       = 2.30D0       ! Parameter for the screening cluster along the z-direction
      RCUTXY      = 2.30D0       ! Parameter for the screening cluster along the x-y plane
      QBOUND      = 1.D-7        ! Convergence parameter for the potential
      MIXING      = 0.001D0      ! Magnitude of the mixing parameter
      TOLRDIF     = 1.0D0        ! Tolerance for r<tolrdif (a.u.) to handle vir. atoms
      SOLVER      = 'BS        ' ! Set the BS-solver as default

   end subroutine init_convergence_variables

   !-------------------------------------------------------------------------
   ! subroutine: init_potential_variables
   !> @brief set default values for the potential variables
   !> @details The idea behind this kind of routine is to separate the initialization
   !> of variables such that one can use this to modularize the code
   !> @author Jonathan Chico
   !> @date 26.12.2017
   !-------------------------------------------------------------------------
   subroutine init_potential_variables()

      implicit none

      KWS         = 1         ! FP/ASA potential ()
      KXC         = 2         ! VWN potential
      INS         = 1         ! FP/ASA calculation ()
      ICST        = 2         ! Number of Born approximation
      IRID        = 350       ! Shape functions parameters in non-spherical part
      NREF        = 1         ! Number of reference potentials
      LPOT        = 3         ! Maximum l for expansion of the potential
      NFUND       = 289       ! Shape functions parameters in non-spherical part
      NGSHD       = 60000     ! Shape functions parameters in non-spherical part
      IPAND       = 50        ! Number of panels in non-spherical part
      NTOTD       = 80        ! IPAND+30
      IFILE       = 13        ! File identifier of the potential file
      LMPOT       = 16        ! (LPOT+1)**2
      INSREF      = 0         ! INS For reference potential
      KNOSPH      = 1         ! Switch for spherical/non-spherical(0/1) program. Same obs. as for KREL applies.
      KSHAPE      = 2         ! FP/ASA calculation (2/0)
      NCELLD      = 1         ! Number of cells (shapes) in non-spherical part
      NSPOTD      = 2         ! Number of potentials for storing non-sph. potentials (2*KREL+(1-KREL)*NSPIND)*NSATYPD
      NSATYPD     = 1         ! (NATYPD-1)*KNOSPH+1
      VCONST      = 0.D0      ! Potential shift
      LAMBDA_XC   = 1.0D0     ! Scale magnetic moment (0 < Lambda_XC < 1, 0=zero moment, 1= full moment)
      I13         = 'potential                               ' ! Default name of potential file
      I19         = 'shapefun                                ' ! Default name of shapefunction file

   end subroutine init_potential_variables

   !-------------------------------------------------------------------------
   ! subroutine: init_angular_momentum_variables
   !> @brief set default values for the angular momentum variables
   !> @details The idea behind this kind of routine is to separate the initialization
   !> of variables such that one can use this to modularize the code
   !> @author Jonathan Chico
   !> @date 26.12.2017
   !-------------------------------------------------------------------------
   subroutine init_angular_momentum_variables()

      implicit none

      LMAX        = 3   ! Maximum l for expansion
      NCLEB       = 784 ! (LMAX*2+1)**2 * (LMAX+1)**2
      LMMAX       = 16  ! (LMAX+1)**2
      LMMAXD      = 32  ! (KREL+KORBIT+1)*(LMAX+1)**2

   end subroutine init_angular_momentum_variables

   !-------------------------------------------------------------------------
   ! subroutine: init_magnetization_variables
   !> @brief set default values for the magnetisation variables for the calculation
   !> @details the idea behind this kind of routine is to separate the initialization
   !> of variables such that one can use this to modularize the code
   !> @author Jonathan Chico
   !> @date 22.12.2017
   !-------------------------------------------------------------------------
   subroutine init_magnetization_variables()

      implicit none

      KMROT    = 0         ! 0: no rotation of the magnetisation; 1: individual rotation of the magnetisation for every site
      KNOCO    = 0         ! Collinear calculation
      NSPIN    = 2         ! Spin polarised calculation
      NSPIND   = 2         ! KREL+(1-KREL)*(NSPIN+1)
      KHFELD   = 0         ! No external magnetic field
      NSPINDD  = 1         ! NSPIND-KORBIT
      HFIELD   = 0.D0      ! External magnetic field, for initial potential shift in spin polarised case
      LINIPOL  = .False.   ! True: Initial spin polarization; false: no initial spin polarization

   end subroutine init_magnetization_variables

   !-------------------------------------------------------------------------
   ! SUBROUTINE: init_CPA_variables
   !> @brief Set default values for the CPA variables for the calculation
   !> @details The idea behind this kind of routine is to separate the initialization
   !> of variables such that one can use this to modularize the code
   !> @author Jonathan Chico
   !> @date 22.12.2017
   !-------------------------------------------------------------------------
   subroutine init_CPA_variables()

      implicit none

      NCPA     = 0      ! NCPA = 0/1 CPA flag
      ITCPAMAX = 0      ! Max. number of CPA iterations
      CPATOL   = 1D-4   ! Convergency tolerance for CPA-cycle

   end subroutine init_CPA_variables

   !-------------------------------------------------------------------------
   ! SUBROUTINE: init_CPA_variables
   !> @brief Set default values for the LDA+U variables for the calculation
   !> @details The idea behind this kind of routine is to separate the initialization
   !> of variables such that one can use this to modularize the code
   !> @author Jonathan Chico
   !> @date 22.12.2017
   !-------------------------------------------------------------------------
   subroutine init_LDAU_variables()

      implicit none

      NTLDAU      = 0   ! number of atoms on which LDA+U is applied
      IDOLDAU     = 0   ! flag to perform LDA+U
      ITRUNLDAU   = 0   ! iteration index
      KREADLDAU   = 0   ! LDA+U arrays available

   end subroutine init_LDAU_variables

   !-------------------------------------------------------------------------
   ! subroutine: init_mesh_variables
   !> @brief set default values for the mesh variables
   !> @details The idea behind this kind of routine is to separate the initialization
   !> of variables such that one can use this to modularize the code
   !> @author Jonathan Chico
   !> @date 26.12.2017
   !-------------------------------------------------------------------------
   subroutine init_mesh_variables()

      implicit none

      NRD          = 20000     ! Number of real space
      IRMD         = 900       ! Number of radial mesh points in (0,...,RWS)
      IRNSD       = 890       ! Number of radial mesh points in (RMT,...,RWS)
      KPOIBZ      = 250000    ! Number of reciprocal space vectors
      INTERVX     = 0         ! Number of intervals in x-direction for k-net in IB of the BZ
      INTERVY     = 0         ! Number of intervals in y-direction for k-net in IB of the BZ
      INTERVZ     = 0         ! Number of intervals in z-direction for k-net in IB of the BZ
      MAXMESH     = 1

   end subroutine init_mesh_variables


   !-------------------------------------------------------------------------
   ! subroutine: init_all_wrapper
   !> @brief wrapper for initialization subroutines and allocation of test/opt arrays
   !> @author Philipp Ruessmann
   !> @date 22.08.2018
   !-------------------------------------------------------------------------
   subroutine init_all_wrapper()
      use mod_wunfiles, only: t_params
      implicit none
      integer :: i_stat

      ! Set the default values for the I/O variables
      call init_io_variables()
      ! Set the defaults values for the CPA variables
      call init_CPA_variables()
      ! Set the default values for the Slab mode variables
      call init_slab_variables()
      ! Set the default values for the LDA+U variables
      call init_LDAU_variables()
      ! Set the default values for misc variables for the calculation
      call init_misc_variables()
      ! Set the default values for the distinct mesh variables
      call init_mesh_variables()
      ! Set the default values for the unit cell variables
      call init_cell_variables()
      ! Set the default variables for the energy related variables
      call init_energy_variables()
      ! Set the default values for the cluster variables
      call init_cluster_variables()
      ! Set the default values for the potential variables
      call init_potential_variables()
      ! Set the default values for the convergence and solver variables
      call init_convergence_variables()
      ! Set the default values for the relativistic variables
      call init_relativistic_variables()
      ! Set the default values for the magnetisation variables
      call init_magnetization_variables()
      ! Set the default values for the angular momentum related variables
      call init_angular_momentum_variables()

      ! allocate and initialize testc and optc in t_params for run and test options
      allocate(t_params%OPTC(32), stat=i_stat) !CHARACTER*8
      call memocc(i_stat,product(shape(t_params%OPTC))*kind(t_params%OPTC),'t_params%OPTC','main0')
      t_params%OPTC(1:32) =  '        '
      allocate(t_params%TESTC(32), stat=i_stat)
      call memocc(i_stat,product(shape(t_params%TESTC))*kind(t_params%TESTC),'t_params%TESTC','main0')
      t_params%TESTC(1:32) = '        '

   end subroutine init_all_wrapper


   subroutine print_versionserial(iunit, version1, version2, version3, version4, serialnr)
      implicit none
      integer, intent(in) :: iunit
      character (len=*), intent(in) :: version1
      character (len=*), intent(in) :: version2
      character (len=*), intent(in) :: version3
      character (len=*), intent(in) :: version4
      character (len=*), intent(in) :: serialnr

      write(iunit,'(1A)') '     Screened Korringa-Kohn-Rostoker Electronic Structure Code'                          
      write(iunit,'(1A)') '                      for Bulk and Interfaces'                            
      write(iunit,'(1A)') '                    Juelich-Munich 2001 - 2018'                         
      write(iunit,'(1A)') ''
      write(iunit,'(2A)') '  Code version: ',trim(version1)                     
      write(iunit,'(6A)') '  Compile options: ', trim(version2), ' ', trim(version3), ' ', trim(version4)                                    
      write(iunit,'(2A)') '  serial number for files: ', serialnr
   end subroutine print_versionserial

end module mod_main0
