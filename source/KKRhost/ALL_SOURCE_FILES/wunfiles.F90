!-------------------------------------------------------------------------------
! MODULE: mod_wunfiles
!> @brief Module responsible for storing the input variables and primary arrays
!> so that they are distributed via MPI processes.
!> @details Previously this routine wrote unformatted files to disk, so that they
!> would be used by the different executables. Since the advent of the single
!> executable mode, this routine creates a copy of most of the variables in the program
!> as special `type` parameters. This are then used in the MPI communication, and
!> in the rest of the variables used in the code.
!> @author Philipp Rüssmann and many others ...
!> @note Jonatan Chico: 02.01.2018 Modifications to ensure compatibility for the removal of
!> the inc.p file. Also added the memory profiling calls to the allocation/deallocation
!> of the arrays.
!-------------------------------------------------------------------------------
module mod_wunfiles

   use mod_Profiling
   use mod_DataTypes

   implicit none


   !define type that replace wunfiles here, later define bcast routine
   type :: type_params

      integer :: Nscalars = 126

      !     .. Scalars
      integer :: I1
      integer :: NR        !< Number of real space vectors rr
      integer :: IRM       !< Maximum number of radial points
      integer :: LLY       !< LLY <> 0 : apply Lloyds formula
      integer :: INS       !< 0 (MT), 1(ASA), 2(Full Potential)
      integer :: ICC       !< Enables the calculation of off-diagonal elements of the GF.(0=SCF/DOS; 1=cluster; -1=custom)
      integer :: IGF       !< Do not print or print (0/1) the KKRFLEX_* files
      integer :: KTE       !< Calculation of the total energy On/Off (1/0)
      integer :: KXC       !< Type of xc-potential 0=vBH 1=MJW 2=VWN 3=PW91
      integer :: NAEZ      !< Number of atoms in unit cell
      integer :: LMAX      !< Maximum l component in wave function expansion
      integer :: NREF      !< Number of diff. ref. potentials
      integer :: LM2D      !< (2*LMAX+1)**2
      integer :: IRID      !< Shape functions parameters in non-spherical part
      integer :: KREL      !< Switch for non-relativistic/relativistic (0/1) program. Attention: several other parameters depend explicitly on KREL, they are set automatically Used for Dirac solver in ASA
      integer :: KPRE
      integer :: NSRA
      integer :: NEMB      !< Number of sites added to the slab in 2D calculations to extend the structure left and right (down and up)
      integer :: NCLS      !< Number of reference clusters
      integer :: IEND      !< Number of nonzero gaunt coefficients
      integer :: NCPA      !< NCPA = 0/1 CPA flag
      integer :: ICST      !< Number of Born approximation
      integer :: IMIX      !< Type of mixing scheme used (0=straight, 4=Broyden 2nd, 5=Anderson)
      integer :: ITAB
      integer :: LPOT      !< Maximum l component in potential expansion
      integer :: NPOL      !< Number of Matsubara Poles (EMESHT)
      integer :: NPNT1     !< number of E points (EMESHT) for the contour integration
      integer :: NPNT2     !< number of E points (EMESHT) for the contour integration
      integer :: NPNT3     !< number of E points (EMESHT) for the contour integration
      integer :: ITSCF
      integer :: IEMXD     !< Dimension for energy-dependent arrays
      integer :: NPOTD     !< (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP)
      integer :: NATYP     !< Number of kinds of atoms in unit cell
      integer :: IPAND     !< Number of panels in non-spherical part
      integer :: NCLEB     !< Number of Clebsch-Gordon coefficients
      integer :: NCLSD     !< Maximum number of different TB-clusters
      integer :: NFUND     !< Shape functions parameters in non-spherical part
      integer :: NGSHD     !< Shape functions parameters in non-spherical part
      integer :: MMAXD     !< 2*LMAX+1
      integer :: NINEQ     !< Number of ineq. positions in unit cell
      integer :: NSPIN     !< Counter for spin directions
      integer :: KMROT     !< 0: no rotation of the magnetisation; 1: individual rotation of the magnetisation for every site
      integer :: ILTMP
      integer :: NCHEB     !< Number of Chebychev pannels for the new solver
      integer :: NTOTD
      integer :: KVMAD
      integer :: IRNSD     !< Number of radial mesh points in (RMT,...,RWS)
      integer :: KNOCO     !< (0/1) Collinear/Non-collinear magnetism (even in non-relativistic non-spin-orbit case)
      integer :: LMPOT     !< (LPOT+1)**2
      integer :: NLEFT     !< Number of repeated basis for left host to get converged electrostatic potentials
      integer :: NRIGHT    !< Number of repeated basis for right host to get converged electrostatic potentials
      integer :: KORBIT    !< Spin-orbit/non-spin-orbit (1/0) added to the Schroedinger or SRA equations. Works with FP. KREL and KORBIT cannot be both non-zero.
      integer :: NTPERD    !< Parameter in broyden subroutines
      integer :: IELAST
      integer :: NRMAXD    !< NTOTD*(NCHEBD+1)
      integer :: ISHIFT
      integer :: KNOSPH    !< Switch for spherical/non-spherical (0/1) program. Same obs. as for KREL applies.
      integer :: KFORCE    !< Calculation of the forces
      integer :: ITDBRY    !< Number of SCF steps to remember for the Broyden mixing
      integer :: KSHAPE    !< Exact treatment of WS cell
      integer :: NOFGIJ    !< number of GF pairs IJ to be calculated as determined from IJTABCALC<>0
      integer :: NSPIND    !< KREL+(1-KREL)*(NSPIN+1)
      integer :: IRMIND    !< IRM-IRNSD
      integer :: NSPOTD    !< Number of potentials for storing non-sph. potentials
      integer :: NEMBD1    !< NEMB+1
      integer :: LMMAXD    !< (KREL+KORBIT+1)(LMAX+1)^2
      integer :: NEMBD2
      integer :: NACLSD    !< Maximum number of atoms in a TB-cluster
      integer :: LMAXD1
      integer :: NSHELD    !< Number of blocks of the GF matrix that need to be calculated (NATYP + off-diagonals in case of impurity)
      integer :: NCELLD    !< Number of cells (shapes) in non-spherical part
      integer :: LMXSPD    !< (2*LPOT+1)**2
      integer :: NSYMAT
      integer :: NPRINC    !< Number of atoms in one principal layer
      integer :: N1SEMI    !< Number of energy points for the semicore contour
      integer :: N2SEMI    !< Number of energy points for the semicore contour
      integer :: N3SEMI    !< Number of energy points for the semicore contour
      integer :: INVMOD    !< Inversion scheme
      integer :: NQCALC
      integer :: NTLDAU    !< number of atoms on which LDA+U is applied
      integer :: KPOIBZ    !< Number of reciprocal space vectors
      integer :: NSATYPD   !< (NATYP-1)*NSPIN+1
      integer :: IDOLDAU   !< flag to perform LDA+U
      integer :: NLAYERD   !< Number of principal layers (NAEZD/NPRINCD) used in the inversion routines (independent on NATYPD)
      integer :: INTERVX   !< Number of intervals in x-direction for k-net in IB of the BZ
      integer :: INTERVY   !< Number of intervals in y-direction for k-net in IB of the BZ
      integer :: INTERVZ   !< Number of intervals in z-direction for k-net in IB of the BZ
      integer :: NLBASIS   !< Number of basis layers of left host (repeated units)
      integer :: NRBASIS   !< Number of basis layers of right host (repeated units)
      integer :: NSYMAXD
      integer :: WLENGTH   !< Word length for direct access files, compiler dependent ifort/others (1/4)
      integer :: NAEZDPD
      integer :: MAXMESH
      integer :: ITMPDIR
      integer :: NSPINDD   !< NSPIND-KORBIT
      integer :: NPAN_EQ   !< Variables for the pannels for the new solver
      integer :: NPAN_LOG  !< Variables for the pannels for the new solver
      integer :: SCFSTEPS  !< number of scf iterations
      integer :: ITCPAMAX  !< Max. number of CPA iterations
      integer :: NATOMIMP  !< Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
      integer :: NMVECMAX
      integer :: NPOLSEMI  !< Number of poles for the semicore contour
      integer :: NATOMIMPD !< Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
      integer :: ITRUNLDAU !< Iteration index for LDA+U
      integer :: IESEMICORE
      real (kind=dp) :: TK        !< Temperature
      real (kind=dp) :: FCM
      real (kind=dp) :: EMIN      !< Energies needed in EMESHT
      real (kind=dp) :: EMAX      !< Energies needed in EMESHT
      real (kind=dp) :: ALAT      !< Lattice constant in a.u.
      real (kind=dp) :: R_LOG
      real (kind=dp) :: EFOLD
      real (kind=dp) :: DENEF
      real (kind=dp) :: EFERMI    !< Fermi energy
      real (kind=dp) :: CPATOL    !< Convergency tolerance for CPA-cycle
      real (kind=dp) :: MIXING    !< Magnitude of the mixing parameter
      real (kind=dp) :: QBOUND    !< Convergence parameter for the potential
      real (kind=dp) :: TKSEMI    !< Temperature of semi-core contour
      real (kind=dp) :: CHRGOLD
      real (kind=dp) :: TOLRDIF   !< Tolerance for r<tolrdif (a.u.) to handle vir. atoms
      real (kind=dp) :: LASTERR
      real (kind=dp) :: EMUSEMI
      real (kind=dp) :: EBOTSEMI
      real (kind=dp) :: FSEMICORE
      real (kind=dp) :: LAMBDA_XC !< Scale magnetic moment (0 < Lambda_XC < 1, 0=zero moment, 1= full moment)
      real (kind=dp) :: CHRGSEMICORE
      complex (kind=dp) :: DELTAE      !< Energy difference for numerical derivative
      logical :: LNC                !< Coupled equations in two spins (switches true if KREL=1 or KORBIT=1 or KNOCO=1)
      logical :: LRHOSYM
      logical :: LINTERFACE         !< If True a matching with semi-inifinite surfaces must be performed
      character(len=10) :: SOLVER   !< Type of solver
      character(len=80) :: TMPDIR

      !     .. Arrays
      complex (kind=dp), dimension(:), allocatable :: EZ
      complex (kind=dp), dimension(:), allocatable :: WEZ
      complex (kind=dp), dimension(:,:), allocatable :: RC        !< NREL REAL spher. harm. > CMPLX. spher. harm. NREL CMPLX. spher. harm. > REAL spher. harm.
      complex (kind=dp), dimension(:,:), allocatable :: CREL      !< Non-relat. CMPLX. spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. CMPLX. spher. harm.
      complex (kind=dp), dimension(:,:), allocatable :: RREL      !< Non-relat. REAL spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. REAL spher. harm.
      complex (kind=dp), dimension(:,:), allocatable :: PHILDAU
      complex (kind=dp), dimension(:,:,:), allocatable :: SRREL
      complex (kind=dp), dimension(:,:,:), allocatable :: DROTQ   !< Rotation matrices to change between LOCAL/GLOBAL frame of reference for magnetisation <> Oz or noncollinearity
      complex (kind=dp), dimension(:,:,:), allocatable :: DSYMLL
      complex (kind=dp), dimension(:,:,:,:,:), allocatable :: LEFTTINVLL
      complex (kind=dp), dimension(:,:,:,:,:), allocatable :: RIGHTTINVLL
      real (kind=dp), dimension(:), allocatable :: A               !< Constants for exponential R mesh
      real (kind=dp), dimension(:), allocatable :: B               !< Constants for exponential R mesh
      real (kind=dp), dimension(:), allocatable :: EU
      real (kind=dp), dimension(:), allocatable :: EDC
      real (kind=dp), dimension(:), allocatable :: VBC             !< Potential constants
      real (kind=dp), dimension(:), allocatable :: ZAT             !< Nuclear charge
      real (kind=dp), dimension(:), allocatable :: RMT             !< Muffin-tin radius of true system
      real (kind=dp), dimension(:), allocatable :: RWS             !< Wigner Seitz radius
      real (kind=dp), dimension(:), allocatable :: GSH
      real (kind=dp), dimension(:), allocatable :: PHI
      real (kind=dp), dimension(:), allocatable :: UEFF            !< input U parameter for each atom
      real (kind=dp), dimension(:), allocatable :: JEFF            !< input J parameter for each atom
      real (kind=dp), dimension(:), allocatable :: VREF
      real (kind=dp), dimension(:), allocatable :: CONC            !< Concentration of a given atom
      real (kind=dp), dimension(:), allocatable :: THETA
      real (kind=dp), dimension(:), allocatable :: VOLBZ
      real (kind=dp), dimension(:), allocatable :: QMTET           !< \f$ \theta\f$ angle of the agnetization with respect to the z-axis
      real (kind=dp), dimension(:), allocatable :: QMPHI           !< \f$ \phi\f$ angle of the agnetization with respect to the z-axis
      real (kind=dp), dimension(:), allocatable :: RMTREF          !< Muffin-tin radius of reference system
      real (kind=dp), dimension(:), allocatable :: RMTNEW          !< Adapted muffin-tin radius
      real (kind=dp), dimension(:), allocatable :: DENEFAT
      real (kind=dp), dimension(:), allocatable :: EREFLDAU        !< the energies of the projector's wave functions (REAL)
      real (kind=dp), dimension(:), allocatable :: SOCSCALE        !< Spin-orbit scaling
      real (kind=dp), dimension(:,:,:), allocatable :: VINS        !< Non-spherical part of the potential
      real (kind=dp), dimension(:,:), allocatable :: RMESH         !< Radial mesh ( in units a Bohr)
      real (kind=dp), dimension(:,:), allocatable :: RR
      real (kind=dp), dimension(:,:), allocatable :: DRDI          !< Derivative dr/di
      real (kind=dp), dimension(:,:), allocatable :: CSCL          !< Speed of light scaling
      real (kind=dp), dimension(:,:), allocatable :: CLEB          !< GAUNT coefficients (GAUNT)
      real (kind=dp), dimension(:,:), allocatable :: RHOC
      real (kind=dp), dimension(:,:), allocatable :: ESPV
      real (kind=dp), dimension(:,:), allocatable :: RNEW
      real (kind=dp), dimension(:,:), allocatable :: VISP          !< Spherical part of the potential
      real (kind=dp), dimension(:,:), allocatable :: VTREL         !< potential (spherical part)
      real (kind=dp), dimension(:,:), allocatable :: BTREL         !< magnetic field
      real (kind=dp), dimension(:,:), allocatable :: ECORE         !< Core energies
      real (kind=dp), dimension(:,:), allocatable :: RMREL         !< radial mesh
      real (kind=dp), dimension(:,:), allocatable :: RATOM
      real (kind=dp), dimension(:,:), allocatable :: RBASIS        !< Position of atoms in the unit cell in units of bravais vectors
      real (kind=dp), dimension(:,:), allocatable :: SOCSCL
      real (kind=dp), dimension(:,:), allocatable :: VOLCUB
      real (kind=dp), dimension(:,:), allocatable :: RHOORB
      real (kind=dp), dimension(:,:), allocatable :: RCLSIMP
      real (kind=dp), dimension(:,:), allocatable :: DRDIREL       !< derivative of radial mesh
      real (kind=dp), dimension(:,:), allocatable :: ECOREREL
      real (kind=dp), dimension(:,:), allocatable :: CMOMHOST      !< Charge moments of each atom of the (left/right) host
      real (kind=dp), dimension(:,:), allocatable :: QMPHITAB
      real (kind=dp), dimension(:,:), allocatable :: QMTETTAB
      real (kind=dp), dimension(:,:), allocatable :: QMGAMTAB
      real (kind=dp), dimension(:,:), allocatable :: R2DRDIREL     !< \f$ r^2 \frac{\partial}{\partial \mathbf{r}}\frac{\partial}{\partial i}\f$ (r**2 * drdi)
      real (kind=dp), dimension(:,:), allocatable :: RPAN_INTERVALL

      real (kind=dp), dimension(:,:,:), allocatable :: RCLS        !< Real space position of atom in cluster
      real (kind=dp), dimension(:,:,:), allocatable :: RROT
      real (kind=dp), dimension(:,:,:), allocatable :: BZKP
      real (kind=dp), dimension(:,:,:), allocatable :: MVEVI
      real (kind=dp), dimension(:,:,:), allocatable :: THETAS      !< shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
      real (kind=dp), dimension(:,:,:), allocatable :: MVEVIEF
      real (kind=dp), dimension(:,:,:), allocatable :: THETASNEW
      real (kind=dp), dimension(:,:,:,:), allocatable :: R2NEF
      real (kind=dp), dimension(:,:,:,:), allocatable :: WLDAU     !< potential matrix
      real (kind=dp), dimension(:,:,:,:), allocatable :: RHO2NS
      real (kind=dp), dimension(:,:,:,:,:), allocatable :: ULDAU   !< calculated Coulomb matrix elements (EREFLDAU)
      integer, dimension(:), allocatable :: CLS       !< Cluster around atomic sites
      integer, dimension(:), allocatable :: NOQ       !< Number of diff. atom types located
      integer, dimension(:), allocatable :: IMT       !< R point at MT radius
      integer, dimension(:), allocatable :: IRC       !< R point for potential cutting
      integer, dimension(:), allocatable :: NFU
      integer, dimension(:), allocatable :: ZREL      !< atomic number (cast integer)
      integer, dimension(:), allocatable :: LOPT      !< angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
      integer, dimension(:), allocatable :: IPAN      !< Number of panels in non-MT-region
      integer, dimension(:), allocatable :: IQAT      !< The site on which an atom is located on a given site
      integer, dimension(:), allocatable :: ICPA      !< ICPA = 0/1 site-dependent CPA flag
      integer, dimension(:), allocatable :: IRNS      !< Position of atoms in the unit cell in units of bravais vectors
      integer, dimension(:), allocatable :: IRWS      !< R point at WS radius
      integer, dimension(:), allocatable :: NSH1      !< Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
      integer, dimension(:), allocatable :: NSH2      !< Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
      integer, dimension(:), allocatable :: IRMIN     !< Max R for spherical treatment
      integer, dimension(:), allocatable :: NCORE     !< Number of core states
      integer, dimension(:), allocatable :: NACLS     !< Number of atoms in cluster
      integer, dimension(:), allocatable :: NOFKS
      integer, dimension(:), allocatable :: LOFLM     !< l of lm=(l,m) (GAUNT)
      integer, dimension(:), allocatable :: KMESH
      integer, dimension(:), allocatable :: ITLDAU    !< integer pointer connecting the NTLDAU atoms to heir corresponding index in the unit cell
      integer, dimension(:), allocatable :: NSHELL    !< Index of atoms/pairs per shell (ij-pairs); nshell(0) = number of shells
      integer, dimension(:), allocatable :: IQCALC
      integer, dimension(:), allocatable :: REFPOT    !< Ref. pot. card  at position
      integer, dimension(:), allocatable :: NTCELL    !< Index for WS cell
      integer, dimension(:), allocatable :: IXIPOL    !< Constraint of spin pol.
      integer, dimension(:), allocatable :: JWSREL    !< index of the WS radius
      integer, dimension(:), allocatable :: IMAXSH
      integer, dimension(:), allocatable :: ATOMIMP
      integer, dimension(:), allocatable :: HOSTIMP
      integer, dimension(:), allocatable :: IRSHIFT   !< shift of the REL radial mesh with respect no NREL
      integer, dimension(:), allocatable :: IJTABSH   !< Linear pointer, assigns pair (i,j) to a shell in the array GS(*,*,*,NSHELD)
      integer, dimension(:), allocatable :: NPAN_TOT
      integer, dimension(:), allocatable :: IJTABSYM  !< Linear pointer, assigns pair (i,j) to the rotation bringing GS into Gij
      integer, dimension(:), allocatable :: IJTABCALC !< Linear pointer, specifying whether the block (i,j) has to be calculated needs set up for ICC=-1, not used for ICC=1
      integer, dimension(:), allocatable :: NPAN_EQ_AT
      integer, dimension(:), allocatable :: NPAN_LOG_AT
      integer, dimension(:), allocatable :: IJTABCALC_I
      integer, dimension(:), allocatable :: qdos_atomselect
      integer, dimension(:,:), allocatable :: ISH
      integer, dimension(:,:), allocatable :: JSH
      integer, dimension(:,:), allocatable :: ILM_MAP
      integer, dimension(:,:), allocatable :: EZOA       !< EZ of atom at site in cluster
      integer, dimension(:,:), allocatable :: ATOM       !< Atom at site in cluster
      integer, dimension(:,:), allocatable :: LMSP       !< 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
      integer, dimension(:,:), allocatable :: ICLEB      !< Pointer array
      integer, dimension(:,:), allocatable :: LCORE      !< Angular momentum of core states
      integer, dimension(:,:), allocatable :: IRCUT      !< R points of panel borders
      integer, dimension(:,:), allocatable :: KAOEZ      !< Kind of atom at site in elem. cell
      integer, dimension(:,:), allocatable :: NRREL
      integer, dimension(:,:), allocatable :: LMSP1
      integer, dimension(:,:), allocatable :: IFUNM
      integer, dimension(:,:), allocatable :: LLMSP      !< lm=(l,m) of 'nfund'th nonvanishing component of non-spherical pot.
      integer, dimension(:,:), allocatable :: ICHECK
      integer, dimension(:,:), allocatable :: IFUNM1
      integer, dimension(:,:), allocatable :: ITITLE
      integer, dimension(:,:), allocatable :: NKCORE
      integer, dimension(:,:), allocatable :: KAPCORE
      integer, dimension(:,:), allocatable :: IPAN_INTERVALL
      integer, dimension(:,:,:), allocatable :: JEND     !< Pointer array for icleb()
      integer, dimension(:,:,:), allocatable :: IRREL
      logical, dimension(:), allocatable :: VACFLAG
      logical, dimension(:), allocatable :: SYMUNITARY   !< unitary/antiunitary symmetry flag
      character(len=8), dimension(:), allocatable :: OPTC
      character(len=8), dimension(:), allocatable :: TESTC
      character(len=124), dimension(:), allocatable :: TXC

   end type type_params

   type (type_params), save :: t_params

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: WUNFILES
   !> @brief This routine takes the read parameters from the inputcard and stores
   !> them in the t_params type to be distributed via MPI
   !> @details This routine was oiginally meant to write unformated files to then
   !> be read by other executables, now it does the same job via storing types instead
   !> reducing I/O and allowing for MPI communication.
   !> @author Philipp Rüssmann and many others ...
   !----------------------------------------------------------------------------
   subroutine WUNFILES(NPOL,NPNT1,NPNT2,NPNT3,IELAST,TK,EMIN,EMAX,EZ,WEZ,EFERMI, &
      NPOLSEMI,N1SEMI,N2SEMI,N3SEMI,IESEMICORE,TKSEMI,EBOTSEMI,EMUSEMI,FSEMICORE,&
      VINS,VISP,VBC,VTREL,BTREL,RMREL,DRDIREL,R2DRDIREL,ZREL,JWSREL,IRSHIFT,     &
      ITSCF,SCFSTEPS,CMOMHOST,ECORE,LCORE,NCORE,QMTET,QMPHI,QMPHITAB,QMTETTAB,   &
      QMGAMTAB,DROTQ,NSRA,INS,NATYP,NAEZ,NINEQ,NREF,NSPIN,NCLS,ICST,IPAN,IRCUT,  &
      ALAT,ZAT,R,DRDI,REFPOT,RMTREF,VREF,IEND,JEND,CLEB,ICLEB,ATOM,CLS,RCLS,     &
      NACLS,LOFLM,SOLVER,SOCSCL,CSCL,ICC,IGF,NLBASIS,NRBASIS,NCPA,ICPA,ITCPAMAX, &
      CPATOL,RBASIS,RR,EZOA,NSHELL,NSH1,NSH2,IJTABCALC,IJTABCALC_I,ISH,JSH,      &
      IJTABSYM,IJTABSH,NOFGIJ,NQCALC,IQCALC,KMROT,KAOEZ,IQAT,NOQ,CONC,KMESH,     &
      MAXMESH,NSYMAT,SYMUNITARY,RROT,DSYMLL,INVMOD,ICHECK,NATOMIMP,RATOM,ATOMIMP,&
      RC,CREL,RREL,SRREL,NRREL,IRREL,LEFTTINVLL,RIGHTTINVLL,VACFLAG,A,B,IFUNM,   &
      IFUNM1,INTERVX,INTERVY,INTERVZ,ITITLE,LMSP1,NTCELL,THETAS,LPOT,LMPOT,      &
      NRIGHT,NLEFT,LINTERFACE,IMIX,MIXING,QBOUND,FCM,ITDBRY,IRNS,KPRE,KSHAPE,KTE,&
      KVMAD,KXC,LAMBDA_XC,TXC,ISHIFT,IXIPOL,LRHOSYM,KFORCE,LMSP,LLMSP,RMT,RMTNEW,&
      RWS,IMT,IRC,IRMIN,IRWS,NFU,HOSTIMP,GSH,ILM_MAP,IMAXSH,IDOLDAU,ITRUNLDAU,NTLDAU,&
      LOPT,ITLDAU,UEFF,JEFF,EREFLDAU,ULDAU,WLDAU,PHILDAU,IEMXD,IRMIND,IRM,NSPOTD,&
      NPOTD,NEMBD1,LMMAXD,IPAND,NEMBD2,LMAX,NCLEB,NACLSD,NCLSD,LM2D,LMAXD1,MMAXD,&
      NR,NSHELD,NSYMAXD,NAEZDPD,NATOMIMPD,NSPIND,IRID,NFUND,NCELLD,LMXSPD,NGSHD, &
      KREL,NTOTD,NCHEB,NPAN_LOG,NPAN_EQ,NPAN_LOG_AT,NPAN_EQ_AT,R_LOG,NPAN_TOT,   &
      RNEW,RPAN_INTERVALL,IPAN_INTERVALL,NSPINDD,THETASNEW,SOCSCALE,TOLRDIF,LLY, &
      DELTAE,RCLSIMP)
      ! **********************************************************************
      ! *                                                                    *
      ! *  This subroutine is part of the MAIN0 program in the tbkkr package *
      ! *  It writes out different unformatted files meant to provide the    *
      ! *  communication between the other parts (MAIN1a, 1b, 1c and 2)      *
      ! *  during an SCF cycle                                               *
      ! *  v.popescu, munich 2004                                            *
      ! *                                                                    *
      ! **********************************************************************

      use mod_types, only: t_inc, t_tgmat, t_lloyd, t_cpa
      use mod_mympi, only: nranks, MPIatom, MPIadapt

      IMPLICIT NONE
      !     ..
      !     .. Scalar arguments
      integer, intent(in) :: NR        !< Number of real space vectors rr
      integer, intent(in) :: IRM       !< Maximum number of radial points
      integer, intent(in) :: INS       !< 0 (MT), 1(ASA), 2(Full Potential)
      integer, intent(in) :: ICC       !< Enables the calculation of off-diagonal elements of the GF.(0=SCF/DOS; 1=cluster; -1=custom)
      integer, intent(in) :: IGF       !< Do not print or print (0/1) the KKRFLEX_* files
      integer, intent(in) :: LLY       !< LLY <> 0 : apply Lloyds formula
      integer, intent(in) :: KTE       !< Calculation of the total energy On/Off (1/0)
      integer, intent(in) :: KXC       !< Type of xc-potential 0=vBH 1=MJW 2=VWN 3=PW91
      integer, intent(in) :: KPRE
      integer, intent(in) :: LPOT      !< Maximum l component in potential expansion
      integer, intent(in) :: IMIX      !< Type of mixing scheme used (0=straight, 4=Broyden 2nd, 5=Anderson)
      integer, intent(in) :: NCLS      !< Number of reference clusters
      integer, intent(in) :: ICST      !< Number of Born approximation
      integer, intent(in) :: IEND      !< Number of nonzero gaunt coefficients
      integer, intent(in) :: LMAX      !< Maximum l component in wave function expansion
      integer, intent(in) :: NREF      !< Number of diff. ref. potentials
      integer, intent(in) :: NAEZ      !< Number of atoms in unit cell
      integer, intent(in) :: IRID      !< Shape functions parameters in non-spherical part
      integer, intent(in) :: KREL      !< Switch for non-relativistic/relativistic (0/1) program. Attention: several other parameters depend explicitly on KREL, they are set automatically Used for Dirac solver in ASA
      integer, intent(in) :: NCPA      !< NCPA = 0/1 CPA flag
      integer, intent(in) :: LM2D      !< (2*LMAX+1)**2
      integer, intent(in) :: NSRA
      integer, intent(in) :: NPOL      !< Number of Matsubara Poles (EMESHT)
      integer, intent(in) :: NPNT1     !< number of E points (EMESHT) for the contour integration
      integer, intent(in) :: NPNT2     !< number of E points (EMESHT) for the contour integration
      integer, intent(in) :: NPNT3     !< number of E points (EMESHT) for the contour integration
      integer, intent(in) :: ITSCF
      integer, intent(in) :: KVMAD
      integer, intent(in) :: LMPOT     !< (LPOT+1)**2
      integer, intent(in) :: NGSHD     !< Shape functions parameters in non-spherical part
      integer, intent(in) :: MMAXD     !< 2*LMAX+1
      integer, intent(in) :: NPOTD     !< (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP)
      integer, intent(in) :: NATYP     !< Number of kinds of atoms in unit cell
      integer, intent(in) :: IPAND     !< Number of panels in non-spherical part
      integer, intent(in) :: NCLSD     !< Maximum number of different TB-clusters
      integer, intent(in) :: NCLEB     !< Number of Clebsch-Gordon coefficients
      integer, intent(in) :: NFUND     !< Shape functions parameters in non-spherical part
      integer, intent(in) :: IEMXD     !< Dimension for energy-dependent arrays
      integer, intent(in) :: NINEQ     !< Number of ineq. positions in unit cell
      integer, intent(in) :: NSPIN     !< Counter for spin directions
      integer, intent(in) :: KMROT     !< 0: no rotation of the magnetisation; 1: individual rotation of the magnetisation for every site
      integer, intent(in) :: NTOTD
      integer, intent(in) :: NCHEB     !< Number of Chebychev pannels for the new solver
      integer, intent(in) :: NLEFT     !< Number of repeated basis for left host to get converged electrostatic potentials
      integer, intent(in) :: NRIGHT    !< Number of repeated basis for right host to get converged electrostatic potentials
      integer, intent(in) :: ITDBRY    !< Number of SCF steps to remember for the Broyden mixing
      integer, intent(in) :: KSHAPE    !< Exact treatment of WS cell
      integer, intent(in) :: ISHIFT
      integer, intent(in) :: IELAST
      integer, intent(in) :: KFORCE    !< Calculation of the forces
      integer, intent(in) :: NCELLD    !< Number of cells (shapes) in non-spherical part
      integer, intent(in) :: LMXSPD    !< (2*LPOT+1)**2
      integer, intent(in) :: IRMIND    !< IRM-IRNSD
      integer, intent(in) :: NSPOTD    !< Number of potentials for storing non-sph. potentials
      integer, intent(in) :: NEMBD1    !< NEMB+1
      integer, intent(in) :: LMMAXD    !< (KREL+KORBIT+1)(LMAX+1)^2
      integer, intent(in) :: NEMBD2
      integer, intent(in) :: NACLSD    !< Maximum number of atoms in a TB-cluster
      integer, intent(in) :: LMAXD1
      integer, intent(in) :: NOFGIJ    !< number of GF pairs IJ to be calculated as determined from IJTABCALC<>0
      integer, intent(in) :: NSPIND    !< KREL+(1-KREL)*(NSPIN+1)
      integer, intent(in) :: NSHELD    !< Number of blocks of the GF matrix that need to be calculated (NATYP + off-diagonals in case of impurity)
      integer, intent(in) :: NSYMAT
      integer, intent(in) :: INVMOD    !< Inversion scheme
      integer, intent(in) :: NTLDAU    !< number of atoms on which LDA+U is applied
      integer, intent(in) :: NQCALC
      integer, intent(in) :: N1SEMI    !< Number of energy points for the semicore contour
      integer, intent(in) :: N2SEMI    !< Number of energy points for the semicore contour
      integer, intent(in) :: N3SEMI    !< Number of energy points for the semicore contour
      integer, intent(in) :: NPAN_EQ   !< Variables for the pannels for the new solver
      integer, intent(in) :: NSYMAXD
      integer, intent(in) :: NAEZDPD
      integer, intent(in) :: NLBASIS      !< Number of basis layers of left host (repeated units)
      integer, intent(in) :: NRBASIS      !< Number of basis layers of right host (repeated units)
      integer, intent(in) :: NSPINDD      !< NSPIND-KORBIT
      integer, intent(in) :: MAXMESH
      integer, intent(in) :: INTERVX      !< Number of intervals in x-direction for k-net in IB of the BZ
      integer, intent(in) :: INTERVY      !< Number of intervals in y-direction for k-net in IB of the BZ
      integer, intent(in) :: INTERVZ      !< Number of intervals in z-direction for k-net in IB of the BZ
      integer, intent(in) :: IDOLDAU      !< flag to perform LDA+U
      integer, intent(in) :: NPOLSEMI     !< Number of poles for the semicore contour
      integer, intent(inout) :: SCFSTEPS  !< number of scf iterations
      integer, intent(in) :: ITCPAMAX     !< Max. number of CPA iterations
      integer, intent(in) :: NATOMIMP     !< Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
      integer, intent(in) :: NPAN_LOG     !< Variables for the pannels for the new solver
      integer, intent(in) :: NATOMIMPD    !< Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
      integer, intent(in) :: ITRUNLDAU    !< Iteration index for LDA+U
      integer, intent(in) :: IESEMICORE
      !     .. nembd2 = NAEZ+NEMB, lmaxd1=lmaxd+1, naezdpd=NAEZ/nprincd)
      real (kind=dp), intent(in) :: TK        !< Temperature
      real (kind=dp), intent(in) :: FCM
      real (kind=dp), intent(in) :: ALAT      !< Lattice constant in a.u.
      real (kind=dp), intent(inout) :: EMIN   !< Energies needed in EMESHT
      real (kind=dp), intent(in) :: EMAX      !< Energies needed in EMESHT
      real (kind=dp), intent(in) :: R_LOG
      real (kind=dp), intent(in) :: EFERMI    !< Fermi energy
      real (kind=dp), intent(in) :: CPATOL    !< Convergency tolerance for CPA-cycle
      real (kind=dp), intent(in) :: MIXING    !< Magnitude of the mixing parameter
      real (kind=dp), intent(in) :: QBOUND    !< Convergence parameter for the potential
      real (kind=dp), intent(in) :: TKSEMI    !< Temperature of semi-core contour
      real (kind=dp), intent(in) :: EMUSEMI
      real (kind=dp), intent(in) :: TOLRDIF   !< Tolerance for r<tolrdif (a.u.) to handle vir. atoms
      real (kind=dp), intent(in) :: EBOTSEMI
      real (kind=dp), intent(in) :: FSEMICORE
      real (kind=dp), intent(in) :: LAMBDA_XC !< Scale magnetic moment (0 < Lambda_XC < 1, 0=zero moment, 1= full moment)
      logical, intent(in) :: LRHOSYM
      logical, intent(in) :: LINTERFACE         !< If True a matching with semi-inifinite surfaces must be performed
      character(len=10), intent(in) :: SOLVER   !< Type of solver
      complex (kind=dp), intent(in) :: DELTAE      !< Energy difference for numerical derivative
      !     ..
      !     .. Array arguments
      complex (kind=dp), dimension(IEMXD), intent(in) :: EZ
      complex (kind=dp), dimension(IEMXD), intent(in) :: WEZ
      complex (kind=dp), dimension(LMMAXD,LMMAXD), intent(in)  :: RC           !< NREL REAL spher. harm. > CMPLX. spher. harm. NREL CMPLX. spher. harm. > REAL spher. harm.
      complex (kind=dp), dimension(LMMAXD,LMMAXD), intent(in)  :: CREL         !< Non-relat. CMPLX. spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. CMPLX. spher. harm.
      complex (kind=dp), dimension(LMMAXD,LMMAXD), intent(in)  :: RREL         !< Non-relat. REAL spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. REAL spher. harm.
      complex (kind=dp), dimension(IRM,NATYP), intent(in)      :: PHILDAU

      complex (kind=dp), dimension(LMMAXD,LMMAXD,NAEZ), intent(in)      :: DROTQ   !< Rotation matrices to change between LOCAL/GLOBAL frame of reference for magnetisation <> Oz or noncollinearity
      complex (kind=dp), dimension(LMMAXD,LMMAXD,NSYMAXD), intent(in)   :: DSYMLL
      complex (kind=dp), dimension(2,2,LMMAXD), intent(in)              :: SRREL
      complex (kind=dp), dimension(LMMAXD,LMMAXD,NEMBD1,NSPINDD,IEMXD), intent(in) :: LEFTTINVLL
      complex (kind=dp), dimension(LMMAXD,LMMAXD,NEMBD1,NSPINDD,IEMXD), intent(in) :: RIGHTTINVLL
      real (kind=dp), dimension(NATYP), intent(in)  :: A        !< Constants for exponential R mesh
      real (kind=dp), dimension(NATYP), intent(in)  :: B        !< Constants for exponential R mesh
      real (kind=dp), dimension(2), intent(in)      :: VBC      !< Potential constants
      real (kind=dp), dimension(NATYP), intent(in)  :: ZAT      !< Nuclear charge
      real (kind=dp), dimension(NATYP), intent(in)  :: RMT      !< Muffin-tin radius of true system
      real (kind=dp), dimension(NATYP), intent(in)  :: RWS      !< Wigner Seitz radius
      real (kind=dp), dimension(NGSHD), intent(in)  :: GSH
      real (kind=dp), dimension(NATYP), intent(in)  :: CONC     !< Concentration of a given atom
      real (kind=dp), dimension(NREF), intent(in)   :: VREF
      real (kind=dp), dimension(NATYP), intent(in)  :: UEFF     !< input U parameter for each atom
      real (kind=dp), dimension(NATYP), intent(in)  :: JEFF     !< input J parameter for each atom
      real (kind=dp), dimension(NAEZ), intent(in)   :: QMTET    !< \f$ \theta\f$ angle of the agnetization with respect to the z-axis
      real (kind=dp), dimension(NAEZ), intent(in)   :: QMPHI    !< \f$ \phi\f$ angle of the agnetization with respect to the z-axis
      real (kind=dp), dimension(NREF), intent(in)   :: RMTREF   !< Muffin-tin radius of reference system
      real (kind=dp), dimension(NATYP), intent(in)  :: RMTNEW   !< Adapted muffin-tin radius
      real (kind=dp), dimension(NATYP), intent(in)  :: EREFLDAU !< the energies of the projector's wave functions (REAL)
      real (kind=dp), dimension(NATYP), intent(in)  :: SOCSCALE !< Spin-orbit scaling
      real (kind=dp), dimension(IRM,NATYP), intent(in)             :: R        !< Radial mesh ( in units a Bohr)
      real (kind=dp), dimension(3,0:NR), intent(in)                :: RR
      real (kind=dp), dimension(IRM,NATYP), intent(in)             :: DRDI     !< Derivative dr/di
      real (kind=dp), dimension(NCLEB,2), intent(in)               :: CLEB     !< GAUNT coefficients (GAUNT)
      !real (kind=dp), dimension(LMAXD1,NATYP), intent(in)          :: CSCL     !< Speed of light scaling
      real (kind=dp), dimension(KREL*LMAX+1,KREL*NATYP+(1-KREL)), intent(inout) :: CSCL      !< Speed of light scaling
      real (kind=dp), dimension(NTOTD*(NCHEB+1),NATYP), intent(in) :: RNEW
      real (kind=dp), dimension(IRM,NPOTD), intent(in)             :: VISP     !< Spherical part of the potential
      !real (kind=dp), dimension(IRM,NATYP), intent(in)             :: VTREL    !< potential (spherical part)
      real (kind=dp), dimension(IRM*KREL+(1-KREL),NATYP), intent(inout)         :: VTREL     !< potential (spherical part)
      !real (kind=dp), dimension(IRM,NATYP), intent(in)             :: BTREL    !< magnetic field
      real (kind=dp), dimension(IRM*KREL+(1-KREL),NATYP), intent(inout)         :: BTREL     !< magnetic field
      !real (kind=dp), dimension(IRM,NATYP), intent(in)             :: RMREL    !< radial mesh
      real (kind=dp), dimension(IRM*KREL+(1-KREL),NATYP), intent(inout)         :: RMREL     !< radial mesh
      real (kind=dp), dimension(3,NSHELD), intent(in)              :: RATOM
      real (kind=dp), dimension(20,NPOTD), intent(in)              :: ECORE    !< Core energies
      real (kind=dp), dimension(3,NEMBD2), intent(in)              :: RBASIS   !< Position of atoms in the unit cell in units of bravais vectors
      !real (kind=dp), dimension(LMAXD1,NATYP), intent(in)          :: SOCSCL
      real (kind=dp), dimension(KREL*LMAX+1,KREL*NATYP+(1-KREL)), intent(inout) :: SOCSCL
      !real (kind=dp), dimension(IRM,NATYP), intent(in)             :: DRDIREL  !< derivative of radial mesh
      real (kind=dp), dimension(IRM*KREL+(1-KREL),NATYP), intent(in) :: DRDIREL
      real (kind=dp), dimension(NAEZ,3), intent(in)                :: QMPHITAB
      real (kind=dp), dimension(NAEZ,3), intent(in)                :: QMTETTAB
      real (kind=dp), dimension(NAEZ,3), intent(in)                :: QMGAMTAB
      real (kind=dp), dimension(3,NATOMIMPD), intent(in)           :: RCLSIMP
      real (kind=dp), dimension(LMPOT,NEMBD1), intent(in)          :: CMOMHOST !< Charge moments of each atom of the (left/right) host
      !real (kind=dp), dimension(IRM,NATYP), intent(in)             :: R2DRDIREL   !< \f$ r^2 \frac{\partial}{\partial \mathbf{r}}\frac{\partial}{\partial i}\f$ (r**2 * drdi)
      real (kind=dp), dimension(IRM*KREL+(1-KREL),NATYP), intent(in) :: R2DRDIREL
      real (kind=dp), dimension(0:NTOTD,NATYP), intent(in)         :: RPAN_INTERVALL
      real (kind=dp), dimension(48,3,NSHELD), intent(in)              :: RROT
      real (kind=dp), dimension(3,NACLSD,NCLSD), intent(in)           :: RCLS   !< Real space position of atom in cluster
      real (kind=dp), dimension(IRMIND:IRM,LMPOT,NSPOTD), intent(in)  :: VINS   !< Non-spherical part of the potential
      real (kind=dp), dimension(IRID,NFUND,NCELLD), intent(in)        :: THETAS !< shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
      real (kind=dp), dimension(NTOTD*(NCHEB+1),NFUND,NCELLD), intent(in)   :: THETASNEW
      real (kind=dp), dimension(MMAXD,MMAXD,NSPIND,NATYP), intent(in)       :: WLDAU  !< potential matrix
      real (kind=dp), dimension(MMAXD,MMAXD,MMAXD,MMAXD,NATYP), intent(in)  :: ULDAU  !< calculated Coulomb matrix elements (EREFLDAU)
      !     ..
      integer, dimension(NAEZ), intent(in)      :: NOQ       !< Number of diff. atom types located
      integer, dimension(NATYP), intent(in)     :: IMT       !< R point at MT radius
      integer, dimension(NATYP), intent(in)     :: IRC       !< R point for potential cutting
      integer, dimension(NATYP), intent(in)     :: NFU
      integer, dimension(NEMBD2), intent(in)    :: CLS       !< Cluster around atomic sites
      integer, dimension(NAEZ), intent(in)      :: ICPA      !< ICPA = 0/1 site-dependent CPA flag
      integer, dimension(NATYP), intent(in)     :: IPAN      !< Number of panels in non-MT-region
      integer, dimension(NATYP), intent(in)     :: ZREL      !< atomic number (cast integer)
      integer, dimension(NATYP), intent(in)     :: IQAT      !< The site on which an atom is located on a given site
      integer, dimension(NATYP), intent(in)     :: LOPT      !< angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
      integer, dimension(NATYP), intent(in)     :: IRNS      !< Position of atoms in the unit cell in units of bravais vectors
      integer, dimension(NSHELD), intent(in)    :: NSH1      !< Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
      integer, dimension(NSHELD), intent(in)    :: NSH2      !< Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
      integer, dimension(NATYP), intent(in)     :: IRWS      !< R point at WS radius
      integer, dimension(LM2D), intent(in)      :: LOFLM     !< l of lm=(l,m) (GAUNT)
      integer, dimension(NCLSD), intent(in)     :: NACLS     !< Number of atoms in cluster
      integer, dimension(IEMXD), intent(in)     :: KMESH
      integer, dimension(NPOTD), intent(in)     :: NCORE     !< Number of core states
      integer, dimension(NAEZ), intent(in)      :: IQCALC
      integer, dimension(NATYP), intent(in)     :: IRMIN     !< Max R for spherical treatment
      integer, dimension(NATYP), intent(in)     :: NTCELL    !< Index for WS cell
      integer, dimension(NATYP), intent(in)     :: IXIPOL    !< Constraint of spin pol.
      integer, dimension(NATYP), intent(in)     :: JWSREL    !< index of the WS radius
      integer, dimension(NATYP), intent(in)     :: ITLDAU    !< integer pointer connecting the NTLDAU atoms to heir corresponding index in the unit cel
      integer, dimension(NEMBD2), intent(in)    :: REFPOT    !< Ref. pot. card  at position
      integer, dimension(0:LMPOT), intent(in)   :: IMAXSH
      integer, dimension(0:NSHELD), intent(in)  :: NSHELL    !< Index of atoms/pairs per shell (ij-pairs); nshell(0) = number of shells
      integer, dimension(0:NATYP), intent(in)   :: HOSTIMP
      integer, dimension(NATYP), intent(in)     :: IRSHIFT   !< shift of the REL radial mesh with respect no NREL
      integer, dimension(NOFGIJ), intent(in)    :: IJTABSH   !< Linear pointer, assigns pair (i,j) to a shell in the array GS(*,*,*,NSHELD)
      integer, dimension(NATOMIMPD), intent(in) :: ATOMIMP
      integer, dimension(NATYP), intent(in)     :: NPAN_TOT
      integer, dimension(NOFGIJ), intent(in)    :: IJTABSYM  !< Linear pointer, assigns pair (i,j) to the rotation bringing GS into Gij
      integer, dimension(NOFGIJ), intent(in)    :: IJTABCALC !< Linear pointer, specifying whether the block (i,j) has to be calculated needs set up for ICC=-1, not used for ICC=1
      integer, dimension(NATYP), intent(in)     :: NPAN_EQ_AT
      integer, dimension(NATYP), intent(in)     :: NPAN_LOG_AT
      integer, dimension(NOFGIJ), intent(in)    :: IJTABCALC_I
      integer, dimension(NGSHD,3), intent(in)         :: ILM_MAP
      integer, dimension(NSHELD,NOFGIJ), intent(in)   :: ISH
      integer, dimension(NSHELD,NOFGIJ), intent(in)   :: JSH
      integer, dimension(NACLSD,NEMBD2), intent(in)   :: ATOM   !< Atom at site in cluster
      integer, dimension(NACLSD,NEMBD2), intent(in)   :: EZOA   !< EZ of atom at site in cluster
      integer, dimension(NATYP,LMXSPD), intent(in)    :: LMSP   !< 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
      integer, dimension(NATYP,NFUND), intent(in)     :: LLMSP  !< lm=(l,m) of 'nfund'th nonvanishing component of non-spherical pot.
      integer, dimension(LMXSPD,NATYP), intent(in)    :: LMSP1
      integer, dimension(NATYP,LMXSPD), intent(in)    :: IFUNM
      integer, dimension(NATYP,NEMBD2), intent(in)    :: KAOEZ  !< Kind of atom at site in elem. cell
      integer, dimension(0:IPAND,NATYP), intent(in)   :: IRCUT  !< R points of panel borders
      integer, dimension(NCLEB,4), intent(in)         :: ICLEB  !< Pointer array
      integer, dimension(20,NPOTD), intent(in)        :: LCORE  !< Angular momentum of core states
      integer, dimension(2,LMMAXD), intent(in)        :: NRREL
      integer, dimension(NAEZDPD,NAEZDPD), intent(in) :: ICHECK
      integer, dimension(LMXSPD,NATYP), intent(in)    :: IFUNM1
      integer, dimension(20,NPOTD), intent(in)        :: ITITLE
      integer, dimension(0:NTOTD,NATYP), intent(in)   :: IPAN_INTERVALL
      integer, dimension(LMPOT,0:LMAX,0:LMAX), intent(in)   :: JEND    !< Pointer array for icleb()
      integer, dimension(2,2,LMMAXD), intent(in)            :: IRREL
      logical, dimension(2), intent(in)         :: VACFLAG
      logical, dimension(NSYMAXD), intent(in)   :: SYMUNITARY          !< unitary/antiunitary symmetry flag
      character(len=124), dimension(6), intent(in) :: TXC
      ! .. Local scalars
      integer :: I1
      integer :: IC, NACLSMIN, NACLSMAX, NQDOS, IRMDNEW  ! variables for t_inc filling
      integer :: ITMPDIR,ILTMP
      real (kind=dp), dimension(NATYP) :: PHI
      real (kind=dp), dimension(NATYP) :: THETA
      character(len=80) :: TMPDIR
      ! .. External Functions
      LOGICAL OPT,TEST
      EXTERNAL OPT,TEST
      !
      ITMPDIR=0
      ILTMP=0

      ! put information about scf steps into t_inc
      t_inc%I_ITERATION = ITSCF
      t_inc%N_ITERATION = SCFSTEPS
      ! put information for save_wavefun also in:
      t_inc%NSRA     =  NSRA
      t_inc%LMMAXSO  =  2*(LMAX+1)**2
      IRMDNEW = 0
      do I1=1,NATYP
         if (NPAN_TOT(I1)*(NCHEB+1)>IRMDNEW) then
            IRMDNEW = NPAN_TOT(I1)*(NCHEB+1)
         end if
      end do
      t_inc%IRMDNEW = IRMDNEW

      !-------------------------------------------------------------------------
      ! itermdir
      !-------------------------------------------------------------------------
      IF (OPT('ITERMDIR')) THEN
         I1 = 0
         EMIN = 0D0
         t_params%QMTET    = QMTET
         t_params%QMPHI    = QMPHI
         t_params%QMPHITAB = QMPHITAB
         t_params%QMTETTAB = QMTETTAB
         t_params%QMGAMTAB = QMGAMTAB
         t_params%I1       = I1
         t_params%EMIN     = EMIN
         t_params%DROTQ    = DROTQ
      END IF
      !-------------------------------------------------------------------------
      ! lda+u data in this file change
      !-------------------------------------------------------------------------
      IF ( IDOLDAU.EQ.1 ) THEN
         OPEN (67,FILE='ldau.unformatted',FORM='unformatted')
         WRITE (67) ITRUNLDAU,WLDAU,ULDAU,PHILDAU
         CLOSE(67)
      END IF
      !-------------------------------------------------------------------------
      ! nonco_angle file
      !-------------------------------------------------------------------------
      call read_angles(t_params,NATYP,THETA,PHI)
      !-------------------------------------------------------------------------
      ! fill t_inc type meant for simpler passing of basic parameter to other routines
      !-------------------------------------------------------------------------

      ! find maximal cluster information (taken from main1b)
      NACLSMAX = 1
      DO IC = 1,NCLS
         IF (NACLS(IC).GT.NACLSMAX) NACLSMAX = NACLS(IC)
      ENDDO

      ! find NQDOS (taken from main1b)
      IF (OPT('qdos    ')) THEN
         OPEN(67,FILE='qvec.dat')
         READ(67,*) NQDOS
         CLOSE(67)
      ELSE
         NQDOS = 1
      END IF

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! t_inc t_inc t_inc t_inc t_inc t_inc t_inc t_inc t_inc t_inc
      !fill t_inc
      t_inc%LMMAXD   = LMMAXD
      t_inc%NSPIN    = NSPIN
      t_inc%IELAST   = IELAST
      t_inc%NQDOS    = NQDOS
      t_inc%NATYP    = NATYP
      t_inc%LMGF0D   = (LMAX+1)**2  ! see main1b
      t_inc%NCLSD    = NCLS
      t_inc%NACLSMAX   = NACLSMAX
      t_inc%NSHELL0  = NSHELL(0)
      IF(OPT('NEWSOSOL')) t_inc%NEWSOSOL = .true.
      IF(OPT('deci-out')) t_inc%deci_out = .true.
      !t_inc t_inc t_inc t_inc t_inc t_inc t_inc t_inc t_inc t_inc
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! writeout flags writeout flags writeout flags writeout flags writeout flags writeout flags writeout flags writeout flags
      !set logical switches in t_tgmat which control if tmat, gmat and gref are written to files or stored in memory
      if(TEST('tmatfile')) t_tgmat%tmat_to_file = .true.
      if(TEST('gmatfile')) t_tgmat%gmat_to_file = .true.
      if(TEST('greffile')) t_tgmat%gref_to_file = .true.
      if(TEST('projfile')) t_cpa%dmatproj_to_file = .true.


      !bug bug bug bug bug
      ! in case of ASA DIRAC solver (KREL==1) then gmat file has to be written out otherwise something is going wrong.
      if(KREL>0) t_tgmat%gmat_to_file = .true.
      !bug bug bug bug bug

      !some special run options:
      if(OPT('KKRFLEX '))  t_tgmat%tmat_to_file = .true.  ! for KKRFLEX option tmat must be written to file
      if(OPT('qdos    '))  t_tgmat%gmat_to_file = .true.  ! for qdos write gmat to file since it scales with NQDOS and can become huge

      !set logical switches in t_lloyd which control if files are written to files or stored in memory
      !       if(TEST('tmatfile').or.TEST('llyfiles'))
      if(TEST('wrtdtmat').or.TEST('llyfiles')) t_lloyd%dtmat_to_file = .true.
      if(TEST('wrttral ').or.TEST('llyfiles')) t_lloyd%tralpha_to_file = .true.
      if(TEST('wrtcdos ').or.TEST('llyfiles')) t_lloyd%cdos_diff_lly_to_file = .true.
      if(TEST('wrtdgref').or.TEST('llyfiles')) t_lloyd%dgref_to_file = .true.
      if(TEST('wrtgotr ').or.TEST('llyfiles')) t_lloyd%g0tr_to_file = .true.

      ! set verbosity level in t_inc%i_write = 0,1,2 for default, verbose1, verbose2
      t_inc%i_write = 0 !default: write only output.000.txt and reset file after each iteration
      if(TEST('verbose1')) t_inc%i_write = 1 !write on all processors but only the latest iteration
      if(TEST('verbose2')) t_inc%i_write = 2 !write everything
      ! and t_inc_i_time for timing writeout
      t_inc%i_time = 1  !default: only timings from master, all iterations
      if(TEST('timings0')) t_inc%i_time = 0  !only timings from master, only the last iteration
      if(TEST('timings2')) t_inc%i_time = 2  !all timing files, all iterations
      ! writeout flags writeout flags writeout flags writeout flags writeout flags writeout flags writeout flags writeout flags
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! MPI communication scheme
      !set switch for MPIatom test option (see mod_types and mod_mympi)
      ! default values for MPIadapt and MPIatom
      if(NATYP<=IELAST) then
         MPIatom = .true.
      else
         MPIatom = .false.
      end if
      MPIadapt = 1 ! 1 means check timings and then reshuffle ranks if necessary
      ! change default behaviour for the corresponding test flags are found
      if(test('MPIatom ')) then
         MPIatom = .true.
         MPIadapt= 0 ! 0 means no change in communication scheme
      end if
      if(test('MPIenerg')) then
         MPIatom = .false.
         MPIadapt= 0
      end if
      if(test('MPIadapt')) then
         MPIadapt= 2 ! 2 means force run with MPIatom then with MPIenerg and then compare to choose optimal
      end if
      ! so far changing does not work yet, so turn this off:
      MPIadapt = 0

      if(opt('FERMIOUT') .or. OPT('WRTGREEN') .or. OPT('GREENIMP') .or. opt('OPERATOR') ) then ! fswrt
         MPIatom = .true.                                               ! fswrt
         if(SCFSTEPS>1) then                                            ! fswrt
            write(*,*) 'Warning: Setting SCFSTEPS=1 for FERMIOUT option'! fswrt
            SCFSTEPS = 1                                                ! fswrt
         end if                                                         ! fswrt
         if (nranks>NATYP) then                                         ! fswrt
            write(*,*) 'FERMIOUT/WRTGREEN/GREENIMP option chosen'       ! fswrt
            write(*,*) 'Nranks>NATYP', nranks, NATYP                    ! fswrt
            stop 'Please choose Nranks<=NATYP'                          ! fswrt
         end if                                                         ! fswrt
         naclsmin = minval(NACLS(1:NCLS))                                ! fswrt
         if(.not.OPT('GREENIMP') .and. naclsmin<150) then               ! fswrt
            write(*,*) ' !!!  WARNING  !!!'                               ! fswrt
            write(*,*) '   FERMIOUT/WRTGREEN option chosen'                ! fswrt
            write(*,*) '   minimal cluster size smaller than 150 atoms!!!' ! fswrt
            write(*,*) '   should be increased to least 200-300 atoms'     ! fswrt
         end if                                                         ! fswrt
         if(test('MPIenerg')) then                                      ! fswrt
            write(*,*) 'FERMIOUT/WRTGREEN/GREENIMP option chosen'       ! fswrt
            write(*,*) 'found unsupported test option MPIenerg'         ! fswrt
            stop 'Please choose MPIatom instead'                        ! fswrt
         end if                                                         ! fswrt
         if(opt('OPERATOR') .and. opt('FERMIOUT')) then                 ! fswrt
            write(*,*) 'OPERATOR and FERMIOUT cannot be used together'  ! fswrt
            stop 'Please chose only one of the two'                     ! fswrt
         end if                                                         ! fswrt
         if(opt('OPERATOR') .and. ielast.ne.1) then                     ! fswrt
            write(*,*) 'OPERATOR option chosen'                         ! fswrt
            write(*,*) 'energy contour should contain a single point'   ! fswrt
            write(*,*) 'on the real axis only'                          ! fswrt
            stop 'Please correct energy contour'                        ! fswrt
         end if                                                         ! fswrt
      end if                                                            ! fswrt
      if(.not.OPT('OPERATOR') .and. TEST('IMP_ONLY')) then              ! fswrt
         write(*,*) 'test option "IMP_ONLY" can only be used'           ! fswrt
         write(*,*) 'in combination with option "OPERATOR"'             ! fswrt
         stop                                                           ! fswrt
      end if                                                            ! fswrt

      if(test('MPIatom ') .and. test('MPIenerg')) then
         stop "[wunfiles] Found test options 'MPIenerg' and 'MPIatom' which do not work together. Please choose only one of these."
      end if
      ! MPI communication scheme
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      ! all parameters are stored in t_params fomr mod_wunfiles
      ! first fill scalar values
      call fill_t_params_scalars(IEMXD,IRMIND,IRM,LMPOT,NSPOTD,NPOTD,NATYP,   &
         NEMBD1,LMMAXD,NAEZ,IPAND,NEMBD2,NREF,LMAX,NCLEB,NACLSD,NCLSD,LM2D,   &
         LMAXD1,NR,NSHELD,NSYMAXD,NAEZDPD,NATOMIMPD,NOFGIJ,NSPIND,NSPINDD,    &
         IRID,NFUND,NCELLD,LMXSPD,NGSHD,KREL,MMAXD,IELAST,NPOL,NPNT1,NPNT2,   &
         NPNT3,ITSCF,SCFSTEPS,LLY,NSRA,INS,NINEQ,NSPIN,NCLS,ICST,IEND,ICC,IGF,&
         NLBASIS,NRBASIS,NCPA,ITCPAMAX,KMROT,MAXMESH,NSYMAT,NATOMIMP,INVMOD,  &
         NQCALC,INTERVX,INTERVY,INTERVZ,LPOT,NRIGHT,NLEFT,IMIX,ITDBRY,KPRE,   &
         KSHAPE,KTE,KVMAD,KXC,ISHIFT,KFORCE,IDOLDAU,ITRUNLDAU,NTLDAU,NPOLSEMI,&
         N1SEMI,N2SEMI,N3SEMI,IESEMICORE,EBOTSEMI,EMUSEMI,TKSEMI,FSEMICORE,   &
         R_LOG,EMIN,EMAX,TK,EFERMI,ALAT,CPATOL,MIXING,QBOUND,FCM,LAMBDA_XC,   &
         TOLRDIF,LINTERFACE,LRHOSYM,SOLVER,TMPDIR,ITMPDIR,ILTMP,NTOTD,NCHEB,  &
         DELTAE,t_params)

      ! initialize allocatable arrays
      call init_t_params(t_params)

      ! now fill arrays that have just been allocated
      call fill_t_params_arrays(t_params,IEMXD,LMMAXD,NAEZ,NSYMAXD,NEMBD1,    &
         NSPINDD,IRMIND,IRM,LMPOT,NSPOTD,NPOTD,NATYP,NR,NEMBD2,NREF,NCLEB,    &
         NCLSD,NACLSD,NSHELD,NGSHD,NFUND,IRID,NCELLD,MMAXD,LM2D,LMXSPD,LMAXD1,&
         NSPIND,NTOTD,NCHEB,IPAND,LMAX,NOFGIJ,NAEZDPD,NATOMIMPD,EZ,WEZ,DROTQ, &
         DSYMLL,LEFTTINVLL,RIGHTTINVLL,CREL,RC,RREL,SRREL,PHILDAU,VINS,VISP,  &
         VBC,VTREL,BTREL,SOCSCALE,DRDIREL,R2DRDIREL,RMREL,CMOMHOST,ECORE,     &
         QMTET,QMPHI,QMPHITAB,QMTETTAB,QMGAMTAB,ZAT,R,DRDI,RMTREF,VREF,CLEB,  &
         RCLS,SOCSCL,CSCL,RBASIS,RR,CONC,RROT,RATOM,A,B,THETAS,RMT,RMTNEW,RWS,&
         GSH,EREFLDAU,UEFF,JEFF,ULDAU,WLDAU,RPAN_INTERVALL,RNEW,THETASNEW,    &
         LOPT,ITLDAU,IRSHIFT,JWSREL,ZREL,LCORE,NCORE,IPAN,IRCUT,JEND,ICLEB,   &
         ATOM,CLS,NACLS,LOFLM,EZOA,KAOEZ,IQAT,ICPA,NOQ,KMESH,NSHELL,NSH1,NSH2,&
         IJTABCALC,IJTABCALC_I,IJTABSYM,IJTABSH,ISH,JSH,IQCALC,ICHECK,ATOMIMP,&
         REFPOT,IRREL,NRREL,IFUNM1,ITITLE,LMSP1,NTCELL,IXIPOL,IRNS,IFUNM,     &
         LLMSP,LMSP,IMT,IRC,IRMIN,IRWS,NFU,HOSTIMP,ILM_MAP,IMAXSH,NPAN_LOG,       &
         NPAN_EQ,NPAN_LOG_AT,NPAN_EQ_AT,NPAN_TOT,IPAN_INTERVALL,SYMUNITARY,   &
         VACFLAG,TXC,RCLSIMP, KREL)

      ! save information about the energy mesh
      call save_emesh(IELAST,EZ,WEZ,EMIN,EMAX,IESEMICORE,FSEMICORE,NPOL,TK,   &
         NPNT1,NPNT2,NPNT3,EBOTSEMI,EMUSEMI,TKSEMI,NPOLSEMI,N1SEMI,N2SEMI,    &
         N3SEMI,IEMXD,t_params)

   end subroutine WUNFILES


   !----------------------------------------------------------------------------
   ! SUBROUTINE: init_t_params
   !> @brief Allocate initial parameters to be broadcasted via mpi
   !> @author Philipp Rüssmann
   !----------------------------------------------------------------------------
   subroutine init_t_params(t_params)
      ! allocate arrays, has to be done after bcast t_params_scalars for myrank<>master
      ! otherwise are the parameters not set
      implicit none

      type(type_params), intent(inout) :: t_params

      integer :: i_stat

      !-------------------------------------------------------------------------
      ! Allocate complex (kind=dp) arrays
      !-------------------------------------------------------------------------
      allocate(t_params%EZ(t_params%IEMXD),stat=i_stat)   !complex (kind=dp)
      call memocc(i_stat,product(shape(t_params%EZ))*kind(t_params%EZ),'t_params%EZ','init_t_params')
      allocate(t_params%WEZ(t_params%IEMXD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%WEZ))*kind(t_params%WEZ),'t_params%WEZ','init_t_params')
      allocate(t_params%DROTQ(t_params%LMMAXD,t_params%LMMAXD,t_params%NAEZ),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%DROTQ))*kind(t_params%DROTQ),'t_params%DROTQ','init_t_params')
      allocate(t_params%DSYMLL(t_params%LMMAXD,t_params%LMMAXD,t_params%NSYMAXD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%DSYMLL))*kind(t_params%DSYMLL),'t_params%DSYMLL','init_t_params')
      allocate(t_params%LEFTTINVLL(t_params%LMMAXD,t_params%LMMAXD,&
         t_params%NEMBD1,t_params%NSPINDD,t_params%IEMXD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%LEFTTINVLL))*kind(t_params%LEFTTINVLL),'t_params%LEFTTINVLL','init_t_params')
      allocate(t_params%RIGHTTINVLL(t_params%LMMAXD,t_params%LMMAXD,&
         t_params%NEMBD1,t_params%NSPINDD,t_params%IEMXD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%RIGHTTINVLL))*kind(t_params%RIGHTTINVLL),'t_params%RIGHTTINVLL','init_t_params')
      allocate(t_params%CREL(t_params%LMMAXD,t_params%LMMAXD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%CREL))*kind(t_params%CREL),'t_params%CREL','init_t_params')
      allocate(t_params%RC(t_params%LMMAXD,t_params%LMMAXD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%RC))*kind(t_params%RC),'t_params%RC','init_t_params')
      allocate(t_params%RREL(t_params%LMMAXD,t_params%LMMAXD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%RREL))*kind(t_params%RREL),'t_params%RREL','init_t_params')
      allocate(t_params%SRREL(2,2,t_params%LMMAXD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%SRREL))*kind(t_params%SRREL),'t_params%SRREL','init_t_params')
      allocate(t_params%PHILDAU(t_params%IRM,t_params%NATYP), stat=i_stat)
      call memocc(i_stat,product(shape(t_params%PHILDAU))*kind(t_params%PHILDAU),'t_params%PHILDAU','init_t_params')
      !-------------------------------------------------------------------------
      ! End of allocation of complex (kind=dp) arrays
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      ! Allocate real (kind=dp) arrays
      !-------------------------------------------------------------------------
      allocate(t_params%VINS(t_params%IRMIND:t_params%IRM,t_params%LMPOT,t_params%NSPOTD),stat=i_stat)!real (kind=dp)
      call memocc(i_stat,product(shape(t_params%VINS))*kind(t_params%VINS),'t_params%VINS','init_t_params')
      allocate(t_params%VISP(t_params%IRM,t_params%NPOTD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%VISP))*kind(t_params%VISP),'t_params%VISP','init_t_params')
      allocate(t_params%VBC(2),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%VBC))*kind(t_params%VBC),'t_params%VBC','init_t_params')
      allocate(t_params%VTREL(t_params%IRM,t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%VTREL))*kind(t_params%VTREL),'t_params%VTREL','init_t_params')
      allocate(t_params%BTREL(t_params%IRM,t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%BTREL))*kind(t_params%BTREL),'t_params%BTREL','init_t_params')
      allocate(t_params%SOCSCALE(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%SOCSCALE))*kind(t_params%SOCSCALE),'t_params%SOCSCALE','init_t_params')
      allocate(t_params%DRDIREL(t_params%IRM,t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%DRDIREL))*kind(t_params%DRDIREL),'t_params%DRDIREL','init_t_params')
      allocate(t_params%R2DRDIREL(t_params%IRM,t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%R2DRDIREL))*kind(t_params%R2DRDIREL),'t_params%R2DRDIREL','init_t_params')
      allocate(t_params%RMREL(t_params%IRM,t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%RMREL))*kind(t_params%RMREL),'t_params%RMREL','init_t_params')
      allocate(t_params%CMOMHOST(t_params%LMPOT,t_params%NEMBD1),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%CMOMHOST))*kind(t_params%CMOMHOST),'t_params%CMOMHOST','init_t_params')
      allocate(t_params%ECORE(20,t_params%NPOTD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%ECORE))*kind(t_params%ECORE),'t_params%ECORE','init_t_params')
      allocate(t_params%QMTET(t_params%NAEZ),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%QMTET))*kind(t_params%QMTET),'t_params%QMTET','init_t_params')
      allocate(t_params%QMPHI(t_params%NAEZ),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%QMPHI))*kind(t_params%QMPHI),'t_params%QMPHI','init_t_params')
      allocate(t_params%QMPHITAB(t_params%NAEZ,3),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%QMPHITAB))*kind(t_params%QMPHITAB),'t_params%QMPHITAB','init_t_params')
      allocate(t_params%QMTETTAB(t_params%NAEZ,3),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%QMTETTAB))*kind(t_params%QMTETTAB),'t_params%QMTETTAB','init_t_params')
      allocate(t_params%QMGAMTAB(t_params%NAEZ,3),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%QMGAMTAB))*kind(t_params%QMGAMTAB),'t_params%QMGAMTAB','init_t_params')
      allocate(t_params%ZAT(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%ZAT))*kind(t_params%ZAT),'t_params%ZAT','init_t_params')
      allocate(t_params%RMESH(t_params%IRM,t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%RMESH))*kind(t_params%RMESH),'t_params%RMESH','init_t_params')
      allocate(t_params%DRDI(t_params%IRM,t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%DRDI))*kind(t_params%DRDI),'t_params%DRDI','init_t_params')
      allocate(t_params%RMTREF(t_params%NREF),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%RMTREF))*kind(t_params%RMTREF),'t_params%RMTREF','init_t_params')
      allocate(t_params%VREF(t_params%NREF),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%VREF))*kind(t_params%VREF),'t_params%VREF','init_t_params')
      allocate(t_params%CLEB(t_params%NCLEB,2),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%CLEB))*kind(t_params%CLEB),'t_params%CLEB','init_t_params')
      allocate(t_params%RCLS(3,t_params%NACLSD,t_params%NCLSD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%RCLS))*kind(t_params%RCLS),'t_params%RCLS','init_t_params')
      allocate(t_params%SOCSCL(t_params%LMAXD1,t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%SOCSCL))*kind(t_params%SOCSCL),'t_params%SOCSCL','init_t_params')
      allocate(t_params%CSCL(t_params%LMAXD1,t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%CSCL))*kind(t_params%CSCL),'t_params%CSCL','init_t_params')
      allocate(t_params%RBASIS(3,t_params%NEMBD2),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%RBASIS))*kind(t_params%RBASIS),'t_params%RBASIS','init_t_params')
      allocate(t_params%RR(3,0:t_params%NR),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%RR))*kind(t_params%RR),'t_params%RR','init_t_params')
      allocate(t_params%CONC(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%CONC))*kind(t_params%CONC),'t_params%CONC','init_t_params')
      allocate(t_params%RROT(48,3,t_params%NSHELD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%RROT))*kind(t_params%RROT),'t_params%RROT','init_t_params')
      allocate(t_params%RATOM(3,t_params%NSHELD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%RATOM))*kind(t_params%RATOM),'t_params%RATOM','init_t_params')
      allocate(t_params%A(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%A))*kind(t_params%A),'t_params%A','init_t_params')
      allocate(t_params%B(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%B))*kind(t_params%B),'t_params%B','init_t_params')
      allocate(t_params%THETAS(t_params%IRID,t_params%NFUND,t_params%NCELLD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%THETAS))*kind(t_params%THETAS),'t_params%THETAS','init_t_params')
      allocate(t_params%RMT(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%RMT))*kind(t_params%RMT),'t_params%RMT','init_t_params')
      allocate(t_params%RMTNEW(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%RMTNEW))*kind(t_params%RMTNEW),'t_params%RMTNEW','init_t_params')
      allocate(t_params%RWS(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%RWS))*kind(t_params%RWS),'t_params%RWS','init_t_params')
      allocate(t_params%GSH(t_params%NGSHD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%GSH))*kind(t_params%GSH),'t_params%GSH','init_t_params')
      allocate(t_params%EREFLDAU(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%EREFLDAU))*kind(t_params%EREFLDAU),'t_params%EREFLDAU','init_t_params')
      allocate(t_params%UEFF(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%UEFF))*kind(t_params%UEFF),'t_params%UEFF','init_t_params')
      allocate(t_params%JEFF(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%JEFF))*kind(t_params%JEFF),'t_params%JEFF','init_t_params')
      allocate(t_params%ULDAU(t_params%MMAXD,t_params%MMAXD,&
         t_params%MMAXD,t_params%MMAXD,t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%ULDAU))*kind(t_params%ULDAU),'t_params%ULDAU','init_t_params')
      allocate(t_params%WLDAU(t_params%MMAXD,t_params%MMAXD,t_params%NSPIND,t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%WLDAU))*kind(t_params%WLDAU),'t_params%WLDAU','init_t_params')
      allocate(t_params%RPAN_INTERVALL(0:t_params%NTOTD,t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%RPAN_INTERVALL))*kind(t_params%RPAN_INTERVALL),'t_params%RPAN_INTERVALL','init_t_params')
      allocate(t_params%RNEW(t_params%NTOTD*(t_params%NCHEB+1),t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%RNEW))*kind(t_params%RNEW),'t_params%RNEW','init_t_params')
      allocate(t_params%MVEVI(t_params%NATYP,3,t_params%NMVECMAX),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%MVEVI))*kind(t_params%MVEVI),'t_params%MVEVI','init_t_params')
      allocate(t_params%MVEVIEF(t_params%NATYP,3,t_params%NMVECMAX),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%MVEVIEF))*kind(t_params%MVEVIEF),'t_params%MVEVIEF','init_t_params')
      allocate(t_params%THETASNEW(t_params%NTOTD*(t_params%NCHEB+1),&
         t_params%NFUND,t_params%NCELLD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%THETASNEW))*kind(t_params%THETASNEW),'t_params%THETASNEW','init_t_params')
      allocate(t_params%RHO2NS(t_params%IRM,t_params%LMPOT,t_params%NATYP,2),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%RHO2NS))*kind(t_params%RHO2NS),'t_params%RHO2NS','init_t_params')
      allocate(t_params%R2NEF(t_params%IRM,t_params%LMPOT,t_params%NATYP,2),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%R2NEF))*kind(t_params%R2NEF),'t_params%R2NEF','init_t_params')
      allocate(t_params%RHOC(t_params%IRM,t_params%NPOTD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%RHOC))*kind(t_params%RHOC),'t_params%RHOC','init_t_params')
      allocate(t_params%DENEFAT(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%DENEFAT))*kind(t_params%DENEFAT),'t_params%DENEFAT','init_t_params')
      allocate(t_params%ESPV(0:t_params%LMAXD1,t_params%NPOTD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%ESPV))*kind(t_params%ESPV),'t_params%ESPV','init_t_params')
      allocate(t_params%EDC(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%EDC))*kind(t_params%EDC),'t_params%EDC','init_t_params')
      allocate(t_params%EU(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%EU))*kind(t_params%EU),'t_params%EU','init_t_params')
      allocate(t_params%RHOORB(t_params%IRM*t_params%KREL+(1-t_params%KREL),t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%RHOORB))*kind(t_params%RHOORB),'t_params%RHOORB','init_t_params')
      allocate(t_params%ECOREREL(t_params%KREL*20+(1-t_params%KREL),t_params%NPOTD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%ECOREREL))*kind(t_params%ECOREREL),'t_params%ECOREREL','init_t_params')
      allocate(t_params%RCLSIMP(3,t_params%NATOMIMPD), stat=i_stat)    !real (kind=dp)
      call memocc(i_stat,product(shape(t_params%RCLSIMP))*kind(t_params%RCLSIMP),'t_params%ECOREREL','init_t_params')
      if (.not.allocated(t_params%theta)) then
         allocate(t_params%THETA(t_params%natyp),stat=i_stat)    !real (kind=dp)
         call memocc(i_stat,product(shape(t_params%THETA))*kind(t_params%THETA),'t_params%THETA','init_t_params')
         allocate(t_params%PHI(t_params%natyp), stat=i_stat)      !real (kind=dp)
         call memocc(i_stat,product(shape(t_params%PHI))*kind(t_params%PHI),'t_params%PHI','init_t_params')
      end if
      !-------------------------------------------------------------------------
      ! End of allocation of real (kind=dp) arrays
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      ! Allocate the integer arrays
      !-------------------------------------------------------------------------
      allocate(t_params%LOPT(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%LOPT))*kind(t_params%LOPT),'t_params%LOPT','init_t_params')
      allocate(t_params%ITLDAU(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%ITLDAU))*kind(t_params%ITLDAU),'t_params%ITLDAU','init_t_params')
      allocate(t_params%IRSHIFT(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%IRSHIFT))*kind(t_params%IRSHIFT),'t_params%IRSHIFT','init_t_params')
      allocate(t_params%JWSREL(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%JWSREL))*kind(t_params%JWSREL),'t_params%JWSREL','init_t_params')
      allocate(t_params%ZREL(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%ZREL))*kind(t_params%ZREL),'t_params%ZREL','init_t_params')
      allocate(t_params%LCORE(20,t_params%NPOTD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%LCORE))*kind(t_params%LCORE),'t_params%LCORE','init_t_params')
      allocate(t_params%NCORE(t_params%NPOTD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%NPOTD))*kind(t_params%NPOTD),'t_params%NPOTD','init_t_params')
      allocate(t_params%IPAN(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%IPAN))*kind(t_params%IPAN),'t_params%IPAN','init_t_params')
      allocate(t_params%IRCUT(0:t_params%IPAND,t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%IRCUT))*kind(t_params%IRCUT),'t_params%IRCUT','init_t_params')
      allocate(t_params%JEND(t_params%LMPOT,0:t_params%LMAX,0:t_params%LMAX),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%JEND))*kind(t_params%JEND),'t_params%JEND','init_t_params')
      allocate(t_params%ICLEB(t_params%NCLEB,4),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%ICLEB))*kind(t_params%ICLEB),'t_params%ICLEB','init_t_params')
      allocate(t_params%ATOM(t_params%NACLSD,t_params%NEMBD2),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%ATOM))*kind(t_params%ATOM),'t_params%ATOM','init_t_params')
      allocate(t_params%CLS(t_params%NEMBD2),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%CLS))*kind(t_params%CLS),'t_params%CLS','init_t_params')
      allocate(t_params%NACLS(t_params%NCLSD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%NACLS))*kind(t_params%NACLS),'t_params%NACLS','init_t_params')
      allocate(t_params%LOFLM(t_params%LM2D),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%LOFLM))*kind(t_params%LOFLM),'t_params%LOFLM','init_t_params')
      allocate(t_params%EZOA(t_params%NACLSD,t_params%NEMBD2),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%EZOA))*kind(t_params%EZOA),'t_params%EZOA','init_t_params')
      allocate(t_params%KAOEZ(t_params%NATYP,t_params%NEMBD2),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%KAOEZ))*kind(t_params%KAOEZ),'t_params%KAOEZ','init_t_params')
      allocate(t_params%IQAT(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%IQAT))*kind(t_params%IQAT),'t_params%IQAT','init_t_params')
      allocate(t_params%ICPA(t_params%NAEZ),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%ICPA))*kind(t_params%ICPA),'t_params%ICPA','init_t_params')
      allocate(t_params%NOQ(t_params%NAEZ),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%NOQ))*kind(t_params%NOQ),'t_params%NOQ','init_t_params')
      allocate(t_params%KMESH(t_params%IEMXD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%KMESH))*kind(t_params%KMESH),'t_params%KMESH','init_t_params')
      allocate(t_params%NSHELL(0:t_params%NSHELD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%NSHELL))*kind(t_params%NSHELL),'t_params%NSHELL','init_t_params')
      allocate(t_params%NSH1(t_params%NSHELD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%NSH1))*kind(t_params%NSH1),'t_params%NSH1','init_t_params')
      allocate(t_params%NSH2(t_params%NSHELD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%NSH2))*kind(t_params%NSH2),'t_params%NSH2','init_t_params')
      allocate(t_params%IJTABCALC(t_params%NOFGIJ),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%IJTABCALC))*kind(t_params%IJTABCALC),'t_params%IJTABCALC','init_t_params')
      allocate(t_params%IJTABCALC_I(t_params%NOFGIJ),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%IJTABCALC_I))*kind(t_params%IJTABCALC_I),'t_params%IJTABCALC_I','init_t_params')
      allocate(t_params%IJTABSYM(t_params%NOFGIJ),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%IJTABSYM))*kind(t_params%IJTABSYM),'t_params%IJTABSYM','init_t_params')
      allocate(t_params%IJTABSH(t_params%NOFGIJ),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%IJTABSH))*kind(t_params%IJTABSH),'t_params%IJTABSH','init_t_params')
      allocate(t_params%ISH(t_params%NSHELD,t_params%NOFGIJ),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%ISH))*kind(t_params%ISH),'t_params%ISH','init_t_params')
      allocate(t_params%JSH(t_params%NSHELD,t_params%NOFGIJ),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%JSH))*kind(t_params%JSH),'t_params%JSH','init_t_params')
      allocate(t_params%IQCALC(t_params%NAEZ),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%IQCALC))*kind(t_params%IQCALC),'t_params%IQCALC','init_t_params')
      allocate(t_params%ICHECK(t_params%NAEZDPD,t_params%NAEZDPD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%ICHECK))*kind(t_params%ICHECK),'t_params%ICHECK','init_t_params')
      allocate(t_params%ATOMIMP(t_params%NATOMIMPD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%ATOMIMP))*kind(t_params%ATOMIMP),'t_params%ATOMIMP','init_t_params')
      allocate(t_params%REFPOT(t_params%NEMBD2),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%REFPOT))*kind(t_params%REFPOT),'t_params%REFPOT','init_t_params')
      allocate(t_params%IRREL(2,2,t_params%LMMAXD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%IRREL))*kind(t_params%IRREL),'t_params%IRREL','init_t_params')
      allocate(t_params%NRREL(2,t_params%LMMAXD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%NRREL))*kind(t_params%NRREL),'t_params%NRREL','init_t_params')
      allocate(t_params%IFUNM1(t_params%LMXSPD,t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%IFUNM1))*kind(t_params%IFUNM1),'t_params%IFUNM1','init_t_params')
      allocate(t_params%ITITLE(20,t_params%NPOTD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%ITITLE))*kind(t_params%ITITLE),'t_params%ITITLE','init_t_params')
      allocate(t_params%LMSP1(t_params%LMXSPD,t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%LMSP1))*kind(t_params%LMSP1),'t_params%LMSP1','init_t_params')
      allocate(t_params%NTCELL(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%NTCELL))*kind(t_params%NTCELL),'t_params%NTCELL','init_t_params')
      allocate(t_params%IXIPOL(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%IXIPOL))*kind(t_params%IXIPOL),'t_params%IXIPOL','init_t_params')
      allocate(t_params%IRNS(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%IRNS))*kind(t_params%IRNS),'t_params%IRNS','init_t_params')
      allocate(t_params%IFUNM(t_params%NATYP,t_params%LMXSPD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%IFUNM))*kind(t_params%IFUNM),'t_params%IFUNM','init_t_params')
      allocate(t_params%LLMSP(t_params%NATYP,t_params%NFUND),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%LLMSP))*kind(t_params%LLMSP),'t_params%LLMSP','init_t_params')
      allocate(t_params%LMSP(t_params%NATYP,t_params%LMXSPD),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%LMSP))*kind(t_params%LMSP),'t_params%LMSP','init_t_params')
      allocate(t_params%IMT(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%IMT))*kind(t_params%IMT),'t_params%IMT','init_t_params')
      allocate(t_params%IRC(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%IRC))*kind(t_params%IRC),'t_params%IRC','init_t_params')
      allocate(t_params%IRMIN(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%IRMIN))*kind(t_params%IRMIN),'t_params%IRMIN','init_t_params')
      allocate(t_params%IRWS(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%IRWS))*kind(t_params%IRWS),'t_params%IRWS','init_t_params')
      allocate(t_params%NFU(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%NFU))*kind(t_params%NFU),'t_params%NFU','init_t_params')
      allocate(t_params%HOSTIMP(0:t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%HOSTIMP))*kind(t_params%HOSTIMP),'t_params%HOSTIMP','init_t_params')
      allocate(t_params%ILM_MAP(t_params%NGSHD,3),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%ILM_MAP))*kind(t_params%ILM_MAP),'t_params%ILM_MAP','init_t_params')
      allocate(t_params%IMAXSH(0:t_params%LMPOT),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%IMAXSH))*kind(t_params%IMAXSH),'t_params%IMAXSH','init_t_params')
      allocate(t_params%NPAN_LOG_AT(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%NPAN_LOG_AT))*kind(t_params%NPAN_LOG_AT),'t_params%NPAN_LOG_AT','init_t_params')
      allocate(t_params%NPAN_EQ_AT(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%NPAN_EQ_AT))*kind(t_params%NPAN_EQ_AT),'t_params%NPAN_EQ_AT','init_t_params')
      allocate(t_params%NPAN_TOT(t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%NPAN_TOT))*kind(t_params%NPAN_TOT),'t_params%NPAN_EQ_AT','init_t_params')
      allocate(t_params%IPAN_INTERVALL(0:t_params%NTOTD,t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%IPAN_INTERVALL))*kind(t_params%IPAN_INTERVALL),'t_params%IPAN_INTERVALL','init_t_params')
      allocate(t_params%NKCORE(20,t_params%NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(t_params%NKCORE))*kind(t_params%NKCORE),'t_params%NKCORE','init_t_params')
      allocate(t_params%KAPCORE(20,t_params%NPOTD),stat=i_stat) !INTEGER
      call memocc(i_stat,product(shape(t_params%KAPCORE))*kind(t_params%KAPCORE),'t_params%KAPCORE','init_t_params')
      if(.not. allocated(t_params%qdos_atomselect)) then
         allocate(t_params%qdos_atomselect(t_params%NATYP), stat=i_stat) !INTEGER
         call memocc(i_stat,product(shape(t_params%qdos_atomselect))*kind(t_params%qdos_atomselect),'t_params%qdos_atomselect','init_t_params')
      end if
      !-------------------------------------------------------------------------
      ! End of allocation of integer arrays
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      ! Allocate the logical arrays
      !-------------------------------------------------------------------------
      allocate(t_params%SYMUNITARY(t_params%NSYMAXD), stat=i_stat)!LOGICALS
      call memocc(i_stat,product(shape(t_params%SYMUNITARY))*kind(t_params%SYMUNITARY),'t_params%SYMUNITARY','init_t_params')
      allocate(t_params%VACFLAG(2), stat=i_stat)
      call memocc(i_stat,product(shape(t_params%VACFLAG))*kind(t_params%VACFLAG),'t_params%VACFLAG','init_t_params')
      !-------------------------------------------------------------------------
      ! End allocation of logical arrays
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      ! Allocate the character arrays
      !-------------------------------------------------------------------------
      allocate(t_params%TXC(6), stat=i_stat) !CHARACTER*124
      call memocc(i_stat,product(shape(t_params%TXC))*kind(t_params%TXC),'t_params%TXC','init_t_params')

      if(.not.allocated(t_params%TESTC)) then
         allocate(t_params%TESTC(32),stat=i_stat)
         call memocc(i_stat,product(shape(t_params%TESTC))*kind(t_params%TESTC),'t_params%TESTC','init_t_params')
         allocate(t_params%OPTC(32), stat=i_stat) !CHARACTER*8
         call memocc(i_stat,product(shape(t_params%OPTC))*kind(t_params%OPTC),'t_params%OPTC','init_t_params')
      end if
      !-------------------------------------------------------------------------
      ! End allocation of character arrays
      !-------------------------------------------------------------------------

      if (.not.allocated(t_params%nofks)) then
         allocate(t_params%BZKP(3,t_params%KPOIBZ,t_params%MAXMESH), stat=i_stat)! real (kind=dp)
         call memocc(i_stat,product(shape(t_params%BZKP))*kind(t_params%BZKP),'t_params%BZKP','init_t_params')
         allocate(t_params%VOLCUB(t_params%KPOIBZ,t_params%MAXMESH), stat=i_stat)! real (kind=dp)
         call memocc(i_stat,product(shape(t_params%VOLCUB))*kind(t_params%VOLCUB),'t_params%VOLCUB','init_t_params')
         allocate(t_params%VOLBZ(t_params%MAXMESH), stat=i_stat)! real (kind=dp)
         call memocc(i_stat,product(shape(t_params%VOLBZ))*kind(t_params%VOLBZ),'t_params%VOLBZ','init_t_params')
         allocate(t_params%NOFKS(t_params%MAXMESH), stat=i_stat) ! integer
         call memocc(i_stat,product(shape(t_params%NOFKS))*kind(t_params%NOFKS),'t_params%NOFKS','init_t_params')
      end if

   end subroutine init_t_params



#ifdef CPP_MPI
   !----------------------------------------------------------------------------
   ! subroutine: bcast_t_params_scalars
   !> @brief Broadcast scalar parameters via MPI
   !> @author Philipp Rüssmann
   !----------------------------------------------------------------------------
   subroutine bcast_t_params_scalars(t_params)
      ! broadcast scalar parameters, deal with arrays later
      use mpi
      use mod_mympi, only: master

      implicit none

      type(type_params), intent(inout) :: t_params
      integer :: ierr
      !integer :: myMPItype1
      !integer, dimension(t_params%Nscalars) :: blocklen1
      !integer, dimension(t_params%Nscalars) :: etype1
      !integer(kind=MPI_ADDRESS_KIND) :: disp1(t_params%Nscalars), base

      !!INTEGER
      !call MPI_Get_address(t_params%IEMXD,      disp1(1), ierr)
      !call MPI_Get_address(t_params%IRMIND,     disp1(2), ierr)
      !call MPI_Get_address(t_params%IRM,        disp1(3), ierr)
      !call MPI_Get_address(t_params%LMPOT,      disp1(4), ierr)
      !call MPI_Get_address(t_params%NSPOTD,     disp1(5), ierr)
      !call MPI_Get_address(t_params%NPOTD,      disp1(6), ierr)
      !call MPI_Get_address(t_params%NATYP,      disp1(7), ierr)
      !call MPI_Get_address(t_params%NEMBD1,     disp1(8), ierr)
      !call MPI_Get_address(t_params%LMMAXD,     disp1(9), ierr)
      !call MPI_Get_address(t_params%NAEZ,       disp1(10), ierr)
      !call MPI_Get_address(t_params%IPAND,      disp1(11), ierr)
      !call MPI_Get_address(t_params%NEMBD2,     disp1(12), ierr)
      !call MPI_Get_address(t_params%NREF,       disp1(13), ierr)
      !call MPI_Get_address(t_params%LMAX,       disp1(14), ierr)
      !call MPI_Get_address(t_params%NCLEB,      disp1(15), ierr)
      !call MPI_Get_address(t_params%NACLSD,     disp1(16), ierr)
      !call MPI_Get_address(t_params%NCLSD,      disp1(17), ierr)
      !call MPI_Get_address(t_params%LM2D,       disp1(18), ierr)
      !call MPI_Get_address(t_params%LMAXD1,     disp1(19), ierr)
      !call MPI_Get_address(t_params%NR,         disp1(20), ierr)
      !call MPI_Get_address(t_params%NSHELD,     disp1(21), ierr)
      !call MPI_Get_address(t_params%NSYMAXD,    disp1(22), ierr)
      !call MPI_Get_address(t_params%NAEZDPD,    disp1(23), ierr)
      !call MPI_Get_address(t_params%NATOMIMPD,  disp1(24), ierr)
      !call MPI_Get_address(t_params%NOFGIJ,     disp1(25), ierr)
      !call MPI_Get_address(t_params%NSPIND,     disp1(26), ierr)
      !call MPI_Get_address(t_params%NSPINDD,    disp1(27), ierr)
      !call MPI_Get_address(t_params%IRID,       disp1(28), ierr)
      !call MPI_Get_address(t_params%NFUND,      disp1(29), ierr)
      !call MPI_Get_address(t_params%NCELLD,     disp1(30), ierr)
      !call MPI_Get_address(t_params%LMXSPD,     disp1(31), ierr)
      !call MPI_Get_address(t_params%NGSHD,      disp1(32), ierr)
      !call MPI_Get_address(t_params%KREL,       disp1(33), ierr)
      !call MPI_Get_address(t_params%MMAXD,      disp1(34), ierr)
      !call MPI_Get_address(t_params%IELAST,     disp1(35), ierr)
      !call MPI_Get_address(t_params%NPOL,       disp1(36), ierr)
      !call MPI_Get_address(t_params%NPNT1,      disp1(37), ierr)
      !call MPI_Get_address(t_params%NPNT2,      disp1(38), ierr)
      !call MPI_Get_address(t_params%NPNT3,      disp1(39), ierr)
      !call MPI_Get_address(t_params%ITSCF,      disp1(40), ierr)
      !call MPI_Get_address(t_params%SCFSTEPS,   disp1(41), ierr)
      !call MPI_Get_address(t_params%LLY,        disp1(42), ierr)
      !call MPI_Get_address(t_params%NSRA,       disp1(43), ierr)
      !call MPI_Get_address(t_params%INS,        disp1(44), ierr)
      !call MPI_Get_address(t_params%KORBIT,     disp1(45), ierr)
      !call MPI_Get_address(t_params%KNOCO,      disp1(46), ierr)
      !call MPI_Get_address(t_params%NINEQ,      disp1(47), ierr)
      !call MPI_Get_address(t_params%KNOSPH,     disp1(48), ierr)
      !call MPI_Get_address(t_params%NSPIN,      disp1(49), ierr)
      !call MPI_Get_address(t_params%IRNSD,      disp1(50), ierr)
      !call MPI_Get_address(t_params%NPRINC,     disp1(51), ierr)
      !call MPI_Get_address(t_params%NCLS,       disp1(52), ierr)
      !call MPI_Get_address(t_params%ICST,       disp1(53), ierr)
      !call MPI_Get_address(t_params%IEND,       disp1(54), ierr)
      !call MPI_Get_address(t_params%ICC,        disp1(55), ierr)
      !call MPI_Get_address(t_params%IGF,        disp1(56), ierr)
      !call MPI_Get_address(t_params%NLBASIS,    disp1(57), ierr)
      !call MPI_Get_address(t_params%NRBASIS,    disp1(58), ierr)
      !call MPI_Get_address(t_params%NCPA,       disp1(59), ierr)
      !call MPI_Get_address(t_params%ITCPAMAX,   disp1(60), ierr)
      !call MPI_Get_address(t_params%KMROT,      disp1(61), ierr)
      !call MPI_Get_address(t_params%MAXMESH,    disp1(62), ierr)
      !call MPI_Get_address(t_params%NSYMAT,     disp1(63), ierr)
      !call MPI_Get_address(t_params%NATOMIMP,   disp1(64), ierr)
      !call MPI_Get_address(t_params%INVMOD,     disp1(65), ierr)
      !call MPI_Get_address(t_params%NQCALC,     disp1(66), ierr)
      !call MPI_Get_address(t_params%INTERVX,    disp1(67), ierr)
      !call MPI_Get_address(t_params%INTERVY,    disp1(68), ierr)
      !call MPI_Get_address(t_params%INTERVZ,    disp1(69), ierr)
      !call MPI_Get_address(t_params%LPOT,       disp1(70), ierr)
      !call MPI_Get_address(t_params%NRIGHT,     disp1(71), ierr)
      !call MPI_Get_address(t_params%NLEFT,      disp1(72), ierr)
      !call MPI_Get_address(t_params%IMIX,       disp1(73), ierr)
      !call MPI_Get_address(t_params%ITDBRY,     disp1(74), ierr)
      !call MPI_Get_address(t_params%KPRE,       disp1(75), ierr)
      !call MPI_Get_address(t_params%KSHAPE,     disp1(76), ierr)
      !call MPI_Get_address(t_params%KTE,        disp1(77), ierr)
      !call MPI_Get_address(t_params%KVMAD,      disp1(78), ierr)
      !call MPI_Get_address(t_params%KXC,        disp1(79), ierr)
      !call MPI_Get_address(t_params%ISHIFT,     disp1(80), ierr)
      !call MPI_Get_address(t_params%KFORCE,     disp1(81), ierr)
      !call MPI_Get_address(t_params%IDOLDAU,    disp1(82), ierr)
      !call MPI_Get_address(t_params%ITRUNLDAU,  disp1(83), ierr)
      !call MPI_Get_address(t_params%NTLDAU,     disp1(84), ierr)
      !call MPI_Get_address(t_params%NPOLSEMI,   disp1(85), ierr)
      !call MPI_Get_address(t_params%N1SEMI,     disp1(86), ierr)
      !call MPI_Get_address(t_params%N2SEMI,     disp1(87), ierr)
      !call MPI_Get_address(t_params%N3SEMI,     disp1(88), ierr)
      !call MPI_Get_address(t_params%IESEMICORE, disp1(89), ierr)
      !call MPI_Get_address(t_params%ITMPDIR,    disp1(90), ierr)
      !call MPI_Get_address(t_params%ILTMP,      disp1(91), ierr)
      !call MPI_Get_address(t_params%NCHEB,      disp1(92), ierr)
      !call MPI_Get_address(t_params%NTOTD,      disp1(93), ierr)
      !call MPI_Get_address(t_params%WLENGTH,    disp1(94), ierr)
      !call MPI_Get_address(t_params%NTPERD,     disp1(95), ierr)
      !!DOUBPLE PRECISION
      !call MPI_Get_address(t_params%EBOTSEMI,   disp1(96), ierr)
      !call MPI_Get_address(t_params%EMUSEMI,    disp1(97), ierr)
      !call MPI_Get_address(t_params%TKSEMI,     disp1(98), ierr)
      !call MPI_Get_address(t_params%FSEMICORE,  disp1(99), ierr)
      !call MPI_Get_address(t_params%R_LOG,      disp1(100), ierr)
      !call MPI_Get_address(t_params%EMIN,       disp1(101), ierr)
      !call MPI_Get_address(t_params%EMAX,       disp1(102), ierr)
      !call MPI_Get_address(t_params%TK,         disp1(103), ierr)
      !call MPI_Get_address(t_params%EFERMI,     disp1(104), ierr)
      !call MPI_Get_address(t_params%ALAT,       disp1(105), ierr)
      !call MPI_Get_address(t_params%CPATOL,     disp1(106), ierr)
      !call MPI_Get_address(t_params%MIXING,     disp1(107), ierr)
      !call MPI_Get_address(t_params%QBOUND,     disp1(108), ierr)
      !call MPI_Get_address(t_params%FCM,        disp1(109), ierr)
      !call MPI_Get_address(t_params%LAMBDA_XC,  disp1(110), ierr)
      !call MPI_Get_address(t_params%TOLRDIF,    disp1(111), ierr)
      !call MPI_Get_address(t_params%EFOLD,      disp1(112), ierr)
      !call MPI_Get_address(t_params%CHRGOLD,    disp1(113), ierr)
      !!!complex (kind=dp)
      !call MPI_Get_address(t_params%DELTAE,     disp1(114), ierr)
      !!LOGICAL
      !call MPI_Get_address(t_params%LINTERFACE, disp1(115), ierr)
      !call MPI_Get_address(t_params%LRHOSYM,    disp1(116), ierr)
      !!CHARACTER*10
      !call MPI_Get_address(t_params%SOLVER,     disp1(117), ierr)
      !!CHARACTER*80
      !call MPI_Get_address(t_params%TMPDIR,     disp1(118), ierr)
      !!INTEGER
      !call MPI_Get_address(t_params%I1,         disp1(119), ierr)
      !call MPI_Get_address(t_params%NMVECMAX,   disp1(120), ierr)
      !call MPI_Get_address(t_params%ITAB,       disp1(121), ierr)
      !real (kind=dp)
      !call MPI_Get_address(t_params%LASTERR,      disp1(122), ierr)
      !call MPI_Get_address(t_params%DENEF,        disp1(123), ierr)
      !call MPI_Get_address(t_params%CHRGSEMICORE, disp1(124), ierr)
      !INTEGER
      !call MPI_Get_address(t_params%NPAN_LOG    , disp1(125), ierr)
      !call MPI_Get_address(t_params%NPAN_EQ     , disp1(126), ierr)

      !base  = disp1(1)
      !disp1 = disp1 - base

      !blocklen1(1:116)=1
      !blocklen1(117)=10
      !blocklen1(118)=80
      !blocklen1(119:126)=1

      !etype1(1:95) = MPI_INTEGER
      !etype1(96:113) = MPI_DOUBLE_PRECISION
      !etype1(114) = MPI_DOUBLE_COMPLEX
      !etype1(115:116) = MPI_LOGICAL
      !etype1(117:118) = MPI_CHARACTER
      !etype1(119:121) = MPI_INTEGER
      !etype1(122:124) = MPI_DOUBLE_PRECISION
      !etype1(125:126) = MPI_INTEGER

      !call MPI_Type_create_struct(t_params%Nscalars, blocklen1, disp1,  &
      !   etype1, myMPItype1, ierr)
      !if(ierr/=MPI_SUCCESS) stop 'Problem in create_mpimask_t_params'

      !call MPI_Type_commit(myMPItype1, ierr)
      !if(ierr/=MPI_SUCCESS) stop 'error comiting create_mpimsk_t_params'

      !call MPI_Bcast(t_params%Nscalars, 1, myMPItype1, master,MPI_COMM_WORLD, ierr)
      !if(ierr/=MPI_SUCCESS) stop 'error brodcasting t_params'

      !call MPI_Type_free(myMPItype1, ierr)

      ! somehow this parameter gets overlooked in the communication, possibly a but somewhere, but for now this workaround does the job
      !call MPI_Bcast(t_params%NCHEB, 1, MPI_INTEGER, master,MPI_COMM_WORLD, ierr)


      ! try bcast for every variables instead:

!integer
      call MPI_Bcast(t_params%IEMXD     , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%IRMIND    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%IRM       , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%LMPOT     , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NSPOTD    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NPOTD     , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NATYP     , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NEMBD1    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%LMMAXD    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NAEZ      , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%IPAND     , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NEMBD2    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NREF      , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%LMAX      , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NCLEB     , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NACLSD    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NCLSD     , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%LM2D      , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%LMAXD1    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NR        , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NSHELD    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NSYMAXD   , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NAEZDPD   , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NATOMIMPD , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NOFGIJ    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NSPIND    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NSPINDD   , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%IRID      , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NFUND     , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NCELLD    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%LMXSPD    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NGSHD     , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%KREL      , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%MMAXD     , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%IELAST    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NPOL      , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NPNT1     , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NPNT2     , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NPNT3     , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%ITSCF     , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%SCFSTEPS  , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%LLY       , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NSRA      , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%INS       , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%KORBIT    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%KNOCO     , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NINEQ     , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%KNOSPH    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NSPIN     , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%IRNSD     , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NPRINC    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NCLS      , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%ICST      , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%IEND      , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%ICC       , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%IGF       , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NLBASIS   , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NRBASIS   , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NCPA      , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%ITCPAMAX  , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%KMROT     , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%MAXMESH   , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NSYMAT    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NATOMIMP  , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%INVMOD    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NQCALC    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%INTERVX   , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%INTERVY   , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%INTERVZ   , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%LPOT      , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NRIGHT    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NLEFT     , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%IMIX      , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%ITDBRY    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%KPRE      , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%KSHAPE    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%KTE       , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%KVMAD     , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%KXC       , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%ISHIFT    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%KFORCE    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%IDOLDAU   , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%ITRUNLDAU , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NTLDAU    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NPOLSEMI  , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%N1SEMI    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%N2SEMI    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%N3SEMI    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%IESEMICORE, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%ITMPDIR   , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%ILTMP     , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NCHEB     , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NTOTD     , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%WLENGTH   , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NTPERD    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
! forgot to bcast kpoibz?
      call MPI_Bcast(t_params%KPOIBZ    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)

!double precision
      call MPI_Bcast(t_params%EBOTSEMI  , 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%EMUSEMI   , 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%TKSEMI    , 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%FSEMICORE , 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%R_LOG     , 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%EMIN      , 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%EMAX      , 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%TK        , 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%EFERMI    , 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%ALAT      , 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%CPATOL    , 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%MIXING    , 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%QBOUND    , 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%FCM       , 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%LAMBDA_XC , 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%TOLRDIF   , 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%EFOLD     , 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%CHRGOLD   , 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)

!double complex
      call MPI_Bcast(t_params%DELTAE    , 1, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, ierr)

!logical
      call MPI_Bcast(t_params%LINTERFACE, 1, MPI_LOGICAL, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%LRHOSYM   , 1, MPI_LOGICAL, master, MPI_COMM_WORLD, ierr)

! character
      call MPI_Bcast(t_params%SOLVER    , 10, MPI_CHARACTER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%TMPDIR    , 80, MPI_CHARACTER, master, MPI_COMM_WORLD, ierr)

! integer
      call MPI_Bcast(t_params%I1        , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NMVECMAX  , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%ITAB      , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)

!double precision
      call MPI_Bcast(t_params%LASTERR, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)  
      call MPI_Bcast(t_params%DENEF, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)     
      call MPI_Bcast(t_params%CHRGSEMICORE, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)

! integer
      call MPI_Bcast(t_params%NPAN_LOG    , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(t_params%NPAN_EQ     , 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)


   end subroutine bcast_t_params_scalars

   !----------------------------------------------------------------------------
   ! SUBROUTINE: bcast_t_params_scalars
   !> @brief Broadcast arrays via MPI
   !> @author Philipp Rüssmann
   !----------------------------------------------------------------------------
   subroutine bcast_t_params_arrays(t_params)
      ! broadcast arrays from t_params
      use mpi
      use mod_mympi, only: master
      implicit none

      type(type_params), intent(inout) :: t_params
      integer :: ierr
 
      !-------------------------------------------------------------------------
      !complex (kind=dp) arrays
      !-------------------------------------------------------------------------
      call MPI_Bcast(t_params%EZ,t_params%IEMXD,&
         MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%WEZ,t_params%IEMXD,&
         MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%DROTQ,t_params%LMMAXD*t_params%LMMAXD*t_params%NAEZ,&
         MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%DSYMLL,t_params%LMMAXD*t_params%LMMAXD*t_params%NSYMAXD,&
         MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%LEFTTINVLL,t_params%LMMAXD*t_params%LMMAXD*&
         t_params%NEMBD1*t_params%NSPINDD*t_params%IEMXD,&
         MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RIGHTTINVLL,t_params%LMMAXD*t_params%LMMAXD*&
         t_params%NEMBD1*t_params%NSPINDD*t_params%IEMXD,&
         MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%CREL,t_params%LMMAXD*t_params%LMMAXD,&
         MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RC,t_params%LMMAXD*t_params%LMMAXD,&
         MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RREL,t_params%LMMAXD*t_params%LMMAXD,&
         MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%SRREL,2*2*t_params%LMMAXD,&
         MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%PHILDAU,t_params%IRM*t_params%NATYP,&
         MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)

      !-------------------------------------------------------------------------
      !real (kind=dp) arrays
      !-------------------------------------------------------------------------
      call MPI_Bcast(t_params%VINS,((t_params%IRM-t_params%IRMIND+1)*&
         t_params%LMPOT*t_params%NSPOTD),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%VISP,(t_params%IRM*t_params%NPOTD),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%VBC,2,&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%VTREL,(t_params%IRM*t_params%NATYP),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%BTREL,(t_params%IRM*t_params%NATYP),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%SOCSCALE,(t_params%NATYP),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%DRDIREL,(t_params%IRM*t_params%NATYP),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%R2DRDIREL,(t_params%IRM*t_params%NATYP),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RMREL,(t_params%IRM*t_params%NATYP),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%CMOMHOST,(t_params%LMPOT*t_params%NEMBD1),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ECORE,(20*t_params%NPOTD),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%QMTET,(t_params%NAEZ),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%QMPHI,(t_params%NAEZ),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%QMPHITAB,(t_params%NAEZ*3),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%QMTETTAB,(t_params%NAEZ*3),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%QMGAMTAB,(t_params%NAEZ*3),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ZAT,(t_params%NATYP),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RMESH,(t_params%IRM*t_params%NATYP),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%DRDI,(t_params%IRM*t_params%NATYP),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RMTREF,(t_params%NREF),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%VREF,(t_params%NREF),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%CLEB,(t_params%NCLEB*2),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RCLS,(3*t_params%NACLSD*t_params%NCLSD),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%SOCSCL,(t_params%LMAXD1*t_params%NATYP),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%CSCL,(t_params%LMAXD1*t_params%NATYP),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RBASIS,(3*t_params%NEMBD2),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RR,(3*(t_params%NR+1)),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%CONC,(t_params%NATYP),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RROT,(48*3*t_params%NSHELD),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RATOM,(3*t_params%NSHELD),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%A,t_params%NATYP,&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%B,t_params%NATYP,&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%THETAS,(t_params%IRID*t_params%NFUND*t_params%NCELLD),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RMT,(t_params%NATYP),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RMTNEW,(t_params%NATYP),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RWS,(t_params%NATYP),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%GSH,(t_params%NGSHD),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%EREFLDAU,(t_params%NATYP),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%UEFF,(t_params%NATYP),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%JEFF,(t_params%NATYP),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ULDAU,(t_params%MMAXD*t_params%MMAXD*&
         t_params%MMAXD*t_params%MMAXD*t_params%NATYP),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%WLDAU,(t_params%MMAXD*t_params%MMAXD*&
         t_params%NSPIND*t_params%NATYP),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RPAN_INTERVALL,((t_params%NTOTD+1)*t_params%NATYP),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RNEW,(t_params%NTOTD*(t_params%NCHEB+1)*t_params%NATYP),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%THETASNEW,(t_params%NTOTD*(t_params%NCHEB+1)*&
         t_params%NFUND*t_params%NCELLD),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%MVEVI,t_params%NATYP*3*t_params%NMVECMAX,&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%MVEVIEF,&
         t_params%NATYP*3*t_params%NMVECMAX,&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RHO2NS,&
         (t_params%IRM*t_params%LMPOT*t_params%NATYP*2),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%R2NEF,&
         (t_params%IRM*t_params%LMPOT*t_params%NATYP*2),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RHOC,(t_params%IRM*t_params%NPOTD),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%DENEFAT,(t_params%NATYP),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ESPV,(t_params%LMAXD1+1)*t_params%NPOTD,&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%EDC,(t_params%NATYP),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%EU,(t_params%NATYP),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RHOORB,&
         (t_params%IRM*t_params%KREL+(1-t_params%KREL)*t_params%NATYP),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ECOREREL,&
         (t_params%KREL*20+(1-t_params%KREL)*t_params%NPOTD),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%theta,t_params%natyp,&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%phi,t_params%natyp,&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RCLSIMP,3*t_params%NATOMIMPD,&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)

      !-------------------------------------------------------------------------
      !INTEGERS arrays
      !-------------------------------------------------------------------------
      call MPI_Bcast(t_params%LOPT,(t_params%NATYP),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ITLDAU,(t_params%NATYP),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IRSHIFT,(t_params%NATYP),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%JWSREL,(t_params%NATYP),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ZREL,(t_params%NATYP),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%LCORE,(20*t_params%NPOTD),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%NCORE,(t_params%NPOTD),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IPAN,(t_params%NATYP),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IRCUT,((t_params%IPAND+1)*t_params%NATYP),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%JEND,(t_params%LMPOT*(t_params%LMAX+1)*&
         (t_params%LMAX+1)),MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ICLEB,(t_params%NCLEB*4),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ATOM,(t_params%NACLSD*t_params%NEMBD2),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%CLS,(t_params%NEMBD2),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%NACLS,(t_params%NCLSD),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%LOFLM,(t_params%LM2D),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%EZOA,(t_params%NACLSD*t_params%NEMBD2),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%KAOEZ,(t_params%NATYP*t_params%NEMBD2),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IQAT,(t_params%NATYP),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ICPA,(t_params%NAEZ),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%NOQ,(t_params%NAEZ),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%KMESH,(t_params%IEMXD),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%NSHELL,((t_params%NSHELD+1)),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%NSH1,(t_params%NSHELD),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%NSH2,(t_params%NSHELD),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IJTABCALC,(t_params%NOFGIJ),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IJTABCALC_I,(t_params%NOFGIJ),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IJTABSYM,(t_params%NOFGIJ),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IJTABSH,(t_params%NOFGIJ),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ISH,(t_params%NSHELD*t_params%NOFGIJ),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%JSH,(t_params%NSHELD*t_params%NOFGIJ),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IQCALC,(t_params%NAEZ),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ICHECK,(t_params%NAEZDPD*t_params%NAEZDPD),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ATOMIMP,(t_params%NATOMIMPD),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%REFPOT,(t_params%NEMBD2),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IRREL,(2*2*t_params%LMMAXD),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%NRREL,(2*t_params%LMMAXD),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IFUNM1,(t_params%LMXSPD*t_params%NATYP),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ITITLE,(20*t_params%NPOTD),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%LMSP1,(t_params%LMXSPD*t_params%NATYP),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%NTCELL,(t_params%NATYP),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IXIPOL,(t_params%NATYP),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IRNS,(t_params%NATYP),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IFUNM,(t_params%NATYP*t_params%LMXSPD),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%LLMSP,(t_params%NATYP*t_params%NFUND),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%LMSP,(t_params%NATYP*t_params%LMXSPD),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IMT,(t_params%NATYP),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IRC,(t_params%NATYP),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IRMIN,(t_params%NATYP),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IRWS,(t_params%NATYP),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%NFU,(t_params%NATYP),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%HOSTIMP,((t_params%NATYP+1)),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ILM_MAP,(t_params%NGSHD*3),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IMAXSH,((t_params%LMPOT+1)),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%NPAN_LOG_AT,(t_params%NATYP),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%NPAN_EQ_AT,(t_params%NATYP),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%NPAN_TOT,(t_params%NATYP),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IPAN_INTERVALL,((t_params%NTOTD+1)*t_params%NATYP),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%NKCORE,(20*t_params%NATYP),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%KAPCORE,(20*t_params%NPOTD),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%qdos_atomselect,(t_params%NATYP),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr)

      !-------------------------------------------------------------------------
      ! LOGICAL arrays
      !-------------------------------------------------------------------------
      call MPI_Bcast(t_params%SYMUNITARY,(t_params%NSYMAXD),&
         MPI_LOGICAL,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%VACFLAG,2,MPI_LOGICAL,master,&
         MPI_COMM_WORLD,ierr)

      !-------------------------------------------------------------------------
      !CHARACTER arrays
      !-------------------------------------------------------------------------
      call MPI_Bcast(t_params%TXC,6*124, &!6 entries of length 124
         MPI_CHARACTER,master,MPI_COMM_WORLD,ierr) !CHARACTER*124
      call MPI_Bcast(t_params%TESTC,32*8, &!32 entries of length 8
         MPI_CHARACTER,master,MPI_COMM_WORLD,ierr) !CHARACTER*8
      call MPI_Bcast(t_params%OPTC,32*8, &!32 entries of length 8
         MPI_CHARACTER,master,MPI_COMM_WORLD,ierr) !CHARACTER*8

      !-------------------------------------------------------------------------
      ! K-points arrays
      !-------------------------------------------------------------------------
      call MPI_Bcast(t_params%BZKP,(3*t_params%KPOIBZ*t_params%MAXMESH),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr) ! real (kind=dp)
      call MPI_Bcast(t_params%VOLCUB,(t_params%KPOIBZ*t_params%MAXMESH),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr) ! real (kind=dp)
      call MPI_Bcast(t_params%VOLBZ,(t_params%MAXMESH),&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr) ! real (kind=dp)
      call MPI_Bcast(t_params%NOFKS,(t_params%MAXMESH),&
         MPI_INTEGER,master,MPI_COMM_WORLD,ierr) ! integer

   end subroutine bcast_t_params_arrays
#endif

   !----------------------------------------------------------------------------
   ! SUBROUTINE: fill_t_params_scalars
   !> @brief Set the values of the t_params scalars with the input values
   !> @author Philipp Rüssmann
   !----------------------------------------------------------------------------
   subroutine fill_t_params_scalars(IEMXD,IRMIND,IRM,LMPOT,NSPOTD,NPOTD,NATYP,   &
      NEMBD1,LMMAXD,NAEZ,IPAND,NEMBD2,NREF,LMAX,NCLEB,NACLSD,NCLSD,LM2D,LMAXD1,  &
      NR,NSHELD,NSYMAXD,NAEZDPD,NATOMIMPD,NOFGIJ,NSPIND,NSPINDD,IRID,NFUND,      &
      NCELLD,LMXSPD,NGSHD,KREL,MMAXD,IELAST,NPOL,NPNT1,NPNT2,NPNT3,ITSCF,        &
      SCFSTEPS,LLY,NSRA,INS,NINEQ,NSPIN,NCLS,ICST,IEND,ICC,IGF,NLBASIS,NRBASIS,  &
      NCPA,ITCPAMAX,KMROT,MAXMESH,NSYMAT,NATOMIMP,INVMOD,NQCALC,INTERVX,INTERVY, &
      INTERVZ,LPOT,NRIGHT,NLEFT,IMIX,ITDBRY,KPRE,KSHAPE,KTE,KVMAD,KXC,ISHIFT,    &
      KFORCE,IDOLDAU,ITRUNLDAU,NTLDAU,NPOLSEMI,N1SEMI,N2SEMI,N3SEMI,IESEMICORE,  &
      EBOTSEMI,EMUSEMI,TKSEMI,FSEMICORE,R_LOG,EMIN,EMAX,TK,EFERMI,ALAT,CPATOL,   &
      MIXING,QBOUND,FCM,LAMBDA_XC,TOLRDIF,LINTERFACE,LRHOSYM,SOLVER,TMPDIR,      &
      ITMPDIR,ILTMP,NTOTD,NCHEB,DELTAE,t_params)
      ! fill scalars into t_params
      implicit none

      type(type_params), intent(inout) :: t_params
      !     ..
      !     .. Scalar arguments
      integer, intent(in) :: NR        !< Number of real space vectors rr
      integer, intent(in) :: IRM       !< Maximum number of radial points
      integer, intent(in) :: INS       !< 0 (MT), 1(ASA), 2(Full Potential)
      integer, intent(in) :: LLY       !< LLY <> 0 : apply Lloyds formula
      integer, intent(in) :: ICC       !< Enables the calculation of off-diagonal elements of the GF.(0=SCF/DOS; 1=cluster; -1=custom)
      integer, intent(in) :: IGF       !< Do not print or print (0/1) the KKRFLEX_* files
      integer, intent(in) :: KTE       !< Calculation of the total energy On/Off (1/0)
      integer, intent(in) :: KXC       !< Type of xc-potential 0=vBH 1=MJW 2=VWN 3=PW91
      integer, intent(in) :: NREF      !< Number of diff. ref. potentials
      integer, intent(in) :: LM2D      !< (2*LMAX+1)**2
      integer, intent(in) :: KREL      !< Switch for non-relativistic/relativistic (0/1) program. Attention: several other parameters depend explicitly on KREL, they are set automatically Used for Dirac solver in ASA
      integer, intent(in) :: IRID      !< Shape functions parameters in non-spherical part
      integer, intent(in) :: LMAX      !< Maximum l component in wave function expansion
      integer, intent(in) :: NCLS      !< Number of reference clusters
      integer, intent(in) :: ICST      !< Number of Born approximation
      integer, intent(in) :: IEND      !< Number of nonzero gaunt coefficients
      integer, intent(in) :: NSRA
      integer, intent(in) :: LPOT      !< Maximum l component in potential expansion
      integer, intent(in) :: KPRE
      integer, intent(in) :: IMIX      !< Type of mixing scheme used (0=straight, 4=Broyden 2nd, 5=Anderson)
      integer, intent(in) :: NCPA      !< NCPA = 0/1 CPA flag
      integer, intent(in) :: NAEZ      !< Number of atoms in unit cell
      integer, intent(in) :: NPOL      !< Number of Matsubara Poles (EMESHT)
      integer, intent(in) :: NPNT1     !< number of E points (EMESHT) for the contour integration
      integer, intent(in) :: NPNT2     !< number of E points (EMESHT) for the contour integration
      integer, intent(in) :: NPNT3     !< number of E points (EMESHT) for the contour integration
      integer, intent(in) :: ITSCF
      integer, intent(in) :: IEMXD     !< Dimension for energy-dependent arrays
      integer, intent(in) :: LMPOT     !< (LPOT+1)**2
      integer, intent(in) :: NPOTD     !< (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP)
      integer, intent(in) :: NATYP     !< Number of kinds of atoms in unit cell
      integer, intent(in) :: IPAND     !< Number of panels in non-spherical part
      integer, intent(in) :: NCLEB     !< Number of Clebsch-Gordon coefficients
      integer, intent(in) :: NCLSD     !< Maximum number of different TB-clusters
      integer, intent(in) :: NFUND     !< Shape functions parameters in non-spherical part
      integer, intent(in) :: NGSHD     !< Shape functions parameters in non-spherical part
      integer, intent(in) :: MMAXD     !< 2*LMAX+1
      integer, intent(in) :: NINEQ     !< Number of ineq. positions in unit cell
      integer, intent(in) :: NSPIN     !< Counter for spin directions
      integer, intent(in) :: KMROT     !< 0: no rotation of the magnetisation; 1: individual rotation of the magnetisation for every site
      integer, intent(in) :: NTOTD
      integer, intent(in) :: NCHEB     !< Number of Chebychev pannels for the new solver
      integer, intent(in) :: KVMAD
      integer, intent(in) :: ILTMP
      integer, intent(in) :: NLEFT     !< Number of repeated basis for left host to get converged electrostatic potentials
      integer, intent(in) :: NRIGHT    !< Number of repeated basis for right host to get converged electrostatic potentials
      integer, intent(in) :: ITDBRY    !< Number of SCF steps to remember for the Broyden mixing
      integer, intent(in) :: KSHAPE    !< Exact treatment of WS cell
      integer, intent(in) :: ISHIFT
      integer, intent(in) :: KFORCE    !< Calculation of the forces
      integer, intent(in) :: IRMIND    !< IRM-IRNSD
      integer, intent(in) :: NSPOTD    !< Number of potentials for storing non-sph. potentials
      integer, intent(in) :: NEMBD1    !< NEMB+1
      integer, intent(in) :: LMMAXD    !< (KREL+KORBIT+1)(LMAX+1)^2
      integer, intent(in) :: NEMBD2
      integer, intent(in) :: NACLSD    !< Maximum number of atoms in a TB-cluster
      integer, intent(in) :: LMAXD1
      integer, intent(in) :: NSHELD    !< Number of blocks of the GF matrix that need to be calculated (NATYPD + off-diagonals in case of impurity)
      integer, intent(in) :: NOFGIJ    !< number of GF pairs IJ to be calculated as determined from IJTABCALC<>0
      integer, intent(in) :: NSPIND    !< KREL+(1-KREL)*(NSPIN+1)
      integer, intent(in) :: NCELLD    !< Number of cells (shapes) in non-spherical part
      integer, intent(in) :: LMXSPD    !< (2*LPOT+1)**2
      integer, intent(in) :: IELAST
      integer, intent(in) :: NSYMAT
      integer, intent(in) :: INVMOD    !< Inversion scheme
      integer, intent(in) :: NQCALC
      integer, intent(in) :: N1SEMI    !< Number of energy points for the semicore contour
      integer, intent(in) :: N2SEMI    !< Number of energy points for the semicore contour
      integer, intent(in) :: N3SEMI    !< Number of energy points for the semicore contour
      integer, intent(in) :: NTLDAU    !< number of atoms on which LDA+U is applied
      integer, intent(in) :: ITMPDIR
      integer, intent(in) :: MAXMESH
      integer, intent(in) :: NSYMAXD
      integer, intent(in) :: NAEZDPD
      integer, intent(in) :: NSPINDD   !< NSPIND-KORBIT
      integer, intent(in) :: NLBASIS   !< Number of basis layers of left host (repeated units)
      integer, intent(in) :: NRBASIS   !< Number of basis layers of right host (repeated units)
      integer, intent(in) :: INTERVX   !< Number of intervals in x-direction for k-net in IB of the BZ
      integer, intent(in) :: INTERVY   !< Number of intervals in y-direction for k-net in IB of the BZ
      integer, intent(in) :: INTERVZ   !< Number of intervals in z-direction for k-net in IB of the BZ
      integer, intent(in) :: IDOLDAU   !< flag to perform LDA+U
      integer, intent(in) :: SCFSTEPS  !< number of scf iterations
      integer, intent(in) :: ITCPAMAX  !< Max. number of CPA iterations
      integer, intent(in) :: NATOMIMP  !< Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
      integer, intent(in) :: NPOLSEMI  !< Number of poles for the semicore contour
      integer, intent(in) :: NATOMIMPD !< Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
      integer, intent(in) :: ITRUNLDAU !< Iteration index for LDA+U
      integer, intent(in) :: IESEMICORE
      real (kind=dp), intent(in) :: TK           !< Temperature
      real (kind=dp), intent(in) :: FCM
      real (kind=dp), intent(in) :: EMIN         !< Energies needed in EMESHT
      real (kind=dp), intent(in) :: EMAX         !< Energies needed in EMESHT
      real (kind=dp), intent(in) :: ALAT         !< Lattice constant in a.u.
      real (kind=dp), intent(in) :: R_LOG
      real (kind=dp), intent(in) :: TKSEMI       !< Temperature of semi-core contour
      real (kind=dp), intent(in) :: EFERMI       !< Fermi energy
      real (kind=dp), intent(in) :: CPATOL       !< Convergency tolerance for CPA-cycle
      real (kind=dp), intent(in) :: MIXING       !< Magnitude of the mixing parameter
      real (kind=dp), intent(in) :: QBOUND       !< Convergence parameter for the potential
      real (kind=dp), intent(in) :: EMUSEMI
      real (kind=dp), intent(in) :: TOLRDIF      !< Tolerance for r<tolrdif (a.u.) to handle vir. atoms
      real (kind=dp), intent(in) :: EBOTSEMI
      real (kind=dp), intent(in) :: FSEMICORE
      real (kind=dp), intent(in) :: LAMBDA_XC    !< Scale magnetic moment (0 < Lambda_XC < 1, 0=zero moment, 1= full moment)
      complex (kind=dp), intent(in) :: DELTAE         !< Energy difference for numerical derivative
      logical, intent(in) :: LRHOSYM
      logical, intent(in) :: LINTERFACE            !< If True a matching with semi-inifinite surfaces must be performed
      character(len=10), intent(in) :: SOLVER      !< Type of solver
      character(len=80), intent(in) :: TMPDIR
      !     ..
      ! fill scalars:
      ! Integer
      t_params%NR          = NR
      t_params%IRM         = IRM
      t_params%LLY         = LLY
      t_params%INS         = INS
      t_params%ICC         = ICC
      t_params%IGF         = IGF
      t_params%KTE         = KTE
      t_params%KXC         = KXC
      t_params%LM2D        = LM2D
      t_params%LPOT        = LPOT
      t_params%IMIX        = IMIX
      t_params%KPRE        = KPRE
      t_params%NSRA        = NSRA
      t_params%NREF        = NREF
      t_params%LMAX        = LMAX
      t_params%NCLS        = NCLS
      t_params%ICST        = ICST
      t_params%IEND        = IEND
      t_params%NCPA        = NCPA
      t_params%KREL        = KREL
      t_params%IRID        = IRID
      t_params%NAEZ        = NAEZ
      t_params%NPOL        = NPOL
      t_params%NPNT1       = NPNT1
      t_params%NPNT2       = NPNT2
      t_params%NPNT3       = NPNT3
      t_params%NPOTD       = NPOTD
      t_params%NATYP       = NATYP
      t_params%ITSCF       = ITSCF
      t_params%NSPIN       = NSPIN
      t_params%NINEQ       = NINEQ
      t_params%ILTMP       = ILTMP
      t_params%NCHEB       = NCHEB
      t_params%NTOTD       = NTOTD
      t_params%NCHEB       = NCHEB
      t_params%LMPOT       = LMPOT
      t_params%KMROT       = KMROT
      t_params%KVMAD       = KVMAD
      t_params%NGSHD       = NGSHD
      t_params%MMAXD       = MMAXD
      t_params%IEMXD       = IEMXD
      t_params%LMPOT       = LMPOT
      t_params%IPAND       = IPAND
      t_params%NCLEB       = NCLEB
      t_params%NCLSD       = NCLSD
      t_params%NFUND       = NFUND
      t_params%NLEFT       = NLEFT
      t_params%NRIGHT      = NRIGHT
      t_params%ITDBRY      = ITDBRY
      t_params%KSHAPE      = KSHAPE
      t_params%ISHIFT      = ISHIFT
      t_params%KFORCE      = KFORCE
      t_params%IRMIND      = IRMIND
      t_params%NSPOTD      = NSPOTD
      t_params%NEMBD1      = NEMBD1
      t_params%LMMAXD      = LMMAXD
      t_params%NEMBD2      = NEMBD2
      t_params%NACLSD      = NACLSD
      t_params%LMAXD1      = LMAXD1
      t_params%NSHELD      = NSHELD
      t_params%NOFGIJ      = NOFGIJ
      t_params%NSPIND      = NSPIND
      t_params%NCELLD      = NCELLD
      t_params%LMXSPD      = LMXSPD
      t_params%IELAST      = IELAST
      t_params%N1SEMI      = N1SEMI
      t_params%N2SEMI      = N2SEMI
      t_params%N3SEMI      = N3SEMI
      t_params%INVMOD      = INVMOD
      t_params%NQCALC      = NQCALC
      t_params%NSYMAT      = NSYMAT
      t_params%NTLDAU      = NTLDAU
      t_params%NLBASIS     = NLBASIS
      t_params%NRBASIS     = NRBASIS
      t_params%MAXMESH     = MAXMESH
      t_params%NSYMAXD     = NSYMAXD
      t_params%NAEZDPD     = NAEZDPD
      t_params%NSPINDD     = NSPINDD
      t_params%INTERVX     = INTERVX
      t_params%INTERVY     = INTERVY
      t_params%INTERVZ     = INTERVZ
      t_params%IDOLDAU     = IDOLDAU
      t_params%ITMPDIR     = ITMPDIR
      t_params%SCFSTEPS    = SCFSTEPS
      t_params%ITCPAMAX    = ITCPAMAX
      t_params%NATOMIMP    = NATOMIMP
      t_params%NPOLSEMI    = NPOLSEMI
      t_params%NATOMIMPD   = NATOMIMPD
      t_params%ITRUNLDAU   = ITRUNLDAU
      t_params%IESEMICORE  = IESEMICORE
      ! Double precision
      t_params%TK          = TK
      t_params%FCM         = FCM
      t_params%EMIN        = EMIN
      t_params%EMAX        = EMAX
      t_params%ALAT        = ALAT
      t_params%R_LOG       = R_LOG
      t_params%EFERMI      = EFERMI
      t_params%CPATOL      = CPATOL
      t_params%MIXING      = MIXING
      t_params%QBOUND      = QBOUND
      t_params%TKSEMI      = TKSEMI
      t_params%EMUSEMI     = EMUSEMI
      t_params%TOLRDIF     = TOLRDIF
      t_params%EBOTSEMI    = EBOTSEMI
      t_params%FSEMICORE   = FSEMICORE
      t_params%LAMBDA_XC   = LAMBDA_XC
      ! Double complex
      t_params%DELTAE      = DELTAE
      ! Logical
      t_params%LRHOSYM     = LRHOSYM
      t_params%LINTERFACE  = LINTERFACE
      ! Character
      t_params%SOLVER      = SOLVER
      t_params%TMPDIR      = TMPDIR

      t_params%NMVECMAX    = 4

   end subroutine fill_t_params_scalars

   !----------------------------------------------------------------------------
   ! SUBROUTINE: fill_t_params_arrays
   !> @brief Set the values of the t_params arrays with the input values of the arrays
   !> @author Philipp Rüssmann
   !----------------------------------------------------------------------------
   subroutine fill_t_params_arrays(t_params,IEMXD,LMMAXD,NAEZ,NSYMAXD,NEMBD1,    &
      NSPINDD,IRMIND,IRM,LMPOT,NSPOTD,NPOTD,NATYP,NR,NEMBD2,NREF,NCLEB,NCLSD,    &
      NACLSD,NSHELD,NGSHD,NFUND,IRID,NCELLD,MMAXD,LM2D,LMXSPD,LMAXD1,NSPIND,     &
      NTOTD,NCHEB,IPAND,LMAX,NOFGIJ,NAEZDPD,NATOMIMPD,EZ,WEZ,DROTQ,DSYMLL,       &
      LEFTTINVLL,RIGHTTINVLL,CREL,RC,RREL,SRREL,PHILDAU,VINS,VISP,VBC,VTREL,     &
      BTREL,SOCSCALE,DRDIREL,R2DRDIREL,RMREL,CMOMHOST,ECORE,QMTET,QMPHI,QMPHITAB,&
      QMTETTAB,QMGAMTAB,ZAT,R,DRDI,RMTREF,VREF,CLEB,RCLS,SOCSCL,CSCL,RBASIS,RR,  &
      CONC,RROT,RATOM,A,B,THETAS,RMT,RMTNEW,RWS,GSH,EREFLDAU,UEFF,JEFF,ULDAU,    &
      WLDAU,RPAN_INTERVALL,RNEW,THETASNEW,LOPT,ITLDAU,IRSHIFT,JWSREL,ZREL,LCORE, &
      NCORE,IPAN,IRCUT,JEND,ICLEB,ATOM,CLS,NACLS,LOFLM,EZOA,KAOEZ,IQAT,ICPA,NOQ, &
      KMESH,NSHELL,NSH1,NSH2,IJTABCALC,IJTABCALC_I,IJTABSYM,IJTABSH,ISH,JSH,     &
      IQCALC,ICHECK,ATOMIMP,REFPOT,IRREL,NRREL,IFUNM1,ITITLE,LMSP1,NTCELL,IXIPOL,&
      IRNS,IFUNM,LLMSP,LMSP,IMT,IRC,IRMIN,IRWS,NFU,HOSTIMP,ILM_MAP,IMAXSH,NPAN_LOG,  &
      NPAN_EQ,NPAN_LOG_AT,NPAN_EQ_AT,NPAN_TOT,IPAN_INTERVALL,SYMUNITARY,VACFLAG, &
      TXC,RCLSIMP, KREL)
      ! fill arrays after they have been allocated in init_t_params
      !     ..
      implicit none

      type(type_params), intent(inout) :: t_params
      !     ..
      !     .. Scalars for array dimensions
      integer, intent(in) :: NR        !< Number of real space vectors rr
      integer, intent(in) :: IRM       !< Maximum number of radial points
      integer, intent(in) :: LMAX      !< Maximum l component in wave function expansion
      integer, intent(in) :: NREF      !< Number of diff. ref. potentials
      integer, intent(in) :: NAEZ      !< Number of atoms in unit cell
      integer, intent(in) :: IRID      !< Shape functions parameters in non-spherical part
      integer, intent(in) :: LM2D      !< (2*LMAX+1)**2
      integer, intent(in) :: NCLEB     !< Number of Clebsch-Gordon coefficients
      integer, intent(in) :: NCLSD     !< Maximum number of different TB-clusters
      integer, intent(in) :: NPOTD     !< (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP)
      integer, intent(in) :: NATYP     !< Number of kinds of atoms in unit cell
      integer, intent(in) :: LMPOT     !< (LPOT+1)**2
      integer, intent(in) :: NGSHD     !< Shape functions parameters in non-spherical part
      integer, intent(in) :: NFUND     !< Shape functions parameters in non-spherical part
      integer, intent(in) :: MMAXD     !< 2*LMAX+1
      integer, intent(in) :: NTOTD
      integer, intent(in) :: NCHEB     !< Number of Chebychev pannels for the new solver
      integer, intent(in) :: IPAND     !< Number of panels in non-spherical part
      integer, intent(in) :: IEMXD     !< Dimension for energy-dependent arrays
      integer, intent(in) :: LMMAXD    !< (KREL+KORBIT+1)(LMAX+1)^2
      integer, intent(in) :: NEMBD1    !< NEMB+1
      integer, intent(in) :: IRMIND    !< IRM-IRNSD
      integer, intent(in) :: NSPOTD    !< Number of potentials for storing non-sph. potentials
      integer, intent(in) :: NEMBD2
      integer, intent(in) :: NACLSD    !< Maximum number of atoms in a TB-cluster
      integer, intent(in) :: NSHELD    !< Number of blocks of the GF matrix that need to be calculated (NATYPD + off-diagonals in case of impurity)
      integer, intent(in) :: NCELLD    !< Number of cells (shapes) in non-spherical part
      integer, intent(in) :: LMXSPD    !< (2*LPOT+1)**2
      integer, intent(in) :: LMAXD1
      integer, intent(in) :: NOFGIJ    !< number of GF pairs IJ to be calculated as determined from IJTABCALC<>0
      integer, intent(in) :: NSPIND    !< KREL+(1-KREL)*(NSPIN+1)
      integer, intent(in) :: NSYMAXD
      integer, intent(in) :: NSPINDD   !< NSPIND-KORBIT
      integer, intent(in) :: NAEZDPD
      integer, intent(in) :: NPAN_EQ
      integer, intent(in) :: NPAN_LOG
      integer, intent(in) :: NATOMIMPD !< Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
      integer, intent(in) :: KREL
      !     .. Array arguments
      complex (kind=dp), dimension(IEMXD), intent(in) :: EZ
      complex (kind=dp), dimension(IEMXD), intent(in) :: WEZ
      complex (kind=dp), dimension(LMMAXD,LMMAXD), intent(in) :: RC     !< NREL REAL spher. harm. > CMPLX. spher. harm. NREL CMPLX. spher. harm. > REAL spher. harm.
      complex (kind=dp), dimension(LMMAXD,LMMAXD), intent(in) :: CREL   !< Non-relat. CMPLX. spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. CMPLX. spher. harm.
      complex (kind=dp), dimension(LMMAXD,LMMAXD), intent(in) :: RREL   !< Non-relat. REAL spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. REAL spher. harm.
      complex (kind=dp), dimension(IRM,NATYP), intent(in) :: PHILDAU
      complex (kind=dp), dimension(LMMAXD,LMMAXD,NAEZ), intent(in) :: DROTQ   !< Rotation matrices to change between LOCAL/GLOBAL frame of reference for magnetisation <> Oz or noncollinearity
      complex (kind=dp), dimension(2,2,LMMAXD), intent(in) :: SRREL
      complex (kind=dp), dimension(LMMAXD,LMMAXD,NSYMAXD), intent(in) :: DSYMLL
      complex (kind=dp), dimension(LMMAXD,LMMAXD,NEMBD1,NSPINDD,IEMXD), intent(in) :: LEFTTINVLL
      complex (kind=dp), dimension(LMMAXD,LMMAXD,NEMBD1,NSPINDD,IEMXD), intent(in) :: RIGHTTINVLL

      real (kind=dp), dimension(NATYP), intent(in)  :: A                 !< Constants for exponential R mesh
      real (kind=dp), dimension(NATYP), intent(in)  :: B                 !< Constants for exponential R mesh
      real (kind=dp), dimension(2), intent(in)      :: VBC               !< Potential constants
      real (kind=dp), dimension(NATYP), intent(in)  :: RMT               !< Muffin-tin radius of true system
      real (kind=dp), dimension(NATYP), intent(in)  :: RWS               !< Wigner Seitz radius
      real (kind=dp), dimension(NGSHD), intent(in)  :: GSH
      real (kind=dp), dimension(NATYP), intent(in)  :: ZAT               !< Nuclear charge
      real (kind=dp), dimension(NATYP), intent(in)  :: UEFF              !< input U parameter for each atom
      real (kind=dp), dimension(NATYP), intent(in)  :: JEFF              !< input J parameter for each atom
      real (kind=dp), dimension(NATYP), intent(in)  :: CONC              !< Concentration of a given atom
      real (kind=dp), dimension(NREF), intent(in)   :: VREF
      real (kind=dp), dimension(NAEZ), intent(in)   :: QMTET             !< \f$ \theta\f$ angle of the agnetization with respect to the z-axis
      real (kind=dp), dimension(NAEZ), intent(in)   :: QMPHI             !< \f$ \phi\f$ angle of the agnetization with respect to the z-axis
      real (kind=dp), dimension(NREF), intent(in)   :: RMTREF            !< Muffin-tin radius of reference system
      real (kind=dp), dimension(NATYP), intent(in)  :: RMTNEW            !< Adapted muffin-tin radius
      real (kind=dp), dimension(NATYP), intent(in)  :: EREFLDAU          !< the energies of the projector's wave functions (REAL)
      real (kind=dp), dimension(NATYP), intent(in)  :: SOCSCALE          !< Spin-orbit scaling

      real (kind=dp), dimension(IRM,NATYP), intent(in)             :: R          !< Radial mesh ( in units a Bohr)
      real (kind=dp), dimension(3,0:NR), intent(in)                :: RR         !< Set of real space vectors (in a.u.)
      real (kind=dp), dimension(NCLEB,2), intent(in)               :: CLEB       !< GAUNT coefficients (GAUNT)
      real (kind=dp), dimension(IRM,NATYP), intent(in)             :: DRDI       !< Derivative dr/di
      real (kind=dp), dimension(IRM,NPOTD), intent(in)             :: VISP       !< Spherical part of the potential
      !real (kind=dp), dimension(LMAXD1,NATYP), intent(in)          :: CSCL       !< Speed of light scaling
      real (kind=dp), dimension(KREL*LMAX+1,KREL*NATYP+(1-KREL)), intent(inout) :: CSCL      !< Speed of light scaling
      real (kind=dp), dimension(NTOTD*(NCHEB+1),NATYP), intent(in) :: RNEW
      !real (kind=dp), dimension(IRM,NATYP), intent(in)             :: VTREL      !< potential (spherical part)
      real (kind=dp), dimension(IRM*KREL+(1-KREL),NATYP), intent(inout)         :: VTREL     !< potential (spherical part)
      !real (kind=dp), dimension(IRM,NATYP), intent(in)             :: BTREL      !< magnetic field
      real (kind=dp), dimension(IRM*KREL+(1-KREL),NATYP), intent(inout)         :: BTREL     !< magnetic field
      !real (kind=dp), dimension(IRM,NATYP), intent(in)             :: RMREL      !< radial mesh
      real (kind=dp), dimension(IRM*KREL+(1-KREL),NATYP), intent(inout)         :: RMREL     !< radial mesh
      real (kind=dp), dimension(20,NPOTD), intent(in)              :: ECORE      !< Core energies
      real (kind=dp), dimension(3,NSHELD), intent(in)              :: RATOM
      real (kind=dp), dimension(3,NEMBD2), intent(in)              :: RBASIS     !< Position of atoms in the unit cell in units of bravais vectors
      !real (kind=dp), dimension(LMAXD1,NATYP), intent(in)          :: SOCSCL
      real (kind=dp), dimension(KREL*LMAX+1,KREL*NATYP+(1-KREL)), intent(inout) :: SOCSCL
      real (kind=dp), dimension(3,NATOMIMPD), intent(in)           :: RCLSIMP
      !real (kind=dp), dimension(IRM,NATYP), intent(in)             :: DRDIREL    !< derivative of radial mesh
      real (kind=dp), dimension(IRM*KREL+(1-KREL),NATYP), intent(in) :: DRDIREL
      real (kind=dp), dimension(NAEZ,3), intent(in)                :: QMPHITAB
      real (kind=dp), dimension(NAEZ,3), intent(in)                :: QMTETTAB
      real (kind=dp), dimension(NAEZ,3), intent(in)                :: QMGAMTAB
      real (kind=dp), dimension(LMPOT,NEMBD1), intent(in)          :: CMOMHOST   !< Charge moments of each atom of the (left/right) host
      !real (kind=dp), dimension(IRM,NATYP), intent(in)             :: R2DRDIREL  !< \f$ r^2 \frac{\partial}{\partial \mathbf{r}}\frac{\partial}{\partial i}\f$ (r**2 * drdi)
      real (kind=dp), dimension(IRM*KREL+(1-KREL),NATYP), intent(in) :: R2DRDIREL
      real (kind=dp), dimension(0:NTOTD,NATYP), intent(in)         :: RPAN_INTERVALL

      real (kind=dp), dimension(3,NACLSD,NCLSD), intent(in)                 :: RCLS   !< Real space position of atom in cluster
      real (kind=dp), dimension(48,3,NSHELD), intent(in)                    :: RROT
      real (kind=dp), dimension(IRMIND:IRM,LMPOT,NSPOTD), intent(in)        :: VINS   !< Non-spherical part of the potential
      real (kind=dp), dimension(IRID,NFUND,NCELLD), intent(in)              :: THETAS !< shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
      real (kind=dp), dimension(NTOTD*(NCHEB+1),NFUND,NCELLD), intent(in)   :: THETASNEW
      real (kind=dp), dimension(MMAXD,MMAXD,NSPIND,NATYP), intent(in)       :: WLDAU  !< potential matrix
      real (kind=dp), dimension(MMAXD,MMAXD,MMAXD,MMAXD,NATYP), intent(in)  :: ULDAU  !< calculated Coulomb matrix elements (EREFLDAU)
      !     ..
      integer, dimension(NAEZ), intent(in)      :: NOQ        !< Number of diff. atom types located
      integer, dimension(NATYP), intent(in)     :: IMT        !< R point at MT radius
      integer, dimension(NATYP), intent(in)     :: IRC        !< R point for potential cutting
      integer, dimension(NATYP), intent(in)     :: NFU
      integer, dimension(NEMBD2), intent(in)    :: CLS        !< Cluster around atomic sites
      integer, dimension(NATYP), intent(in)     :: IRWS       !< R point at WS radius
      integer, dimension(NATYP), intent(in)     :: IRNS       !< Position of atoms in the unit cell in units of bravais vectors
      integer, dimension(NATYP), intent(in)     :: ZREL       !< atomic number (cast integer)
      integer, dimension(NATYP), intent(in)     :: IQAT       !< The site on which an atom is located on a given site
      integer, dimension(NAEZ), intent(in)      :: ICPA       !< ICPA = 0/1 site-dependent CPA flag
      integer, dimension(NATYP), intent(in)     :: IPAN       !< Number of panels in non-MT-region
      integer, dimension(NATYP), intent(in)     :: LOPT       !< angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
      integer, dimension(NSHELD), intent(in)    :: NSH1       !< Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
      integer, dimension(NSHELD), intent(in)    :: NSH2       !< Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
      integer, dimension(NCLSD), intent(in)     :: NACLS      !< Number of atoms in cluster
      integer, dimension(IEMXD), intent(in)     :: KMESH
      integer, dimension(NATYP), intent(in)     :: IRMIN      !< Max R for spherical treatment
      integer, dimension(LM2D), intent(in)      :: LOFLM      !< l of lm=(l,m) (GAUNT)
      integer, dimension(NPOTD), intent(in)     :: NCORE      !< Number of core states
      integer, dimension(NAEZ), intent(in)      :: IQCALC
      integer, dimension(NATYP), intent(in)     :: ITLDAU     !< integer pointer connecting the NTLDAU atoms to heir corresponding index in the unit cell
      integer, dimension(NATYP), intent(in)     :: JWSREL     !< index of the WS radius
      integer, dimension(NATYP), intent(in)     :: NTCELL     !< Index for WS cell
      integer, dimension(NATYP), intent(in)     :: IXIPOL     !< Constraint of spin pol.
      integer, dimension(NEMBD2), intent(in)    :: REFPOT     !< Ref. pot. card  at position
      integer, dimension(0:NSHELD), intent(in)  :: NSHELL     !< Index of atoms/pairs per shell (ij-pairs); nshell(0) = number of shells
      integer, dimension(0:LMPOT), intent(in)   :: IMAXSH
      integer, dimension(NATYP), intent(in)     :: IRSHIFT
      integer, dimension(NATOMIMPD), intent(in) :: ATOMIMP
      integer, dimension(NOFGIJ), intent(in)    :: IJTABSH    !< Linear pointer, assigns pair (i,j) to a shell in the array GS(*,*,*,NSHELD)
      integer, dimension(0:NATYP), intent(in)   :: HOSTIMP
      integer, dimension(NOFGIJ), intent(in)    :: IJTABSYM   !< Linear pointer, assigns pair (i,j) to the rotation bringing GS into Gij
      integer, dimension(NATYP), intent(in)     :: NPAN_TOT
      integer, dimension(NOFGIJ), intent(in)    :: IJTABCALC  !< Linear pointer, specifying whether the block (i,j) has to be calculated needs set up for ICC=-1, not used for ICC=1
      integer, dimension(NATYP), intent(in)     :: NPAN_EQ_AT
      integer, dimension(NATYP), intent(in)     :: NPAN_LOG_AT
      integer, dimension(NOFGIJ), intent(in)    :: IJTABCALC_I
      integer, dimension(NGSHD,3), intent(in)         :: ILM_MAP
      integer, dimension(NSHELD,NOFGIJ), intent(in)   :: ISH
      integer, dimension(NSHELD,NOFGIJ), intent(in)   :: JSH
      integer, dimension(NATYP,LMXSPD), intent(in)    :: LMSP     !< 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
      integer, dimension(NACLSD,NEMBD2), intent(in)   :: ATOM     !< Atom at site in cluster
      integer, dimension(NACLSD,NEMBD2), intent(in)   :: EZOA     !< EZ of atom at site in cluster
      integer, dimension(NCLEB,4), intent(in)         :: ICLEB    !< Pointer array
      integer, dimension(20,NPOTD), intent(in)        :: LCORE    !< Angular momentum of core states
      integer, dimension(2,LMMAXD), intent(in)        :: NRREL
      integer, dimension(NATYP,NFUND), intent(in)     :: LLMSP    !< lm=(l,m) of 'nfund'th nonvanishing component of non-spherical pot.
      integer, dimension(LMXSPD,NATYP), intent(in)    :: LMSP1
      integer, dimension(NATYP,LMXSPD), intent(in)    :: IFUNM
      integer, dimension(0:IPAND,NATYP), intent(in)   :: IRCUT    !< R points of panel borders
      integer, dimension(NATYP,NEMBD2), intent(in)    :: KAOEZ    !< Kind of atom at site in elem. cell
      integer, dimension(20,NPOTD), intent(in)        :: ITITLE
      integer, dimension(LMXSPD,NATYP), intent(in)    :: IFUNM1
      integer, dimension(NAEZDPD,NAEZDPD), intent(in) :: ICHECK
      integer, dimension(0:NTOTD,NATYP), intent(in)   :: IPAN_INTERVALL
      integer, dimension(LMPOT,0:LMAX,0:LMAX), intent(in) :: JEND !< Pointer array for icleb()
      integer, dimension(2,2,LMMAXD), intent(in)          :: IRREL

      logical, dimension(2), intent(in)       :: VACFLAG
      logical, dimension(NSYMAXD), intent(in) :: SYMUNITARY       !< unitary/antiunitary symmetry flag
      character(len=124), dimension(6), intent(in) :: TXC
      !     ..
      !fill arrays:

      t_params%EZ             = EZ
      t_params%WEZ            = WEZ
      t_params%DROTQ          = DROTQ
      t_params%DSYMLL         = DSYMLL
      t_params%LEFTTINVLL     = LEFTTINVLL
      t_params%RIGHTTINVLL    = RIGHTTINVLL
      t_params%CREL           = CREL
      t_params%RC             = RC
      t_params%RREL           = RREL
      t_params%SRREL          = SRREL
      t_params%PHILDAU        = PHILDAU
      t_params%VINS           = VINS
      t_params%VISP           = VISP
      t_params%VBC            = VBC
      if(t_params%KREL.gt.0) t_params%VTREL          = VTREL
      if(t_params%KREL.gt.0) t_params%BTREL          = BTREL
      t_params%SOCSCALE       = SOCSCALE
      if(t_params%KREL.gt.0) t_params%DRDIREL        = DRDIREL
      if(t_params%KREL.gt.0) t_params%R2DRDIREL      = R2DRDIREL
      if(t_params%KREL.gt.0) t_params%RMREL          = RMREL
      t_params%CMOMHOST       = CMOMHOST
      t_params%ECORE          = ECORE
      t_params%QMTET          = QMTET
      t_params%QMPHI          = QMPHI
      t_params%QMPHITAB       = QMPHITAB
      t_params%QMTETTAB       = QMTETTAB
      t_params%QMGAMTAB       = QMGAMTAB
      t_params%ZAT            = ZAT
      t_params%RMESH          = R
      t_params%DRDI           = DRDI
      t_params%RMTREF         = RMTREF
      t_params%VREF           = VREF
      t_params%CLEB           = CLEB
      t_params%RCLS           = RCLS
      t_params%SOCSCL         = SOCSCL
      t_params%CSCL           = CSCL
      t_params%RBASIS         = RBASIS
      t_params%RR             = RR
      t_params%CONC           = CONC
      t_params%RROT           = RROT
      t_params%RATOM          = RATOM
      t_params%A              = A
      t_params%B              = B
      t_params%THETAS         = THETAS
      t_params%RMT            = RMT
      t_params%RMTNEW         = RMTNEW
      t_params%RWS            = RWS
      t_params%GSH            = GSH
      t_params%EREFLDAU       = EREFLDAU
      t_params%UEFF           = UEFF
      t_params%JEFF           = JEFF
      if(t_params%IDOLDAU==1) t_params%ULDAU          = ULDAU
      t_params%WLDAU          = WLDAU
      t_params%RCLSIMP        = RCLSIMP
      t_params%RPAN_INTERVALL = RPAN_INTERVALL
      t_params%RNEW           = RNEW
      t_params%THETASNEW      = THETASNEW
      t_params%LOPT           = LOPT
      t_params%ITLDAU         = ITLDAU
      t_params%IRSHIFT        = IRSHIFT
      t_params%JWSREL         = JWSREL
      t_params%ZREL           = ZREL
      t_params%LCORE          = LCORE
      t_params%NCORE          = NCORE
      t_params%IPAN           = IPAN
      t_params%IRCUT          = IRCUT
      t_params%JEND           = JEND
      t_params%ICLEB          = ICLEB
      t_params%ATOM           = ATOM
      t_params%CLS            = CLS
      t_params%NACLS          = NACLS
      t_params%LOFLM          = LOFLM
      t_params%EZOA           = EZOA
      t_params%KAOEZ          = KAOEZ
      t_params%IQAT           = IQAT
      t_params%ICPA           = ICPA
      t_params%NOQ            = NOQ
      t_params%KMESH          = KMESH
      t_params%NSHELL         = NSHELL
      t_params%NSH1           = NSH1
      t_params%NSH2           = NSH2
      t_params%IJTABCALC      = IJTABCALC
      t_params%IJTABCALC_I    = IJTABCALC_I
      t_params%IJTABSYM       = IJTABSYM
      t_params%IJTABSH        = IJTABSH
      t_params%ISH            = ISH
      t_params%JSH            = JSH
      t_params%IQCALC         = IQCALC
      t_params%ICHECK         = ICHECK
      t_params%ATOMIMP        = ATOMIMP
      t_params%REFPOT         = REFPOT
      t_params%IRREL          = IRREL
      t_params%NRREL          = NRREL
      t_params%IFUNM1         = IFUNM1
      t_params%ITITLE         = ITITLE
      t_params%LMSP1          = LMSP1
      t_params%NTCELL         = NTCELL
      t_params%IXIPOL         = IXIPOL
      t_params%IRNS           = IRNS
      t_params%IFUNM          = IFUNM
      t_params%LLMSP          = LLMSP
      t_params%LMSP           = LMSP
      t_params%IMT            = IMT
      t_params%IRC            = IRC
      t_params%IRMIN          = IRMIN
      t_params%IRWS           = IRWS
      t_params%NFU            = NFU
      t_params%HOSTIMP        = HOSTIMP
      t_params%ILM_MAP            = ILM_MAP
      t_params%IMAXSH         = IMAXSH
      t_params%NPAN_LOG       = NPAN_LOG
      t_params%NPAN_EQ        = NPAN_EQ
      t_params%NPAN_LOG_AT    = NPAN_LOG_AT
      t_params%NPAN_EQ_AT     = NPAN_EQ_AT
      t_params%NPAN_TOT       = NPAN_TOT
      t_params%IPAN_INTERVALL = IPAN_INTERVALL
      t_params%SYMUNITARY     = SYMUNITARY
      t_params%VACFLAG        = VACFLAG
      t_params%TXC            = TXC
      t_params%NTOTD          = NTOTD

   end subroutine fill_t_params_arrays

   !----------------------------------------------------------------------------
   ! SUBROUTINE: get_params_1a
   !> @brief Set the values of the local variables according to the stored t_params
   !> so that they can be passed between different control modules, specifically for main1a
   !> @author Philipp Rüssmann
   !> @note JC: NPAN_EQ seems to have been passed here as an array, while in the
   !> rest of the routines it is an scalar. Why?
   !----------------------------------------------------------------------------
   subroutine get_params_1a(t_params,IPAND,NATYPD,IRMD,NACLSD,IELAST,NCLSD,NREFD,   &
      NCLEB,NEMBD,NAEZD,LM2D,NSRA,INS,NSPIN,ICST,IPAN,IRCUT,LMAX,NCLS,NINEQ,       &
      IDOLDAU,LLY,KREL,ATOM,CLS,ICLEB,LOFLM,NACLS,REFPOT,IRWS,IEND,EZ,VINS,IRMIN,&
      ITMPDIR,ILTMP,ALAT,DRDI,RMESH,ZAT,RCLS,IEMXD,VISP,RMTREF,VREF,CLEB,CSCL,   &
      SOCSCALE,SOCSCL,EREFLDAU,UEFF,JEFF,SOLVER,TMPDIR,DELTAE,TOLRDIF,NPAN_LOG,  &
      NPAN_EQ,NCHEB,NPAN_TOT,IPAN_INTERVALL,RPAN_INTERVALL,RNEW,NTOTD,NRMAXD,    &
      R_LOG,NTLDAU,ITLDAU,LOPT,VTREL,BTREL,DRDIREL,R2DRDIREL,RMREL,IRMIND,LMPOT, &
      NSPOTD,NPOTD,JWSREL,ZREL,ITSCF,NATOMIMPD,NATOMIMP,ATOMIMP,IQAT, NAEZ, NATYP, NREF)
      ! get relevant parameters from t_params
      !     ..
#ifdef CPP_MPI
      use mpi
#endif
      implicit none

      type(type_params), intent(in) :: t_params
      integer, intent(in) :: IRMD          !< Maximum number of radial points
      integer, intent(in) :: KREL         !< Switch for non-relativistic/relativistic (0/1) program. Attention: several other parameters depend explicitly on KREL, they are set automatically Used for Dirac solver in ASA
      integer, intent(in) :: NEMBD         !< Number of 'embedding' positions
      integer, intent(in) :: LM2D         !< (2*LMAX+1)**2
      integer, intent(in) :: NCLSD        !< Maximum number of different TB-clusters
      integer, intent(in) :: NCLEB        !< Number of Clebsch-Gordon coefficients
      integer, intent(in) :: NTOTD
      integer, intent(in) :: IPAND        !< Number of panels in non-spherical part
      integer, intent(in) :: IEMXD        !< Dimension for energy-dependent arrays
      integer, intent(in) :: LMPOT        !< (LPOT+1)**2
      integer, intent(in) :: NPOTD        !< (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP)
      integer, intent(in) :: NACLSD       !< Maximum number of atoms in a TB-cluster
      integer, intent(in) :: NRMAXD
      integer, intent(in) :: IRMIND       !< IRM-IRNSD
      integer, intent(in) :: NSPOTD       !< Number of potentials for storing non-sph. potentials
      integer, intent(in) :: NATOMIMPD    !< Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
      integer, intent(inout) :: LLY       !< LLY <> 0 : apply Lloyds formula
      integer, intent(inout) :: INS       !< 0 (MT), 1(ASA), 2(Full Potential)
      integer, intent(inout) :: NSRA
      integer, intent(inout) :: LMAX      !< Maximum l component in wave function expansion
      integer, intent(in) :: NREFD      !< Number of diff. ref. potentials
      integer, intent(out) :: NREF      !< Number of diff. ref. potentials
      integer, intent(inout) :: ICST      !< Number of Born approximation
      integer, intent(inout) :: NCLS      !< Number of reference clusters
      integer, intent(inout) :: IEND      !< Number of nonzero gaunt coefficients
      integer, intent(in) :: NAEZD      !< Number of atoms in unit cell
      integer, intent(out) :: NAEZ      !< Number of atoms in unit cell
      integer, intent(in) :: NATYPD     !< Number of kinds of atoms in unit cell
      integer, intent(out) :: NATYP     !< Number of kinds of atoms in unit cell
      integer, intent(inout) :: NSPIN     !< Counter for spin directions
      integer, intent(inout) :: NINEQ     !< Number of ineq. positions in unit cell
      integer, intent(inout) :: ILTMP
      integer, intent(inout) :: ITSCF
      integer, intent(inout) :: NCHEB     !< Number of Chebychev pannels for the new solver
      integer, intent(inout) :: NTLDAU    !< number of atoms on which LDA+U is applied
      integer, intent(inout) :: IELAST
      integer, intent(inout) :: ITMPDIR
      integer, intent(inout) :: IDOLDAU   !< flag to perform LDA+U
      integer, intent(inout) :: NATOMIMP  !< Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
      integer, dimension(NAEZD+NEMBD), intent(inout) :: CLS         !< Cluster around atomic sites
      integer, dimension(NATYPD), intent(inout)     :: LOPT        !< angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
      integer, dimension(NATYPD), intent(inout)     :: IRWS        !< R point at WS radius
      integer, dimension(NATYPD), intent(inout)     :: IPAN        !< Number of panels in non-MT-region
      integer, dimension(NATYPD), intent(inout)     :: ZREL        !< atomic number (cast integer)
      integer, dimension(NATYPD), intent(inout)     :: IQAT        !< The site on which an atom is located on a given site
      integer, dimension(LM2D), intent(inout)      :: LOFLM       !< l of lm=(l,m) (GAUNT)
      integer, dimension(NATYPD), intent(inout)     :: IRMIN       !< Max R for spherical treatment
      integer, dimension(NCLSD), intent(inout)     :: NACLS       !< Number of atoms in cluster
      integer, dimension(NAEZD+NEMBD), intent(inout) :: REFPOT      !< Ref. pot. card  at position
      integer, dimension(NATYPD), intent(inout)     :: ITLDAU      !< integer pointer connecting the NTLDAU atoms to heir corresponding index in the unit cell
      integer, dimension(NATYPD), intent(inout)     :: JWSREL      !< index of the WS radius
      integer, dimension(NATYPD), intent(inout)     :: NPAN_EQ
      integer, dimension(NATYPD), intent(inout)     :: NPAN_LOG
      integer, dimension(NATYPD), intent(inout)     :: NPAN_TOT
      integer, dimension(NATOMIMPD), intent(inout) :: ATOMIMP
      integer, dimension(NACLSD,NAEZD+NEMBD), intent(inout) :: ATOM    !< Atom at site in cluster
      integer, dimension(NCLEB,4), intent(inout)          :: ICLEB   !< Pointer array
      integer, dimension(0:IPAND,NATYPD), intent(inout)    :: IRCUT   !< R points of panel borders
      integer, dimension(0:NTOTD,NATYPD), intent(inout)    :: IPAN_INTERVALL
      real (kind=dp), intent(inout) :: ALAT                        !< Lattice constant in a.u.
      real (kind=dp), intent(inout) :: R_LOG
      real (kind=dp), intent(inout) :: TOLRDIF                     !< Tolerance for r<tolrdif (a.u.) to handle vir. atoms
      real (kind=dp), dimension(NATYPD), intent(inout) :: ZAT       !< Nuclear charge
      real (kind=dp), dimension(NREFD), intent(inout)  :: VREF
      real (kind=dp), dimension(NATYPD), intent(inout) :: UEFF      !< input U parameter for each atom
      real (kind=dp), dimension(NATYPD), intent(inout) :: JEFF      !< input J parameter for each atom
      real (kind=dp), dimension(NREFD), intent(inout)  :: RMTREF    !< Muffin-tin radius of reference system
      real (kind=dp), dimension(NATYPD), intent(inout) :: EREFLDAU  !< the energies of the projector's wave functions (REAL)
      real (kind=dp), dimension(NATYPD), intent(inout) :: SOCSCALE  !< Spin-orbit scaling
      real (kind=dp), dimension(NCLEB,2), intent(inout)                         :: CLEB      !< GAUNT coefficients (GAUNT)
      real (kind=dp), dimension(IRMD,NPOTD), intent(inout)                       :: VISP      !< Spherical part of the potential
      real (kind=dp), dimension(IRMD,NATYPD), intent(inout)                       :: DRDI      !< Derivative dr/di
      real (kind=dp), dimension(NRMAXD,NATYPD), intent(inout)                    :: RNEW
      real (kind=dp), dimension(KREL*LMAX+1,KREL*NATYPD+(1-KREL)), intent(inout) :: CSCL      !< Speed of light scaling
      real (kind=dp), dimension(IRMD,NATYPD), intent(inout)                       :: RMESH
      real (kind=dp), dimension(IRMD*KREL+(1-KREL),NATYPD), intent(inout)         :: VTREL     !< potential (spherical part)
      real (kind=dp), dimension(IRMD*KREL+(1-KREL),NATYPD), intent(inout)         :: BTREL     !< magnetic field
      real (kind=dp), dimension(IRMD*KREL+(1-KREL),NATYPD), intent(inout)         :: RMREL     !< radial mesh
      real (kind=dp), dimension(KREL*LMAX+1,KREL*NATYPD+(1-KREL)), intent(inout) :: SOCSCL
      real (kind=dp), dimension(IRMD*KREL+(1-KREL),NATYPD), intent(inout)         :: DRDIREL   !< derivative of radial mesh
      real (kind=dp), dimension(IRMD*KREL+(1-KREL),NATYPD), intent(inout)         :: R2DRDIREL !< \f$ r^2 \frac{\partial}{\partial \mathbf{r}}\frac{\partial}{\partial i}\f$ (r**2 * drdi)
      real (kind=dp), dimension(0:NTOTD,NATYPD), intent(inout)                   :: RPAN_INTERVALL
      real (kind=dp), dimension(3,NACLSD,NCLSD), intent(inout)          :: RCLS  !< Real space position of atom in cluster
      real (kind=dp), dimension(IRMIND:IRMD,LMPOT,NSPOTD), intent(inout) :: VINS  !< Non-spherical part of the potential

      complex (kind=dp), intent(inout) :: DELTAE      !< Energy difference for numerical derivative
      complex (kind=dp), dimension(IEMXD), intent(inout) :: EZ
      character(len=10), intent(inout) :: SOLVER   !< Type of solver
      character(len=80), intent(inout) :: TMPDIR


      NSRA    = t_params%NSRA
      INS     = t_params%INS
      NAEZ    = t_params%NAEZ
      NATYP   = t_params%NATYP
      NSPIN   = t_params%NSPIN
      ICST    = t_params%ICST
      IPAN    = t_params%IPAN
      IRCUT   = t_params%IRCUT
      LMAX    = t_params%LMAX
      NCLS    = t_params%NCLS
      NINEQ   = t_params%NINEQ
      NREF    = t_params%NREF
      IDOLDAU = t_params%IDOLDAU
      LLY     = t_params%LLY

      !-------------------------------------------------------------------------
      ! Consistency check
      !-------------------------------------------------------------------------
      if( (KREL.EQ.1) .AND. (INS.NE.0) ) then
         write(6,*) ' FULL-POTENTIAL RELATIVISTIC mode not implemented '
         stop ' set INS = 0 in the input'
      endif
      !
      if( NSRA.LE.2 ) then
         if( KREL.EQ.1 ) stop ' KVREL <= 1 in input, but relativistic program used'
      else
         if( KREL.EQ.0 ) stop ' KVREL > 1 in input, but non-relativistic program used'
      endif
      !-------------------------------------------------------------------------
      ! End of consistency check
      !-------------------------------------------------------------------------

      ALAT           = t_params%ALAT
      ZAT            = t_params%ZAT
      DRDI           = t_params%DRDI
      RMESH          = t_params%RMESH
      RMTREF         = t_params%RMTREF
      VREF           = t_params%VREF
      IEND           = t_params%IEND
      CLEB           = t_params%CLEB
      RCLS           = t_params%RCLS
      ATOM           = t_params%ATOM
      CLS            = t_params%CLS
      ICLEB          = t_params%ICLEB
      LOFLM          = t_params%LOFLM
      NACLS          = t_params%NACLS
      REFPOT         = t_params%REFPOT
      IRWS           = t_params%IRWS
      IRMIN          = t_params%IRMIN
      TOLRDIF        = t_params%TOLRDIF
      DELTAE         = t_params%DELTAE
      SOCSCALE       = t_params%SOCSCALE
      TMPDIR         = t_params%TMPDIR
      ITMPDIR        = t_params%ITMPDIR
      ILTMP          = t_params%ILTMP
      NPAN_LOG       = t_params%NPAN_LOG_AT
      NPAN_EQ        = t_params%NPAN_EQ_AT
      NCHEB          = t_params%NCHEB
      R_LOG          = t_params%R_LOG
      NPAN_TOT       = t_params%NPAN_TOT
      RNEW           = t_params%RNEW
      RPAN_INTERVALL = t_params%RPAN_INTERVALL
      IPAN_INTERVALL = t_params%IPAN_INTERVALL

      if(KREL.eq.1) then
         SOLVER = t_params%SOLVER
         SOCSCL = t_params%SOCSCL
         CSCL   = t_params%CSCL
      endif

      if( IDOLDAU.eq.1 ) then
         NTLDAU   = t_params%NTLDAU
         ITLDAU   = t_params%ITLDAU
         LOPT     = t_params%LOPT
         UEFF     = t_params%UEFF
         JEFF     = t_params%JEFF
         EREFLDAU = t_params%EREFLDAU
      endif
      !-------------------------------------------------------------------------
      ! Energy_mesh
      !-------------------------------------------------------------------------
      IELAST = t_params%IELAST
      EZ     = t_params%EZ
      !-------------------------------------------------------------------------
      ! Input_potential
      !-------------------------------------------------------------------------
      VINS  = t_params%VINS
      VISP  = t_params%VISP
      IF (KREL.EQ.1) THEN
         RMREL     = t_params%RMREL
         DRDIREL   = t_params%DRDIREL
         R2DRDIREL = t_params%R2DRDIREL
         ZREL      = t_params%ZREL
         JWSREL    = t_params%JWSREL
         VTREL     = t_params%VTREL
         BTREL     = t_params%BTREL
      END IF
      ITSCF    = t_params%ITSCF
      !-------------------------------------------------------------------------
      ! Cluster atoms
      !-------------------------------------------------------------------------
      NATOMIMP   = t_params%NATOMIMP
      ATOMIMP    = t_params%ATOMIMP
      IQAT       = t_params%IQAT

   end subroutine get_params_1a

   !----------------------------------------------------------------------------
   ! SUBROUTINE: get_params_1b
   !> @brief Set the values of the local variables according to the stored t_params
   !> so that they can be passed between different control modules, specifically for main1b
   !> @author Philipp Rüssmann
   !----------------------------------------------------------------------------
   subroutine get_params_1b(t_params,NATYPD, NAEZD, NATYP,NACLSD,IELAST,NPOL,NCLSD,NREFD, NREF,NEMBD,   &
      NAEZ,NSRA,INS,NSPIN,LMAX,NCLS,LLY,KREL,ATOM,CLS,NACLS,REFPOT,EZ,ITMPDIR,   &
      ILTMP,ALAT,RCLS,IEMXD,RMTREF,VREF,TMPDIR,NSHELD,NPRINCD,KPOIBZ,ATOMIMP,    &
      NATOMIMPD,ICC,IGF,NLBASIS,NRBASIS,NCPA,ICPA,ITCPAMAX,CPATOL,NRD,IDECI,      &
      RBASIS,RR,EZOA,NSHELL,KMROT,KAOEZ,ISH,JSH,NSH1,NSH2,NOQ,IQAT,NOFGIJ,       &
      NATOMIMP,CONC,KMESH,MAXMESH,NSYMAT,NQCALC,RATOM,RROT,DROTQ,IJTABCALC,      &
      IJTABCALC_I,IJTABSYM,IJTABSH,IQCALC,DSYMLL,INVMOD,ICHECK,SYMUNITARY,RC,    &
      CREL,RREL,SRREL,NRREL,IRREL,LEFTTINVLL,RIGHTTINVLL,VACFLAG,NOFKS,VOLBZ,    &
      BZKP,VOLCUB,WEZ,NEMBD1,LMMAXD,NSYMAXD,NSPINDD,MAXMSHD,RCLSIMP)
      ! get relevant parameters from t_params
      !     ..
      implicit none

      type(type_params), intent(in) :: t_params

      integer, intent(in) :: NATYPD
      integer, intent(in) :: NAEZD
      integer, intent(in) :: NEMBD
      integer, intent(in) :: NRD
      integer, intent(in) :: NREFD
      integer, intent(in) :: KREL
      integer, intent(in) :: IEMXD
      integer, intent(in) :: NCLSD
      integer, intent(in) :: NACLSD
      integer, intent(in) :: NEMBD1
      integer, intent(in) :: LMMAXD
      integer, intent(in) :: NSHELD
      integer, intent(in) :: KPOIBZ
      integer, intent(in) :: NSYMAXD
      integer, intent(in) :: NSPINDD
      integer, intent(in) :: MAXMSHD
      integer, intent(in) :: NPRINCD
      integer, intent(in) :: NATOMIMPD
      integer, intent(inout) :: LLY
      integer, intent(inout) :: ICC
      integer, intent(inout) :: IGF
      integer, intent(inout) :: INS
      integer, intent(inout) :: NPOL
      integer, intent(inout) :: NSRA
      integer, intent(inout) :: LMAX
      integer, intent(inout) :: NCLS
      integer, intent(out) :: NREF
      integer, intent(inout) :: NCPA
      integer, intent(out) :: NAEZ
      integer, intent(out) :: NATYP
      integer, intent(inout) :: NSPIN
      integer, intent(inout) :: ILTMP
      integer, intent(inout) :: IDECI
      integer, intent(inout) :: KMROT
      integer, intent(inout) :: INVMOD
      integer, intent(inout) :: NSYMAT
      integer, intent(inout) :: IELAST
      integer, intent(inout) :: NQCALC
      integer, intent(inout) :: NOFGIJ
      integer, intent(inout) :: MAXMESH
      integer, intent(inout) :: ITMPDIR
      integer, intent(inout) :: NLBASIS
      integer, intent(inout) :: NRBASIS
      integer, intent(inout) :: NATOMIMP
      integer, intent(inout) :: ITCPAMAX
      integer, dimension(NAEZD+NEMBD), intent(inout) :: CLS
      integer, dimension(NAEZD), intent(inout)      :: NOQ
      integer, dimension(NAEZD), intent(inout)      :: ICPA
      integer, dimension(NATYPD), intent(inout)     :: IQAT
      integer, dimension(NSHELD), intent(inout)    :: NSH1
      integer, dimension(NSHELD), intent(inout)    :: NSH2
      integer, dimension(IEMXD), intent(inout)     :: KMESH
      integer, dimension(NCLSD), intent(inout)     :: NACLS
      integer, dimension(MAXMSHD), intent(inout)   :: NOFKS
      integer, dimension(0:NSHELD), intent(inout)  :: NSHELL
      integer, dimension(NAEZD), intent(inout)      :: IQCALC
      integer, dimension(NAEZD+NEMBD), intent(inout) :: REFPOT
      integer, dimension(NATOMIMPD), intent(inout) :: ATOMIMP
      integer, dimension(NOFGIJ), intent(inout)    :: IJTABSH
      integer, dimension(NOFGIJ), intent(inout)    :: IJTABSYM
      integer, dimension(NOFGIJ), intent(inout)    :: IJTABCALC
      integer, dimension(NOFGIJ), intent(inout)    :: IJTABCALC_I
      integer, dimension(NSHELD,NOFGIJ), intent(inout)      :: ISH
      integer, dimension(NSHELD,NOFGIJ), intent(inout)      :: JSH
      integer, dimension(NACLSD,NAEZD+NEMBD), intent(inout)   :: EZOA
      integer, dimension(NACLSD,NAEZD+NEMBD), intent(inout)   :: ATOM
      integer, dimension(2,LMMAXD), intent(inout)           :: NRREL
      integer, dimension(NATYPD,NAEZD+NEMBD), intent(inout)    :: KAOEZ
      integer, dimension(NAEZD/NPRINCD,NAEZD/NPRINCD), intent(inout) :: ICHECK
      integer, dimension(2,2,LMMAXD), intent(inout) :: IRREL
      real (kind=dp), intent(inout) :: ALAT
      real (kind=dp), intent(inout) :: CPATOL
      real (kind=dp), dimension(NREFD), intent(inout)      :: VREF
      real (kind=dp), dimension(NATYPD), intent(inout)     :: CONC
      real (kind=dp), dimension(NREFD), intent(inout)      :: RMTREF
      real (kind=dp), dimension(MAXMSHD), intent(inout)   :: VOLBZ
      real (kind=dp), dimension(3,0:NRD), intent(inout)          :: RR
      real (kind=dp), dimension(3,NSHELD), intent(inout)        :: RATOM
      real (kind=dp), dimension(3,NAEZD+NEMBD), intent(inout)     :: RBASIS
      real (kind=dp), dimension(KPOIBZ,MAXMSHD), intent(inout)  :: VOLCUB
      real (kind=dp), dimension(3,NATOMIMPD), intent(inout)     :: RCLSIMP
      real (kind=dp), dimension(48,3,NSHELD), intent(inout)        :: RROT
      real (kind=dp), dimension(3,NACLSD,NCLSD), intent(inout)     :: RCLS
      real (kind=dp), dimension(3,KPOIBZ,MAXMSHD), intent(inout)   :: BZKP
      complex (kind=dp), dimension(IEMXD), intent(inout) :: EZ
      complex (kind=dp), dimension(IEMXD), intent(inout) :: WEZ
      complex (kind=dp), dimension(LMMAXD,LMMAXD), intent(inout) :: RC
      complex (kind=dp), dimension(LMMAXD,LMMAXD), intent(inout) :: CREL
      complex (kind=dp), dimension(LMMAXD,LMMAXD), intent(inout) :: RREL
      complex (kind=dp), dimension(LMMAXD,LMMAXD,NAEZD), intent(inout)      :: DROTQ
      complex (kind=dp), dimension(LMMAXD,LMMAXD,NSYMAXD), intent(inout)   :: DSYMLL
      complex (kind=dp), dimension(2,2,LMMAXD), intent(inout)              :: SRREL
      complex (kind=dp), dimension(LMMAXD,LMMAXD,NEMBD1,NSPINDD,IEMXD), intent(inout) :: LEFTTINVLL
      complex (kind=dp), dimension(LMMAXD,LMMAXD,NEMBD1,NSPINDD,IEMXD), intent(inout) :: RIGHTTINVLL
      character(len=80), intent(inout) :: TMPDIR
      logical, dimension(2), intent(inout)         :: VACFLAG
      logical, dimension(NSYMAXD), intent(inout)   :: SYMUNITARY

      integer :: i,l,id

      !     .. External Functions ..
      LOGICAL OPT,TEST
      EXTERNAL OPT,TEST

      NSRA     = t_params%NSRA
      INS      = t_params%INS
      NATYP    = t_params%NATYP
      NAEZ     = t_params%NAEZ
      NSPIN    = t_params%NSPIN
      LMAX     = t_params%LMAX
      NREF     = t_params%NREF
      ICC      = t_params%ICC
      IGF      = t_params%IGF
      NLBASIS  = t_params%NLBASIS
      NRBASIS  = t_params%NRBASIS
      NCPA     = t_params%NCPA
      ICPA     = t_params%ICPA
      ITCPAMAX = t_params%ITCPAMAX
      CPATOL   = t_params%CPATOL
      !-------------------------------------------------------------------------
      ! Consistency check
      !-------------------------------------------------------------------------
      !
      IF ( (KREL.EQ.1) .AND. (INS.NE.0) ) THEN
         WRITE(6,*) ' FULL-POTENTIAL RELATIVISTIC mode not implemented '
         STOP ' set INS = 0 in the input'
      END IF
      !
      IF ( NSRA.LE.2 ) THEN
         IF ( KREL.EQ.1 ) STOP ' KVREL <= 1 in input, but relativistic program used'
      ELSE
         IF ( KREL.EQ.0 ) STOP ' KVREL > 1 in input, but non-relativistic program used'
      END IF
      IDECI = 0
      !-------------------------------------------------------------------------
      ALAT        = t_params%ALAT
      RBASIS      = t_params%RBASIS
      REFPOT      = t_params%REFPOT
      RMTREF      = t_params%RMTREF
      VREF        = t_params%VREF
      RCLS        = t_params%RCLS
      RR          = t_params%RR
      ATOM        = t_params%ATOM
      CLS         = t_params%CLS
      NCLS        = t_params%NCLS
      EZOA        = t_params%EZOA
      NACLS       = t_params%NACLS
      NSHELL      = t_params%NSHELL
      KMROT       = t_params%KMROT
      KAOEZ       = t_params%KAOEZ
      IQAT        = t_params%IQAT
      NOQ         = t_params%NOQ
      CONC        = t_params%CONC
      KMESH       = t_params%KMESH
      MAXMESH     = t_params%MAXMESH
      LLY         = t_params%LLY
      NSYMAT      = t_params%NSYMAT
      NATOMIMP    = t_params%NATOMIMP
      NQCALC      = t_params%NQCALC
      RATOM       = t_params%RATOM
      RROT        = t_params%RROT
      NSH1        = t_params%NSH1
      NSH2        = t_params%NSH2
      DROTQ       = t_params%DROTQ
      IJTABCALC   = t_params%IJTABCALC
      IJTABCALC_I = t_params%IJTABCALC_I
      do I=1,NSHELL(0)
         do L=1,NOFGIJ
            ISH(I,L) = t_params%ISH(I,L)
            JSH(I,L) = t_params%JSH(I,L)
         end do
      end do
      IJTABSYM   = t_params%IJTABSYM
      IJTABSH    = t_params%IJTABSH
      IQCALC     = t_params%IQCALC
      DSYMLL     = t_params%DSYMLL
      INVMOD     = t_params%INVMOD
      ICHECK     = t_params%ICHECK
      ATOMIMP    = t_params%ATOMIMP
      SYMUNITARY = t_params%SYMUNITARY
      TMPDIR     = t_params%TMPDIR
      ITMPDIR    = t_params%ITMPDIR
      ILTMP      = t_params%ILTMP
      RCLSIMP    = t_params%RCLSIMP
      IF ( KREL.EQ.1 ) then
         RC    = t_params%RC
         CREL  = t_params%CREL
         RREL  = t_params%RREL
         SRREL = t_params%SRREL
         NRREL = t_params%NRREL
         IRREL = t_params%IRREL
      endif
      IF ( OPT('DECIMATE') ) THEN
         LEFTTINVLL  = t_params%LEFTTINVLL
         RIGHTTINVLL = t_params%RIGHTTINVLL
         VACFLAG     = t_params%VACFLAG
         IDECI = 1
      END IF
      !-------------------------------------------------------------------------
      ! K-points
      !-------------------------------------------------------------------------
      IF(TEST('kptsfile')) THEN
         OPEN (52,FILE='kpoints',FORM='formatted')
         REWIND (52)
         write(1337,*) 'kpoints read from kpoints file due to test option "kptsfile"'
      ENDIF
      DO L = 1,MAXMESH
         IF(TEST('kptsfile')) THEN
            READ (52,FMT='(I8,f15.10)') NOFKS(L),VOLBZ(L)
            write(1337,*) 'kpts:',NOFKS(L),VOLBZ(L),t_params%NOFKS(L),t_params%VOLBZ(L)
         ELSE
            NOFKS(L) = t_params%NOFKS(L)
            VOLBZ(L) = t_params%VOLBZ(L)
         ENDIF
         DO I=1,NOFKS(L)
            IF(TEST('kptsfile')) THEN
               READ (52,FMT=*) (BZKP(ID,I,L),ID=1,3),VOLCUB(I,L)
            ELSE
               DO ID=1,3
                  !                   write(*,*)'bzkp', BZKP(ID,I,L), t_params%BZKP(ID,I,L)
                  BZKP(ID,I,L) = t_params%BZKP(ID,I,L)
               ENDDO
               !                write(*,*)'volcub', VOLCUB(I,L), t_params%VOLCUB(I,L)
               VOLCUB(I,L) = t_params%VOLCUB(I,L)
            ENDIF
            !             write(*,'(A,4F21.17,I5,F21.17)') 'bzkmesh input:',
            !      &             (BZKP(ID,I,L),ID=1,3),VOLCUB(I,L),NOFKS(L),VOLBZ(L)
         END DO
      END DO
      IF(TEST('kptsfile')) CLOSE (52)
      !-------------------------------------------------------------------------
      ! Energy_mesh
      !-------------------------------------------------------------------------
      IELAST = t_params%IELAST
      EZ     = t_params%EZ
      WEZ    = t_params%WEZ
      NPOL   = t_params%NPOL
      !-------------------------------------------------------------------------
      ! Itermdir
      !-------------------------------------------------------------------------
      IF (OPT('ITERMDIR')) THEN
         DROTQ = t_params%DROTQ
         IF ( KMROT.EQ.0 ) KMROT = 1
      END IF

   end subroutine get_params_1b


   !----------------------------------------------------------------------------
   ! SUBROUTINE: get_params_1c
   !> @brief Set the values of the local variables according to the stored t_params
   !> so that they can be passed between different control modules, specifically for main1c
   !> @author Philipp Rüssmann
   !----------------------------------------------------------------------------
   subroutine get_params_1c(t_params,KREL,NAEZD,NATYPD,NCLEB,LM2D,NCHEB,IPAND,     &
      LMPOTD,LMAXD,LMXSPD,NFUND,NPOTD,NTOTD,MMAXD,IEMXD,IRMD,NSRA,INS,NSPIN,NACLS1, &
      ICST,KMROT,IQAT,IDOLDAU,IRWS,IPAN,IRCUT,IEND,ICLEB,LOFLM,JEND,IFUNM1,LMSP1,&
      NFU,LLMSP,LCORE,NCORE,NTCELL,IRMIN,ITITLE,INTERVX,INTERVY,INTERVZ,LLY,     &
      ITMPDIR,ILTMP,NPAN_EQ,IPAN_INTERVALL,NPAN_LOG,NPAN_TOT,NTLDAU,LOPT,ITLDAU, &
      IELAST,IESEMICORE,NPOL,IRSHIFT,JWSREL,ZREL,ITRUNLDAU,QMTET,QMPHI,CONC,ALAT,&
      ZAT,DRDI,RMESH,A,B,CLEB,THETAS,SOCSCALE,RPAN_INTERVALL,CSCL,RNEW,SOCSCL,   &
      THETASNEW,EFERMI,EREFLDAU,UEFF,JEFF,EMIN,EMAX,TK,VINS,VISP,ECORE,DRDIREL,  &
      R2DRDIREL,RMREL,VTREL,BTREL,WLDAU,ULDAU,EZ,WEZ,PHILDAU,TMPDIR,SOLVER,      &
      NSPIND,NSPOTD,IRMIND,LMAXD1,NCELLD,IRID,R_LOG, NAEZ, NATYP, LMAX)
      ! get relevant parameters from t_params
      !     ..
      implicit none

      type(type_params), intent(in) :: t_params

      integer, intent(in) :: IRMD
      integer, intent(in) :: KREL
      integer, intent(in) :: IRID
      integer, intent(in) :: LM2D
      integer, intent(in) :: NCLEB
      integer, intent(in) :: NTOTD
      integer, intent(in) :: IPAND
      integer, intent(in) :: LMPOTD
      integer, intent(in) :: NFUND
      integer, intent(in) :: NPOTD
      integer, intent(in) :: MMAXD
      integer, intent(in) :: IEMXD
      integer, intent(in) :: LMXSPD
      integer, intent(in) :: NSPIND
      integer, intent(in) :: NSPOTD
      integer, intent(in) :: IRMIND
      integer, intent(in) :: LMAXD1
      integer, intent(in) :: NCELLD
      integer, intent(inout) :: INS
      integer, intent(inout) :: LLY
      integer, intent(inout) :: NSRA
      integer, intent(inout) :: ICST
      integer, intent(in) :: NAEZD
      integer, intent(out) :: NAEZ
      integer, intent(in) :: LMAXD
      integer, intent(out) :: LMAX
      integer, intent(inout) :: NPOL
      integer, intent(inout) :: IEND
      integer, intent(inout) :: NCHEB
      integer, intent(in) :: NATYPD
      integer, intent(out) :: NATYP
      integer, intent(inout) :: NSPIN
      integer, intent(inout) :: KMROT
      integer, intent(inout) :: ILTMP
      integer, intent(inout) :: NTLDAU
      integer, intent(inout) :: IELAST
      integer, intent(inout) :: NACLS1
      integer, intent(inout) :: IDOLDAU
      integer, intent(inout) :: INTERVX
      integer, intent(inout) :: INTERVY
      integer, intent(inout) :: INTERVZ
      integer, intent(inout) :: ITMPDIR
      integer, intent(inout) :: ITRUNLDAU
      integer, intent(inout) :: IESEMICORE
      integer, dimension(NATYPD), intent(inout) :: NFU
      integer, dimension(NATYPD), intent(inout) :: LOPT
      integer, dimension(NATYPD), intent(inout) :: IQAT
      integer, dimension(NATYPD), intent(inout) :: IRWS
      integer, dimension(NATYPD), intent(inout) :: IPAN
      integer, dimension(NATYPD), intent(inout) :: ZREL
      integer, dimension(LM2D), intent(inout) :: LOFLM
      integer, dimension(NPOTD), intent(inout) :: NCORE
      integer, dimension(NATYPD), intent(inout) :: IRMIN
      integer, dimension(NATYPD), intent(inout) :: NTCELL
      integer, dimension(NATYPD), intent(inout) :: ITLDAU
      integer, dimension(NATYPD), intent(inout) :: JWSREL
      integer, dimension(NATYPD), intent(inout) :: NPAN_EQ
      integer, dimension(NATYPD), intent(inout) :: IRSHIFT
      integer, dimension(NATYPD), intent(inout) :: NPAN_LOG
      integer, dimension(NATYPD), intent(inout) :: NPAN_TOT
      integer, dimension(LMPOTD,0:LMAXD,0:LMAXD), intent(inout)   :: JEND
      integer, dimension(0:IPAND,NATYPD), intent(inout)         :: IRCUT
      integer, dimension(NCLEB,4), intent(inout)               :: ICLEB
      integer, dimension(LMXSPD,NATYPD), intent(inout)          :: LMSP1
      integer, dimension(NATYPD,NFUND), intent(inout)           :: LLMSP
      integer, dimension(20,NPOTD), intent(inout)              :: LCORE
      integer, dimension(LMXSPD,NATYPD), intent(inout)          :: IFUNM1
      integer, dimension(20,NPOTD), intent(inout)              :: ITITLE
      integer, dimension(0:NTOTD,NATYPD), intent(inout)         :: IPAN_INTERVALL

      real (kind=dp), intent(inout) :: TK
      real (kind=dp), intent(inout) :: EMIN
      real (kind=dp), intent(inout) :: EMAX
      real (kind=dp), intent(inout) :: ALAT
      real (kind=dp), intent(inout) :: R_LOG
      real (kind=dp), intent(inout) :: EFERMI

      real (kind=dp), dimension(NATYPD), intent(inout)  :: A
      real (kind=dp), dimension(NATYPD), intent(inout)  :: B
      real (kind=dp), dimension(NATYPD), intent(inout)  :: ZAT
      real (kind=dp), dimension(NATYPD), intent(inout)  :: CONC
      real (kind=dp), dimension(NATYPD), intent(inout)  :: UEFF
      real (kind=dp), dimension(NATYPD), intent(inout)  :: JEFF
      real (kind=dp), dimension(NAEZD), intent(inout)   :: QMPHI
      real (kind=dp), dimension(NAEZD), intent(inout)   :: QMTET
      real (kind=dp), dimension(NATYPD), intent(inout)  :: SOCSCALE
      real (kind=dp), dimension(NATYPD), intent(inout)  :: EREFLDAU
      real (kind=dp), dimension(IRMD,NATYPD), intent(inout)             :: DRDI
      real (kind=dp), dimension(IRMD,NATYPD), intent(inout)             :: RMESH
      real (kind=dp), dimension(NCLEB,2), intent(inout)               :: CLEB
      !real (kind=dp), dimension(LMAXD1,NATYPD), intent(inout)          :: CSCL
      real (kind=dp), dimension(KREL*LMAXD+1,KREL*NATYPD+(1-KREL)), intent(inout) :: CSCL      !< Speed of light scaling
      real (kind=dp), dimension(IRMD,NPOTD), intent(inout)             :: VISP
      real (kind=dp), dimension(NTOTD*(NCHEB+1),NATYPD), intent(inout) :: RNEW
      !real (kind=dp), dimension(IRMD,NATYPD), intent(inout)             :: RMREL
      real (kind=dp), dimension(IRMD*KREL+(1-KREL),NATYPD), intent(inout)         :: RMREL     !< radial mesh
      !real (kind=dp), dimension(IRMD,NATYPD), intent(inout)             :: VTREL
      real (kind=dp), dimension(IRMD*KREL+(1-KREL),NATYPD), intent(inout)         :: VTREL     !< potential (spherical part)
      !real (kind=dp), dimension(IRMD,NATYPD), intent(inout)             :: BTREL
      real (kind=dp), dimension(IRMD*KREL+(1-KREL),NATYPD), intent(inout)         :: BTREL     !< magnetic field
      real (kind=dp), dimension(20,NPOTD), intent(inout)              :: ECORE
      !real (kind=dp), dimension(LMAXD1,NATYPD), intent(inout)          :: SOCSCL
      real (kind=dp), dimension(KREL*LMAXD+1,KREL*NATYPD+(1-KREL)), intent(inout) :: SOCSCL
      !real (kind=dp), dimension(IRMD,NATYPD), intent(inout)             :: DRDIREL
      real (kind=dp), dimension(IRMD*KREL+(1-KREL),NATYPD), intent(inout)         :: DRDIREL
      !real (kind=dp), dimension(IRMD,NATYPD), intent(inout)             :: R2DRDIREL
      real (kind=dp), dimension(IRMD*KREL+(1-KREL),NATYPD), intent(inout)         :: R2DRDIREL !< \f$ r^2 \frac{\partial}{\partial \mathbf{r}}\frac{\partial}{\partial i}\f$ (r**2 * drdi)
      real (kind=dp), dimension(0:NTOTD,NATYPD), intent(inout)         :: RPAN_INTERVALL
      real (kind=dp), dimension(IRMIND:IRMD,LMPOTD,NSPOTD), intent(inout)        :: VINS
      real (kind=dp), dimension(IRID,NFUND,NCELLD), intent(inout)              :: THETAS
      real (kind=dp), dimension(NTOTD*(NCHEB+1),NFUND,NCELLD), intent(inout)   :: THETASNEW
      real (kind=dp), dimension(MMAXD,MMAXD,NSPIND,NATYPD), intent(inout) :: WLDAU
      real (kind=dp), dimension(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD), intent(inout) :: ULDAU
      complex (kind=dp), dimension(IEMXD), intent(inout) :: EZ
      complex (kind=dp), dimension(IEMXD), intent(inout) :: WEZ
      complex (kind=dp), dimension(IRMD,NATYPD), intent(inout) :: PHILDAU
      character(len=10), intent(inout) :: SOLVER
      character(len=80), intent(inout) :: TMPDIR
      !     .. External Functions ..
      LOGICAL OPT,TEST
      EXTERNAL OPT,TEST

      !
      NSRA    = t_params%NSRA
      INS     = t_params%INS
      NATYP   = t_params%NATYP
      NAEZ    = t_params%NAEZ
      NSPIN   = t_params%NSPIN
      ICST    = t_params%ICST
      IPAN    = t_params%IPAN
      IRCUT   = t_params%IRCUT
      KMROT   = t_params%KMROT
      IQAT    = t_params%IQAT
      CONC    = t_params%CONC
      QMTET   = t_params%QMTET
      QMPHI   = t_params%QMPHI
      IDOLDAU = t_params%IDOLDAU
      LMAX    = t_params%LMAX
      IRWS    = t_params%IRWS
      !-------------------------------------------------------------------------
      ! Consistency check
      !-------------------------------------------------------------------------
      IF ( ( KREL.EQ.1 ) .AND. ( INS.NE.0 ) ) THEN
         WRITE(6,*) ' FULL-POTENTIAL RELATIVISTIC mode not implemented '
         STOP ' set INS = 0 in the input'
      END IF
      !
      IF ( NSRA.LE.2 ) THEN
         IF ( KREL.EQ.1 ) STOP ' KVREL <= 1 in input, but relativistic program used'
      ELSE
         IF ( KREL.EQ.0 ) STOP ' KVREL > 1 in input, but non-relativistic program used'
      END IF
      !-------------------------------------------------------------------------
      ALAT           = t_params%ALAT
      ZAT            = t_params%ZAT
      DRDI           = t_params%DRDI
      RMESH          = t_params%RMESH
      A              = t_params%A
      B              = t_params%B
      IEND           = t_params%IEND
      CLEB           = t_params%CLEB
      ICLEB          = t_params%ICLEB
      LOFLM          = t_params%LOFLM
      JEND           = t_params%JEND
      THETAS         = t_params%THETAS
      IFUNM1         = t_params%IFUNM1
      LMSP1          = t_params%LMSP1
      NFU            = t_params%NFU
      LLMSP          = t_params%LLMSP
      LCORE          = t_params%LCORE
      NCORE          = t_params%NCORE
      NTCELL         = t_params%NTCELL
      IRMIN          = t_params%IRMIN
      ITITLE         = t_params%ITITLE
      INTERVX        = t_params%INTERVX
      INTERVY        = t_params%INTERVY
      INTERVZ        = t_params%INTERVZ
      NACLS1         = t_params%NACLS(1)
      LLY            = t_params%LLY
      SOCSCALE       = t_params%SOCSCALE
      TMPDIR         = t_params%TMPDIR
      ITMPDIR        = t_params%ITMPDIR
      ILTMP          = t_params%ILTMP
      NPAN_LOG       = t_params%NPAN_LOG_AT
      NPAN_EQ        = t_params%NPAN_EQ_AT
      NCHEB          = t_params%NCHEB
      R_LOG          = t_params%R_LOG
      NPAN_TOT       = t_params%NPAN_TOT
      RNEW           = t_params%RNEW
      RPAN_INTERVALL = t_params%RPAN_INTERVALL
      IPAN_INTERVALL = t_params%IPAN_INTERVALL
      THETASNEW      = t_params%THETASNEW
      IF (KREL.EQ.1) THEN
         SOLVER = t_params%SOLVER
         SOCSCL = t_params%SOCSCL
         CSCL   = t_params%CSCL
      ENDIF
      IF ( IDOLDAU.EQ.1 ) THEN
         NTLDAU   = t_params%NTLDAU
         ITLDAU   = t_params%ITLDAU
         LOPT     = t_params%LOPT
         UEFF     = t_params%UEFF
         JEFF     = t_params%JEFF
         EREFLDAU = t_params%EREFLDAU
      ENDIF
      !-------------------------------------------------------------------------
      ! Energy_mesh
      !-------------------------------------------------------------------------
      IELAST     = t_params%IELAST
      EZ         = t_params%EZ
      WEZ        = t_params%WEZ
      EMIN       = t_params%EMIN
      EMAX       = t_params%EMAX
      IESEMICORE = t_params%IESEMICORE
      NPOL       = t_params%NPOL
      TK         = t_params%TK
      IF ( NPOL.EQ.0 ) EFERMI = t_params%EFERMI
      !-------------------------------------------------------------------------
      ! Input_potential
      !-------------------------------------------------------------------------
      VINS  = t_params%VINS
      VISP  = t_params%VISP
      ECORE = t_params%ECORE
      IF (KREL.EQ.1) THEN
         RMREL     = t_params%RMREL
         DRDIREL   = t_params%DRDIREL
         R2DRDIREL = t_params%R2DRDIREL
         ZREL      = t_params%ZREL
         JWSREL    = t_params%JWSREL
         IRSHIFT   = t_params%IRSHIFT
         VTREL     = t_params%VTREL
         BTREL     = t_params%BTREL
      END IF
      IF ( TEST('Vspher  ') ) VINS(IRMIND:IRMD,2:LMPOTD,1:NSPOTD) = 0.D0
      !-------------------------------------------------------------------------
      ! Itermdir
      !-------------------------------------------------------------------------
      IF ( OPT('ITERMDIR') ) THEN
         QMTET = t_params%QMTET
         QMPHI = t_params%QMPHI
      END IF
      !-------------------------------------------------------------------------
      ! LDA+U
      !-------------------------------------------------------------------------
      IF ( IDOLDAU.EQ.1 ) THEN
         ITRUNLDAU = t_params%ITRUNLDAU
         WLDAU     = t_params%WLDAU
         ULDAU     = t_params%ULDAU
         PHILDAU   = t_params%PHILDAU
      END IF

   end subroutine get_params_1c

   !----------------------------------------------------------------------------
   ! SUBROUTINE: get_params_2
   !> @brief Set the values of the local variables according to the stored t_params
   !> so that they can be passed between different control modules, specifically for main2
   !> @author Philipp Rüssmann
   !----------------------------------------------------------------------------
   subroutine get_params_2(t_params,KREL,NATYP,IPAND,NPOTD,NATOMIMPD,LMXSPD,     &
      NFUND,LMPOT,NCELLD,IRMD,NEMBD1,NEMBD,IRMIND,NSRA,INS,NSPIN,IPAN,IRCUT,LCORE, &
      NCORE,LMAX,NTCELL,LPOT,NLBASIS,NRBASIS,NRIGHT,NLEFT,NATOMIMP,ATOMIMP,IMIX, &
      QBOUND,FCM,ITDBRY,IRNS,KPRE,KSHAPE,KTE,KVMAD,KXC,ICC,ISHIFT,IXIPOL,KFORCE, &
      IFUNM,LMSP,IMT,IRC,IRMIN,IRWS,LLMSP,ITITLE,NFU,HOSTIMP,ILM_MAP,IMAXSH,IELAST,  &
      NPOL,NPNT1,NPNT2,NPNT3,ITSCF,SCFSTEPS,IESEMICORE,KAOEZ,IQAT,NOQ,LLY,       &
      NPOLSEMI,N1SEMI,N2SEMI,N3SEMI,ZREL,JWSREL,IRSHIFT,MIXING,LAMBDA_XC,A,B,    &
      THETAS,DRDI,R,ZAT,RMT,RMTNEW,RWS,EMIN,EMAX,TK,ALAT,EFOLD,CHRGOLD,CMOMHOST, &
      CONC,GSH,EBOTSEMI,EMUSEMI,TKSEMI,VINS,VISP,RMREL,DRDIREL,VBC,FSOLD,        &
      R2DRDIREL,ECORE,EZ,WEZ,TXC,LINTERFACE,LRHOSYM,NGSHD,NAEZ,IRID,NSPOTD,IEMXD)
      ! get relevant parameters from t_params
      !     ..
      IMPLICIT NONE
      type(type_params), intent(in) :: t_params

      integer, intent(in) :: IRMD
      integer, intent(in) :: IRID
      integer, intent(in) :: NEMBD
      integer, intent(in) :: KREL
      integer, intent(in) :: IPAND
      integer, intent(in) :: NPOTD
      integer, intent(in) :: NFUND
      integer, intent(in) :: NGSHD
      integer, intent(in) :: IEMXD
      integer, intent(in) :: LMXSPD
      integer, intent(in) :: NCELLD
      integer, intent(in) :: NEMBD1
      integer, intent(in) :: IRMIND
      integer, intent(in) :: NSPOTD
      integer, intent(in) :: NATOMIMPD
      integer, intent(inout) :: KTE
      integer, intent(inout) :: KXC
      integer, intent(inout) :: ICC
      integer, intent(inout) :: LLY
      integer, intent(inout) :: INS
      integer, intent(inout) :: NSRA
      integer, intent(inout) :: IMIX
      integer, intent(inout) :: LPOT
      integer, intent(inout) :: KPRE
      integer, intent(inout) :: LMAX
      integer, intent(inout) :: NAEZ
      integer, intent(inout) :: NPOL
      integer, intent(inout) :: NPNT1
      integer, intent(inout) :: NPNT2
      integer, intent(inout) :: NPNT3
      integer, intent(inout) :: ITSCF
      integer, intent(inout) :: NSPIN
      integer, intent(inout) :: NATYP
      integer, intent(inout) :: KVMAD
      integer, intent(inout) :: LMPOT
      integer, intent(inout) :: NLEFT
      integer, intent(inout) :: NRIGHT
      integer, intent(inout) :: ISHIFT
      integer, intent(inout) :: KFORCE
      integer, intent(inout) :: N1SEMI
      integer, intent(inout) :: N2SEMI
      integer, intent(inout) :: N3SEMI
      integer, intent(inout) :: IELAST
      integer, intent(inout) :: KSHAPE
      integer, intent(inout) :: ITDBRY
      integer, intent(inout) :: NLBASIS
      integer, intent(inout) :: NRBASIS
      integer, intent(inout) :: NATOMIMP
      integer, intent(inout) :: SCFSTEPS
      integer, intent(inout) :: NPOLSEMI
      integer, intent(inout) :: IESEMICORE
      integer, dimension(NAEZ), intent(inout)      :: NOQ
      integer, dimension(NATYP), intent(inout)     :: IMT
      integer, dimension(NATYP), intent(inout)     :: IRC
      integer, dimension(NATYP), intent(inout)     :: NFU
      integer, dimension(NATYP), intent(inout)     :: IQAT
      integer, dimension(NATYP), intent(inout)     :: ZREL
      integer, dimension(NATYP), intent(inout)     :: IRWS
      integer, dimension(NATYP), intent(inout)     :: IPAN
      integer, dimension(NATYP), intent(inout)     :: IRNS
      integer, dimension(NPOTD), intent(inout)     :: NCORE
      integer, dimension(NATYP), intent(inout)     :: IRMIN
      integer, dimension(0:LMPOT), intent(inout)   :: IMAXSH
      integer, dimension(NATYP), intent(inout)     :: JWSREL
      integer, dimension(NATYP), intent(inout)     :: IXIPOL
      integer, dimension(NATYP), intent(inout)     :: NTCELL
      integer, dimension(0:NATYP), intent(inout)   :: HOSTIMP
      integer, dimension(NATYP), intent(inout)     :: IRSHIFT
      integer, dimension(NATOMIMPD), intent(inout) :: ATOMIMP
      integer, dimension(NGSHD,3), intent(inout)         :: ILM_MAP
      integer, dimension(NATYP,LMXSPD), intent(inout)    :: LMSP
      integer, dimension(20,NPOTD), intent(inout)        :: LCORE
      integer, dimension(NATYP,NFUND), intent(inout)     :: LLMSP
      integer, dimension(NATYP,LMXSPD), intent(inout)    :: IFUNM
      integer, dimension(NATYP,NAEZ+NEMBD), intent(inout) :: KAOEZ
      integer, dimension(0:IPAND,NATYP), intent(inout)   :: IRCUT
      integer, dimension(20,NPOTD), intent(inout)        :: ITITLE
      real (kind=dp), intent(inout) :: TK
      real (kind=dp), intent(inout) :: FCM
      real (kind=dp), intent(inout) :: EMIN
      real (kind=dp), intent(inout) :: EMAX
      real (kind=dp), intent(inout) :: ALAT
      real (kind=dp), intent(inout) :: EFOLD
      real (kind=dp), intent(inout) :: FSOLD
      real (kind=dp), intent(inout) :: QBOUND
      real (kind=dp), intent(inout) :: TKSEMI
      real (kind=dp), intent(inout) :: MIXING
      real (kind=dp), intent(inout) :: CHRGOLD
      real (kind=dp), intent(inout) :: EMUSEMI
      real (kind=dp), intent(inout) :: EBOTSEMI
      real (kind=dp), intent(inout) :: LAMBDA_XC
      real (kind=dp), dimension(NATYP), intent(inout)  :: A
      real (kind=dp), dimension(NATYP), intent(inout)  :: B
      real (kind=dp), dimension(2), intent(inout)      :: VBC
      real (kind=dp), dimension(NATYP), intent(inout)  :: RWS
      real (kind=dp), dimension(NGSHD), intent(inout)  :: GSH
      real (kind=dp), dimension(NATYP), intent(inout)  :: ZAT
      real (kind=dp), dimension(NATYP), intent(inout)  :: RMT
      real (kind=dp), dimension(NATYP), intent(inout)  :: CONC
      real (kind=dp), dimension(NATYP), intent(inout)  :: RMTNEW
      real (kind=dp), dimension(IRMD,NATYP), intent(inout)                :: R
      real (kind=dp), dimension(IRMD,NATYP), intent(inout)                :: DRDI
      real (kind=dp), dimension(IRMD,NPOTD), intent(inout)                :: VISP
      real (kind=dp), dimension(20,NPOTD), intent(inout)                 :: ECORE
      real (kind=dp), dimension(IRMD*KREL+(1-KREL),NATYP), intent(inout)  :: RMREL
      real (kind=dp), dimension(IRMD*KREL+(1-KREL),NATYP), intent(inout)  :: DRDIREL
      real (kind=dp), dimension(IRMD*KREL+(1-KREL),NATYP), intent(inout)  :: R2DRDIREL
      real (kind=dp), dimension(LMPOT,NEMBD1), intent(inout)             :: CMOMHOST
      real (kind=dp), dimension(IRMIND:IRMD,LMPOT,NSPOTD), intent(inout)  :: VINS
      real (kind=dp), dimension(IRID,NFUND,NCELLD), intent(inout)        :: THETAS
      complex (kind=dp), dimension(IEMXD), intent(inout) :: EZ
      complex (kind=dp), dimension(IEMXD), intent(inout) :: WEZ
      character(len=124), dimension(6), intent(inout) :: TXC
      logical, intent(inout) :: LRHOSYM
      logical, intent(inout) :: LINTERFACE
      !     .. External Functions ..
      LOGICAL OPT,TEST
      EXTERNAL OPT,TEST

      NSRA       = t_params%NSRA
      INS        = t_params%INS
      NATYP      = t_params%NATYP
      NAEZ       = t_params%NAEZ
      NSPIN      = t_params%NSPIN
      IPAN       = t_params%IPAN
      IRCUT      = t_params%IRCUT
      LCORE      = t_params%LCORE
      NCORE      = t_params%NCORE
      NTCELL     = t_params%NTCELL
      LMAX       = t_params%LMAX
      LPOT       = t_params%LPOT
      LMPOT      = t_params%LMPOT
      NLBASIS    = t_params%NLBASIS
      NRBASIS    = t_params%NRBASIS
      NRIGHT     = t_params%NRIGHT
      NLEFT      = t_params%NLEFT
      LINTERFACE = t_params%LINTERFACE
      ATOMIMP    = t_params%ATOMIMP
      NATOMIMP   = t_params%NATOMIMP

      !-------------------------------------------------------------------------
      ! Consistency check
      !-------------------------------------------------------------------------
      IF ( (KREL.EQ.1) .AND. (INS.NE.0) ) THEN
         WRITE(6,*) ' FULL-POTENTIAL RELATIVISTIC mode not implemented '
         STOP ' set INS = 0 in the input'
      END IF
      IF ( NSRA.LE.2 ) THEN
         IF ( KREL.EQ.1 ) STOP ' KVREL <= 1 in input, but relativistic program used'
      ELSE
         IF ( KREL.EQ.0 ) STOP ' KVREL > 1 in input, but non-relativistic program used'
      END IF
      !-------------------------------------------------------------------------
      IMIX      = t_params%IMIX
      MIXING    = t_params%MIXING
      QBOUND    = t_params%QBOUND
      FCM       = t_params%FCM
      ITDBRY    = t_params%ITDBRY
      IRNS      = t_params%IRNS
      KPRE      = t_params%KPRE
      KSHAPE    = t_params%KSHAPE
      KTE       = t_params%KTE
      KVMAD     = t_params%KVMAD
      KXC       = t_params%KXC
      LAMBDA_XC = t_params%LAMBDA_XC
      TXC       = t_params%TXC
      ICC       = t_params%ICC
      ISHIFT    = t_params%ISHIFT
      IXIPOL    = t_params%IXIPOL
      LRHOSYM   = t_params%LRHOSYM
      KFORCE    = t_params%KFORCE
      A         = t_params%A
      B         = t_params%B
      DRDI      = t_params%DRDI
      R         = t_params%RMESH
      THETAS    = t_params%THETAS
      ZAT       = t_params%ZAT
      IFUNM     = t_params%IFUNM
      LMSP      = t_params%LMSP
      RMT       = t_params%RMT
      RMTNEW    = t_params%RMTNEW
      RWS       = t_params%RWS
      IMT       = t_params%IMT
      IRC       = t_params%IRC
      IRMIN     = t_params%IRMIN
      IRWS      = t_params%IRWS
      ITITLE    = t_params%ITITLE
      LLMSP     = t_params%LLMSP
      NFU       = t_params%NFU
      HOSTIMP   = t_params%HOSTIMP
      ALAT      = t_params%ALAT
      KAOEZ     = t_params%KAOEZ
      IQAT      = t_params%IQAT
      NOQ       = t_params%NOQ
      CONC      = t_params%CONC
      GSH       = t_params%GSH
      ILM_MAP       = t_params%ILM_MAP
      IMAXSH    = t_params%IMAXSH
      LLY       = t_params%LLY
      !-------------------------------------------------------------------------
      ! Energy_mesh
      !-------------------------------------------------------------------------
      IELAST     = t_params%IELAST
      EZ         = t_params%EZ
      WEZ        = t_params%WEZ
      EMIN       = t_params%EMIN
      EMAX       = t_params%EMAX
      IESEMICORE = t_params%IESEMICORE
      FSOLD      = t_params%FSEMICORE
      NPOL       = t_params%NPOL
      TK         = t_params%TK
      NPNT1      = t_params%NPNT1
      NPNT2      = t_params%NPNT2
      NPNT3      = t_params%NPNT3
      EBOTSEMI   = t_params%EBOTSEMI
      EMUSEMI    = t_params%EMUSEMI
      TKSEMI     = t_params%TKSEMI
      NPOLSEMI   = t_params%NPOLSEMI
      N1SEMI     = t_params%N1SEMI
      N2SEMI     = t_params%N2SEMI
      N3SEMI     = t_params%N3SEMI
      !-------------------------------------------------------------------------
      ! Input_potential
      !-------------------------------------------------------------------------
      VINS  = t_params%VINS
      VISP  = t_params%VISP
      ECORE = t_params%ECORE
      VBC   = t_params%VBC
      IF (KREL.EQ.1) THEN
         RMREL     = t_params%RMREL
         DRDIREL   = t_params%DRDIREL
         R2DRDIREL = t_params%R2DRDIREL
         ZREL      = t_params%ZREL
         JWSREL    = t_params%JWSREL
         IRSHIFT   = t_params%IRSHIFT
      END IF
      ITSCF    = t_params%ITSCF
      SCFSTEPS = t_params%SCFSTEPS
      EFOLD    = t_params%EFOLD
      CHRGOLD  = t_params%CHRGOLD
      CMOMHOST = t_params%CMOMHOST
      IF ( TEST('Vspher  ') ) VINS(IRMIND:IRMD,2:LMPOT,1:NSPOTD) = 0.D0

   end subroutine get_params_2

   !----------------------------------------------------------------------------
   ! SUBROUTINE: save_emesh
   !> @brief Store the values of the local variables related to the energy mesh,
   !> in the t_params data types
   !> @author Philipp Rüssmann
   !----------------------------------------------------------------------------
   subroutine save_emesh(IELAST,EZ,WEZ,EMIN,EMAX,IESEMICORE,FSEMICORE,NPOL,TK,   &
      NPNT1,NPNT2,NPNT3,EBOTSEMI,EMUSEMI,TKSEMI,NPOLSEMI,N1SEMI,N2SEMI,N3SEMI,   &
      IEMXD,t_params)
      ! save information of energy mesh in t_params
      implicit none

      type(type_params), intent(inout) :: t_params

      integer, intent(in) :: NPOL
      integer, intent(in) :: NPNT1
      integer, intent(in) :: NPNT2
      integer, intent(in) :: NPNT3
      integer, intent(in) :: IEMXD
      integer, intent(in) :: IELAST
      integer, intent(in) :: N1SEMI
      integer, intent(in) :: N2SEMI
      integer, intent(in) :: N3SEMI
      integer, intent(in) :: NPOLSEMI
      integer, intent(in) :: IESEMICORE
      real (kind=dp), intent(in) :: TK
      real (kind=dp), intent(in) :: EMIN
      real (kind=dp), intent(in) :: EMAX
      real (kind=dp), intent(in) :: TKSEMI
      real (kind=dp), intent(in) :: EMUSEMI
      real (kind=dp), intent(in) :: EBOTSEMI
      real (kind=dp), intent(in) :: FSEMICORE
      complex (kind=dp), dimension(IEMXD), intent(in) :: EZ
      complex (kind=dp), dimension(IEMXD), intent(in) :: WEZ

      t_params%IELAST     = IELAST
      t_params%EZ         = EZ
      t_params%WEZ        = WEZ
      t_params%EMIN       = EMIN
      t_params%EMAX       = EMAX
      t_params%IESEMICORE = IESEMICORE
      t_params%FSEMICORE  = FSEMICORE
      t_params%NPOL       = NPOL
      t_params%TK         = TK
      t_params%NPNT1      = NPNT1
      t_params%NPNT2      = NPNT2
      t_params%NPNT3      = NPNT3
      t_params%EBOTSEMI   = EBOTSEMI
      t_params%EMUSEMI    = EMUSEMI
      t_params%TKSEMI     = TKSEMI
      t_params%NPOLSEMI   = NPOLSEMI
      t_params%N1SEMI     = N1SEMI
      t_params%N2SEMI     = N2SEMI
      t_params%N3SEMI     = N3SEMI

   end subroutine save_emesh

   !----------------------------------------------------------------------------
   ! SUBROUTINE: save_scfinfo
   !> @brief Store the values of the local variables related to the SCF parameters
   !> in the t_params data types
   !> @author Philipp Rüssmann
   !----------------------------------------------------------------------------
   subroutine save_scfinfo(t_params,VINS,VISP,ECORE,VBC,RMREL,DRDIREL,R2DRDIREL, &
      ZREL,JWSREL,IRSHIFT,VTREL,BTREL,ITSCF,SCFSTEPS,EFOLD,CHRGOLD,CMOMHOST,KREL,&
      IRMIND,IRM,LMPOT,NSPOTD,NATYP,NPOTD,NEMBD1)
      ! save information that is needed in next iteration and that is changeing, i.e. potential etc.
      implicit none

      type(type_params) , intent(inout) :: t_params

      integer, intent(in) :: IRM
      integer, intent(in) :: KREL
      integer, intent(in) :: LMPOT
      integer, intent(in) :: NATYP
      integer, intent(in) :: NPOTD
      integer, intent(in) :: ITSCF
      integer, intent(in) :: NEMBD1
      integer, intent(in) :: NSPOTD
      integer, intent(in) :: IRMIND
      integer, intent(in) :: SCFSTEPS
      integer, dimension(NATYP), intent(in) :: ZREL
      integer, dimension(NATYP), intent(in) :: JWSREL
      integer, dimension(NATYP), intent(in) :: IRSHIFT
      real (kind=dp), intent(in) :: EFOLD
      real (kind=dp), intent(in) :: CHRGOLD
      real (kind=dp), dimension(2), intent(in) :: VBC
      real (kind=dp), dimension(IRM,NPOTD), intent(in) :: VISP
      real (kind=dp), dimension(IRM*KREL+(1-KREL),NATYP), intent(in) :: VTREL
      real (kind=dp), dimension(IRM*KREL+(1-KREL),NATYP), intent(in) :: BTREL
      real (kind=dp), dimension(20,NPOTD), intent(in) :: ECORE
      real (kind=dp), dimension(IRM*KREL+(1-KREL),NATYP), intent(in) :: RMREL
      real (kind=dp), dimension(IRM*KREL+(1-KREL),NATYP), intent(in) :: DRDIREL
      real (kind=dp), dimension(LMPOT,NEMBD1), intent(in) :: CMOMHOST
      real (kind=dp), dimension(IRM*KREL+(1-KREL),NATYP), intent(in) :: R2DRDIREL
      real (kind=dp), dimension(IRMIND:IRM,LMPOT,NSPOTD), intent(in) :: VINS

      t_params%VINS  = VINS
      t_params%VISP  = VISP
      t_params%ECORE = ECORE
      t_params%VBC   = VBC
      IF (KREL.EQ.1) THEN
         t_params%RMREL     = RMREL
         t_params%DRDIREL   = DRDIREL
         t_params%R2DRDIREL = R2DRDIREL
         t_params%ZREL      = ZREL
         t_params%JWSREL    = JWSREL
         t_params%IRSHIFT   = IRSHIFT
         t_params%VTREL     = VTREL
         t_params%BTREL     = BTREL
      END IF
      t_params%ITSCF    = ITSCF
      t_params%SCFSTEPS = SCFSTEPS
      t_params%EFOLD    = EFOLD
      t_params%CHRGOLD  = CHRGOLD
      t_params%CMOMHOST = CMOMHOST
   end subroutine save_scfinfo

   !----------------------------------------------------------------------------
   ! SUBROUTINE: save_density
   !> @brief Store the values of the local variables related to the electronic density
   !> in the t_params data types
   !> @author Philipp Rüssmann
   !----------------------------------------------------------------------------
   subroutine save_density(t_params,RHO2NS,R2NEF,RHOC,DENEF,DENEFAT,ESPV,ECORE,  &
      IDOLDAU,LOPT,EU,EDC,CHRGSEMICORE,RHOORB,ECOREREL,NKCORE,KAPCORE,KREL,NATYP,&
      NPOTD,IRM,LMPOT,LMAXD1)
      ! save density after it has been calculated in main1c, is further processed in main2
      implicit none

      type(type_params) , intent(inout) :: t_params

      integer, intent(in) :: IRM
      integer, intent(in) :: KREL
      integer, intent(in) :: NATYP
      integer, intent(in) :: NPOTD
      integer, intent(in) :: LMPOT
      integer, intent(in) :: LMAXD1
      integer, intent(in) :: IDOLDAU
      integer, dimension(NATYP), intent(in) :: LOPT
      integer, dimension(20,NATYP), intent(in) :: NKCORE
      integer, dimension(20,NPOTD), intent(in) :: KAPCORE
      real (kind=dp), intent(in) :: DENEF
      real (kind=dp), intent(in) :: CHRGSEMICORE
      real (kind=dp), intent(in) :: EU(NATYP)
      real (kind=dp), intent(in) :: EDC(NATYP)
      real (kind=dp), intent(in) :: DENEFAT(NATYP)
      real (kind=dp), intent(in) :: RHOC(IRM,NPOTD)
      real (kind=dp), intent(in) :: ESPV(0:LMAXD1,NPOTD)
      real (kind=dp), intent(in) :: ECORE(20,NPOTD)
      real (kind=dp), intent(in) :: RHOORB(IRM*KREL+(1-KREL),NATYP)
      real (kind=dp), intent(in) :: ECOREREL(KREL*20+(1-KREL),NPOTD)
      real (kind=dp), intent(in) :: R2NEF(IRM,LMPOT,NATYP,2)
      real (kind=dp), intent(in) :: RHO2NS(IRM,LMPOT,NATYP,2)

      t_params%RHO2NS       = RHO2NS
      t_params%R2NEF        = R2NEF
      t_params%RHOC         = RHOC
      t_params%DENEF        = DENEF
      t_params%DENEFAT      = DENEFAT
      t_params%ESPV         = ESPV
      t_params%ECORE        = ECORE
      t_params%IDOLDAU      = IDOLDAU
      t_params%LOPT         = LOPT
      t_params%EU           = EU
      t_params%EDC          = EDC
      t_params%CHRGSEMICORE = CHRGSEMICORE
      IF (KREL.EQ.1) THEN
         t_params%RHOORB   = RHOORB
         t_params%ECOREREL = ECOREREL
         t_params%NKCORE   = NKCORE
         t_params%KAPCORE  = KAPCORE
      ENDIF

   end subroutine save_density

   !----------------------------------------------------------------------------
   ! SUBROUTINE: read_density
   !> @brief Store the values of the t_params data types related to the electronic
   !> density in local variables
   !> @author Philipp Rüssmann
   !----------------------------------------------------------------------------
   subroutine read_density(t_params,RHO2NS,R2NEF,RHOC,DENEF,DENEFAT,ESPV,ECORE,  &
      IDOLDAU,LOPT,EU,EDC,CHRGSEMICORE,RHOORB,ECOREREL,NKCORE,KAPCORE,KREL,NATYP,&
      NPOTD,IRM,LMPOT,LMAXD1)
      ! read density in main2
      implicit none

      type(type_params) , intent(inout) :: t_params

      integer, intent(in) :: IRM
      integer, intent(in) :: KREL
      integer, intent(in) :: NATYP
      integer, intent(in) :: NPOTD
      integer, intent(in) :: LMPOT
      integer, intent(in) :: LMAXD1
      integer, intent(out) :: IDOLDAU
      integer, dimension(NATYP), intent(out) :: LOPT
      integer, dimension(20,NATYP), intent(out) :: NKCORE
      integer, dimension(20,NPOTD), intent(out) :: KAPCORE
      real (kind=dp), intent(out) :: RHO2NS(IRM,LMPOT,NATYP,2)
      real (kind=dp), intent(out) :: DENEF
      real (kind=dp), intent(out) :: CHRGSEMICORE
      real (kind=dp), dimension(NATYP), intent(out) :: EU
      real (kind=dp), dimension(NATYP), intent(out) :: EDC
      real (kind=dp), dimension(NATYP), intent(out) :: DENEFAT
      real (kind=dp), dimension(IRM,NPOTD), intent(out) :: RHOC
      real (kind=dp), dimension(0:LMAXD1,NPOTD), intent(out) :: ESPV
      real (kind=dp), dimension(20,NPOTD), intent(out) :: ECORE
      real (kind=dp), dimension(IRM*KREL+(1-KREL),NATYP), intent(out) :: RHOORB
      real (kind=dp), dimension(KREL*20+(1-KREL),NPOTD), intent(out) :: ECOREREL
      real (kind=dp), dimension(IRM,LMPOT,NATYP,2), intent(out) :: R2NEF

      RHO2NS       = t_params%RHO2NS
      R2NEF        = t_params%R2NEF
      RHOC         = t_params%RHOC
      DENEF        = t_params%DENEF
      DENEFAT      = t_params%DENEFAT
      ESPV         = t_params%ESPV
      ECORE        = t_params%ECORE
      IDOLDAU      = t_params%IDOLDAU
      LOPT         = t_params%LOPT
      EU           = t_params%EU
      EDC          = t_params%EDC
      CHRGSEMICORE = t_params%CHRGSEMICORE

      IF (KREL.EQ.1) THEN
         RHOORB   = t_params%RHOORB
         ECOREREL = t_params%ECOREREL
         NKCORE   = t_params%NKCORE
         KAPCORE  = t_params%KAPCORE
      ENDIF

   end subroutine read_density

   !----------------------------------------------------------------------------
   ! SUBROUTINE: read_angles
   !> @brief Read the angles variables associated with the angles of magnetic
   !> moments in a non-collinear calcula
   !> @author Philipp Rüssmann
   !----------------------------------------------------------------------------
   subroutine read_angles(t_params,NATYP,THETA,PHI)
      ! read nonco_angles
      use mod_types, only: t_inc
      use mod_mympi, only: myrank, master
      use mod_version_info

      implicit none

      type(type_params), intent(inout) :: t_params

      integer, intent(in) :: NATYP
      real (kind=dp), dimension(NATYP), intent(out) :: THETA
      real (kind=dp), dimension(NATYP), intent(out) :: PHI

      logical :: LREAD, LCHECKANGLES
      integer :: I1,i_stat
      real (kind=dp) :: TH1, PH1
      real (kind=dp), parameter :: PI=4.d0*datan(1.d0), eps=1d-5

      ! if executed first in wunfiles theta is not allocated, thus read angles from file
      if(.not.allocated(t_params%THETA)) then

         THETA(:) = 0.D0
         PHI(:) = 0.D0
         LREAD = .FALSE.
         LCHECKANGLES = .false.
         INQUIRE(file='nonco_angle.dat',EXIST=LREAD)
         IF (LREAD) THEN
            OPEN(UNIT=10,FILE='nonco_angle.dat',FORM='FORMATTED')
            call version_check_header(10)
            DO I1 = 1,NATYP
               read(10,*) TH1, PH1
               if((ABS(TH1).LT.(PI+eps)   .AND. ABS(TH1).GT.eps).OR. (ABS(PH1).LT.(2*PI+eps) .AND. ABS(PH1).GT.eps)) then
                  LCHECKANGLES = .true.
               endif
               THETA(I1)=TH1*(PI/180.0D0)
               PHI(I1)  =PH1*(PI/180.0D0)
            ENDDO
            CLOSE(10)
            if ( LCHECKANGLES.and.((t_inc%i_write>0).or.(myrank==master)) ) then
               write(1337,*) 'WARNING: Check if your nonco_angels file is correct! Found only values that are smaller than pi for theta and 2pi for phi, respectively. But angles are given in degree (0...360)'
            endif
            write(1337,'(A)') '      I1  THETA[deg]  PHI[deg]'
            do I1=1,NATYP
               write(1337,'(I8,2F12.6)') I1,THETA(I1)*180D0/PI,PHI(I1)*180D0/PI
            end do!i1
         ENDIF ! LREAD

         !now save this also to t_params
         allocate(t_params%THETA(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(t_params%THETA))*kind(t_params%THETA),'t_params%THETA','read_angles')

         allocate(t_params%PHI(NATYP),stat=i_stat)
         call memocc(i_stat,product(shape(t_params%PHI))*kind(t_params%PHI),'t_params%PHI','read_angles')

         t_params%THETA = THETA
         t_params%PHI   = PHI

      else ! not first run: information saved in t_params

         THETA = t_params%THETA
         PHI   = t_params%PHI

      endif

   end subroutine read_angles

end module mod_wunfiles
