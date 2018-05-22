module rinput

   implicit none

contains
   !-------------------------------------------------------------------------------
   ! SUBROUTINE: RINPUT13
   !> @brief Routine to read the information from the input file
   !> @author Bernd Zimmermann
   !> @note Jonathan Chico: Added calls for the different parameters options,
   !> that are being changed from the inc.p to now the reader. Of this way if there
   !> is a problem the user only has to add a line or change a value of the inc.p
   !> not recompile the code
   !-------------------------------------------------------------------------------
   subroutine RINPUT13(NR,KTE,IGF,IRM,KXC,LLY,ICC,INS,KWS,IPE,IPF,IPFE,ICST,LM2D,   &
      IMIX,LPOT,NAEZ,NEMB,NREF,NCLS,NPOL,LMAX,KCOR,KEFG,KHYP,KPRE,NCLSD,MMAXD,NPOTD,&
      KVMAD,LMMAX,LMPOT,NCHEB,NLEFT,IFILE,KVREL,NSPIN,NATYP,NINEQ,NPNT1,NPNT2,NPNT3,&
      LMXSPD,LMMAXD,LMGF0D,LASSLD,NEMBD1,IRMIND,NOFGIJ,KFROZN,ISHIFT,N1SEMI,N2SEMI, &
      N3SEMI,NSTEPS,INSREF,KSHAPE,ITDBRY,NRIGHT,KFORCE,NSPINDD,NSATYPD,IVSHIFT,     &
      KHFIELD,NLBASIS,NRBASIS,INTERVX,INTERVY,INTERVZ,NPAN_EQ,NPAN_LOG,NPOLSEMI,TK, &
      FCM,EMIN,EMAX,RMAX,GMAX,ALAT,R_LOG,RCUTZ,RCUTXY,ESHIFT,QBOUND,HFIELD,MIXING,  &
      ABASIS,BBASIS,CBASIS,VCONST,TKSEMI,TOLRDIF,EMUSEMI,EBOTSEMI,FSEMICORE,        &
      LAMBDA_XC,DELTAE,LRHOSYM,LINIPOL,LCARTESIAN,LINTERFACE,IMT,CLS,LMXC,IRNS,IRWS,&
      NTCELL,REFPOT,INIPOL,IXIPOL,HOSTIMP,KFG,VBC,ZPERLEFT,ZPERIGHT,BRAVAIS,RMT,ZAT,&
      RWS,MTFAC,RMTREF,RMTNEW,RMTREFAT,FPRADIUS,TLEFT,TRIGHT,RBASIS,SOCSCALE,CSCL,  &
      SOCSCL,SOLVER,I12,I13,I19,I25,I40,TXC,DROTQ,NCPA,ITCPAMAX,CPATOL,NOQ,IQAT,    &
      ICPA,KAOEZ,CONC,KMROT,QMTET,QMPHI,KREADLDAU,LOPT,UEFF,JEFF,EREFLDAU,NTOTD)

      use Profiling
      use Constants
      use mod_wunfiles, only: t_params
      use memoryhandling
      use global_variables
      use mod_types, only: t_inc
      use mod_save_wavefun, only: t_wavefunctions
      use mod_version_info

      implicit none
      !     ..
      !> @note VP : there should be some crosscheck of competing options
      !>            e.g., XCPL and CONDUCT cannot be done simultaneously
      !>            neither SOC1 and SOC2 manipulation etc.
      !     ..
      !     .. External Subroutines ..
      EXTERNAL RCSTOP
      !     ..
      !     .. Intrinsic Functions ..
      INTRINSIC MIN
      !     ..
      !     .. Scalar Arguments ..
      integer, intent(inout) :: NR        !< Number of real space vectors rr
      integer, intent(inout) :: KTE       !< Calculation of the total energy On/Off (1/0)
      integer, intent(inout) :: IGF       !< Do not print or print (0/1) the KKRFLEX_* files
      integer, intent(inout) :: IRM       !< Maximum number of radial points
      integer, intent(inout) :: KXC       !< Type of xc-potential 0=vBH 1=MJW 2=VWN 3=PW91
      integer, intent(inout) :: LLY       !< LLY <> 0 : apply Lloyds formula
      integer, intent(inout) :: ICC       !< Enables the calculation of off-diagonal elements of the GF.(0=SCF/DOS; 1=cluster; -1=custom)
      integer, intent(inout) :: INS       !< 0 (MT), 1(ASA), 2(Full Potential)
      integer, intent(inout) :: KWS       !< 0 (MT), 1(ASA)
      integer, intent(inout) :: IPE       !< Not real used, IPFE should be 0
      integer, intent(inout) :: IPF       !< Not real used, IPFE should be 0
      integer, intent(inout) :: IPFE      !< Not real used, IPFE should be 0
      integer, intent(inout) :: ICST      !< Number of Born approximation
      integer, intent(inout) :: LM2D
      integer, intent(inout) :: IMIX      !< Type of mixing scheme used (0=straight, 4=Broyden 2nd, 5=Anderson)
      integer, intent(inout) :: LPOT      !< Maximum l component in potential expansion
      integer, intent(inout) :: NAEZ      !< Number of atoms in unit cell
      integer, intent(inout) :: NEMB      !< Number of 'embedding' positions
      integer, intent(inout) :: NREF      !< Number of diff. ref. potentials
      integer, intent(inout) :: NCLS      !< Number of reference clusters
      integer, intent(inout) :: NPOL      !< Number of Matsubara Pols (EMESHT)
      integer, intent(inout) :: LMAX      !< Maximum l component in wave function expansion
      integer, intent(inout) :: KCOR
      integer, intent(inout) :: KEFG
      integer, intent(inout) :: KHYP
      integer, intent(inout) :: KPRE
      integer, intent(inout) :: NCLSD     !< Maximum number of different TB-clusters
      integer, intent(inout) :: MMAXD
      integer, intent(inout) :: NPOTD
      integer, intent(inout) :: KVMAD
      integer, intent(inout) :: LMMAX
      integer, intent(inout) :: LMPOT
      integer, intent(inout) :: NTOTD     !< NTOTD = IPAND+30
      integer, intent(inout) :: NCHEB     !< Number of Chebychev pannels for the new solver
      integer, intent(inout) :: NLEFT     !< Number of repeated basis for left host to get converged  electrostatic potentials
      integer, intent(inout) :: IFILE     !< Unit specifier for potential card
      integer, intent(inout) :: KVREL     !< 0,1 : non / scalar relat. calculation
      integer, intent(inout) :: NSPIN     !< Counter for spin directions
      integer, intent(inout) :: NATYP     !< Number of kinds of atoms in unit cell
      integer, intent(inout) :: NINEQ     !< Number of ineq. positions in unit cell
      integer, intent(inout) :: NPNT1     !< number of E points (EMESHT) for the contour integration
      integer, intent(inout) :: NPNT2     !< number of E points (EMESHT) for the contour integration
      integer, intent(inout) :: NPNT3     !< number of E points (EMESHT) for the contour integration
      integer, intent(inout) :: LMXSPD    !< (2*LPOT+1)**2
      integer, intent(inout) :: LMMAXD
      integer, intent(inout) :: LMGF0D
      integer, intent(inout) :: LASSLD
      integer, intent(inout) :: NEMBD1
      integer, intent(inout) :: IRMIND
      integer, intent(inout) :: NOFGIJ
      integer, intent(inout) :: KFROZN
      integer, intent(inout) :: ISHIFT
      integer, intent(inout) :: N1SEMI    !< Number of energy points for the semicore contour
      integer, intent(inout) :: N2SEMI    !< Number of energy points for the semicore contour
      integer, intent(inout) :: N3SEMI    !< Number of energy points for the semicore contour
      integer, intent(inout) :: NSTEPS    !< number of iterations
      integer, intent(inout) :: INSREF    !< INS for reference pot. (usual 0)
      integer, intent(inout) :: KSHAPE    !< Exact treatment of WS cell
      integer, intent(inout) :: ITDBRY    !< Number of SCF steps to remember for the Broyden mixing
      integer, intent(inout) :: NRIGHT    !< Number of repeated basis for right host to get converged  electrostatic potentials
      integer, intent(inout) :: KFORCE    !< Calculation of the forces
      integer, intent(inout) :: NSPINDD   !< !< NSPIND-KORBIT
      integer, intent(inout) :: NSATYPD   !< (NATYPD-1)*KNOSPH+1
      integer, intent(inout) :: IVSHIFT
      integer, intent(inout) :: KHFIELD   !< 0,1: no / yes external magnetic field
      integer, intent(inout) :: NLBASIS   !< Number of basis layers of left host (repeated units)
      integer, intent(inout) :: NRBASIS   !< Number of basis layers of right host (repeated units)
      integer, intent(inout) :: INTERVX   !< Number of intervals in x-direction for k-net in IB of the BZ
      integer, intent(inout) :: INTERVY   !< Number of intervals in y-direction for k-net in IB of the BZ
      integer, intent(inout) :: INTERVZ   !< Number of intervals in z-direction for k-net in IB of the BZ
      integer, intent(inout) :: NPAN_EQ   !< Number of intervals from [R_LOG] to muffin-tin radius Used in conjunction with runopt NEWSOSOL
      integer, intent(inout) :: NPAN_LOG  !< Number of intervals from nucleus to [R_LOG] Used in conjunction with runopt NEWSOSOL
      integer, intent(inout) :: NPOLSEMI  !< Number of poles for the semicore contour
      double precision, intent(inout) :: TK      !< Temperature
      double precision, intent(inout) :: FCM     !< Factor for increased linear mixing of magnetic part of potential compared to non-magnetic part.
      double precision, intent(inout) :: EMIN    !< Lower value (in Ryd) for the energy contour
      double precision, intent(inout) :: EMAX    !< Maximum value (in Ryd) for the DOS calculation Controls also [NPT2] in some cases
      double precision, intent(inout) :: RMAX    !< Ewald summation cutoff parameter for real space summation
      double precision, intent(inout) :: GMAX    !< Ewald summation cutoff parameter for reciprocal space summation
      double precision, intent(inout) :: ALAT    !< Lattice constant (in a.u.)
      double precision, intent(inout) :: R_LOG   !< Radius up to which log-rule is used for interval width. Used in conjunction with runopt NEWSOSOL
      double precision, intent(inout) :: RCUTZ   !< Parameter for the screening cluster along the z-direction
      double precision, intent(inout) :: RCUTXY  !< Parameter for the screening cluster along the x-y plane
      double precision, intent(inout) :: ESHIFT
      double precision, intent(inout) :: QBOUND  !< Convergence parameter for the potential
      double precision, intent(inout) :: HFIELD  !< External magnetic field, for initial potential shift in spin polarised case
      double precision, intent(inout) :: MIXING  !< Magnitude of the mixing parameter
      double precision, intent(inout) :: ABASIS  !< Scaling factors for rbasis
      double precision, intent(inout) :: BBASIS  !< Scaling factors for rbasis
      double precision, intent(inout) :: CBASIS  !< Scaling factors for rbasis
      double precision, intent(inout) :: VCONST  !< Potential shift in the first iteration
      double precision, intent(inout) :: TKSEMI  !< Temperature for semi-core contour
      double precision, intent(inout) :: TOLRDIF !< For distance between scattering-centers smaller than [<TOLRDIF>], free GF is set to zero. Units are Bohr radii.
      double precision, intent(inout) :: EMUSEMI   !< Top of semicore contour in Ryd.
      double precision, intent(inout) :: EBOTSEMI  !< Bottom of semicore contour in Ryd
      double precision, intent(inout) :: FSEMICORE !< Initial normalization factor for semicore states (approx. 1.)
      double precision, intent(inout) :: LAMBDA_XC !< Scale magnetic moment (0 < Lambda_XC < 1,0=zero moment, 1= full moment)
      double complex, intent(inout) :: DELTAE      !< LLY Energy difference for numerical derivative
      logical, intent(inout) :: LRHOSYM
      logical, intent(inout) :: LINIPOL    !< True: Initial spin polarization; false: no initial spin polarization
      logical, intent(inout) :: LCARTESIAN !< True: Basis in cartesian coords; false: in internal coords
      logical, intent(inout) :: LINTERFACE !< If True a matching with semi-inifinite surfaces must be performed
      !     .. Array Arguments ..
      integer, dimension(:), allocatable, intent(out) :: IMT    !< R point at MT radius
      integer, dimension(:), allocatable, intent(out) :: CLS    !< Cluster around atomic sites
      integer, dimension(:), allocatable, intent(out) :: LMXC
      integer, dimension(:), allocatable, intent(out) :: IRNS   !< Position of atoms in the unit cell in units of bravais vectors
      integer, dimension(:), allocatable, intent(out) :: IRWS   !< R point at WS radius
      integer, dimension(:), allocatable, intent(out) :: NTCELL !< Index for WS cell
      integer, dimension(:), allocatable, intent(out) :: REFPOT !< Ref. pot. card  at position
      integer, dimension(:), allocatable, intent(out) :: INIPOL !< Initial spin polarisation
      integer, dimension(:), allocatable, intent(out) :: IXIPOL !< Constraint of spin pol.
      integer, dimension(:), allocatable, intent(out) :: HOSTIMP
      integer, dimension(:,:), allocatable, intent(out) :: KFG
      double precision, dimension(2), intent(inout) :: VBC        !< Potential constants
      double precision, dimension(3), intent(inout) :: ZPERLEFT   !< Vector to define how to repeat the basis of the left host
      double precision, dimension(3), intent(inout) :: ZPERIGHT   !< Vector to define how to repeat the basis of the right host
      double precision, dimension(3,3), intent(inout) :: BRAVAIS  !< Bravais lattice vectors
      double precision, dimension(:), allocatable, intent(out) :: RMT      !< Muffin-tin radius of true system
      double precision, dimension(:), allocatable, intent(out) :: ZAT      !< Nuclear charge
      double precision, dimension(:), allocatable, intent(out) :: RWS      !< Wigner Seitz radius
      double precision, dimension(:), allocatable, intent(out) :: MTFAC    !< Scaling factor for radius MT
      double precision, dimension(:), allocatable, intent(out) :: RMTREF   !< Muffin-tin radius of reference system
      double precision, dimension(:), allocatable, intent(out) :: RMTNEW   !< Adapted muffin-tin radius
      double precision, dimension(:), allocatable, intent(out) :: RMTREFAT
      double precision, dimension(:), allocatable, intent(out) :: FPRADIUS !< R point at which full-potential treatment starts
      double precision, dimension(:,:), allocatable, intent(out) :: TLEFT  !< Vectors of the basis for the left host
      double precision, dimension(:,:), allocatable, intent(out) :: TRIGHT !< vectors of the basis for the right host
      double precision, dimension(:,:), allocatable, intent(out) :: RBASIS !< Position of atoms in the unit cell in units of bravais vectors
      !     variables for spin-orbit/speed of light scaling
      double precision, dimension(:), allocatable, intent(out) :: SOCSCALE !< Spin-orbit scaling
      double precision, dimension(:,:), allocatable, intent(out) :: CSCL   !< Speed of light scaling
      double precision, dimension(:,:), allocatable, intent(out) :: SOCSCL
      character(len=10), intent(inout) :: SOLVER !< Type of solver
      character(len=40), intent(inout) :: I12 !< File identifiers
      character(len=40), intent(inout) :: I13 !< Potential file name
      character(len=40), intent(inout) :: I19 !< Shape function file name
      character(len=40), intent(inout) :: I25 !< Scoef file name
      character(len=40), intent(inout) :: I40 !< File identifiers
      character(len=124), dimension(6), intent(inout) :: TXC
      double complex, dimension(:,:,:), allocatable, intent(out) :: DROTQ !< Rotation matrices to change between LOCAL/GLOBAL frame of reference for magnetisation <> Oz or noncollinearity
      !----------------------------------------------------------------------------
      !> @note CPA variables. Routine has been modified to look for
      !>     the token ATOMINFOC and only afterwards, if not found, for the
      !>     old token ATOMINFO. The only necessary extra information
      !>     required is the site IQAT(IATOM) on which the atom IATOM
      !>     is located and the occupancy (concentration) CONC(IATOM).
      !>     The rest of CPA variables are deduced from these two.
      !>     The tolerance for the CPA-cycle and the number of CPA iterations
      !>     can be modified adding the token <CPAINFO> in the input file.
      !----------------------------------------------------------------------------
      integer, intent(inout) :: NCPA             !< ncpa = 0/1 CPA flag
      integer, intent(inout) :: ITCPAMAX         !< max. number of CPA iterations
      double precision, intent(inout)  :: CPATOL !< convergency tolerance for CPA-cycle
      integer, dimension(:), allocatable, intent(out) :: NOQ  !< number of diff. atom types located
      integer, dimension(:), allocatable, intent(out) :: IQAT !< the site on which an atom is located on a given site
      integer, dimension(:), allocatable, intent(out) :: ICPA !< icpa = 0/1 site-dependent CPA flag
      integer, dimension(:,:), allocatable, intent(out) :: KAOEZ !< atom types located at a given site
      double precision, dimension(:), allocatable, intent(out) :: CONC !< concentration of a given atom

      !----------------------------------------------------------------------------
      !> @note Variables storing the magnetization direction information.
      !>     QMTET/QMPHI(NAEZ) give the angles to which the magnetic moment
      !>     on a given site is rotated against the z-axis. Default values
      !>     0.0 and 0.0, i.e., magnetic moment parallel to the z-axis.
      !>     The angles are read in after the token RBASISANG is found
      !>     (sought in input file prior to RBASIS token)
      !>
      !>   *  KMROT                                                           *
      !>   *  0: no rotation of the magnetisation                             *
      !>   *  1: individual rotation of the magnetisation for every site      *
      !>   ( see also the routine < FINDGROUP > and ff)
      !----------------------------------------------------------------------------
      integer, intent(inout) :: KMROT !< 0: no rotation of the magnetisation; 1: individual rotation of the magnetisation for every site
      double precision, dimension(:), allocatable, intent(out) :: QMTET !< \f$ \theta\f$ angle of the agnetization with respect to the z-axis
      double precision, dimension(:), allocatable, intent(out) :: QMPHI !< \f$ \phi\f$ angle of the agnetization with respect to the z-axis
      ! ---------------------------------------------------------------------------
      ! LDA+U
      integer, intent(inout) :: KREADLDAU !< LDA+U arrays available
      integer, dimension(:), allocatable, intent(inout) :: LOPT !< angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
      double precision, dimension(:), allocatable, intent(out) :: UEFF !< input U parameter for each atom
      double precision, dimension(:), allocatable, intent(out) :: JEFF !< input J parameter for each atom
      double precision, dimension(:), allocatable, intent(out) :: EREFLDAU !< the energies of the projector's wave functions (REAL)
      ! LDA+U
      ! ---------------------------------------------------------------------------

      !----------------------------------------------------------------------------
      ! Local variables
      !----------------------------------------------------------------------------
      ! for OPERATOR option
      logical :: lexist, operator_imp
      ! IVSHIFT test option
      logical :: TEST,OPT
      EXTERNAL TEST,OPT
      !     ..
      !     .. Local Scalars ..
      integer :: NDIM !< Dimension for the Bravais lattice for slab or bulk (2/3)
      integer :: NASOC
      integer :: I,IL,J,IER,IER2,I1,II,IR,IDOSEMICORE,i_stat,i_all
      double precision :: SOSCALE,CTLSCALE
      double precision :: BRYMIX,STRMIX,TX,TY,TZ
      character(len=43) :: TSHAPE
      character(len=256) :: UIO  ! NCOLIO=256
      logical :: LNEW !< Logical variable for old/new treatment of left and right host
      logical :: MANSOC
      logical :: MANCTL
      logical :: LATOMINFO !< Logical variable for old/new treatment of the ATOMINFO
      !.. Local CPA variables
      integer :: IO,IA,IQ,IPRINT
      double precision :: SUM
      character(len=3), dimension(0:1) :: CPAFLAG

      !     .. Local Arrays ..
      integer, dimension(:), allocatable :: ISP
      integer, dimension(:), allocatable :: IMANSOC
      double precision, dimension(10) :: DVEC
      character(len=4), dimension(3) :: TSPIN
      character(len=8), dimension(3) :: TKWS
      character(len=2), dimension(-2:-1) :: SOCII
      character(len=43), dimension(0:3) :: TINS
      character(len=43), dimension(0:3) :: TKCOR
      character(len=43), dimension(0:2) :: TVREL
      !     ..
      !     .. Data statements ..
      DATA TSPIN/'non-','    ','    '/
      DATA TSHAPE/' exact cell treatment (shape correction)  '/
      DATA TVREL/&
      &     ' non relativistic calculation              ',&
      &     ' s.r.a. calculation                        ',&
      &     ' fully relativistic calculation            '/
      DATA TKCOR/&
      &     ' frozen core approximation                 ',&
      &     ' core relaxation s.r.a.                    ',&
      &     ' core relaxation nonsra                    ',&
      &     ' core relaxation                           '/
      DATA TINS/' spherical averaged input potential        ',&
      &     ' non spherical input potential for cluster ',&
      &     ' non spherical input potential for cluster ',&
      &     ' non spherical input potential             '/
      DATA TKWS/' full mt','   ws   ',' full ws'/
      !
      DATA CPAFLAG/' NO','YES'/
      DATA SOCII/'xy','zz'/
      !     ..
      !
      !------------ array set up and definition of input parameter -----------
      !
      ! concatenate name & serial number
      TXC(1) = ' Morruzi,Janak,Williams  #serial: ' // serialnr
      TXC(2) = ' von Barth,Hedin         #serial: ' // serialnr
      TXC(3) = ' Vosko,Wilk,Nusair       #serial: ' // serialnr
      TXC(4) = ' GGA PW91                #serial: ' // serialnr
      TXC(5) = ' GGA PBE                 #serial: ' // serialnr
      TXC(6) = ' GGA PBEsol              #serial: ' // serialnr

      IPRINT = 0

      open(111,FILE='inputcard_generated.txt') ! Write out found or assumed values
      call version_print_header(111)

      NEMB = 0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Read RUNNING options
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call IoInput('RUNOPT          ',UIO,1,7,IER)
      if (IER.NE.0) then
         write(111,*) 'RUNOPT not found'
      else
         read (UNIT=UIO,FMT=980)(t_params%OPTC(I),I=1,8)
         write(111,FMT='(A6)') 'RUNOPT'
         write(111,FMT=980)  (t_params%OPTC(I),I=1,8)
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Read TEST options
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call IoInput('TESTOPT         ',UIO,1,7,IER)
      if (IER.NE.0) then
         write(111,*) 'TESTOPT not found'
      else
         read(UNIT=UIO,FMT=980)(t_params%TESTC(i),i=1,8)
         call IoInput('TESTOPT         ',UIO,2,7,IER)
         read(UNIT=UIO,FMT=980)(t_params%TESTC(8+i),i=1,8)
         write(111,FMT='(A7)') 'TESTOPT'
         write(111,FMT=980)(t_params%TESTC(i),i=1,8)
         write(111,FMT=980)(t_params%TESTC(8+i),i=1,8)
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Begin lattice structure definition
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call IoInput('ALATBASIS       ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) ALAT
         write(111,*) 'ALATBASIS=',ALAT
      else
         write(111,*) 'ALATBASIS not found in inputcard'
         write(*,*) 'rinput13: ALATBASIS not found in inputcard'
         stop 'rinput13: ALATBASIS not found in inputcard'
      endif

      ! Set 2-d or 3-d geometry
      LINTERFACE = .FALSE.
      call IoInput('INTERFACE       ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) LINTERFACE
         write(111,*) 'INTERFACE=',LINTERFACE
      else
         write(111,*) 'Default INTERFACE= ',LINTERFACE
      endif

      NDIM = 3
      if (LINTERFACE) NDIM = 2
      if (.NOT.LINTERFACE.AND..NOT.OPT('SUPRCELL')) then
         write(1337,*) '3D-calculation, adding run-option "full inv" for full inversion.'
         call ADDOPT('full inv')
      endif

      if (OPT('WRTGREEN')) then
         write(1337,*) 'WRTGREEN option found'
         write(1337,*) 'adding run-opt "full inv" for full inversion.'
         write(1337,*) 'adding run-opt "fix mesh"'
         call ADDOPT('full inv')
         call ADDOPT('fix mesh')
      end if

      write(111,*) 'Bravais vectors in units of ALAT'
      BRAVAIS(1:3,1:3) = 0D0
      do I = 1,NDIM
         call IOINPUT('BRAVAIS         ',UIO,I,7,IER)
         if (IER.NE.0) stop 'RINPUT: BRAVAIS NOT FOUND'
         read (UNIT=UIO,FMT=*) (BRAVAIS(J,I),J=1,NDIM)
      end do
      write(111,FMT='(A7)') 'BRAVAIS'
      do I = 1,NDIM
         write(111,*) (BRAVAIS(J,I),J=1,NDIM)
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Read the number of atoms in the unit cell
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call IoInput('NAEZ            ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) NAEZ
         write(111,*) 'NAEZ=',NAEZ
      else
         write(111,*) 'NAEZ not found'
         stop 'NAEZ not found in <RINPUT13>'
      endif
      !if (NAEZ.GT.NAEZD) then
      !   write(6,*) ' set NAEZD to at least ',NAEZ
      !   stop ' in < RINPUT13 > '
      !end if

      NREF=NAEZ
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Read the atom types, if no CPA NATYP=NAEZ
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      NATYP = NAEZ
      call IoInput('NATYP           ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) NATYP
         write(111,*) 'NATYP=',NATYP
      else
         write(111,*) 'Default NATYP= ',NAEZ
      endif
      !if (NATYP.GT.NATYPD) then
      !   write(6,*) 'RINPUT13: NATYP > NATYPD',NATYP,NATYPD
      !   stop ' IN < RINPUT13 > '
      !end if
      if (NATYP.LT.NAEZ) then
         write(6,*) 'RINPUT13: NATYP < NAEZ ',NATYP,NAEZ
         stop ' IN < RINPUT13 > '
      end if

      allocate(ISP(NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(ISP))*kind(ISP),'ISP','rinput13')
      ISP=0

      LCARTESIAN = .FALSE.
      IER = 0
      call IOINPUT('CARTESIAN       ',UIO,1,7,IER)
      if ( IER.EQ.0 ) then
         read (UNIT=UIO,FMT=*) LCARTESIAN
         write(111,*) 'CARTESIAN= ',LCARTESIAN
      else
         write(111,*) 'Default CARTESIAN= ',LCARTESIAN
      endif

      ! Jonathan Chico: This call needs to be done before the rest as one needs to
      ! find out the value of NEMB to be able to allocate several arrays
      if (LINTERFACE) then
         write(1337,9410)

         NRIGHT = 10
         call IoInput('NRIGHTHO        ',UIO,1,7,IER)
         if (IER.EQ.0) then
            read (UNIT=UIO,FMT=*) NRIGHT
            write(111,*) 'NRIGHTHO=',NRIGHT
         else
            write(111,*) 'Default NRIGHTHO=',NRIGHT
         endif

         NLEFT = 10
         call IoInput('NLEFTHOS        ',UIO,1,7,IER)
         if (IER.EQ.0) then
            read (UNIT=UIO,FMT=*) NLEFT
            write(111,*) 'NLEFTHOS=',NLEFT
         else
            write(111,*) 'Default NLEFTHOS=',NLEFT
         endif

         call IoInput('<NLBASIS>       ',UIO,1,7,IER)
         if (IER.NE.0) then
            write(1337,*) 'rinput13: <NLBASIS> not found in inputcard'
            IER = 0
            call IoInput('NLBASIS         ',UIO,1,7,IER)
            if (IER.NE.0) then
               write(*,*) 'rinput13: NLBASIS also not found in inputcard'
               stop 'rinput13: NLBASIS not found in inputcard'
            endif
         endif
         if (IER.EQ.0) then
            read (UNIT=UIO,FMT=*) NLBASIS
            write(111,*) '<NLBASIS>=',NLBASIS
         endif

         call IoInput('<NRBASIS>       ',UIO,1,7,IER)
         if (IER.NE.0) then
            write(1337,*) 'rinput13: <NRBASIS> not found in inputcard'
            IER = 0
            call IoInput('NRBASIS         ',UIO,1,7,IER)
            if (IER.NE.0) then
               write(*,*) 'rinput13: NRBASIS also not found in inputcard'
               stop 'rinput13: NRBASIS not found in inputcard'
            endif
         endif
         if (IER.EQ.0) then
            read (UNIT=UIO,FMT=*) NRBASIS
            write(111,*) '<NRBASIS>=',NRBASIS
         endif

         NEMB = NLBASIS + NRBASIS
         write(1337,*) 'Number of embedded atoms NEMB=NLBASIS + NRBASIS=',NEMB
         !if(NEMB.GT.NEMBD) then
         !   write(6,*) 'Please, increase the parameter nembd (',nembd,') in inc.p to',nemb
         !   stop 'ERROR in NEMBD.'
         !endif

         IER = 0
         ! Check if the keywords exist for old/new treatment of left and right host
         call IoInput('LEFTBASIS       ',UIO,1,7,IER)
         if (IER.EQ.0) then
            LNEW = .FALSE.
         else
            LNEW = .TRUE.
            IER = 0
            call IoInput('<RBLEFT>        ',UIO,1,7,IER)
         endif
         if (IER.NE.0) then
            write(*,*) 'rinput13: LEFTBASIS or <RBLEFT> not found in inputcard'
            stop 'rinput13: LEFTBASIS or <RBLEFT> not found in inputcard'
         endif
         IER = 0
         call IoInput('RIGHBASIS       ',UIO,1,7,IER)
         if (IER.EQ.0) then
            LNEW = .FALSE.
         else
            LNEW = .TRUE.
            IER = 0
            call IoInput('<RBRIGHT>       ',UIO,1,7,IER)
         endif
         if (IER.NE.0) then
            write(*,*) 'rinput13: RIGHBASIS or <RBRIGHT> not found in inputcard'
            stop 'rinput13: RIGHBASIS or <RBRIGHT> not found in inputcard'
         endif
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Allocate the unit cell arrays
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call allocate_cell(1,NAEZ,NEMB,NATYP,CLS,IMT,IRWS,IRNS,NTCELL,REFPOT,&
      KFG,KAOEZ,RMT,ZAT,RWS,MTFAC,RMTREF,RMTREFAT,RMTNEW,RBASIS,LMXC)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End of allocation of the unit cell arrays
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Allocate the right and left hosts for slab calculation
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call allocate_semi_inf_host(1,NEMB,TLEFT,TRIGHT)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End of allocation of the right and left hosts for slab calculation
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Basis atoms
      write(111,FMT='(A16)') '<RBASIS>        '
      do I=1,NAEZ
         call IoInput('<RBASIS>        ',UIO,I,7,IER)
         if (IER.EQ.0) then
            read (UNIT=UIO,FMT=*) (RBASIS(J,I), J=1,3)
            write(111,FMT='(3E24.12)') (RBASIS(J,I), J=1,3)
         else
            IER=0
            call IoInput('RBASIS          ',UIO,I,7,IER)
            if (IER.EQ.0) then
               read (UNIT=UIO,FMT=*) (RBASIS(J,I), J=1,3)
               write(111,FMT='(3E24.12)') (RBASIS(J,I), J=1,3)
            else
               write(*,*) 'RINPUT13: Keyword <RBASIS> or RBASIS not found. Stopping.'
               stop 'RINPUT13: RBASIS'
            endif
         endif
      enddo                         ! I=1,NAEZ
      call IDREALS(RBASIS(1,1),3*NAEZ,IPRINT)

      DVEC(1:3) = 1.D0
      call IoInput('BASISCALE       ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) (DVEC(I),I=1,3)
         write(111,FMT='(A10,3E12.4)') 'BASISCALE=',DVEC(1:3)
      else
         write(111,FMT='(A18,3E12.4)') 'Default BASISCALE=',DVEC(1:3)
      endif

      call IDREALS(DVEC(1),3,IPRINT)
      ABASIS = DVEC(1)
      BBASIS = DVEC(2)
      CBASIS = DVEC(3)

      write(1337,2019) ABASIS,BBASIS,CBASIS
      write(1337,2107)
      write(1337,2014) ALAT

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Begin read left- and right-host information in 2d-case.
      ! Set up the embeding positions
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (LINTERFACE) then
         !-------------------------------------------------------------------------
         !> @note In leftbasis and rightbasis, kaoez is used only in decimation case.
         !> Then it indicates the correspondence of the atom-coordinate given
         !> by leftbasis and rightbasis to the left- and right-host t-matrix read in
         !> by decimaread. For the slab case, kaoez is not used in the embedded positions.
         !-------------------------------------------------------------------------
         if (LNEW) then

            write(111,FMT='(A82)') '<RBLEFT>                                                      <RMTREFL>   <KAOEZL>'
            do I=1,NLBASIS
               call IoInput('<RBLEFT>        ',UIO,I,7,IER)
               read (UNIT=UIO,FMT=*) (TLEFT(I1,I),I1=1,3)
               KAOEZ(1,NAEZ+I) = I            ! Default
               call IoInput('<KAOEZL>        ',UIO,I,7,IER)
               if (IER.EQ.0) read (UNIT=UIO,FMT=*) KAOEZ(1,NAEZ+I)
               call IoInput('<RMTREFL>       ',UIO,I,7,IER)
               if (IER.EQ.0) read (UNIT=UIO,FMT=*) RMTREFAT(NAEZ+I)
               write (111,FMT='(3E20.12,3X,F9.6,3X,I5)') (TLEFT(I1,I),I1=1,3),RMTREFAT(NAEZ+I),KAOEZ(1,NAEZ+I)
            enddo
            write(111,FMT='(A82)') '<RBRIGHT>                                                     <RMTREFR>   <KAOEZL>'
            do I=1,NRBASIS
               call IoInput('<RBRIGHT>       ',UIO,I,7,IER)
               read (UNIT=UIO,FMT=*) (TRIGHT(I1,I),I1=1,3)
               KAOEZ(1,NAEZ+NLBASIS+I) = I     ! Default
               call IoInput('<KAOEZR>        ',UIO,I,7,IER)
               if (IER.EQ.0) read (UNIT=UIO,FMT=*) KAOEZ(1,NAEZ+NLBASIS+I)
               call IoInput('<RMTREFR>       ',UIO,I,7,IER)
               if (IER.EQ.0) read (UNIT=UIO,FMT=*) RMTREFAT(NAEZ+NLBASIS+I)
               write (111,FMT='(3E20.12,3X,F9.6,3X,I5)') (TRIGHT(I1,I),I1=1,3),RMTREFAT(NAEZ+NLBASIS+I),KAOEZ(1,NAEZ+NLBASIS+I)
            enddo

         else ! (LNEW) now old-style input

            do I=1,NLBASIS
               call IoInput('LEFTBASIS       ',UIO,I,7,IER)
               read (UNIT=UIO,FMT=*) (TLEFT(I1,I),I1=1,3),II,IR
               KAOEZ(1,NAEZ+I) = II    ! changed 1.11.99
               REFPOT(NAEZ+I) = IR
            end do
            do I=1,NRBASIS
               call IoInput('RIGHBASIS       ',UIO,I,7,IER)
               read (UNIT=UIO,FMT=*) (TRIGHT(I1,I),I1=1,3),II,IR
               KAOEZ(1,NAEZ+NLBASIS+I) = II  ! changed 1.11.99
               REFPOT(NAEZ+NLBASIS+I) = IR
            end do
         endif

         call IDREALS(TLEFT,3*(NEMB+1),IPRINT)
         call IDREALS(TRIGHT,3*(NEMB+1),IPRINT)


         ! Put The additional atoms in the "embeding" positions

         do I=1,NLBASIS
            RBASIS(1:3,NAEZ+I) = TLEFT(1:3,I)
         end do
         do I=1,NRBASIS
            RBASIS(1:3,NAEZ+NLBASIS+I) = TRIGHT(1:3,I)
         end do
         !-------------------------------------------------------------------------
         ! In RBASIS we have first the basis atoms or the interface
         ! atoms then the left host then the right host the host
         ! goes in the NEMB positions
         !
         ! IN CASE OF CPA the host is treated as an effective
         ! CPA medium, that is, there is only one kind of atom
         ! occupying a crystallographic site.
         !
         !-------------------------------------------------------------------------
         call IoInput('ZPERIODL        ',UIO,1,7,IER)
         if (IER.NE.0) then
            write(*,*) 'rimput13: ZPERIODL not found in inputcard'
            stop 'rimput13: ZPERIODL not found in inputcard'
         else
            read (UNIT=UIO,FMT=*) (ZPERLEFT(I1),I1=1,3)
            write(111,FMT='(A9,3E20.12)') 'ZPERIODL=',(ZPERLEFT(I1),I1=1,3)
         endif
         call IDREALS(ZPERLEFT(1),3,IPRINT)

         call IoInput('ZPERIODR        ',UIO,1,7,IER)
         if (IER.NE.0) then
            write(*,*) 'rimput13: ZPERIODR not found in inputcard'
            stop 'rimput13: ZPERIODR not found in inputcard'
         else
            read (UNIT=UIO,FMT=*) (ZPERIGHT(I1),I1=1,3)
            write(111,FMT='(A9,3E20.12)') 'ZPERIODR=',(ZPERIGHT(I1),I1=1,3)
         endif
         call IDREALS(ZPERIGHT(1),3,IPRINT)

         write(1337,9430) NLEFT,NLBASIS
         write(1337,9440) NRIGHT,NRBASIS
         write(1337,9450) (ZPERLEFT(i1),I1=1,3)
         write(1337,9460) (ZPERIGHT(i1),I1=1,3)
         write(1337,9465)
         write(1337,9470)
         do I=NLEFT,1,-1
            do I1=NLBASIS,1,-1
               tx = TLEFT(1,i1) + (I-1)*ZPERLEFT(1)
               ty = TLEFT(2,i1) + (I-1)*ZPERLEFT(2)
               tz = TLEFT(3,i1) + (I-1)*ZPERLEFT(3)
               write(1337,9420) (I-1)*NLBASIS+i1, tx,ty,tz,KAOEZ(1,I1)
            end do
         end do
         write(1337,9475)
         do I=1,NAEZ
            write(1337,9420) I, (RBASIS(I1,I),I1=1,3)
         end do
         write(1337,9480)
         do I=1,NRIGHT
            do I1=1,NRBASIS
               tx = TRIGHT(1,i1) + (I-1)*ZPERIGHT(1)
               ty = TRIGHT(2,i1) + (I-1)*ZPERIGHT(2)
               tz = TRIGHT(3,i1) + (I-1)*ZPERIGHT(3)
               write(1337,9420) (I-1)*NRBASIS+i1,tx,ty,tz,KAOEZ(1,I1)
            end do
         end do

      end if ! LINTERFACE
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End read left- and right-host information in 2d-case.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Although NSPIN is fixed to 1 in REL mode,
      ! NSPIN should be used as 1 or 2 at this stage
      ! to indicate a non- or spin-polarised potential
      ! that has to be read in. NSPIN is set to 1 before
      ! being passed to the subsequent programs.
      ! < TESTDIM > has been accordingly modified
      call IoInput('NSPIN           ',UIO,1,7,IER)
      if (IER.NE.0) then
         write(111,*) 'NSPIN not found'
         stop 'NSPIN not found'
      else
         read (UNIT=UIO,FMT=*) NSPIN
         write(111,*) 'NSPIN=',NSPIN
      endif

      write(1337,2010) NSPIN
      write(1337,2104)

      ! Atomic number
      call IoInput('<ZATOM>         ',UIO,1,7,IER)
      if (IER.EQ.0) then
         write(111,'(A10)') '<ZATOM>   '
         do I = 1,NATYP
            call IoInput('<ZATOM>         ',UIO,I,7,IER)
            if (IER.EQ.0) then
               read (UNIT=UIO,FMT=*) ZAT(I)
               write(111,FMT='(F6.3)') ZAT(I)
            endif
         enddo
      else
         write(111,*) 'zatom will be read in from pot-file'
      endif

      ! Angular momentum cutoff
      call IoInput('LMAX            ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) LMAX
         write(111,*) 'LMAX=',LMAX
      else
         stop 'LMAX not found'
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Allocation of CPA arrays
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call allocate_cpa(1,NAEZ,NEMB,NATYP,NOQ,ICPA,IQAT,HOSTIMP,CONC)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End of allocation of CPA arrays
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do I=1,NAEZ
         ICPA(I) = 0
         NOQ(I) = 1
      end do
      NCPA = 0

      do I = 1,NAEZ
         KAOEZ(1,I) = I       ! default
         IQAT(I) = I          ! Basis-Site of atom I
      enddo
      if (NATYP.EQ.NAEZ) CONC(1:NATYP) = 1.D0

      ! CPA calculation, read concentrations
      if (NATYP.GT.NAEZ) then

         NCPA = 1
         NOQ(1:NAEZ) = 0 ! re-initialize

         IER = 0
         IER2 = 0
         call IoInput('<SITE>          ',UIO,1,7,IER)
         call IoInput('<CPA-CONC>      ',UIO,1,7,IER2)
         if (IER.NE.0.OR.IER2.NE.0) then
            write(1337,*) '<SITE> or <CPA-CONC> not found, will search for ATOMINFOC'
         else

            write(111,FMT='(A18)') '<SITE>  <CPA-CONC>'
            do I = 1,NATYP
               call IoInput('<SITE>          ',UIO,I,7,IER)
               read (UNIT=UIO,FMT=*) IQAT(I)
               call IoInput('<CPA-CONC>      ',UIO,I,7,IER)
               read (UNIT=UIO,FMT=*) CONC(I)
               write(111,FMT='(I5,4X,E16.8)') IQAT(I),CONC(I)
            enddo

            do I = 1,NATYP
               IQ = IQAT(I)
               NOQ(IQ) = NOQ(IQ) + 1
               if ( NOQ(IQ) .GT. 1 ) ICPA(IQ) = 1
               KAOEZ(NOQ(IQ),IQ) = I
            enddo

            do IQ=1,NAEZ
               SUM = 0D0
               if (NOQ(IQ).LT.1) then
                  write(6,*) 'RINPUT13: CPA: SITE',IQ,'HAS NO ASSIGNED ATOM'
                  stop 'RINPUT13: CPA'
               endif
               do IO=1,NOQ(IQ)
                  SUM = SUM + CONC(KAOEZ(IO,IQ))
               end do
               if ( ABS(SUM-1.D0).GT.1D-6) then
                  write(6,*) ' SITE ', IQ, ' CONCENTRATION <> 1.0 !'
                  write(6,*) ' CHECK YOUR <ATOMINFO-CPA> INPUT '
                  stop       ' IN <RINPUT99>'
               end if
            end do

         endif
      endif ! (NATYP.GT.NAEZ)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End atom type information
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Begin relativistic treatment information
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      KCOR = 2
      !      call IoInput('KCOR      ',UIO,1,7,IER)
      !                      read (UNIT=UIO,FMT=*) kcor

      KVREL = 1   ! 0=Schroedinger / 1=SRA / 2=Dirac
      call IoInput('KVREL           ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) kvrel
         write(111,*) 'KVREL= ',KVREL
      else
         write(111,*) 'Default KVREL= ',KVREL
      endif

      ! Set KREL value depending on KVREL
      if (KVREL.ne.0) then
         KREL=1
      else
         KREL=0
      endif

      call IoInput('KORBIT          ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) KORBIT
         write(111,*) 'KORBIT= ',KORBIT
      else
         write(111,*) 'Default KORBIT= ',KORBIT
      endif

      !----------------------------------------------------------------------------
      ! Start of the reading of variables that used to be in the inc.p
      !----------------------------------------------------------------------------
      !> @note JC: Read the IRM value from the inputcard. This in principle can be determined from
      !> the potential file, hence maybe it is best to do it that way instead
      call IoInput('IRM             ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) IRM
         write(111,*) 'IRM= ',IRM
      else
         write(111,*) 'Default IRM= ',IRM
      endif

      call IoInput('IRNSD           ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) IRNSD
         write(111,*) 'IRNSD= ',IRNSD
      else
         write(111,*) 'Default IRNSD= ',IRNSD
      endif

      call IoInput('NSHELD          ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) NSHELD
         write(111,*) 'NSHELD= ',NSHELD
      else
         write(111,*) 'Default NSHELD= ',NSHELD
      endif

      call IoInput('KNOSPH          ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) KNOSPH
         write(111,*) 'KNOSPH= ',KNOSPH
      else
         write(111,*) 'Default KNOSPH= ',KNOSPH
      endif

      call IoInput('IEMXD           ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) IEMXD
         write(111,*) 'IEMXD= ',IEMXD
      else
         write(111,*) 'Default IEMXD= ',IEMXD
      endif

      call IoInput('NR              ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) NR
         write(111,*) 'NR= ',NR
      else
         write(111,*) 'Default NR= ',NR
      endif

      KVREL = 1   ! 0=Schroedinger / 1=SRA / 2=Dirac
      CALL IoInput('KVREL           ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) kvrel
         WRITE(111,*) 'KVREL= ',KVREL
      ELSE
         WRITE(111,*) 'Default KVREL= ',KVREL
      ENDIF
      ! store KVREL to be used later on
      t_inc%KVREL = KVREL

      call IoInput('KPOIBZ          ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) KPOIBZ
         write(111,*) 'KPOIBZ= ',KPOIBZ
      else
         write(111,*) 'Default KPOIBZ= ',KPOIBZ
      endif

      call IoInput('NMAXD           ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) NMAXD
         write(111,*) 'NMAXD= ',NMAXD
      else
         write(111,*) 'Default NMAXD= ',NMAXD
      endif

      call IoInput('ISHLD           ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) ISHLD
         write(111,*) 'ISHLD= ',ISHLD
      else
         write(111,*) 'Default ISHLD= ',ISHLD
      endif

      call IoInput('KNOCO           ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) KNOCO
         write(111,*) 'KNOCO= ',KNOCO
      else
         write(111,*) 'Default KNOCO= ',KNOCO
      endif

      call IoInput('NTREFD          ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) NTREFD
         write(111,*) 'NTREFD= ',NTREFD
      else
         write(111,*) 'Default NTREFD= ',NTREFD
      endif

      call IoInput('NATOMIMPD       ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) NATOMIMPD
         write(111,*) 'NATOMIMPD= ',NATOMIMPD
      else
         write(111,*) 'Default NATOMIMPD= ',NATOMIMPD
      endif

      call IoInput('NPRINCD         ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) NPRINCD
         write(111,*) 'NPRINCD= ',NPRINCD
      else
         write(111,*) 'Default NPRINCD= ',NPRINCD
      endif

      call IoInput('IPAND           ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) IPAND
         write(111,*) 'IPAND= ',IPAND
      else
         write(111,*) 'Default IPAND= ',IPAND
      endif

      call IoInput('NFUND           ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) NFUND
         write(111,*) 'NFUND= ',NFUND
      else
         write(111,*) 'Default NFUND= ',NFUND
      endif

      call IoInput('IRID            ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) IRID
         write(111,*) 'IRID= ',IRID
      else
         write(111,*) 'Default IRID= ',IRID
      endif

      call IoInput('NGSHD           ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) NGSHD
         write(111,*) 'NGHSD= ',NGSHD
      else
         write(111,*) 'Default NGSHD= ',NGSHD
      endif

      call IoInput('WLENGTH         ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) WLENGTH
         write(111,*) 'WLENGTH= ',WLENGTH
      else
         write(111,*) 'Default WLENGTH= ',WLENGTH
      endif

      call IoInput('NACLSD          ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) NACLSD
         write(111,*) 'NACLSD= ',NACLSD
      else
         write(111,*) 'Default NACLSD= ',NACLSD
      endif

      call IoInput('NTOTD           ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) NTOTD
         write(111,*) 'NTOTD= ',NTOTD
      else
         write(111,*) 'Default NTOTD= ',NTOTD
      endif
      !----------------------------------------------------------------------------
      ! End of variables that used to be in the inc.
      !----------------------------------------------------------------------------
      !----------------------------------------------------------------------------
      ! Calculate derived parameters
      !----------------------------------------------------------------------------
      LM2D     = (2*LMAX+1)**2
      NCLSD    = NAEZ + NEMB
      MMAXD    = 2*LMAX+1
      NCLEB    = (LMAX*2+1)**2 * (LMAX+1)**2
      NSPIND   = KREL+(1-KREL)*(NSPIN+1)
      NPOTD    = (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP
      LMMAXD   = (KREL+KORBIT+1)*(LMAX+1)**2
      LMGF0D   = (LMAX+1)**2
      LASSLD   = 4*LMAX
      NEMBD1   = NEMB+1
      IRMIND   = IRM-IRNSD
      NOFGIJ   = NATOMIMPD**2+1
      NTPERD   = NATYP-NTREFD
      NSPINDD  = NSPIND-KORBIT
      NSATYPD  = (NATYP-1)*KNOSPH+1
      NSPOTD   = (2*KREL + (1-KREL)*NSPIND) * NSATYPD
      if (KREL.NE.0.OR.KORBIT.NE.0.OR.KNOCO.NE.0) then
         LNC=.true.
      else
         LNC=.false.
      endif
      !----------------------------------------------------------------------------
      ! End of calculation of the derived parameters
      !----------------------------------------------------------------------------

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Allocation of SOC arrays
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call allocate_SOC(1,KREL,NATYP,LMAX,SOCSCALE,CSCL,SOCSCL)
      allocate(IMANSOC(NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(IMANSOC))*kind(IMANSOC),'IMANSOC','rinput13')
      IMANSOC = 0
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End of allocation of SOC arrays
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (OPT('NEWSOSOL')) then ! Spin-orbit
         if ( OPT('NEWSOSOL') .AND. (NSPIN.NE.2) ) stop ' set NSPIN = 2 for SOC solver in inputcard'
         NPAN_LOG = 30
         NPAN_EQ = 30
         NCHEB = 10
         R_LOG = 0.1D0
         call IoInput('NPAN_LOG        ',UIO,1,7,IER)
         if (IER.EQ.0) read (UNIT=UIO,FMT=*) NPAN_LOG
         call IoInput('NPAN_EQ         ',UIO,1,7,IER)
         if (IER.EQ.0) read (UNIT=UIO,FMT=*) NPAN_EQ
         call IoInput('NCHEB           ',UIO,1,7,IER)
         if (IER.EQ.0) read (UNIT=UIO,FMT=*) NCHEB
         call IoInput('R_LOG           ',UIO,1,7,IER)
         if (IER.EQ.0) read (UNIT=UIO,FMT=*) R_LOG
         write(111,*) 'NPAN_LOG= ',NPAN_LOG
         write(111,*) 'NPAN_EQ= ',NPAN_EQ
         write(111,*) 'NCHEB= ',NCHEB
         write(111,*) 'R_LOG= ',R_LOG
      endif

      call IoInput('<SOCSCL>        ',UIO,1,7,IER)
      if (IER.EQ.0) then
         write(111,'(A10)') '<SOCSCL>  '
         do I = 1,NATYP
            call IoInput('<SOCSCL>        ',UIO,I,7,IER)
            if (IER.EQ.0) then
               read (UNIT=UIO,FMT=*) SOCSCALE(I)
               write(111,FMT='(F6.3)') SOCSCALE(I)
            endif
         enddo
         !        read (UNIT=UIO,FMT=*) (SOCSCALE(I1),I1=1,NATYP)                       !Bernd - old way
         !        write(111,FMT='(A10,50E10.2)') '<SOCSCL>= ',(SOCSCALE(I1),I1=1,NATYP) !Bernd - old way
      else
         write(111,FMT='(A18,50E10.2)') 'Default <SOCSCL>= ',(SOCSCALE(I1),I1=1,NATYP)
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End relativistic treatment information
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Begin cell control
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call IoInput('<FPRADIUS>      ',UIO,1,7,IER)
      if (IER.EQ.0) then
         write(111,'(A10)') '<FPRADIUS>'
         do I = 1,NATYP
            call IoInput('<FPRADIUS>      ',UIO,I,7,IER)
            if (IER.EQ.0) then
               read (UNIT=UIO,FMT=*) FPRADIUS(I)
            endif
            write(111,FMT='(F6.3)') FPRADIUS(I)
         enddo
      else
         write(111,*) 'fpradius will be read in from pot-file'
      endif

      !
      INS = 1
      call IoInput('INS             ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) ins
         write(111,*) 'INS= ',INS
      else
         write(111,*) 'Default INS= ',INS
      endif

      KSHAPE = 2
      if (INS.EQ.0) KSHAPE = 0
      call IoInput('KSHAPE          ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) KSHAPE
         write(111,*) 'KSHAPE= ',KSHAPE
      else
         write(111,*) 'Default KSHAPE= ',KSHAPE
      endif

      if ( (KREL.EQ.1).AND.(KSHAPE.NE.0) ) then
         write(1337,*) ' WARNING : KSHAPE set to ZERO for REL case'
         write(111,*) ' WARNING : kshape set to ZERO for REL case'
         KSHAPE = 0
      end if

      ! Read cell information
      write(1337,*) 'Cell information <SHAPE>:'
      write(111,FMT='(A16)') '<SHAPE>         '
      do I = 1,NATYP
         NTCELL(I) = IQAT(I) ! Default: Different shape function per atom
         call IoInput('<SHAPE>         ',UIO,I,7,IER)
         if (IER.EQ.0) then
            read (UNIT=UIO,FMT=*) NTCELL(I)
            write(111,FMT='(I6)') NTCELL(I)
         endif
      enddo
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End cell control
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Begin exchange correlation treatment information
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      KXC = 2 ! 0=vBH 1=MJW 2=VWN 3=PW91
      call IoInput('KEXCOR          ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) kxc
         write(111,*) 'KEXCOR= ',KXC
      else
         write(111,*) 'Default KEXCOR= ',KXC
      endif

      ! Scale magnetic moment (0 < Lambda_XC < 1,  0=zero moment, 1= full moment)
      LAMBDA_XC = 1.D0
      call IoInput('LAMBDA_XC       ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) LAMBDA_XC
         write(111,*) 'LAMBDA_XC= ',LAMBDA_XC
      else
         write(111,*) 'Default LAMBDA_XC= ',LAMBDA_XC
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! LDA+U treatment
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Allocate the LDA+U arrays
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call allocate_ldau(1,NATYP,LOPT,UEFF,JEFF,EREFLDAU)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End of LDA+U array allocation
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (OPT('LDA+U   ')) then

         !Check for LDA+U consistency -- if INS=0 suppress it
         if ((INS.EQ.0) ) then
            write (1337,*)
            write (1337,*)&
            &        ' WARNING: LDA+U should be used only in NON-SPHERICAL',&
            &        ' case (INS=1) '
            write (1337,*) ' Running option LDA+U will be ignored'
            write (1337,*)
            do I=1,32
               if (t_params%OPTC(I)(1:8).EQ.'LDA+U   ') t_params%OPTC(I)='        '
            end do
         end if

         ! -> get number of atoms for lda+u:

         IER = 0
         call IoInput('NAT_LDAU        ',UIO,1,7,IER)
         if ( IER.NE.0 ) then
            NASOC = NATYP
         else
            read (UNIT=UIO,FMT=*) NASOC
            if ( NASOC.GT.NATYP ) stop ' main0: NAT_LDAU > NATYP'
         end if

         ! -> read in UEFF,JEFF,LOPT,EREFLDAU for the desired atoms

         IL = 0
         do I=1,NASOC
            IER = 0
            call IoInput('LDAU_PARA       ',UIO,I,7,IER)
            if ( IER.EQ.0 ) then
               read (UNIT=UIO,FMT=*) I1,LOPT(I1),UEFF(I1),JEFF(I1),EREFLDAU(I1)
               IL = IL + 1
            end if
         enddo
         if ( IL.NE.NASOC ) then
            write(6,*) ' ERROR: LDA+U invoked for ',NASOC,' atoms'
            write(6,*) '        Some (all) parameters are missing in the input-file'
            stop
         end if
         KREADLDAU = 0
         IER = 0
         call IoInput('KREADLDAU       ',UIO,1,7,IER)
         if ( IER.EQ.0 ) read (UNIT=UIO,FMT=*) KREADLDAU

      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End exchange correlation treatment information
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Begin external field control
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      KHFIELD = 0
      HFIELD = 0.D0
      call IoInput('HFIELD          ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) HFIELD
         if (HFIELD.NE.0.D0) then
            KHFIELD = 1
            write(*,*) 'WARNING: HFIELD>0.0 found, set KHFIELD to 1'
            write(1337,*) 'WARNING: HFIELD>0.0 found, set KHFIELD to 1'
         end if
         write(111,*) 'HFIELD= ',HFIELD
      else
         write(111,*) 'Default HFIELD= ',HFIELD
      endif

      VCONST = 0.D0
      call IoInput('VCONST          ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) vconst
         write(111,*) 'VCONST= ',VCONST
      else
         write(111,*) 'Default VCONST= ',VCONST
      endif


      if (TEST('atptshft')) then
         write(1337,*) 'read IN IVSHIFT'
         call IoInput('IVSHIFT         ',UIO,1,7,IER)
         read (UNIT=UIO,FMT=*) ivshift
      endif

      ! Initial polarization
      LINIPOL = .FALSE.
      call IoInput('LINIPOL         ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) LINIPOL
         write(111,*) 'LINIPOL= ',LINIPOL
      else
         write(111,*) 'Default: LINIPOL= ',LINIPOL
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Allocate magnetization arrays
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call allocate_magnetization(1,NAEZ,NATYP,LMMAXD,INIPOL,IXIPOL,QMTET,&
      QMPHI,DROTQ)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End of allocation of magnetization arrays
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (LINIPOL) then
         INIPOL(1:NATYP) = 1
         call IoInput('XINIPOL         ',UIO,1,7,IER)
         if (IER.EQ.0) then
            read (UNIT=UIO,FMT=*) (inipol(I),I=1,natyp)
            write(111,FMT='(A10,80I2)') 'XINIPOL=  ',(INIPOL(I),I=1,NATYP)
         else
            write(111,FMT='(A18,80I2)') 'Default XINIPOL=  ',(INIPOL(I),I=1,NATYP)
         endif
      endif


      write (1337,2021) (INIPOL(I),I=1,NATYP)
      write(1337,2103)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End external field control
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Begin Green function calculation control (diag./non-diag)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IGF = 0
      call IoInput('IGREENFUN       ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) igf
         write(111,*) 'IGREENFUN= ',IGF
      else
         write(111,*) 'Default IGREENFUN= ',IGF
      endif
      
      IF (OPT('OPERATOR')) THEN
        ! check if impurity files are present (otherwise no imp.
        ! wavefunctions can be calculated)
        operator_imp = .true.
        inquire(file='potential_imp', exist=lexist)
        if (.not.lexist) operator_imp = .false.
        inquire(file='shapefun_imp', exist=lexist)
        if (.not.lexist) operator_imp = .false.
        inquire(file='scoef', exist=lexist)
        if (.not.lexist) operator_imp = .false.
      ELSE
        operator_imp = .false.
      END IF
      if (OPT('KKRFLEX ') .or. OPT('WRTGREEN') .or. OPT('GREENIMP') .or. operator_imp) then
         write(1337,*) 'Setting IGREENFUN=1 for KKRFLEX/WRTGREEN/GREENIMP/OPERATOR options'
         IGF = 1
      end if

      ICC = 0
      call IoInput('ICC             ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) ICC
         write(111,*) 'ICC= ',ICC
      else
         write(111,*) 'Default ICC= ',ICC
      endif
      if (OPT('KKRFLEX ') .or. OPT('WRTGREEN') .or. OPT('GREENIMP') .or. operator_imp) then
         write(1337,*) 'Setting ICC=1 for KKRFLEX/WRTGREEN/GREENIMP/OPERATOR  options'
         ICC = 1
      end if
      if ( ( OPT('XCPL    ') ).OR.( OPT('CONDUCT ') ) ) ICC = -1

      if ( ICC.NE.0 .AND. IGF.EQ.0 ) IGF = 1
      if ( ICC.EQ.0 .AND. IGF.NE.0 ) ICC = -1
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End Green function calculation control (diag./non-diag)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Begin accuracy parameters
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Brilloun zone mesh
      INTERVX = 10
      INTERVY = 10
      INTERVZ = 10
      call IoInput('BZDIVIDE        ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) INTERVX,INTERVY,INTERVZ
         write(111,FMT='(A9,3I5)') 'BZDIVIDE=',INTERVX,INTERVY,INTERVZ
      else
         write(111,FMT='(A17,3I5)') 'Default BZDIVIDE=',INTERVX,INTERVY,INTERVZ
      endif
      write(1337,2104)
      write(1337,2015) INTERVX,INTERVY,INTERVZ
      write(1337,2102)

      if(OPT('GREENIMP')) then
         write(*,*) 'WARNING! Found option GREENIMP: resetting BZDIVIDE to 1,1,1'
         write(1337,*) 'WARNING! Found option GREENIMP: resetting BZDIVIDE to 1,1,1'
         INTERVX = 1
         INTERVY = 1
         INTERVZ = 1
      endif

      ! Energy contour
      NPOL = 7
      !      if (OPT('dos     ').OR.OPT('DOS     ')) NPOL = 0
      call IoInput('NPOL            ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) NPOL
         write(111,*) 'NPOL=',NPOL
      else
         write(111,*) 'Default NPOL=',NPOL
      endif

      call IoInput('EMIN            ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) EMIN
         write(111,*) 'EMIN= ',EMIN
      else if (NPOL.EQ.0) then
         EMIN = -1.D0
         write(111,*) 'Default for DOS: EMIN= ',EMIN
      else
         write(1337,*) 'Error in rinput13: EMIN not found'
         write(111,*) 'Error in rinput13: EMIN not found'
         stop 'Error in rinput13: EMIN not found'
      endif

      EMAX = 1.D0
      call IoInput('EMAX            ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) EMAX
         write(111,*) ' EMAX=',EMAX
      else
         write(111,*) 'Default  EMAX=',EMAX
      endif

      TK = 800.D0
      call IoInput('TEMPR           ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) TK
         write(111,*) 'TEMPR=',TK
      else
         write(111,*) 'Default TEMPR=',TK
      endif

      NPNT1 = 3
      if (NPOL.EQ.0) NPNT1 = 0
      call IoInput('NPT1            ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) NPNT1
         write(111,*) ' NPT1=',NPNT1
      else
         write(111,*) 'Default  NPT1=',NPNT1
      endif

      NPNT2 = NINT((EMAX-EMIN)*20.D0) ! 20 pts/Ryd
      if (NPOL.EQ.0) NPNT2 = NINT((EMAX-EMIN)*100.D0) ! For dos, 100 pts/Ryd
      call IoInput('NPT2            ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) NPNT2
         write(111,*) ' NPT2=',NPNT2
      else
         write(111,*) 'Default  NPT2=',NPNT2
      endif

      NPNT3 = 3
      if (NPOL.EQ.0) NPNT3 = 0
      call IoInput('NPT3            ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) NPNT3
         write(111,*) ' NPT3=',NPNT3
      else
         write(111,*) 'Default  NPT3=',NPNT3
      endif

      ! -> semicore
      ! initialise variables
      IDOSEMICORE = 0
      EBOTSEMI = EMIN
      EMUSEMI = EBOTSEMI
      NPOLSEMI = 0
      N1SEMI = 0
      N2SEMI = 0
      N3SEMI = 0
      FSEMICORE = 1.D0

      IER = 0
      if ( OPT('SEMICORE') ) then
         call IoInput('EBOTSEMI        ',UIO,1,7,IER)
         if ( IER.NE.0 ) GOTO 99800
         read (UNIT=UIO,FMT=*) EBOTSEMI
         call IoInput('EMUSEMI         ',UIO,1,7,IER)
         if ( IER.NE.0 ) GOTO 99800
         read (UNIT=UIO,FMT=*) EMUSEMI

         ! -> EMUSEMI < EBOT
         if ( EMUSEMI.GE.EMIN ) GOTO 99800
         call IoInput('TKSEMI          ',UIO,1,7,IER)
         if ( IER.NE.0 ) GOTO 99800
         read (UNIT=UIO,FMT=*) TKSEMI

         call IoInput('NPOLSEMI        ',UIO,1,7,IER)
         if ( IER.NE.0 ) GOTO 99800
         read (UNIT=UIO,FMT=*) NPOLSEMI
         call IoInput('N1SEMI          ',UIO,1,7,IER)
         if ( IER.NE.0 ) GOTO 99800
         read (UNIT=UIO,FMT=*) N1SEMI
         call IoInput('N2SEMI          ',UIO,1,7,IER)
         if ( IER.NE.0 ) GOTO 99800
         read (UNIT=UIO,FMT=*) N2SEMI
         call IoInput('N3SEMI          ',UIO,1,7,IER)
         if ( IER.NE.0 ) GOTO 99800
         read (UNIT=UIO,FMT=*) N3SEMI
         call IoInput('FSEMICORE       ',UIO,1,7,IER)
         if ( IER.NE.0 ) GOTO 99800
         read (UNIT=UIO,FMT=*) FSEMICORE
         IDOSEMICORE = 1
         99800    continue
         if ( IDOSEMICORE.EQ.0 ) then
            write (1337,*)
            write (1337,*) ' WARNING: SEMICORE used',&
            &           ' with incomplete/incorrect contour description'
            write (1337,*) ' Running option SEMICORE will be ignored'
            write (111,*)
            write (111,*) ' WARNING: SEMICORE used',&
            &           ' with incomplete/incorrect contour description'
            write (111,*) ' Running option SEMICORE will be ignored'
            do I=1,32
               if (t_params%OPTC(I)(1:8).EQ.'SEMICORE') t_params%OPTC(I)='        '
            end do
         end if
      end if

      ! CPA convergence parameters
      CPATOL = 1D-4
      ITCPAMAX = 20
      call IOINPUT('CPAINFO         ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) CPATOL,ITCPAMAX
      else
         write(111,*) 'Default cpainfo:'
      endif
      write(111,FMT='(A7)') 'CPAINFO'
      write(111,FMT='(E12.4,I5)') CPATOL,ITCPAMAX

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Begin screening cluster information
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      RCUTZ = 11.D0/ALAT  ! Default 11 Bohr radii
      call IoInput('RCLUSTZ         ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) RCUTZ
         write(111,*) 'RCLUSTZ=',RCUTZ
      else
         write(111,*) 'Default RCLUSTZ=',RCUTZ
      endif

      RCUTXY = RCUTZ
      call IoInput('RCLUSTXY        ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) RCUTXY
         write(111,*) 'RCLUSTXY=',RCUTXY
      else
         write(111,*) 'Default RCLUSTXY=',RCUTXY
      endif

      write(1337,*) 'Parameters used for the cluster calculation'
      if (abs(rcutz-rcutxy).lt.1.d-4) then
         write(1337,*) 'Clusters inside spheres with radius R = ',rcutz
      else
         write(1337,*) 'Clusters inside cylinders with '
         write(1337,*) 'Rz = ',rcutz,' Rxy = ',rcutxy
      end if
      write(1337,2104)
      write(1337,2018)                 ! rbasis
      write(1337,2101)
      do i=1,naez
         write(1337,2025) i,(rbasis(j,i),j=1,3),&
         QMTET(I),QMPHI(I),ICPA(I),NOQ(I),(KAOEZ(J,I),J=1,NOQ(I))
      enddo

      do I=1,NAEZ
         call IoInput('<RMTREF>        ',UIO,I,7,IER)
         if (IER.EQ.0) then
            read (UNIT=UIO,FMT=*) RMTREFAT(I)
         endif
      enddo
      if (IER.EQ.0) then
         write(111,FMT='(A18)') '        <RMTREF>  '
      else
         write(111,FMT='(A18)') 'Default <RMTREF>  '
      endif
      do I=1,NAEZ
         write(111,FMT='(9X,F9.6)') RMTREFAT(I)
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End screening cluster information
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Number of Born iterations
      ICST = 2
      call IoInput('ICST            ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) ICST
         write(111,*) 'ICST=',ICST
      else
         write(111,*) 'Default ICST=',ICST
      endif

      ! Usage of Lloyd's formula
      LLY = 0 ! LLY Default=0 : do not apply Lloyds formula
      if (OPT('LLOYD   ').OR.OPT('Lloyd   ').OR.OPT('lloyd   ')) LLY = 1
      call IoInput('<LLOYD>         ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) LLY
         write(111,*) '<LLOYD>=',LLY
      else
         write(111,*) 'Default <LLOYD>=',LLY
      endif
      if (LLY.NE.0) write(1337,*) 'Applying Lloyds formula, LLY=',LLY

      DELTAE = (1.D-5,0.D0) ! Difference for numer. derivative in Lloyds formula
      call IoInput('<DELTAE>        ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) DELTAE
         write(111,*) '<DELTAE>=',DELTAE
      else
         write(111,*) 'Default <DELTAE>=',DELTAE
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End accuracy parameters
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Begin old-type of ATOMINFO
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      LATOMINFO = .FALSE.
      ! Initialize all clusters to 1
      CLS(1:NAEZ+NEMB) = 1
      write(1337,*) 'ATOMINFOC or ATOMINFO:'
      do I=1,NATYP
         call IoInput('ATOMINFOC       ',UIO,I+1,7,IER)
         IA = 1
         if ( IER.EQ.0 ) then
            LATOMINFO = .TRUE.
            read (UNIT=UIO,FMT=*)   ZAT(I),           &
            LMXC(I),          &
            (KFG(J,I),J=1,4), &
            J,                &
            IER,              &
            NTCELL(I),        &
            MTFAC(I),         &
            IRNS(I),          &
            RMTREF(IER),IQAT(I),CONC(I)
            IQ = IQAT(I)
            REFPOT(IQ) = IER
            RMTREFAT(I) = RMTREF(IER)
            CLS(IQ) = J
            NOQ(IQ) = NOQ(IQ) + 1
            if ( NOQ(IQ) .GT. 1 ) then
               ICPA(IQ) = 1
               NCPA = 1
            end if
            KAOEZ(NOQ(IQ),IQ) = I
         else
            IER = 0
            call IoInput('ATOMINFO        ',UIO,I+1,7,IER)
            if (IER.EQ.0) then
               LATOMINFO = .TRUE.
               read (UNIT=UIO,FMT=*)   ZAT(I),           &
               LMXC(I),          &
               (KFG(J,I),J=1,4), &
               J,                &
               REFPOT(I),        &
               NTCELL(I),        &
               MTFAC(I),         &
               IRNS(I),          &
               RMTREF(REFPOT(I))
               IQAT(I) = I
               RMTREFAT(I) = RMTREF(REFPOT(I))
               CLS(I) = J
               CONC(I) = 1.D0
               NOQ(I) = 1
               KAOEZ(1,I) = I
            endif
         end if
      end do

      ! If old-style ATOMINFO is present, and if a 2-dim calculation is performed,
      ! and also if the RMTREF of the "outside region" is not read in explicitly
      ! (LNEW is false) then assign the RMTREF of the outside region according to
      ! the already-read-in REFPOT under LEFTBASIS  and RIGHBASIS.
      if (LATOMINFO.AND.LINTERFACE.AND..NOT.LNEW) then
         do I = NAEZ + 1,NAEZ + NEMB
            RMTREFAT(I) = RMTREF(REFPOT(I))
         enddo
      endif

      NCLS = 0
      NREF = 0

      ! Determine total number of clusters
      do I=1,NATYP
         NCLS = MAX(NCLS,CLS(IQAT(I)))
      enddo

      ! Determine total number of different reference potentials
      do I=1,NAEZ + NEMB
         NREF = MAX(NREF,REFPOT(I))
      enddo

      !in line 1792  this is done: NINEQ = NAEZ, so here NINEQ is still undefinded
      !so we move this writeout back
      !
      !write(6,2016) NCLS,NREF,NINEQ
      !write(6,2110)
      !write(6,2103)

      do IQ=1,NAEZ
         SUM = 0D0
         if (NOQ(IQ).LT.1) then
            write(6,*) 'RINPUT13: CPA: SITE',IQ,'HAS NO ASSIGNED ATOM'
            stop 'RINPUT13: CPA'
         endif
         do IO=1,NOQ(IQ)
            SUM = SUM + CONC(KAOEZ(IO,IQ))
         end do
         if ( ABS(SUM-1.D0).GT.1D-6) then
            write(6,*) ' SITE ', IQ, ' CONCENTRATION <> 1.0 !'
            write(6,*) ' CHECK YOUR <ATOMINFO-CPA> INPUT '
            stop       ' IN <RINPUT99>'
         end if
      end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End old-type of ATOMINFO
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Write out atominfo
      write(1337,2028) NATYP
      write(1337,2104)
      write(1337,1029) (                     &
      ZAT(I),           &
      LMXC(I),          &
      (KFG(J,I),J=1,4), &
      CLS(IQAT(I)),     &
      REFPOT(IQAT(I)),  &
      NTCELL(I),        &
      MTFAC(I),         &
      IRNS(I),          &
      IQAT(I),CONC(I),I=1,NATYP)
      write(1337,2108)
      write(1337,2104)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Begin SCF convergence control
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      NSTEPS = 1
      call IoInput('NSTEPS          ',UIO,1,7,IER)
      if (IER.NE.0) then
         write(111,*) 'Default NSTEPS=',NSTEPS
      else
         read (UNIT=UIO,FMT=*) NSTEPS
      endif
      if (NPOL.EQ.0) then
         NSTEPS = 1
         write(1337,*) 'NPOL=0, setting NSTEPS to 1'
      endif
      if (IGF.NE.0) then
         NSTEPS = 1
         write(1337,*) 'IGF.NE.0, setting NSTEPS to 1'
      endif
      if (ICC.NE.0) then
         NSTEPS = 1
         write(1337,*) 'ICC.NE.0, setting NSTEPS to 1'
      endif
      if (OPT('XCPL    ')) then
         NSTEPS = 1
         write(1337,*) 'RUNOPT XCPL used, setting NSTEPS to 1'
      endif
      if (OPT('KKRFLEX ')) then
         NSTEPS = 1
         write(1337,*) 'RUNOPT KKRFLEX used, setting NSTEPS to 1'
      endif

      write(1337,2011) NSTEPS
      write(1337,2104)

      IMIX = 0
      call IoInput('IMIX            ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) imix
         write(111,*) 'IMIX= ',IMIX
      else
         write(111,*) 'Default IMIX= ',IMIX
      endif
      if (NPOL.EQ.0) then
         write(1337,*) 'NPOL=0, setting IMIX= 0'
         IMIX = 0
      endif

      STRMIX = 0.01D0
      call IoInput('STRMIX          ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) STRMIX
         write(111,*) 'STRMIX= ',STRMIX
      else
         write(111,*) 'Default STRMIX= ',STRMIX
      endif
      if (NPOL.EQ.0) then
         write(1337,*) 'NPOL=0, setting STRMIX= 0.'
         STRMIX = 0
      endif

      ITDBRY = 40
      call IoInput('ITDBRY          ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) itdbry
         write(111,*) 'ITDBRY= ',ITDBRY
      else
         write(111,*) 'Default ITDBRY= ',ITDBRY
      endif

      FCM = 20.D0
      call IoInput('FCM             ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) fcm
         write(111,*) 'FCM= ',FCM
      else
         write(111,*) 'Default FCM= ',FCM
      endif

      QBOUND = 1.D-7
      call IoInput('QBOUND          ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) qbound
         write(111,*) 'QBOUND= ',QBOUND
      else
         write(111,*) 'Default QBOUND= ',QBOUND
      endif

      BRYMIX = 0.01D0
      call IoInput('BRYMIX          ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) brymix
         write(111,*) 'BRYMIX= ',BRYMIX
      else
         write(111,*) 'Default BRYMIX= ',BRYMIX
      endif

      call IOINPUT('RMAX            ',UIO,1,7,IER)
      if ( IER.NE.0 ) stop 'rinput13: RMAX not in the inputcard'
      read (UNIT=UIO,FMT=*) RMAX
      write(111,*) 'RMAX= ',RMAX

      call IOINPUT('GMAX            ',UIO,1,7,IER)
      if ( IER.NE.0 ) stop 'rinput13: GMAX not in the inputcard'
      read (UNIT=UIO,FMT=*) GMAX
      write(111,*) 'GMAX= ',GMAX
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End SCF convergence control
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Begin file name definitions
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IL=1
      call IoInput('FILES           ',UIO,IL,7,IER)
      if (IER.EQ.0) then
         call IoInput('FILES           ',UIO,IL,7,IER)
         read (UNIT=UIO,FMT='(A40)')  I12
         call IoInput('FILES           ',UIO,IL+1,7,IER)
         read (UNIT=UIO,FMT='(A40)')  I13
         call IoInput('FILES           ',UIO,IL+2,7,IER)
         read (UNIT=UIO,FMT='(A40)')  I40
         call IoInput('FILES           ',UIO,IL+3,7,IER)
         read (UNIT=UIO,FMT='(A40)')  I19
         call IoInput('FILES           ',UIO,IL+4,7,IER)
         read (UNIT=UIO,FMT='(A40)')  I25
      else
         I13 = 'potential                               ' ! 40 chars
         I19 = 'shapefun                                ' ! 40 chars
         I25 = 'scoef                                   ' ! 40 chars
         I12 = '                                        ' ! 40 chars (not used)
         I40 = '                                        ' ! 40 chars (not used)
      endif

      write(1337,*) 'I12="',I12,'"'
      write(1337,*) 'I13="',I13,'"'
      write(1337,*) 'I40="',I40,'"'
      write(1337,*) 'I19="',I19,'"'
      write(1337,*) 'I25="',I25,'"'

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End file name definitions
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IFILE = 13
      call IoInput('<IFILE>         ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) ifile
         write(111,*) '<IFILE>= ',IFILE
      else
         write(111,*) 'Default <IFILE>= ',IFILE
      endif

      IPE = 1    ! Used to print out in calrmt
      ISHIFT = 0
      call IoInput('ISHIFT          ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) ishift
         write(111,*) 'ISHIFT= ',ISHIFT
      else
         write(111,*) 'Default ISHIFT= ',ISHIFT
      endif
      if (  OPT('rigid-ef').OR. OPT('DECIMATE') ) then
         ISHIFT = 2
         write(1337,*) ' Rigid Fermi Energy, ISHIFT is set to ',ISHIFT
         write(111,*) ' Rigid Fermi Energy, ishift is set to ',ISHIFT
      end if
      if ( TEST('no-neutr').OR.OPT('no-neutr') )   then
         ISHIFT = 1
         write(1337,*) 'No charge neutrality required, ISHIFT is set to',ISHIFT
         write(111,*) 'No charge neutrality required, ISHIFT is set to',ISHIFT
      endif

      ESHIFT = 0.D0
      INSREF = 0
      KWS = 2
      KHYP = 0

      TOLRDIF = 0.5D0 ! Set free GF to zero for r<tolrdif (a.u.)(vir. atoms)
      call IoInput('<TOLRDIF>       ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) TOLRDIF
         write(111,*) '<TOLRDIF>=',TOLRDIF
      else
         write(111,*) 'Default <TOLRDIF>=',TOLRDIF
      endif

      ! -------------------------------------------------
      KTE = 1
      !     call IoInput('KTE       ',UIO,1,7,IER)
      !                      read (UNIT=UIO,FMT=*) kte

      KPRE = 1
      !     call IoInput('KPRE      ',UIO,1,7,IER)
      !                      read (UNIT=UIO,FMT=*) kpre

      KEFG = 0
      !     call IoInput('KEFG      ',UIO,1,7,IER)
      !                      read (UNIT=UIO,FMT=*) kefg

      KVMAD = 0
      call IoInput('KVMAD           ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) kvmad
         write(111,*) 'KVMAD= ',KVMAD
      else
         write(111,*) 'Default KVMAD= ',KVMAD
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Determination of properties at Fermi level
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ( OPT('GF-EF   ') ) then
         IGF = 1
         if (NPOL.GT.0) NPOL = 0
         if (NPOL.LT.0) then
            NPNT1 = 0
            NPNT3 = 0
         end if
         NPNT2 = 1
      end if

      if (OPT('DOS-EF  ')) then
         NPOL = 0
         NPNT2 = 1
      end if
      ! ----------------------------------------------------------------------
      ! ---------------------------------------------------------------------
      !
      KFORCE = 0
      if (INS.GT.0) then
         call IoInput('KFORCE          ',UIO,1,7,IER)
         if (IER.EQ.0) then
            read (UNIT=UIO,FMT=*) KFORCE
            write(111,*) 'KFORCE= ',KFORCE
         else
            write(111,*) 'Default KFORCE= ',KFORCE
         endif
      end if

      KFROZN = KCOR
      if (KCOR.EQ.0) KCOR = 2

      ! ------------------------------------------------------------------------
      write (1337,9210) LMAX
      write (1337,9301)
      write (1337,9220) EMIN,EMAX,TK
      write (1337,9302)
      write (1337,9230) NPOL,NPNT1,NPNT2,NPNT3
      write (1337,9304)
      write (1337,9303)
      write (1337,9250) IFILE,IPE,ISHIFT,ESHIFT
      write (1337,9305)
      write (1337,9260) KSHAPE,IRM,INS,ICST,INSREF
      write (1337,9309)
      write (1337,9270) KCOR,KVREL,KWS,KHYP,KHFIELD,KXC
      write (1337,9306)
      write (1337,9330) KTE,KPRE,KEFG,KVMAD
      write (1337,9309)
      write (1337,9290) IMIX,IGF,ICC
      write (1337,9304)
      write (1337,9300) ITDBRY
      write (1337,9307)
      write (1337,9310) STRMIX,FCM,QBOUND
      write (1337,9302)
      write (1337,9320) BRYMIX
      write (1337,9308)
      write (1337,9280) HFIELD,VCONST
      ! ------------------------------------------------------------------------

      ! ------------------------------------------------------------------------
      IPF = 1337
      IPFE = IPF + 3

      if (OPT('SEARCHEF')) then
         IMIX=0
         MIXING=0.0d0
         STRMIX=MIXING
         ITDBRY=1
         QBOUND=1.0d-10
         write(1337,'(1X,A)') 'Option SEARCHEF used overriding INPUT for'
         write(1337,'(1X,A)') 'IMIX,MIX,QBOUND,ITDBRY: 0, 0.0, 1E-10, 1'
         write(1337,*)
      endif

      if (IMIX.GT.2) then
         FCM = 1.0D0
         MIXING = BRYMIX
      else
         MIXING = STRMIX
      end if

      if (IMIX.GE.6) write (1337,FMT=9110) (IMIX-5),ITDBRY - 1

      write (1337,FMT=9090) MIXING,QBOUND
      !--------------------------------------------------------
      write (1337,FMT=9091) CPAFLAG(NCPA)
      if (NCPA.NE.0) write(1337,9092) ITCPAMAX,CPATOL
      !--------------------------------------------------------

      LMMAX = (LMAX+1)**2
      LPOT  = 2*LMAX
      LMPOT = (LPOT+1)*(LPOT+1)
      LMXSPD = (2*LPOT+1)**2

      write (1337,FMT=9020) LMAX,LMAX,NATYP,NATYP,IRM,IRM,NSPIN,NSPIND

      if (INS.GT.0) then
         write (1337,FMT=9130)
         write (1337,FMT=9140)
         do I = 1,NATYP
            write (1337,FMT=9150) I,IRNS(I),IRNSD
            if (IRNS(I).GT.IRNSD) call RCSTOP('19      ')
         enddo
      end if

      write (1337,FMT=9130)

      if (KHFIELD.EQ.1) write (1337,FMT=9030) HFIELD
      if (KVREL.LE.1 ) then
         write (1337,FMT=9050) TSPIN(NSPIN)
      else
         write (1337,FMT=9050) TSPIN(NSPIN+1)
      end if
      write (1337,FMT=9170) TVREL(KVREL)
      write (1337,FMT=9170) TKCOR(KFROZN)
      if (KSHAPE.EQ.0) then
         write (1337,FMT=9070) TKWS(KWS+1)
      else
         write (1337,FMT=9170) TSHAPE
      end if

      write (1337,FMT=9100) TXC(KXC+1)
      if (INS.GT.0) write (1337,FMT=9160) TINS(INS),ICST
      write (1337,FMT=9080)

      VBC(1) = VCONST
      VBC(2) = VBC(1)

      LRHOSYM = .FALSE.
      call IoInput('LRHOSYM         ',UIO,1,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) lrhosym
         write(111,*) 'LRHOSYM= ',LRHOSYM
      else
         write(111,*) 'Default LRHOSYM= ',LRHOSYM
      endif

      if ( (NCPA.NE.0).AND.LRHOSYM ) then
         write(1337,*) ' WARNING : CHARGE SYMMETRISATION NOT ALLOWED FOR CPA '
         write(1337,*) '        YOUR SETTING IN INPUT FILE IS OVERRIDDEN'
         write(111,*) ' WARNING : CHARGE SYMMETRISATION NOT ALLOWED FOR CPA '
         write(111,*) '    YOUR SETTING IN INPUT FILE IS OVERRIDDEN'
         LRHOSYM = .FALSE.
      end if

      if (LRHOSYM) then
         call IoInput('IXIPOL          ',UIO,1,7,IER)
         read (UNIT=UIO,FMT=*) (ixipol(I),I=1,natyp)
         write (1337,2022) (ixipol(i),i=1,natyp)
         write (1337,2103)
         do I=1,NATYP
            if ( IXIPOL(I).NE.0 .AND. ABS(IXIPOL(ABS(IXIPOL(I)))).NE.I) then
               write(6,*) 'Error in IXIPOL at atom ',I,'.'
               stop 'IXIPOL'
            end if
         end do
      else
         do I=1,NATYP
            IXIPOL(I) = 0
         end do
         write (1337,2022) (ixipol(i),i=1,natyp)
         write (1337,2103)
      end if
      write(1337,2023) NAEZ,NEMB
      write(1337,2110)

      NINEQ = NAEZ
      write(1337,2016) NCLS,NREF,NINEQ
      write(1337,2110)
      write(1337,2103)

      !----------------------------------------------------------------------------
      KMROT = 0

      do I=1,NAEZ
         !-------------------------------------------------------------------------
         ! Atoms equivalent by inversional symmetry
         !-------------------------------------------------------------------------
         QMTET(I)=0D0
         QMPHI(I)=0D0
         IER = 0
         call IoInput('RBASISANG       ',UIO,I,7,IER)

         if( IER.EQ.0 ) then
            read (UNIT=UIO,FMT=*) (RBASIS(J,I), J=1,3),QMTET(I),QMPHI(I)
            if( ABS(QMTET(I)) .GT. 1D-6 ) KMROT = 1
            if( ABS(QMPHI(I)) .GT. 1D-6 ) KMROT = 1
         endif
      enddo                         ! I=1,NAEZ
      call IDREALS(RBASIS(1,1),3*NAEZ,IPRINT)
      !----------------------------------------------------------------------------

      if (nemb.gt.0) write(1337,*)
      write(1337,2031) ((rbasis(j,i),j=1,3),i,refpot(i),i=naez+1,naez+nemb)

      ! ------------------------------------------------------------------------
      if ( .not. OPT('VIRATOMS') ) then
         do I=1,NAEZ
            do IO=1,NOQ(I)
               if (KAOEZ(IO,I).LT.1) stop 'Error in KAOEZ'
            end do
         enddo
      end if
      ! ------------------------------------------------------------------------
      write(1337,2111)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Check for DECIMATE consistency
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (OPT('DECIMATE')) then
         if ( MOD(NPRINCD,NLBASIS).NE.0 ) then
            write(6,*) ' Decimation cannot continue '
            write(6,*) 'NPRINCD=',NPRINCD,' NLBASIS=',NLBASIS
            stop
         end if
         if ( MOD(NPRINCD,NRBASIS).NE.0 )  then
            write(6,*) ' Decimation cannot continue '
            write(6,*) 'NPRINCD=',NPRINCD,' NRBASIS=',NRBASIS
            stop
         end if
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Check for ITERMDIR consistency -- if KMROT=0 suppress it
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ( (OPT('ITERMDIR')).AND.(KMROT.EQ.0) ) then
         write (1337,*)
         write (1337,*)' WARNING: ITERMDIR running option used with collinear/',&
         &        'parallel Oz starting'
         write (1337,*)'          system (KMROT = 0 ). Please check token',&
         &        ' RBASISANG in your input'
         write (1337,*) ' Running option ITERMDIR will be ignored'
         write (1337,*)
         do I=1,32
            if (t_params%OPTC(I)(1:8).EQ.'ITERMDIR') t_params%OPTC(I)='        '
         end do
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Check for XCPL consistency
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      MANCTL = ( KMROT.EQ.0 ).AND.( KREL.EQ.0 ).AND.( NSPIN.GT.1 )
      if ( (OPT('XCPL    ') ).AND.( .NOT.MANCTL ) ) then
         write (1337,*)
         write (1337,*)' WARNING: XCPL running option requires collinear ',&
         &        'magnetic systems'
         write (1337,*)' in a NON/SCALAR/SCALAR+SOC relativistic mode (KREL=0)'
         write (1337,*) ' Running option XCPL will be ignored'
         write (1337,*)
         do I=1,32
            if (t_params%OPTC(I)(1:8).EQ.'XCPL    ') t_params%OPTC(I)='        '
         end do
      end if

      write(1337,62) (t_params%OPTC(I),I=1,8)
      62   format(79('-')/' EXECUTION OPTIONS:'/1X,A8,7('//',A8)/79('-'))
      write(1337,52) (t_params%TESTC(I),I=1,16)
      52   format(79('-')/' TEST OPTIONS:'/2(1X,A8,7('//',A8)/)/79('-'))
      980  format(8A8)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Initialise SOLVER, SOC and CTL parameters in REL case
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CSCL(1:LMAX+1,1:NATYP) = CVLIGHT
      MANSOC=.FALSE.
      MANCTL=.FALSE.

      if (KREL.EQ.1) then
         SOLVER='BS        '
         call IoInput('SOLVER          ',UIO,0,7,IER)
         if (IER.EQ.0) then
            read (UNIT=UIO,FMT=*) SOLVER
            if ( SOLVER(1:2) .EQ. 'BS' ) then
               SOLVER = 'BS        '
            else
               if( SOLVER .NE. 'ABM-OP    ' ) SOLVER='ABM-OP    '
            end if
         end if

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! SOC-MAN
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! For Dirac-ASA
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if (OPT('SOC     ')) then
            call IOInput('SOSCALE         ',UIO,0,7,IER)
            if (IER.EQ.0) then
               read (UNIT=UIO,FMT=*) SOSCALE
               if (SOSCALE.GT.-2.5D0) then
                  if (SOSCALE.GE.0.0D0) then           ! SOC-I
                     SOLVER='ABM-SOC   '
                     MANSOC=.TRUE.
                  else                                  ! SOC-II
                     SOLVER       = 'ABM-SOC-II'
                     MANSOC=.TRUE.
                     do I=1,NATYP
                        SOCSCL(1:LMAX+1,I) = SOSCALE
                     end do
                     write(1337,99010) SOCII(NINT(SOSCALE))
                  end if
               else
                  write(1337,99001) '< SOC >'
                  write(1337,99003)
               end if
            else
               write(1337,99002) '< SOC >'
               write(1337,99003)
            end if

            if ( MANSOC .AND. (SOSCALE.GE.0D0) ) then
               IMANSOC(1:NATYP) = 1
               !-------------------------------------------------------------------
               ! Now look for a possible include/exclude list (SOCLIST= +/- NASOC)
               ! if SOCLIST is not found, ALL the atoms will have SOC modified with
               ! SOSCALE (+NASOC=only NASOC atoms, -NASOC=all but these NASOC atoms)
               ! Note that this is allowed only for SOC-I manipulation
               !-------------------------------------------------------------------
               call IOInput('SOCLIST         ',UIO,0,7,IER)
               if (IER.EQ.0) then
                  read(UNIT=UIO,FMT=*) NASOC,(ISP(I),I=1,ABS(NASOC))

                  if (NASOC.NE.0) then
                     if (NASOC.LT.0) then ! exclude this atoms
                        do I=1,-NASOC
                           IMANSOC(ISP(I)) = 0
                        end do
                     else
                        IMANSOC(1:NATYP) = 0
                        do I=1,NASOC
                           IMANSOC(ISP(I)) = 1
                        end do
                     end if
                  end if
               end if

               write(1337,2100)
               do I=1,NATYP
                  if (IMANSOC(I).EQ.1) then
                     SOCSCL(1:LMAX+1,I)=SOSCALE
                  end if
               end do
               write(1337,99004)
               if (NASOC.EQ.0) write(1337,99005)
               if (NASOC.GT.0) then
                  write(1337,99006)
                  write(1337,99008) (ISP(I),I=1,NASOC)
               end if
               if (NASOC.LT.0) then
                  write(1337,99007)
                  write(1337,99008) (ISP(I),I=1,ABS(NASOC))
               end if
               write(1337,99009) SOSCALE
               write(1337,2100)
            end if
         end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! SOC-MAN
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         write(1337,'('' SOLVER used for the DIRAC equation : '',2X,A)') SOLVER
         write(1337,2100)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! CTL-MAN
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (OPT('CSCALE  ')) then
            call IOInput('CTLSCALE        ',UIO,0,7,IER)
            if (IER.EQ.0) then
               read (UNIT=UIO,FMT=*) CTLSCALE
               if (CTLSCALE.GE.1D-12) then
                  MANCTL=.TRUE.
               else
                  write(1337,99001) '< CSCALE >'
                  write(1337,99011)
               end if
            else
               write(1337,99002) '< CSCALE >'
               write(1337,99011)
            end if

            if (MANCTL) then
               CSCL(1:LMAX+1,1:NATYP) = CSCL(1:LMAX+1,1:NATYP)/DSQRT(CTLSCALE)
               write(1337,99012)
               write(1337,99005)
               write(1337,99009) 1.D0/DSQRT(CTLSCALE)
            end if
            write(1337,2100)
         end if

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! CTL-MAN
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! LDA+U
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(OPT('qdos    ')) then
         allocate(t_params%qdos_atomselect(NATYP), stat=i_stat) !INTEGER
         call memocc(i_stat,product(shape(t_params%qdos_atomselect))*kind(t_params%qdos_atomselect),'t_params%qdos_atomselect','rinput13')

         t_params%qdos_atomselect(1:NATYP) = 1
         !for now this is not used. Later this should be used to speed up the qdos calculations if not all atoms are supposed to be calculated Then if fullinv was not chosen then tmatrix is only needed for the principle layer of the atom of interest and the calculation of G(k) can be done only on that subblock.
         !          call IoInput('qdosatoms       ',UIO,1,7,IER)
         !          if (IER.EQ.0) then
         !            read (UNIT=UIO,FMT=*) (t_params%qdos_atomselect(I),I=1,NATYP)
         !            write(111,FMT='(A10,80I2)') 'qdosatoms=  ', (t_params%qdos_atomselect(I),I=1,NATYP)
         !          else
         !            write(111,FMT='(A18,80I2)') 'Default qdosatoms=  ', (t_params%qdos_atomselect(I),I=1,NATYP)
         !          endif
         !
         !          write (1337,'(A)') 'atom selective writeout for qdos:'
         !          write (1337,'(A,1000I5)') 'qdosatoms=',  (t_params%qdos_atomselect(I),I=1,NATYP)

      end if

      ! =============================================================         ! fswrt
      ! check and correct some settings automatically for FERMIOUT writeout   ! fswrt
           if(OPT('FERMIOUT').or.OPT('OPERATOR')) then                       ! fswrt
             if (NSTEPS/=1) then                                             ! fswrt
               write(6,2012)                                                 ! fswrt
               NSTEPS = 1                                                    ! fswrt
             end if                                                          ! fswrt
             if (.not.TEST('STOP1B  ')) CALL ADDTEST('STOP1B  ')             ! fswrt
             if (.not.TEST('STOP1B  ')) stop 'addtest failed for STOP1B'     ! fswrt
           end if                                                            ! fswrt
      ! =============================================================         ! fswrt

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! WF_SAVE
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call IOInput('MEMWFSAVE       ',UIO,0,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) t_wavefunctions%maxmem_number
         write(1337,*) '< MEMWFSAVE >', t_wavefunctions%maxmem_number
         write(111,*) 'MEMWFSAVE=',t_wavefunctions%maxmem_number
      else
         t_wavefunctions%maxmem_number = 0
         write(1337,*) '< MEMWFSAVE >, use default:', t_wavefunctions%maxmem_number
         write(111,*) 'Default MEMWFSAVE= ',t_wavefunctions%maxmem_number
      end if
      call IOInput('UNITMEMWFSAVE   ',UIO,0,7,IER)
      if (IER.EQ.0) then
         read (UNIT=UIO,FMT=*) t_wavefunctions%maxmem_units
         write(1337,*) '< UNITMEMWFSAVE >', t_wavefunctions%maxmem_units, ' (max memory= UNITMEMWFSAVE*1024**MEMWFSAVE)'
         write(111,*) 'UNITMEMWFSAVE=',t_wavefunctions%maxmem_units
      else
         t_wavefunctions%maxmem_units = 2
         write(1337,*) '< UNITMEMWFSAVE >, use default:', t_wavefunctions%maxmem_units, '(MB) (max memory= MEMWFSAVE*1024**UNITMEMWFSAVE)'
         write(111,*) 'Default UNITMEMWFSAVE= ',t_wavefunctions%maxmem_units, '(MB)'
      end if

      !default flags: save only rll from main1a>tmatnewsolver since left solutions can be calculated always in main1c>rhovalnew and sll is not used
      t_wavefunctions%save_rll     = .true.
      t_wavefunctions%save_sll     = .false.
      t_wavefunctions%save_rllleft = .false.
      t_wavefunctions%save_sllleft = .false.

      if(opt('OPERATOR')) then
         write(1337,*) 'Found option "OPERATOR"'
         write(1337,*) 'Overwrite MEMWFSAVE with big numbers'
         t_wavefunctions%maxmem_number = 5
         t_wavefunctions%maxmem_units = 3
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! WF_SAVE
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      write(1337,2100)
      write(1337,2040) KMROT
      write(1337,2110)
      write(1337,*) ' >>>>>>>>> RINPUT13 EXITS NOW <<<<<<<<<< '

      close(111) ! Close file inputcard_generated.txt

      if (allocated(IMANSOC)) then
         i_all=-product(shape(IMANSOC))*kind(IMANSOC)
         deallocate(IMANSOC,stat=i_stat)
         call memocc(i_stat,i_all,'IMANSOC','rinput13')
      endif

      return
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! INPUT END
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      1029 format((F4.0,I4,4x,4I1,3I4,F8.4,I4,I5,1x,f8.5))
      ! ------------------------------------------------------------------------
      2010 format(' NSPIN '/I4)
      2011 format(' NSTEPS'/I4)
      2012 format(' WARINING: Setting NSTEPS to 1 for runoption FERMOUT')
      2014 format('          ALAT = ',F15.8)
      2015 format('   INTERVX   INTERVY   INTERVZ'/3I10)
      2016 format('    NCLS    NREF   NINEQ'/,3I8)
      2018 format(' RBASIS'/,&
      &     'SITE                BASIS VECTORS                 ',&
      &     'THETA   PHI CPA OCC KAOEZ')
      2019 format('         ABASIS         BBASIS         CBASIS'/3F15.8)
      2021 format(' INIPOL'/,(10I4))
      2022 format(' IXIPOL'/,(10I4))
      2023 format('    NAEZ    NEMB  '/,2I8)
      2025 format((i4,3F15.8,2F6.1,2(1x,I3),4I3))
      2028 format(' NATYP '/,I4/,&
      &     '   Z lmx     KFG cls pot ntc  MTFAC irns SITE  CONC')
      2031 format((3F15.8,2I6))
      2032 format(' NTCELLR'/,(10I4))
      2040 format(' KMROT'/,4I8)
      ! ------------------------------------------------------------------------
      2100 format(79(1H-))
      2101 format(   3(1H-),1H+  , 3(14(1H-),1H+),  30(1H-))
      2102 format( 3(9(1H-),1H+) ,49(1H-))
      2103 format(10(3(1H-),1H+) ,39(1H-))
      2104 format(   3(1H-),1H+  ,75(1H-))
      2107 format( 3(14(1H-),1H+),34(1H-))
      2108 format( 2(3(1H-),1H+),  7(1H-),1H+,      3(3(1H-),1H+),&
      &          7(1H-),1H+,   3(1H-),1H+,      39(1H-))
      2110 format( 3(7(1H-),1H+) ,55(1H-))
      2111 format( 7(7(1H-),1H+) ,23(1H-))
      9020 format (/,33x,'check of dimension-data consistency',/,33x,&
      &       35 ('-'),/,40x,'lmax   : (',i6,',',i6,')',/,40x,&
      &       'natyp  : (',i6,',',i6,')',/,40x,'irm    : (',i6,',',i6,&
      &       ')',/,40x,'nspin  : (',i6,',',i6,')',/)
      9030 format (1x,10 ('*'),' external magnetic field applied hfield=',&
      &       f8.5)
      9050 format (20x,a4,'spin polarized calculation')
      9070 format (1x,20x,' calculation with',a8,'-potential')
      9080 format (1x,79 ('*'))
      9090 format (' mixing factor used           :',f15.6,/,&
      &        ' convergence quality required :',1p,d15.2)
      9091 format (' make use of CPA algorithm    :',1x,a14)
      9092 format ('         max. iterations      :',i15,/,&
      &        '         req. CPA convergency :',1p,d15.2)
      9100 format (1x,20x,a24,'exchange-correlation potential')
      9110 format (/,20x,'broyden"s method # :',i3,&
      &       ' is used up to iteration-      ',/,20x,'depth :',i3,&
      &       '  then jacobian is fixed and potential      ',/,20x,&
      &       'is updated using that jacobian')
      9120 format (13x,' in case of calculating non - spherical wavefcts ',&
      &       'the parameter lmaxd has to be set equal lmax ')
      9130 format (/)
      9140 format (20x,'full potential calculation ',&
      &       '- cut off of non spherical potential',/,' >',/)
      9150 format (31x,'representive atom no.',i3,' irns :',i5,' irnsd :',i5)
      9160 format (21x,a43,/,21x,' using',i3,'-th. born approximation ')
      9170 format (21x,a43)
      9210 format (' lmax'/,i4)
      9220 format ('          EMIN        EMAX        TK'/,3f12.6)
      9230 format ('   NPOL  NPNT1  NPNT2  NPNT3'/,4i7)
      9250 format ('  IFILE    IPE ISHIFT ESHIFT'/,3i7,f12.6)
      9260 format (' KSHAPE    IRM    INS   ICST INSREF'/,5i7)
      9270 format ('   KCOR  KVREL    KWS   KHYP KHFIELD   KXC'/,6i7)
      9280 format (' external magnetic hfield     :',f15.4/,&
      &        ' VCONST                       :',f15.6)
      9290 format ('   IMIX    IGF    ICC'/,3i7)
      9300 format (' ITDBRY'/,i7)
      9310 format ('      STRMIX        FCM       QBOUND'/,3f12.6)
      9320 format ('      BRYMIX'/,f12.6)
      9330 format ('    KTE   KPRE   KEFG  KVMAD '/,5i7)
      9301 format(   3(1H-),1H+  ,75(1H-))
      9302 format( 3(11(1H-),1H+),43(1H-))
      9303 format(3(6(1H-),1H+) ,58(1H-))
      9304 format(4(6(1H-),1H+) ,51(1H-))
      9305 format(3(6(1H-),1H+),11(1H-),1H+ ,46(1H-))
      9306 format(6(6(1H-),1H+) ,37(1H-))
      9307 format(6(1H-),1H+,72(1H-))
      9308 format(11(1H-),1H+,67(1H-))
      9309 format(5(6(1H-),1H+) ,44(1H-))
      9410 format('*** SLAB - INTERFACE CALCULATION ***'/)
      9420 format(I5,3F14.8,I5)
      9430 format('Number of LEFT  Host Layers : ',I5,' with ',I5,' basis')
      9440 format('Number of RIGHT Host Layers : ',I5,' with ',I5,' basis')
      9450 format('Left  side periodicity : ',3F10.5)
      9460 format('Right side periodicity : ',3F10.5)
      9465 format('    Geommetry used : '/,&
      &       ' ATOM       TX          TY          TZ ')
      9470 format('--------------- Left  Host -------------- ')
      9475 format('---------------   S L A B  -------------- ')
      9480 format('--------------- Right Host -------------- ')
      99001 format(/,1X,&
      &     "WARNING: Option ",A," used with an INVALID ",&
      &     "scaling parameter.")
      99002 format(/,1X,&
      &     "WARNING: Option ",A," found but NO value given for the",&
      &     " scaling parameter.")
      99003 format(15X,'++++++++++   SOC option will be IGNORED   ++++++++++',&
      &     /,1X,'Please use SOCSCALE= XXX (real>-2.5) in the inputcard',&
      &     ' to make your option valid ',/)
      99004 format(1X,'The SOC will be SCALED',$)
      99005 format(' for ALL the atoms in the unit cell.')
      99006 format(' for the FOLLOWING atoms in the unit cell :')
      99007 format(' for all the atoms in the unit cell EXCLUDING :')
      99008 format(1X,6(2X,I3))
      99009 format(1X,'Scaling factor = ',1P,D9.2)
      99010 format(1X,'The SOC is manipulated',' -- part of the SOC kept: ',A)
      99011 format(15X,'+++++++++  CSCALE option will be IGNORED  ++++++++++',&
      &     /,1X,'Please use CTLSCALE= X (real>=1D-12) in the inputcard',&
      &     ' to make your option valid ',/)
      99012 format(1X,'The CLIGHT will be SCALED',$)

   end subroutine RINPUT13
   !---------------------------------------------------------------------
   !---------------------------------------------------------------------

   subroutine ADDOPT(STRING)
      use mod_wunfiles, only: t_params

      implicit none

      integer, parameter :: NOPTD=32
      character(len=8) :: STRING
      integer :: II
      logical, external :: OPT

      if (.not.OPT('        ')) then
         write(*,*) 'Error in ADDOPT for ',STRING,' : No free slots in array OPTC.'
         stop 'Error in ADDOPT: No free slots in array OPTC.'
      endif

      if (.not.OPT(STRING)) then
         II = 1
         do while (II.le.NOPTD)
            if (t_params%OPTC(II).eq.'        ') then
               t_params%OPTC(II) = STRING
               II = NOPTD + 1
            endif
            II = II + 1
         enddo
      endif

   end subroutine ADDOPT

   subroutine ADDTEST(STRING)
      use mod_types, only: t_inc
      use mod_wunfiles, only: t_params
      implicit none
      integer, parameter :: NTSTD=64
      character(len=8) :: STRING
      integer :: II
      logical, external :: TEST
       
      if(t_inc%i_write) write(1337, *) 'in ADDTEST: adding option ', STRING
    
      if (.not.TEST('        ')) then
        write(*,*) 'Error in ADDTEST for ',STRING,' : No free slots in array TESTC.'
        stop 'Error in ADDTEST: No free slots in array TESTC.'
      endif
       
      if (.not.TEST(STRING)) then
        II = 1
        do while (II.le.NTSTD)
           if (t_params%TESTC(II).eq.'        ') then
              t_params%TESTC(II) = STRING
              II = NTSTD + 1
           endif
           II = II + 1
        enddo
      endif
     
   end subroutine ADDTEST

end module rinput
