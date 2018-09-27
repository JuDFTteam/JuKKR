!-------------------------------------------------------------------------------
!> Summary: 
!> Author: 
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!>
!> 
!-------------------------------------------------------------------------------
! -------------------------------------------------------------------------------
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
! -------------------------------------------------------------------------------
#ifdef CPP_HYBRID
#define CPP_OMPSTUFF
#endif
#ifdef CPP_OMP
#define CPP_OMPSTUFF
#endif

module mod_main0


  use :: mod_datatypes, only: dp

#ifdef CPP_TIMING
  use :: mod_timing
#endif
  use :: mod_wunfiles
  use :: mod_types, only: t_imp
  use :: mod_constants
  use :: memoryhandling
  use :: global_variables
  use :: mod_create_newmesh
  use :: mod_rhoqtools, only: rhoq_save_rmesh
  use :: rinput
  use :: mod_addvirtual14
  use :: mod_bzkint0
  use :: mod_calcrotmat
  use :: mod_changerep
  use :: mod_cinit
  use :: mod_clsgen_tb
  use :: mod_convol
  use :: mod_deciopt
  use :: mod_drvbastrans
  use :: mod_epathtb
  use :: mod_gaunt2
  use :: mod_gaunt
  use :: mod_generalpot
  use :: mod_getbr3
  use :: mod_gfmask
  use :: mod_lattix99
  use :: mod_madelung2d
  use :: mod_madelung3d
  use :: mod_outpothost
  use :: mod_outtmathost
  use :: mod_readimppot
  use :: mod_relpotcvt
  use :: mod_rinit
  use :: mod_scalevec
  use :: mod_setgijtab
  use :: mod_shape_corr
  use :: mod_startb1
  use :: mod_startldau
  use :: mod_testdim
  use :: mod_write_tbkkr_files
  use :: mod_writehoststructure

  implicit none

  private
  public :: main0, bshift_ns
  ! ------------- > scalars > ------------- 
  !integers
  public :: kte, kws, kxc, igf, icc, ins, irm, ipe, ipf, ipfe, kcor, kefg, khyp, kpre, nprinc, nsra, lpot, imix, iend
  public :: icst, naez, nemb, lmax, ncls, nref, npol, npnt1, npnt2, npnt3, lmmax, nvirt, lmpot, kvmad, itscf, ncheb, nineq
  public :: natyp, ifile, kvrel, nspin, nleft, nright, invmod, khfeld, itdbry, insref, kshape, ielast, ishift, kfrozn, nsymat
  public :: nqcalc, kforce, n1semi, n2semi, n3semi, nlayer, nlbasis, nrbasis, intervx, intervy, intervz, maxmesh, npan_eq
  public :: npan_log, npolsemi, scfsteps, natomimp, iesemicore, idosemicore
  !real(kind=dp)
  public :: tk, fcm, e2in, emin, emax, alat, rmax, gmax, r_log, rcutz, rcutxy, qbound, vconst, hfield, mixing, abasis, bbasis
  public :: cbasis, efermi, eshift, tksemi, tolrdif, alatnew, volume0, emusemi, ebotsemi, fsemicore, lambda_xc
  !character
  public :: solver, i12, i13, i19, i25, i40
  !logicals
  public :: lrhosym, linipol, lcartesian
  ! ------------- < scalars < ------------- 
  ! ------------- > arrays > ------------- 
  !integer
  public :: isymindex, cls, irc, imt, nfu, nsh1, nsh2, lmxc, ipan, irns, irws, kmesh, irmin, loflm, nacls, ncore, imaxsh, nshell
  public :: inipol, ixipol, refpot, ntcell, iqcalc, iofgij, jofgij, atomimp, ijtabsh, ijtabsym, npan_tot, ijtabcalc, npan_eq_at
  public :: npan_log_at, ijtabcalc_i, ish, jsh, ilm_map, kfg, atom, ezoa, lmsp, lcore, icleb, ircut, llmsp, lmsp1, kaoez, ifunm
  public :: ifunm1, ititle, icheck, ipan_intervall, jend, kmrot, ncpa, itcpamax, noq, iqat, icpa, hostimp, zrel, jwsrel, irshift
  public :: nrrel, ntldau, idoldau, itrunldau, kreadldau, lopt, itldau, lly, ivshift, irrel
  !real
  public :: vbc, zperight, zperleft, recbv, bravais, rsymat, a, b, wg, gsh, zat, rmt, rws, vref, vref_temp, mtfac, rmtnew, rmtref
  public :: rmtref_temp, rmtrefat, fpradius, socscale, rmesh, s, rr, drdi, dror, cleb, visp, cscl, rnew, ratom, ecore, tleft, tright
  public :: socscl, rbasis, rclsimp, cmomhost, rpan_intervall, rs, yrg, vins, rcls, rrot, qmtet, qmphi, qmgam, qmgamtab, qmphitab, qmtettab
  public :: cpatol, conc, fact, vtrel, btrel, rmrel, drdirel, r2drdirel, thesme, thetas, thetasnew, ueff, jeff, erefldau, wldau, uldau
  !complex
  public :: ez, dez, wez, dsymll, dsymll1, lefttinvll, righttinvll, rc, crel, rrel, srrel, drotq, phildau, deltae
  !character
  public :: txc
  !logical
  public :: vacflag, para, symunitary, emeshfile
  ! ------------- < arrays < ------------- 

  
  ! definition of common variables

  integer :: kte                   !! Calculation of the total energy On/Off (1/0)
  integer :: kws                   !! 0 (MT), 1(ASA)
  integer :: kxc                   !! Type of xc-potential 0=vBH 1=MJW 2=VWN 3=PW91
  integer :: igf                   !! Do not print or print (0/1) the KKRFLEX_* files
  integer :: icc                   !! Enables the calculation of off-diagonal elements of the GF.(0=SCF/DOS; 1=cluster; -1=custom)
  integer :: ins                   !! 0 (MT), 1(ASA), 2(Full Potential)
  integer :: irm                   !! Maximum number of radial points
  integer :: ipe                   !! Not real used, IPFE should be 0
  integer :: ipf                   !! Not real used, IPFE should be 0
  integer :: ipfe                  !! Not real used, IPFE should be 0
  integer :: kcor
  integer :: kefg
  integer :: khyp
  integer :: kpre
  integer :: nprinc
  integer :: nsra
  integer :: lpot                  !! Maximum l component in potential expansion
  integer :: imix                  !! Type of mixing scheme used (0=straight, 4=Broyden 2nd, 5=Anderson)
  integer :: iend                  !! Number of nonzero gaunt coefficients
  integer :: icst                  !! Number of Born approximation
  integer :: naez                  !! Number of atoms in unit cell
  integer :: nemb                  !! Number of 'embedding' positions
  integer :: lmax                  !! Maximum l component in wave function expansion
  integer :: ncls                  !! Number of reference clusters
  integer :: nref                  !! Number of diff. ref. potentials
  integer :: npol                  !! Number of Matsubara Poles (EMESHT)
  integer :: npnt1                 !! number of E points (EMESHT) for the contour integration
  integer :: npnt2                 !! number of E points (EMESHT) for the contour integration
  integer :: npnt3                 !! number of E points (EMESHT) for the contour integration
  integer :: lmmax                 !! (LMAX+1)^2
  integer :: nvirt
  integer :: lmpot                 !! (LPOT+1)**2
  integer :: kvmad
  integer :: itscf
  integer :: ncheb                 !! Number of Chebychev pannels for the new solver
  integer :: nineq                 !! Number of ineq. positions in unit cell
  integer :: natyp                 !! Number of kinds of atoms in unit cell
  integer :: ifile                 !! Unit specifier for potential card
  integer :: kvrel                 !! 0,1,2 : non / scalar relat. / full Dirac calculation
  integer :: nspin                 !! Counter for spin directions
  integer :: nleft                 !! Number of repeated basis for left host to get converged electrostatic potentials
  integer :: nright                !! Number of repeated basis for right host to get converged electrostatic potentials
  integer :: invmod                !! Inversion scheme
  integer :: khfeld                !! 0,1: no / yes external magnetic field
  integer :: itdbry                !! Number of SCF steps to remember for the Broyden mixing
  integer :: insref                !! INS for reference pot. (usual 0)
  integer :: kshape                !! Exact treatment of WS cell
  integer :: ielast
  integer :: ishift
  integer :: kfrozn
  integer :: nsymat
  integer :: nqcalc
  integer :: kforce                !! Calculation of the forces
  integer :: n1semi                !! Number of energy points for the semicore contour
  integer :: n2semi                !! Number of energy points for the semicore contour
  integer :: n3semi                !! Number of energy points for the semicore contour
  integer :: nlayer                !! Number of principal layer
  integer :: nlbasis               !! Number of basis layers of left host (repeated units)
  integer :: nrbasis               !! Number of basis layers of right host (repeated units)
  integer :: intervx               !! Number of intervals in x-direction for k-net in IB of the BZ
  integer :: intervy               !! Number of intervals in y-direction for k-net in IB of the BZ
  integer :: intervz               !! Number of intervals in z-direction for k-net in IB of the BZ
  integer :: maxmesh
  integer :: npan_eq               !! Number of intervals from [R_LOG] to muffin-tin radius Used in conjunction with runopt NEWSOSOL
  integer :: npan_log              !! Number of intervals from nucleus to [R_LOG] Used in conjunction with runopt NEWSOSOL
  integer :: npolsemi              !! Number of poles for the semicore contour
  integer :: scfsteps              !! number of scf iterations
  integer :: natomimp              !! Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
  integer :: iesemicore
  integer :: idosemicore
  real (kind=dp) :: tk             !! Temperature
  real (kind=dp) :: fcm            !! Factor for increased linear mixing of magnetic part of potential compared to non-magnetic part.
  real (kind=dp) :: e2in
  real (kind=dp) :: emin           !! Lower value (in Ryd) for the energy contour
  real (kind=dp) :: emax           !! Maximum value (in Ryd) for the DOS calculation Controls also [NPT2] in some cases
  real (kind=dp) :: alat           !! Lattice constant in a.u.
  real (kind=dp) :: rmax           !! Ewald summation cutoff parameter for real space summation
  real (kind=dp) :: gmax           !! Ewald summation cutoff parameter for reciprocal space summation
  real (kind=dp) :: r_log          !! Radius up to which log-rule is used for interval width. Used in conjunction with runopt NEWSOSOL
  real (kind=dp) :: rcutz          !! Parameter for the screening cluster along the z-direction
  real (kind=dp) :: rcutxy         !! Parameter for the screening cluster along the x-y plane
  real (kind=dp) :: qbound         !! Convergence parameter for the potential
  real (kind=dp) :: vconst         !! Potential shift in the first iteration
  real (kind=dp) :: hfield         !! External magnetic field, for initial potential shift in spin polarised case
  real (kind=dp) :: mixing         !! Magnitude of the mixing parameter
  real (kind=dp) :: abasis         !! Scaling factors for rbasis
  real (kind=dp) :: bbasis         !! Scaling factors for rbasis
  real (kind=dp) :: cbasis         !! Scaling factors for rbasis
  real (kind=dp) :: efermi         !! Fermi energy
  real (kind=dp) :: eshift
  real (kind=dp) :: tksemi         !! Temperature of semi-core contour
  real (kind=dp) :: tolrdif        !! For distance between scattering-centers smaller than [<TOLRDIF>], free GF is set to zero. Units are Bohr radii.
  real (kind=dp) :: alatnew
  real (kind=dp) :: volume0
  real (kind=dp) :: emusemi        !! Top of semicore contour in Ryd.
  real (kind=dp) :: ebotsemi       !! Bottom of semicore contour in Ryd
  real (kind=dp) :: fsemicore      !! Initial normalization factor for semicore states (approx. 1.)
  real (kind=dp) :: lambda_xc      !! Scale magnetic moment (0 < Lambda_XC < 1, 0=zero moment, 1= full moment)
  character (len=10) :: solver                               !! Type of solver
  character (len=40) :: i12                               !! File identifiers
  character (len=40) :: i13                               !! Potential file name
  character (len=40) :: i19                               !! Shape function file name
  character (len=40) :: i25                               !! Scoef file name
  character (len=40) :: i40                               !! File identifiers
  logical :: lrhosym
  logical :: linipol               !! True: Initial spin polarization; false: no initial spin polarization
  logical :: lcartesian            !! True: Basis in cartesian coords; false: in internal coords

  ! ..
  ! .. Arrays ..
  integer, dimension (nsymaxd) :: isymindex
  integer, dimension (:), allocatable :: cls !! Cluster around atomic sites
  integer, dimension (:), allocatable :: irc !! R point for potential cutting
  integer, dimension (:), allocatable :: imt !! R point at MT radius
  integer, dimension (:), allocatable :: nfu !! number of shape function components in cell 'icell'
  integer, dimension (:), allocatable :: nsh1 !! Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
  integer, dimension (:), allocatable :: nsh2 !! Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
  integer, dimension (:), allocatable :: lmxc
  integer, dimension (:), allocatable :: ipan !! Number of panels in non-MT-region
  integer, dimension (:), allocatable :: irns !! Position of atoms in the unit cell in units of bravais vectors
  integer, dimension (:), allocatable :: irws !! R point at WS radius
  integer, dimension (:), allocatable :: kmesh
  integer, dimension (:), allocatable :: irmin !! Max R for spherical treatment
  integer, dimension (:), allocatable :: loflm !! l of lm=(l,m) (GAUNT)
  integer, dimension (:), allocatable :: nacls !! Number of atoms in cluster
  integer, dimension (:), allocatable :: ncore !! Number of core states
  integer, dimension (:), allocatable :: imaxsh
  integer, dimension (:), allocatable :: nshell !! Index of atoms/pairs per shell (ij-pairs); nshell(0) = number of shells
  integer, dimension (:), allocatable :: inipol !! Initial spin polarisation
  integer, dimension (:), allocatable :: ixipol !! Constraint of spin pol.
  integer, dimension (:), allocatable :: refpot !! Ref. pot. card  at position
  integer, dimension (:), allocatable :: ntcell !! Index for WS cell
  integer, dimension (:), allocatable :: iqcalc
  integer, dimension (:), allocatable :: iofgij !! Linear pointers, similar to NSH1/NSH2 but giving the actual index of sites I,J = 1,NATOMIMP in the cluster
  integer, dimension (:), allocatable :: jofgij !! Linear pointers, similar to NSH1/NSH2 but giving the actual index of sites I,J = 1,NATOMIMP in the cluster
  integer, dimension (:), allocatable :: atomimp
  integer, dimension (:), allocatable :: ijtabsh !! Linear pointer, assigns pair (i,j) to a shell in the array GS(*,*,*,NSHELD)
  integer, dimension (:), allocatable :: ijtabsym !! Linear pointer, assigns pair (i,j) to the rotation bringing GS into Gij
  integer, dimension (:), allocatable :: npan_tot
  integer, dimension (:), allocatable :: ijtabcalc !! Linear pointer, specifying whether the block (i,j) has to be calculated needs set up for ICC=-1, not used for ICC=1
  integer, dimension (:), allocatable :: npan_eq_at
  integer, dimension (:), allocatable :: npan_log_at
  integer, dimension (:), allocatable :: ijtabcalc_i
  integer, dimension (:, :), allocatable :: ish
  integer, dimension (:, :), allocatable :: jsh
  integer, dimension (:, :), allocatable :: ilm_map
  integer, dimension (:, :), allocatable :: kfg
  integer, dimension (:, :), allocatable :: atom !! Atom at site in cluster
  integer, dimension (:, :), allocatable :: ezoa !! EZ of atom at site in cluster
  integer, dimension (:, :), allocatable :: lmsp !! 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
  integer, dimension (:, :), allocatable :: lcore !! Angular momentum of core states
  integer, dimension (:, :), allocatable :: icleb !! Pointer array
  integer, dimension (:, :), allocatable :: ircut !! R points of panel borders
  integer, dimension (:, :), allocatable :: llmsp !! lm=(l,m) of 'nfund'th nonvanishing component of non-spherical pot.
  integer, dimension (:, :), allocatable :: lmsp1
  integer, dimension (:, :), allocatable :: kaoez !! Kind of atom at site in elem. cell
  integer, dimension (:, :), allocatable :: ifunm
  integer, dimension (:, :), allocatable :: ifunm1
  integer, dimension (:, :), allocatable :: ititle
  integer, dimension (:, :), allocatable :: icheck
  integer, dimension (:, :), allocatable :: ipan_intervall
  integer, dimension (:, :, :), allocatable :: jend !! Pointer array for icleb()
  real (kind=dp), dimension (2) :: vbc !! Potential constants
  real (kind=dp), dimension (3) :: zperight !! Vector to define how to repeat the basis of the right host
  real (kind=dp), dimension (3) :: zperleft !! Vector to define how to repeat the basis of the left host
  real (kind=dp), dimension (3, 3) :: recbv !! Reciprocal basis vectors
  real (kind=dp), dimension (3, 3) :: bravais !! Bravais lattice vectors
  real (kind=dp), dimension (64, 3, 3) :: rsymat
  real (kind=dp), dimension (:), allocatable :: a !! Constants for exponential R mesh
  real (kind=dp), dimension (:), allocatable :: b !! Constants for exponential R mesh
  real (kind=dp), dimension (:), allocatable :: wg !! Integr. weights for Legendre polynomials
  real (kind=dp), dimension (:), allocatable :: gsh
  real (kind=dp), dimension (:), allocatable :: zat !! Nuclear charge
  real (kind=dp), dimension (:), allocatable :: rmt !! Muffin-tin radius of true system
  real (kind=dp), dimension (:), allocatable :: rws !! Wigner Seitz radius
  real (kind=dp), dimension (:), allocatable :: vref
  real (kind=dp), dimension (:), allocatable :: vref_temp
  real (kind=dp), dimension (:), allocatable :: mtfac !! Scaling factor for radius MT
  real (kind=dp), dimension (:), allocatable :: rmtnew !! Adapted muffin-tin radius
  real (kind=dp), dimension (:), allocatable :: rmtref !! Muffin-tin radius of reference system
  real (kind=dp), dimension (:), allocatable :: rmtref_temp !! Muffin-tin radius of reference system
  real (kind=dp), dimension (:), allocatable :: rmtrefat
  real (kind=dp), dimension (:), allocatable :: fpradius !! R point at which full-potential treatment starts
  real (kind=dp), dimension (:), allocatable :: socscale !! Spin-orbit scaling
  real (kind=dp), dimension (:, :), allocatable :: rmesh !! Radial mesh ( in units a Bohr)
  real (kind=dp), dimension (:, :), allocatable :: s
  real (kind=dp), dimension (:, :), allocatable :: rr !! Set of real space vectors (in a.u.)
  real (kind=dp), dimension (:, :), allocatable :: drdi !! Derivative dr/di
  real (kind=dp), dimension (:, :), allocatable :: dror
  real (kind=dp), dimension (:, :), allocatable :: cleb !! GAUNT coefficients (GAUNT)
  real (kind=dp), dimension (:, :), allocatable :: visp !! Spherical part of the potential
  real (kind=dp), dimension (:, :), allocatable :: cscl !! Speed of light scaling
  real (kind=dp), dimension (:, :), allocatable :: rnew
  real (kind=dp), dimension (:, :), allocatable :: ratom
  real (kind=dp), dimension (:, :), allocatable :: ecore !! Core energies
  real (kind=dp), dimension (:, :), allocatable :: tleft !! Vectors of the basis for the left host
  real (kind=dp), dimension (:, :), allocatable :: tright !! Vectors of the basis for the right host
  real (kind=dp), dimension (:, :), allocatable :: socscl
  real (kind=dp), dimension (:, :), allocatable :: rbasis !! Position of atoms in the unit cell in units of bravais vectors
  real (kind=dp), dimension (:, :), allocatable :: rclsimp
  real (kind=dp), dimension (:, :), allocatable :: cmomhost !! Charge moments of each atom of the (left/right) host
  real (kind=dp), dimension (:, :), allocatable :: rpan_intervall
  real (kind=dp), dimension (:, :, :), allocatable :: rs
  real (kind=dp), dimension (:, :, :), allocatable :: yrg !! Spherical harmonics (GAUNT2)
  real (kind=dp), dimension (:, :, :), allocatable :: vins !! Non-spherical part of the potential
  real (kind=dp), dimension (:, :, :), allocatable :: rcls !! Real space position of atom in cluster
  real (kind=dp), dimension (:, :, :), allocatable :: rrot
  complex (kind=dp), dimension (:), allocatable :: ez
  complex (kind=dp), dimension (:), allocatable :: dez
  complex (kind=dp), dimension (:), allocatable :: wez
  complex (kind=dp), dimension (:, :, :), allocatable :: dsymll
  complex (kind=dp), dimension (:, :, :), allocatable :: dsymll1
  complex (kind=dp), dimension (:, :, :, :, :), allocatable :: lefttinvll
  complex (kind=dp), dimension (:, :, :, :, :), allocatable :: righttinvll
  character (len=124), dimension (6) :: txc
  logical, dimension (2) :: vacflag

  ! -------------------------------------------------------------------------
  ! Magnetisation angles -- description see RINPUT13
  ! -------------------------------------------------------------------------
  integer :: kmrot                 !! 0: no rotation of the magnetisation; 1: individual rotation of the magnetisation for every site
  real (kind=dp), dimension (:), allocatable :: qmtet !! \f$ \theta\f$ angle of the agnetization with respect to the z-axis
  real (kind=dp), dimension (:), allocatable :: qmphi !! \f$ \phi\f$ angle of the agnetization with respect to the z-axis
  ! -------------------------------------------------------------------------
  ! CPA variables
  ! -------------------------------------------------------------------------
  integer :: ncpa                  !! NCPA = 0/1 CPA flag
  integer :: itcpamax              !! Max. number of CPA iterations
  integer, dimension (:), allocatable :: noq !! Number of diff. atom types located
  integer, dimension (:), allocatable :: iqat !! The site on which an atom is located on a given site
  integer, dimension (:), allocatable :: icpa !! ICPA = 0/1 site-dependent CPA flag

  ! -------------------------------------------------------------------------
  !> @note ITERMDIR running option introduced Apr 2003 -- Munich
  !>              (H. Ebert + V. Popescu) allows a self-consistent
  !>              determination of the magnetic configuration in REL mode
  ! -------------------------------------------------------------------------
  real (kind=dp), dimension (:), allocatable :: qmgam
  real (kind=dp), dimension (:, :), allocatable :: qmgamtab
  real (kind=dp), dimension (:, :), allocatable :: qmphitab
  real (kind=dp), dimension (:, :), allocatable :: qmtettab
  ! -------------------------------------------------------------------------
  !> @note changes for impurity 20/02/2004 -- v.popescu according to
  !>                                          n.papanikolaou VINS()
  ! -------------------------------------------------------------------------
  integer, dimension (:), allocatable :: hostimp
  real (kind=dp) :: cpatol         !! Convergency tolerance for CPA-cycle
  real (kind=dp), dimension (:), allocatable :: conc !! Concentration of a given atom
  ! -------------------------------------------------------------------------------
  complex (kind=dp), dimension (:, :), allocatable :: rc !! NREL REAL spher. harm. > CMPLX. spher. harm. NREL CMPLX. spher. harm. > REAL spher. harm.
  complex (kind=dp), dimension (:, :), allocatable :: crel !! Non-relat. CMPLX. spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. CMPLX. spher. harm.
  complex (kind=dp), dimension (:, :), allocatable :: rrel !! Non-relat. REAL spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. REAL spher. harm.
  complex (kind=dp), dimension (:, :, :), allocatable :: srrel
  complex (kind=dp), dimension (:, :, :), allocatable :: drotq !! Rotation matrices to change between LOCAL/GLOBAL frame of reference for magnetisation <> Oz or noncollinearity
  integer, dimension (:), allocatable :: zrel !! atomic number (cast integer)
  integer, dimension (:), allocatable :: jwsrel !! index of the WS radius
  integer, dimension (:), allocatable :: irshift !! shift of the REL radial mesh with respect no NREL
  integer, dimension (:, :), allocatable :: nrrel
  integer, dimension (:, :, :), allocatable :: irrel
  real (kind=dp), dimension (0:100) :: fact
  real (kind=dp), dimension (:, :), allocatable :: vtrel !! potential (spherical part)
  real (kind=dp), dimension (:, :), allocatable :: btrel !! magnetic field
  real (kind=dp), dimension (:, :), allocatable :: rmrel !! radial mesh
  real (kind=dp), dimension (:, :), allocatable :: drdirel !! derivative of radial mesh
  real (kind=dp), dimension (:, :), allocatable :: r2drdirel !! \f$ r^2 \frac{\partial}{\partial \mathbf{r}}\frac{\partial}{\partial i}\f$ (r**2 * drdi)
  real (kind=dp), dimension (:, :, :), allocatable :: thesme
  logical :: para
  logical, dimension (nsymaxd) :: symunitary !! unitary/antiunitary symmetry flag

  ! -------------------------------------------------------------------------
  ! LDA+U LDA+U LDA+U
  ! -------------------------------------------------------------------------
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

  ! -------------------------------------------------------------------------
  integer :: ntldau                !! number of atoms on which LDA+U is applied
  integer :: idoldau               !! flag to perform LDA+U
  integer :: itrunldau             !! Iteration index for LDA+U
  integer :: kreadldau             !! LDA+U arrays available
  integer, dimension (:), allocatable :: lopt !! angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
  integer, dimension (:), allocatable :: itldau !! integer pointer connecting the NTLDAU atoms to heir corresponding index in the unit cell
  real (kind=dp), dimension (:), allocatable :: ueff !! input U parameter for each atom
  real (kind=dp), dimension (:), allocatable :: jeff !! input J parameter for each atom
  real (kind=dp), dimension (:), allocatable :: erefldau !! the energies of the projector's wave functions (REAL)
  ! ..
  ! .. distinguish between spin-dependent and spin-independent
  ! .. quantities
  real (kind=dp), dimension (:, :, :, :), allocatable :: wldau !! potential matrix
  real (kind=dp), dimension (:, :, :, :, :), allocatable :: uldau !! calculated Coulomb matrix elements (EREFLDAU)
  complex (kind=dp), dimension (:, :), allocatable :: phildau
  ! -------------------------------------------------------------------------
  ! LDA+U LDA+U LDA+U
  ! -------------------------------------------------------------------------
  ! Lloyds formula
  integer :: lly                   !! LLY <> 0 : apply Lloyds formula
  complex (kind=dp) :: deltae      !! Energy difference for numerical derivative

  ! SUSC (BEGIN: modifications by Manuel and Benedikt)             ! susc
  ! LOGICAL THAT CHECKS WHETHER ENERGY MESH FILE EXISTS            ! susc
  logical :: emeshfile             ! susc
  ! SUSC (END:   modifications by Manuel and Benedikt)             ! susc

  ! ruess: IVSHIFT test option
  integer :: ivshift

  ! allocations:
  real (kind=dp), dimension (:, :, :), allocatable :: thetas !! shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
  real (kind=dp), dimension (:, :, :), allocatable :: thetasnew


contains


  !-------------------------------------------------------------------------------
  !> Summary: 
  !> Author: 
  !> Category: KKRhost, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> 
  !-------------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------
  ! SUBROUTINE: main0
  !> @brief Main wrapper to handle input reading, allocation of arrays, and
  !> preparation of all the necessary data structures for a calculation.
  !> @author Philipp Rüssmann, Bernd Zimmermann, Phivos Mavropoulos, R. Zeller,
  !> and many others ...
  !> @note Jonathan Chico: Re-wrote the module from Fixed format Fortran to
  !> fortran90 Free Form. Also performed modifications to get rid of the inc.p
  !> file. 22.12.2017
  ! ----------------------------------------------------------------------------
  subroutine main0()

    use :: mod_types
#ifdef CPP_OMPSTUFF
    use :: omp_lib                 ! necessary for omp functions
#endif
#ifdef CPP_MPI
    use :: mpi
#endif
    use :: mod_mympi, only: nranks
    use :: mod_version
    use :: mod_version_info
    use :: mod_md5sums

    implicit none
    ! .. Local Scalars ..
    integer :: i
    integer :: j
    integer :: i1
    integer :: ie
    integer :: lm
    integer :: ns
    integer :: isvatom, nvatom
    integer :: i_stat
    integer :: irec
    integer :: lrecabmad
    real (kind=dp) :: zattemp
    integer :: ierr
    real (kind=dp), allocatable :: tmp_rr(:, :)
    ! for OPERATOR option
    logical :: lexist, operator_imp

#ifdef CPP_OMPSTUFF
    ! .. OMP ..
    integer :: nth, ith            ! total number of threads, thread number
#endif
    ! ..
    ! .. External Functions ..
    logical :: opt, test
    external :: opt, test
    ! -------------------------------------------------------------------------
    ! Write version info:
    ! -------------------------------------------------------------------------
    call print_versionserial(6, version1, version2, version3, version4, serialnr)
    call print_versionserial(1337, version1, version2, version3, version4, serialnr)

#ifdef CPP_OMPSTUFF
    ! $omp parallel shared(nth) private(ith)
    ith = omp_get_thread_num()
    if (ith==0) then
      nth = omp_get_num_threads()
      write (*, '(/79("*")//1X,A,I5//79("*")/)') 'Number of OpenMP threads used:', nth
      write (1337, '(1X,A,I5)') 'Number of OpenMP threads used:', nth
    end if
    ! $omp end parallel
#endif

#ifdef CPP_MPI
    write (*, '(1X,A,I5//79("*")/)') 'Number of MPI ranks used:', nranks
    write (1337, '(1X,A,I5//79("*")/)') 'Number of MPI ranks used:', nranks
#endif
    ! -------------------------------------------------------------------------
    ! End write version info
    ! -------------------------------------------------------------------------


    ! allocate and initialize default values
    call init_all_wrapper()


    ! -------------------------------------------------------------------------
    ! Reading of the inputcard, and allocation of several arrays
    !> @note JC: have added reading calls for the parameters that used to be in
    !> the inc.p and can now be modified via the inputcard directly
    ! -------------------------------------------------------------------------
    call rinput13(kte, igf, kxc, lly, icc, ins, kws, ipe, ipf, ipfe, icst, imix, lpot, naez, nemb, nref, ncls, npol, lmax, kcor, kefg, khyp, kpre, kvmad, lmmax, lmpot, ncheb, &
      nleft, ifile, kvrel, nspin, natyp, nineq, npnt1, npnt2, npnt3, kfrozn, ishift, n1semi, n2semi, n3semi, scfsteps, insref, kshape, itdbry, nright, kforce, ivshift, khfeld, &
      nlbasis, nrbasis, intervx, intervy, intervz, npan_eq, npan_log, npolsemi, tk, fcm, emin, emax, rmax, gmax, alat, r_log, rcutz, rcutxy, eshift, qbound, hfield, mixing, abasis, &
      bbasis, cbasis, vconst, tksemi, tolrdif, emusemi, ebotsemi, fsemicore, lambda_xc, deltae, lrhosym, linipol, lcartesian, imt, cls, lmxc, irns, irws, ntcell, refpot, inipol, &
      ixipol, hostimp, kfg, vbc, zperleft, zperight, bravais, rmt, zat, rws, mtfac, rmtref, rmtnew, rmtrefat, fpradius, tleft, tright, rbasis, socscale, cscl, socscl, solver, i12, &
      i13, i19, i25, i40, txc, drotq, ncpa, itcpamax, cpatol, noq, iqat, icpa, kaoez, conc, kmrot, qmtet, qmphi, kreadldau, lopt, ueff, jeff, erefldau)

    ! Some consistency checks
    if ((krel<0) .or. (krel>1)) stop ' set KREL=0/1 (non/fully) relativistic mode in the inputcard'
    if ((krel==1) .and. (nspind==2)) stop ' set NSPIND = 1 for KREL = 1 in the inputcard'

    ! Set the calculation of several parameters
    nrefd = naez                   ! +nemb ! can be changed later on when it is determined in clsgen_tb

    irm = irmd

    ntotd = ipand + 30
    ncelld = naez
    nrmaxd = ntotd*(ncheb+1)
    nchebd = ncheb

    natypd = natyp
    naezd = naez

    lmaxd = lmax
    alm = naezd*lmmaxd
    almgf0 = naezd*lmgf0d
    ndim_slabinv = nprincd*lmmaxd
    nembd = nembd1 - 1
    nembd2 = nembd + naez
    irmind = irmd - irnsd
    nofgij = natomimpd*natomimpd + 1
    lpotd = lpot
    lmpotd = (lpot+1)**2
    lmmaxso = lmmaxd               ! lmmaxd already doubled in size! (KREL+1)*LMMAXD

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocation calls
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> @note Jonathan Chico: The main idea here is to allocate all the needed arrays so that the inc.p
    !> file becomes irrelevant. In principle the philosophy would be to modularize
    !> the code such that each module has its own global variables and allocation routine
    !> e.g. a module called CPA_control could have defined all the needed CPA variables
    !> as well as the allocation calls, this module would be used in the needed routines
    !> and the arrays would only be allocated if a CPA calculation is actually performed
    !> in the current way ALL arrays are allocated which could cause an unnecessary memory
    !> consumption

    ! Call to allocate the arrays associated with the potential
    call allocate_potential(1, irmd, natypd, npotd, ipand, nfund, lmxspd, lmpotd, irmind, nspotd, nfu, irc, ncore, irmin, lmsp, lmsp1, ircut, lcore, llmsp, ititle, fpradius, visp, &
      ecore, vins)
    ! Call to allocate the arrays associated with the LDA+U potential
    call allocate_ldau_potential(1, irmd, natypd, mmaxd, nspind, itldau, wldau, uldau, phildau)
    ! Call to allocate the arrays associated with the energy
    call allocate_energies(1, iemxd, ez, dez, wez)
    ! Call to allocate the arrays associated with the relativistic corrections
    call allocate_relativistic(1, krel, irmd, naezd, natypd, zrel, jwsrel, irshift, vtrel, btrel, rmrel, drdirel, r2drdirel, qmgam, qmgamtab, qmphitab, qmtettab)
    ! Call to allocate the arrays associated with the relativistic transformations
    call allocate_rel_transformations(1, lmmaxd, nrrel, irrel, rc, crel, rrel, srrel)
    ! Call to allocate the arrays associated with the clusters
    call allocate_clusters(1, naezd, lmaxd, ncleb, nclsd, nembd1, nsheld, naclsd, lmpotd, natomimpd, nsh1, nsh2, nacls, nshell, atomimp, atom, ezoa, icleb, jend, ratom, rclsimp, &
      cmomhost, rcls)
    ! Call to allocate the arrays associated with the expansion of the Green function
    call allocate_expansion(1, lm2d, irid, nfund, ntotd, ncleb, lassld, ncelld, nchebd, loflm, wg, cleb, yrg, thetas, thetasnew)
    ! Call to allocate the arrays associated with the integration mesh
    call allocate_mesh(1, irmd, natypd, a, b, rmesh, drdi)
    ! Call to allocate the arrays associated with the pannels for the new solver
    call allocate_pannels(1, natypd, ntotd, ipan, npan_tot, npan_eq_at, npan_log_at, ipan_intervall, rpan_intervall)
    ! Call to allocate misc arrays
    call allocate_misc(1, nrd, irmd, irid, lmaxd, naezd, natypd, nfund, nrefd, iemxd, ntotd, nsheld, lmmaxd, nembd1, nchebd, ncelld, lmxspd, nspindd, nsymaxd, nprincd, ifunm, &
      ifunm1, icheck, vref, s, rr, dror, rnew, rs, rrot, thesme, dsymll, dsymll1, lefttinvll, righttinvll)
    ! Call to allocate the arrays associated with the Green function
    call allocate_green(1, naezd, iemxd, ngshd, nsheld, lmpotd, nofgij, ish, jsh, kmesh, imaxsh, iqcalc, iofgij, jofgij, ijtabsh, ijtabsym, ijtabcalc, ijtabcalc_i, ilm_map, gsh)

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End of allocation calls
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Deal with the lattice
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call lattix99(linterface, alat, natyp, naez, conc, rws, bravais, recbv, volume0, rr, nrd, natyp)

    ! rr has changed, fix allocation of array to new nrd size
    allocate (tmp_rr(3,0:nrd), stat=ierr)
    if (ierr/=0) stop 'error allocating tmp_rr in main0'
    tmp_rr(:, :) = rr(1:3, 0:nrd)
    deallocate (rr, stat=ierr)
    if (ierr/=0) stop 'error deallocating rr in main0'
    allocate (rr(3,0:nrd), stat=ierr)
    if (ierr/=0) stop 'error reallocating rr in main0'
    rr(:, :) = tmp_rr(:, :)
    deallocate (tmp_rr, stat=ierr)
    if (ierr/=0) stop 'error allocating tmp_rr in main0'

    call scalevec(lcartesian, rbasis, abasis, bbasis, cbasis, nlbasis, nrbasis, nleft, nright, zperleft, zperight, tleft, tright, linterface, naez, nemb, bravais, kaoez, noq, naez, &
      natyp, nemb)
    ! After SCALEVEC all basis positions are in cartesian coords.

    nvirt = 0
    if (opt('VIRATOMS')) then
      write (1337, *) 'Calling ADDVIRATOMS'
      call addviratoms14(linterface, nvirt, naez, naez, natyp, nemb, nemb, rbasis, .true., bravais, ncls, nineq, refpot, kaoez, noq, nref, rmtrefat, i25)
    end if

    call clsgen_tb(naez, nemb, nvirt, rr, rbasis, kaoez, zat, cls, ncls, nacls, atom, ezoa, nlbasis, nrbasis, nleft, nright, zperleft, zperight, tleft, tright, rmtref, rmtrefat, &
      vref, refpot, nref, rcls, rcutz, rcutxy, alat, natyp, nclsd, nrd, naclsd, nrefd, nembd, linterface, nprinc)
    ! overwrite nrefd and change allocations accordingly
    allocate (rmtref_temp(nrefd), vref_temp(nrefd))
    rmtref_temp(:) = rmtref(:)
    vref_temp(:) = vref(:)
    nrefd = nref
    deallocate (rmtref, vref)
    allocate (rmtref(nrefd), vref(nrefd), stat=i_stat)
    rmtref(:) = rmtref_temp(1:nref)
    vref(:) = vref_temp(1:nref)
    deallocate (rmtref_temp, vref_temp)

    ! overwrite nprincd if chosen too small
    if (nprincd<nprinc) then
      nlayer = naez/nprinc
      if (nlayer*nprinc/=naez) nprinc = naez
      write (*, *) 'Automatically overwriting nprincd with ', nprinc
      write (1337, *) 'Automatically overwriting nprincd with ', nprinc
      nprincd = nprinc
      ! update parameter that depend on nprincd
      ndim_slabinv = nprincd*lmmaxd
      ! change allocations of arrays that have nprincd
      deallocate (icheck, stat=i_stat)
      call memocc(i_stat, product(shape(icheck))*kind(icheck), 'ICHECK', 'main0')
      allocate (icheck(naez/nprincd,naez/nprincd), stat=i_stat)
      call memocc(i_stat, product(shape(icheck))*kind(icheck), 'ICHECK', 'main0')
      icheck = 0
      ! do not forget to update nlayerd as well!
      nlayerd = nlayer
    end if


    ! Now the clusters, reference potentials and muffin-tin radii have been set.
    ! -------------------------------------------------------------------------
    ! Consistency check
    ! -------------------------------------------------------------------------
    if ((krel==1) .and. (ins/=0)) then
      write (6, *) ' FULL-POTENTIAL RELATIVISTIC mode not implemented '
      stop ' set INS = 0 in the input'
    end if

    if (kvrel<=1) then
      if (krel==1) stop ' KVREL <= 1 in input, but relativistic program used'
    else
      if (krel==0) stop ' KVREL > 1 in input, but non-relativistic program used'
    end if
    ! -------------------------------------------------------------------------

    e2in = emax
    nsra = 1
    if (kvrel>=1) nsra = 2
    if (kvrel>=2) nsra = 3

    call testdim(nspin, naez, nemb, natyp, ins, insref, nref, irns, nlayer, krel, nspind, nprincd, knosph, irnsd, korbit)

    if (ins>0) open (19, file=i19, status='old', form='formatted')
    if (ifile==13) open (ifile, file=i13, status='old', form='formatted')
    if (icc>0) open (25, file=i25, status='unknown', form='formatted')

    call startb1(ifile, 1337, 1337, ipe, krel, kws, lmax, 1, natyp, alatnew, rmtnew, rmt, ititle, imt, irc, vconst, ins, irns, fpradius, nspin, vins, irmin, kshape, ntcell, ircut, &
      ipan, thetas, ifunm, nfu, llmsp, lmsp, e2in, vbc, dror, rs, s, visp, rws, ecore, lcore, ncore, drdi, rmesh, zat, a, b, irws, 1, lmpot, irmind, irmd, lmxspd, ipand, irid, &
      irnsd, natyp, ncelld, nfund, nspotd, ivshift, npotd)


    ! find md5sums for potential and shapefunction
    call get_md5sums(ins, i13, i19)
    write (1337, '(A,A)') 'Doing calculation with potential: ', md5sum_potential
    if (ins>0) then
      write (1337, '(A,A)') 'Doing calculation with shapefun: ', md5sum_shapefun
    end if

    if (test('rhoqtest')) then
      call rhoq_save_rmesh(natyp, irm, ipand, irmin, irws, ipan, rmesh, ntcell, ircut, r_log, npan_log, npan_eq)

      ! check consistency
      if (opt('GREENIMP') .or. opt('WRTGREEN')) then
        write (*, *) 'warning! rhoqtest cannot be used together with '
        write (*, *) '''GREENIMP'' or ''WRTGREEN'' options'
        stop
      end if

      ! ! enforce MPIatom here
      ! if (TEST('MPIenerg')) stop 'rhoqtest assumes MPIatom'
      ! CALL ADDTEST('MPIatom')

      if (nranks>1) then
        write (*, *) 'at the moment rhoqtest does not work with MPI.'
        write (*, *) 'compile hybrid version and use OMP level only.'
        stop
      end if

    end if                         ! TEST('rhoqtest')

    if (test('Vspher  ')) then
      write (1337, *) 'TEST OPTION Vspher,', 'keeping only spherical component of potential.'
      vins(irmind:irmd, 2:lmpot, 1:nspotd) = 0.d0
    end if

    if (opt('zeropot ') .or. test('zeropot ')) then
      write (1337, *) 'Using OPT zeropot, setting potential to zero.'
      write (1337, *) 'Using OPT zeropot, setting nuclear charge to zero.'
      vins(irmind:irmd, 1:lmpot, 1:nspotd) = 0.d0
      visp(1:irmd, 1:npotd) = 0.d0
      zat(1:natyp) = 0.d0
    end if

    do i1 = 1, natyp
      do lm = 1, lmxspd
        ifunm1(lm, i1) = ifunm(i1, lm)
        lmsp1(lm, i1) = lmsp(i1, lm)
      end do
    end do
    ! -------------------------------------------------------------------------
    ! update Fermi energy, adjust energy window according to running options
    ! -------------------------------------------------------------------------
    if (npol==0) efermi = e2in
    if (opt('GF-EF   ') .or. opt('DOS-EF  ')) then
      emin = e2in
      if (opt('GF-EF   ')) then
        write (1337, fmt=130)
      else
        write (1337, fmt=140)
      end if
    end if

    if (abs(e2in-emax)>1d-10 .and. npol/=0) emax = e2in
    ! -------------------------------------------------------------------------
    if (opt('GENPOT  ')) then
      rewind (3)
      call generalpot(3, 1, natyp, nspin, zat, alat, rmt, rmtnew, rws, rmesh, drdi, visp, irws, a, b, ins, irns, lpot, vins, qbound, irc, kshape, e2in, vbc, ecore, lcore, ncore, &
        lmpot, irmd, irmind)
      close (3)
    end if
    ! -------------------------------------------------------------------------
    ! --> Apply external magnetic field
    ! -------------------------------------------------------------------------
    ! from startb1 moved here
    if (khfeld==1) then
      ! ---> maybe apply a magnetic field
      call bshift_ns(irmd, irid, ipand, lmpot, npotd, natyp, nspin, ngshd, nfund, ncelld, irmind, lmxspd, kshape, irc, irmin, inipol, ntcell, imaxsh, ilm_map, lmsp, ifunm, ircut, &
        hfield, gsh, rmesh, thesme, thetas, visp, vins)
    end if
    if (test('vpotout ')) then     ! ruess
      open (unit=54633163, file='test_vpotout_bshift')
      do i1 = 1, natyp*nspin
        write (54633163, *) '# visp of atom ', i1
        write (54633163, '(50000E14.7)') visp(:, i1)
      end do                       ! iatom
      do i1 = 1, natyp*nspin
        write (54633163, *) '# vins of atom ', i1
        write (54633163, '(50000E14.7)') vins(:, :, i1)
      end do                       ! iatom
      close (54633163)
    end if
    ! -------------------------------------------------------------------------
    ! Deal with the potential in the RELATIVISTIC CASE
    ! -------------------------------------------------------------------------
    para = .true.
    if (krel+korbit==1) then
      ! ----------------------------------------------------------------------
      if (nspin==1) then
        ! -------------------------------------------------------------------
        ! for paramagnetic (NSPIN=1) input potential fill up also the
        ! V(DOWN), ECORE(DOWN), LCORE(DOWN), NCORE(DOWN) and ITITLE(DOWN)
        ! arrays (needed)
        ! -------------------------------------------------------------------
        do i = natyp, 1, -1
          j = 2*i - 1
          call dcopy(irmd, visp(1,i), 1, visp(1,j), 1)
          call dcopy(irmd, visp(1,j), 1, visp(1,j+1), 1)

          call dcopy(20, ecore(1,i), 1, ecore(1,j), 1)
          call dcopy(20, ecore(1,j), 1, ecore(1,j+1), 1)

          ncore(j) = ncore(i)
          ncore(j+1) = ncore(j)

          do i1 = 1, 20
            lcore(i1, j) = lcore(i1, i)
            lcore(i1, j+1) = lcore(i1, j)
            ititle(i1, j) = ititle(i1, i)
            ititle(i1, j+1) = ititle(i1, j)
          end do
        end do
        ! -------------------------------------------------------------------
      else                         ! NSPIN.eq.1
        ! -------------------------------------------------------------------
        ! --> check whether, although NSPIN=2 at input, the system is
        ! paramagnetic (useful for symmetry cosiderations)
        ! -------------------------------------------------------------------
        do i = 1, 2*natyp - 1, 2
          do j = 1, irmd
            if (abs(visp(j,i)-visp(j,i+1))>1d-5) para = .false.
          end do
        end do
        if (para) then
          do i = 1, 2*natyp - 1, 2
            call dcopy(irmd, visp(1,i), 1, visp(1,i+1), 1)
          end do
        end if

      end if                       ! NSPIN.eq.1
      ! ----------------------------------------------------------------------
      ! finally, convert input potential to the internal relativistic
      ! form VTREL,BTREL. Set up auxiliary arrays (needed in the REL
      ! routines) ZREL, JWSREL, RMREL, DRDIREL, R2DRDIREL, IRSHIFT
      ! ----------------------------------------------------------------------
      if (krel==1) then
        ! call this only if relativisitic solver is used
        call relpotcvt(1, visp, zat, rmesh, drdi, ircut, vtrel, btrel, zrel, rmrel, jwsrel, drdirel, r2drdirel, irshift, ipand, irmd, npotd, natyp)
      end if
    end if                         ! KREL+KORBIT.EQ.1
    ! -------------------------------------------------------------------------
    ! set up energy contour
    ! -------------------------------------------------------------------------
    idosemicore = 0
    if (opt('SEMICORE')) idosemicore = 1

    call epathtb(ez, dez, e2in, ielast, iesemicore, idosemicore, emin, emax, tk, npol, npnt1, npnt2, npnt3, ebotsemi, emusemi, tksemi, npolsemi, n1semi, n2semi, n3semi, iemxd)

    ! SUSC (BEGIN: modifications by Manuel and Benedikt)             ! susc
    ! ! susc
    if (opt('EMESH   ')) then      ! susc
      ! write out the energy mesh and the corresponding                ! susc
      ! weights to a file called 'emesh.scf'                           ! susc
      write (*, '("main0: Runflag emesh is set.")') ! susc
      write (*, '("       File emesh.scf will be written!")') ! susc
      write (*, *) 'writing emesh.scf file...' ! susc
      open (file='emesh.scf', unit=12111984, status='replace') ! susc
      write (12111984, '(5x,i0)') ielast ! susc
      do ie = 1, ielast            ! susc
        write (12111984, '(4es16.8)') ez(ie), dez(ie) ! susc
      end do                       ! susc
      close (12111984)             ! susc
      write (*, '("       Finished writing file emesh.scf.")') ! susc
    end if                         ! susc
    ! ! susc
    ! ! susc
    if (opt('KKRSUSC ')) then      ! susc
      ! read in 'emesh.dat' with new energy mesh-points                ! susc
      inquire (file='emesh.dat', exist=emeshfile) ! susc
      write (*, '("main0: Runflag KKRSUSC is set.")') ! susc
      if (emeshfile) then          ! susc
        write (*, '("main0: File emesh.dat exists and will ")', advance='no') ! susc
        write (*, '("be read in.")') ! susc
        write (*, '("       Energy contour from inputcard ")', advance='no') ! susc
        write (*, '("will be overwritten!")') ! susc
        open (file='emesh.dat', unit=50) ! susc
        read (50, *) ielast        ! susc
        if (ielast>iemxd) stop 'ielast > iemxd!' ! susc
        do ie = 1, ielast          ! susc
          read (50, '(4es16.8)') ez(ie), dez(ie) ! susc
          write (*, '(i8,4es16.8)') ie, ez(ie), dez(ie) ! susc
        end do                     ! susc
        close (50)                 ! susc
        write (*, '("       Finished reading in file emesh.dat.")') ! susc
      else                         ! susc
        stop 'main0: Runflag KKRSUSC but cannot find file emesh.dat!' ! susc
      end if                       ! susc
    end if                         ! susc
    ! ! susc
    ! ! susc
    ! still missing: check here whether scfsteps is > 1              ! susc
    ! if scfsteps>1 --> option a) stop program here                ! susc
    ! option b) set it to 1 and continue         ! susc
    if (opt('KKRSUSC ') .and. scfsteps>1) then ! susc
      write (*, '("main0: Runflag KKRSUSC is set ")') ! susc
      write (*, '("but scfsteps = ",i0)') scfsteps ! susc
      write (*, '("       Here we enforce scfsteps = 1")') ! susc
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! susc
      scfsteps = 1                 ! susc
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! susc
    end if                         ! susc
    ! ! susc
    ! SUSC (END:   modifications by Manuel and Benedikt)             ! susc

    do ie = 1, ielast
      wez(ie) = -2.0_dp/pi*dez(ie)
      if (ie<=iesemicore) wez(ie) = wez(ie)*fsemicore
    end do
    ! -------------------------------------------------------------------------
    ! update energy contour for Fermi-surface generation                       ! fswrt=fermi-surface write
    ! -------------------------------------------------------------------------
    if (opt('FERMIOUT')) then      ! fswrt
      if (aimag(ez(1))>0d0) stop 'E has imaginary part' ! fswrt
      ielast = 3                   ! fswrt
      ez(2) = ez(1) + cmplx(1.0d-03, 0.0d0, kind=dp) ! fswrt
      ez(3) = ez(1) - cmplx(1.0d-03, 0.0d0, kind=dp) ! fswrt
    end if                         ! fswrt
    ! -------------------------------------------------------------------------
    ! update the value of NSPIN to be consistent with REL mode
    ! -------------------------------------------------------------------------
    if (krel==1) nspin = 1

    call gaunt2(wg, yrg, 4*lmax)
    call gaunt(lmax, lpot, wg, yrg, cleb, loflm, icleb, iend, jend, ncleb, lmax, lmgf0d, lmpot)

    ! -------------------------------------------------------------------------
    ! set up of GAUNT coefficients C(l,m;l',m';l'',m'') for all
    ! nonvanishing (l'',m'')-components of the shape functions THETAS
    ! -------------------------------------------------------------------------
    if (kshape/=0) then
      call shape_corr(lpot, natyp, gsh, ilm_map, imaxsh, lmsp, ntcell, wg, yrg, lassld, lmpot, natyp, ngshd)
    end if
    ! -------------------------------------------------------------------------
    ! calculate Madelung constants (needed only for SCF calculations)
    ! -------------------------------------------------------------------------
    ! fivos      IF ( SCFSTEPS.GT.1 .OR. ICC.GT.0 ) THEN
    if (npol/=0 .or. opt('DECIMATE')) then ! No madelung calculation in case of DOS., needed for demination nevertheless
      ! OPEN(99,FILE='madelinfo.txt')

      !> @note Use option 'ewald2d' if the madelung summation is to be carried out in
      !> single-slab mode, otherwise it is carried out in repeated (periodic)
      !> slab mode.
      !> Reason: the 2d-mode gives wrong results sometimes [e.g. in diamond
      !> structure (110)].
      if (linterface .and. (opt('ewald2d ') .or. opt('DECIMATE'))) then ! ewald2d
        write (*, *) 'Calling MADELUNG2D'
        ! -------------------------------------------------------------------
        ! 2D case
        ! -------------------------------------------------------------------
        call madelung2d(lpot, yrg, wg, naez, alat, volume0, bravais, recbv, rbasis, rmax, gmax, nlbasis, nleft, zperleft, tleft, nrbasis, nright, zperight, tright, lmxspd, lassld, &
          lpot, lmpot, nmaxd, ishld, nembd1, wlength)
        write (*, *) 'Exited MADELUNG2D'
      else
        ! -------------------------------------------------------------------
        ! 3D case
        ! -------------------------------------------------------------------
        if (linterface) then
          call getbr3(nembd1, nlbasis, alat, tleft, nrbasis, tright, bravais, recbv, volume0)
        end if

        write (*, *) 'Calling MADELUNG3D'
        call madelung3d(lpot, yrg, wg, naez, alat, volume0, bravais, recbv, rbasis, rmax, gmax, naez, lmxspd, lassld, lpot, lmpot, nmaxd, ishld, nemb, wlength)
        write (*, *) 'Exited MADELUNG3D'
      end if

      ! CLOSE(99)
    else                           ! NPOL==0
      ! write dummy files

      ! real (kind=dp) AVMAD(LMPOT,LMPOT),BVMAD(LMPOT)
      lrecabmad = wlength*2*lmpot*lmpot + wlength*2*lmpot
      open (69, access='direct', recl=lrecabmad, file='abvmad.unformatted', form='unformatted')
      do i = 1, naez
        do j = 1, naez
          irec = j + naez*(i-1)
          write (69, rec=irec) 0.0d0, 0.0d0 ! AVMAD,BVMAD
        end do
      end do
      close (69)

    end if                         ! npol==0
    ! -------------------------------------------------------------------------
    ! fivos      END IF
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! Set up I,J pairs for ICC = -1
    ! -------------------------------------------------------------------------
    if (icc<0) then
      call setgijtab(linterface, icc, naez, iqat, rbasis, bravais, natomimp, atomimp, rclsimp, ijtabcalc, iofgij, jofgij, nqcalc, iqcalc, natomimpd, ijtabcalc_i)
    end if

    ! -------------------------------------------------------------------------

    dsymll = (0d0, 0d0)
    dsymll1 = (0d0, 0d0)

    call bzkint0(nshell, naez, natyp, noq, rbasis, kaoez, icc, bravais, recbv, atomimp, rsymat, isymindex, nsymat, i25, natomimp, nsh1, nsh2, rclsimp, ratom, ijtabsym, ijtabsh, &
      ijtabcalc, iofgij, jofgij, nofgij, ish, jsh, rrot, dsymll1, para, qmtet, qmphi, symunitary, hostimp, intervx, intervy, intervz, ielast, ez, kmesh, maxmesh, maxmshd, nsymaxd, &
      krel+korbit, lmax, lmmaxd, kpoibz, naez, natyp, natomimpd, nsheld, nemb)

    ! -------------------------------------------------------------------------

    if (opt('KKRFLEX ')) then

      call writehoststructure(bravais, natyp, rbasis, naez, nemb)

      open (58, file='kkrflex_atominfo', form='FORMATTED')
      call version_print_header(58, '; '//md5sum_potential//'; '//md5sum_shapefun)
      nvatom = 0
      do i = 1, natomimp
        if (kaoez(1,atomimp(i))==-1) nvatom = nvatom + 1
      end do
      write (58, '(500A)') '#NATOM   NTOTATOM'
      write (58, *) natomimp, natomimp - nvatom
      write (58, '(500A)') '#Impurity positions x,y,z|Core Charge|Virtual Atom?|Remove Atom?|LMAX'
      do i = 1, natomimp
        if (kaoez(1,atomimp(i))==-1) then
          zattemp = 0.d0
          isvatom = 1
          nvatom = nvatom + 1
        else
          isvatom = 0
          zattemp = zat(kaoez(1,atomimp(i)))
        end if
        write (58, '(3F25.16,F6.2,3I5)')(rclsimp(j,i), j=1, 3), zattemp, isvatom, 0, lmax
      end do
      close (58)
    end if

    ! -------------------------------------------------------------------------
    ! fivos: write out nshell and nsh1,nsh2 into standard output and in file shells.dat
    ! -------------------------------------------------------------------------
    if (icc/=0 .and. .not. opt('KKRFLEX ')) then
      open (58, file='shells.dat')
      write (1337, *) 'Writing out shells (also in shells.dat):' ! fivos
      write (1337, *) 'itype,jtype,iat,jat,r(iat),r(jat)' ! fivos
      write (1337, *) nshell(0), 'NSHELL(0)' ! fivos
      write (58, *) nshell(0), 'NSHELL(0)' ! fivos
      do i1 = 1, nshell(0)         ! fivos
        write (1337, *) i1, nshell(i1), 'No. of shell, No. of atoms in shell' ! fivos
        write (58, *) i1, nshell(i1), 'No. of shell, No. of atoms in shell' ! fivos
        do lm = 1, nshell(i1)      ! fivos
          write (1337, *) 'ish(i1,lm)', ish(i1, lm)
          if (ish(i1,lm)>0 .and. jsh(i1,lm)>0) then ! fix bernd
            write (1337, 100) nsh1(i1), nsh2(i1), ish(i1, lm), jsh(i1, lm) & ! fivos
              , (rclsimp(i,ish(i1,lm)), i=1, 3), (rclsimp(i,jsh(i1,lm)), i=1, 3) ! fivos
            write (58, 100) nsh1(i1), nsh2(i1), ish(i1, lm), jsh(i1, lm), & ! fivos
              (rclsimp(i,ish(i1,lm)), i=1, 3), (rclsimp(i,jsh(i1,lm)), i=1, 3) ! fivos
          else                     ! fix bernd
            write (1337, 110) nsh1(i1), nsh2(i1), ish(i1, lm), jsh(i1, lm) ! fix bernd
            write (58, 110) nsh1(i1), nsh2(i1), ish(i1, lm), jsh(i1, lm) ! fix bernd
          end if                   ! fix bernd
100       format (4i5, 6f16.6)     ! fivos
110       format (4i5)             ! fix bernd
        end do                     ! fivos
      end do                       ! fivos
      write (1337, *) '###################'
      close (58)
    end if
    ! -------------------------------------------------------------------------
    ! end fivos
    ! -------------------------------------------------------------------------

    call gfmask(linterface, icheck, icc, invmod, nsh1, nsh2, naez, nshell(0), naez, nprincd)
    ! -------------------------------------------------------------------------
    ! set up transformation matrices between REL/NREL representations
    ! -------------------------------------------------------------------------
    if ((krel+korbit)==1) then
      call drvbastrans(rc, crel, rrel, srrel, nrrel, irrel, lmax+1, lmmaxd, 2*(lmax+1), lmmaxd+2*(lmax+1), mmaxd, 2*(lmax+1)*mmaxd)
    end if
    if (opt('NEWSOSOL')) then
      do ns = 1, nsymat
        call changerep(dsymll1(1,1,ns), 'REL>RLM', dsymll(1,1,ns), lmmaxd, lmmaxd, rc, crel, rrel, 'DSYMLL', 0)
      end do
      ! DSYMLL(:,:,:)=DSYMLL1(:,:,:)
    else
      dsymll(:, :, :) = dsymll1(:, :, :)
    end if

    ! -------------------------------------------------------------------------
    ! for the case that the magnetisation is rotated with respect to
    ! the (001)-direction (KMROT<>0) calculate the rotation matrices
    ! to switch between the CRYSTAL and LOCAL frames of reference
    ! -------------------------------------------------------------------------
    call cinit(lmmaxd*lmmaxd*naez, drotq)

    if (kmrot/=0) then
      fact(0) = 1.0d0
      do i = 1, 100
        fact(i) = fact(i-1)*dble(i)
      end do

      do i1 = 1, naez
        call calcrotmat(mmaxd, (krel+korbit)*3, qmphi(i1), qmtet(i1), 0.0_dp, drotq(1,1,i1), fact, lmmaxd)
      end do
    end if
    ! -------------------------------------------------------------------------
    ! Treat decimation I/O cases
    ! -------------------------------------------------------------------------
    if (opt('deci-pot')) then
      call outpothost(alat, ins, krel+korbit, kmrot, nspin, naez, natyp, e2in, bravais, rbasis, qmtet, qmphi, noq, kaoez, iqat, zat, conc, ipan, ircut, solver, socscl, cscl, irws, &
        rmtnew, rws, rmesh, drdi, visp, irshift, rmrel, drdirel, vtrel, btrel, lmax, natyp, naez, ipand, irmd)
    end if
    if (opt('deci-out')) then
      if (nranks>1) stop 'ERROR: deci-out does not work with MPI!'
      call outtmathost(alat, ins, krel+korbit, kmrot, nspin, naez, lmmax, bravais, rbasis, qmtet, qmphi, e2in, tk, npol, npnt1, npnt2, npnt3)
    end if
    if (opt('DECIMATE')) then
      call deciopt(alat, ins, krel+korbit, kvrel, kmrot, nspin, naez, lmmax, bravais, tk, npol, npnt1, npnt2, npnt3, ez, ielast, kaoez, lefttinvll, righttinvll, vacflag, nlbasis, &
        nrbasis, cmomhost, vref, rmtref, nref, refpot(naez), lmax, lmgf0d, lmmaxd, lm2d, nembd1, iemxd, nspindd, lmpot, natyp, irmd, ipand)
    end if
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! ITERMDIR  -- initialise
    ! -------------------------------------------------------------------------
    if (opt('ITERMDIR')) then
      write (1337, *)
      write (1337, *) 'Angle mixing scheme will be applied '
      write (1337, *)
      do i1 = 1, naez
        qmphitab(i1, 1) = qmphi(i1)
        qmtettab(i1, 1) = qmtet(i1)
        qmgamtab(i1, 1) = qmgam(i1)
        do i = 2, 3
          qmphitab(i1, i) = 0d0
          qmtettab(i1, i) = 0d0
          qmgamtab(i1, i) = 0d0
        end do
      end do
    end if

    ! -------------------------------------------------------------------------
    ! LDA+U -- initialise
    ! -------------------------------------------------------------------------

    if (opt('LDA+U   ')) then
      call startldau(itrunldau, idoldau, kreadldau, lopt, ueff, jeff, erefldau, natyp, nspin, wldau, uldau, phildau, irws, ntldau, itldau, irmd, natyp, nspind, mmaxd)
    end if
    ! -------------------------------------------------------------------------
    ! LDA+U
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! write out information for the other program parts           =
    ! -------------------------------------------------------------------------

    ! new solver for full-potential, spin-orbit, initialise
    if (opt('NEWSOSOL')) then
      call create_newmesh(natyp, irmd, ipand, irid, ntotd, nfund, ncheb, ntotd*(ncheb+1), nspin, rmesh, irmin, ipan, ircut, r_log, npan_log, npan_eq, npan_log_at, npan_eq_at, &
        npan_tot, rnew, rpan_intervall, ipan_intervall, ncelld, ntcell, thetas, thetasnew)
    end if

    call wunfiles(npol, npnt1, npnt2, npnt3, ielast, tk, emin, emax, ez, wez, efermi, npolsemi, n1semi, n2semi, n3semi, iesemicore, tksemi, ebotsemi, emusemi, fsemicore, vins, &
      visp, vbc, vtrel, btrel, rmrel, drdirel, r2drdirel, zrel, jwsrel, irshift, itscf, scfsteps, cmomhost, ecore, lcore, ncore, qmtet, qmphi, qmphitab, qmtettab, qmgamtab, drotq, &
      nsra, ins, natypd, naezd, nineq, nref, nspin, ncls, icst, ipan, ircut, alat, zat, rmesh, drdi, refpot, rmtref, vref, iend, jend, cleb, icleb, atom, cls, rcls, nacls, loflm, &
      solver, socscl, cscl, icc, igf, nlbasis, nrbasis, ncpa, icpa, itcpamax, cpatol, rbasis, rr, ezoa, nshell, nsh1, nsh2, ijtabcalc, ijtabcalc_i, ish, jsh, ijtabsym, ijtabsh, &
      nofgij, nqcalc, iqcalc, kmrot, kaoez, iqat, noq, conc, kmesh, maxmesh, nsymat, symunitary, rrot, dsymll, invmod, icheck, natomimp, ratom, atomimp, rc, crel, rrel, srrel, &
      nrrel, irrel, lefttinvll, righttinvll, vacflag, a, b, ifunm, ifunm1, intervx, intervy, intervz, ititle, lmsp1, ntcell, thetas, lpotd, lmpotd, nright, nleft, linterface, imix, &
      mixing, qbound, fcm, itdbry, irns, kpre, kshape, kte, kvmad, kxc, lambda_xc, txc, ishift, ixipol, lrhosym, kforce, lmsp, llmsp, rmt, rmtnew, rws, imt, irc, irmin, irws, nfu, &
      hostimp, gsh, ilm_map, imaxsh, idoldau, itrunldau, ntldau, lopt, itldau, ueff, jeff, erefldau, uldau, wldau, phildau, iemxd, irmind, irmd, nspotd, npotd, nembd1, lmmaxd, &
      ipand, nembd2, lmax, ncleb, naclsd, nclsd, lm2d, lmax+1, mmaxd, nrd, nsheld, nsymaxd, naez/nprincd, natomimpd, nspind, irid, nfund, ncelld, lmxspd, ngshd, krel, ntotd, ncheb, &
      npan_log, npan_eq, npan_log_at, npan_eq_at, r_log, npan_tot, rnew, rpan_intervall, ipan_intervall, nspindd, thetasnew, socscale, tolrdif, lly, deltae, rclsimp)

    if (opt('FERMIOUT')) then      ! fswrt
      call write_tbkkr_files(lmax, nemb, ncls, natyp, naez, ielast, ins, alat, & ! fswrt
        bravais, recbv, rbasis, cls, nacls, rcls, ezoa, atom, rr, nspin, nrd, korbit, & ! fswrt
        nclsd, naclsd)             ! fswrt
    end if                         ! fswrt

    if (opt('OPERATOR')) then
      ! check if impurity files are present (otherwise no imp.
      ! wavefunctions can be calculated)
      operator_imp = .true.
      inquire (file='potential_imp', exist=lexist)
      if (.not. lexist) operator_imp = .false.
      inquire (file='shapefun_imp', exist=lexist)
      if (.not. lexist) operator_imp = .false.
      inquire (file='scoef', exist=lexist)
      if (.not. lexist) operator_imp = .false.
    else
      operator_imp = .false.
    end if

    if (opt('GREENIMP') .or. operator_imp) then ! GREENIMP
      ! fill array dimensions and allocate arrays in t_imp          ! GREENIMP
      call init_params_t_imp(t_imp, ipand, natyp, irmd, irid, nfund, & ! GREENIMP
        nspin, irmind, lmpot)      ! GREENIMP
      call init_t_imp(t_inc, t_imp) ! GREENIMP
      ! GREENIMP
      ! next read impurity potential and shapefunction              ! GREENIMP
      call readimppot(natomimp, ins, 1337, 0, 0, 2, nspin, lpot, & ! GREENIMP
        t_imp%ipanimp, t_imp%thetasimp, t_imp%ircutimp & ! GREENIMP
        , t_imp%irwsimp, khfeld, hfield, t_imp%vinsimp, & ! GREENIMP
        t_imp%vispimp, t_imp%irminimp, & ! GREENIMP
        t_imp%rimp, t_imp%zimp, irmd, irnsd, irid, nfund, ipand) ! GREENIMP
    end if                         ! GREENIMP


    if (ishift==2) then            ! fxf
      open (67, file='vmtzero', form='formatted') ! fxf
      write (67, 120) vbc(1)       ! fxf
      close (67)                   ! fxf
120   format (d20.12)              ! fxf
    end if                         ! fxf

    ! Check for inputcard consistency in case of qdos option
    if (opt('qdos    ')) then
      write (1337, *)
      write (1337, *) '     < QDOS > : consistency check '
      if ((npol/=0) .and. (npnt1==0) .and. (npnt3==0)) then
        stop 'For qdos calculation change enery contour to dos path'
      end if
      if (tk>50.d0) write (*, *) &
        'WARNING:  high energy smearing due to high value of TEMPR for energy contour integration could not be of advantage. Consider changeing ''TEMPR'' to lower value'
      if (tk>50.d0) write (1337, *) &
        'WARNING:  high energy smearing due to high value of TEMPR for energy contour integration could not be of advantage. Consider changeing ''TEMPR'' to lower value'
      write (1337, *) '       QDOS: consistecy check complete'
    end if

    ! -------------------------------------------------------------------------

    write (1337, '(79("="),/,31X,"< KKR0 finished >",/,79("="),/)')
130 format (5x, 'INFO:  Output of cluster Green function at E Fermi')
140 format (5x, 'INFO:  Determination of DOS at E Fermi')

  end subroutine main0


  !-------------------------------------------------------------------------------
  !> Summary: 
  !> Author: 
  !> Category: KKRhost, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> 
  !-------------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------
  ! SUBROUTINE: BSHIFT_NS
  !> @brief Adds a constant (=VSHIFT) to the potentials of atoms
  ! ----------------------------------------------------------------------------
  subroutine bshift_ns(irm, irid, ipand, lmpot, npotd, natyp, nspin, ngshd, nfund, ncelld, irmind, lmxspd, kshape, irc, irmin, inipol, ntcell, imaxsh, ilm_map, lmsp, ifunm, ircut, &
    hfield, gsh, rmesh, thesme, thetas, visp, vins)

    use :: mod_datatypes

    implicit none

    ! Adds a constant (=VSHIFT) to the potentials of atoms

    ! Parameters:
    ! Input
    integer, intent (in) :: irm
    integer, intent (in) :: irid
    integer, intent (in) :: ipand
    integer, intent (in) :: lmpot
    integer, intent (in) :: npotd
    integer, intent (in) :: natyp  !! Number of kinds of atoms in unit cell
    integer, intent (in) :: nspin  !! Counter for spin directions
    integer, intent (in) :: ngshd
    integer, intent (in) :: nfund
    integer, intent (in) :: ncelld
    integer, intent (in) :: irmind
    integer, intent (in) :: lmxspd
    integer, intent (in) :: kshape !! exact treatment of WS cell
    integer, dimension (natyp), intent (in) :: irc !! r point for potential cutting
    integer, dimension (natyp), intent (in) :: irmin !! max r for spherical treatment
    integer, dimension (natyp), intent (in) :: inipol !! initial spin polarisation
    integer, dimension (natyp), intent (in) :: ntcell !! index for WS cell
    integer, dimension (0:lmpot), intent (in) :: imaxsh
    integer, dimension (ngshd, 3), intent (in) :: ilm_map
    integer, dimension (natyp, lmxspd), intent (in) :: lmsp !! 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
    integer, dimension (natyp, lmxspd), intent (in) :: ifunm
    integer, dimension (0:ipand, natyp), intent (in) :: ircut !! r points of panel borders
    real (kind=dp), intent (in) :: hfield !! External magnetic field, for initial potential shift in spin polarised case
    real (kind=dp), dimension (ngshd), intent (in) :: gsh
    real (kind=dp), dimension (irm, natyp), intent (in) :: rmesh
    real (kind=dp), dimension (irid, nfund, ncelld), intent (in) :: thesme
    real (kind=dp), dimension (irid, nfund, ncelld), intent (in) :: thetas !! shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion

    ! Input/Output:
    real (kind=dp), dimension (irm, npotd), intent (inout) :: visp !! Spherical part of the potential
    real (kind=dp), dimension (irmind:irm, lmpot, nspotd), intent (inout) :: vins !! Non-spherical part of the potential

    ! Inside
    integer :: ispin, ih, ipot, ir, lm, imt1, irc1, irmin1
    real (kind=dp) :: rfpi, vshift
    real (kind=dp), dimension (irm) :: pshiftr
    real (kind=dp), dimension (irm, lmpot) :: pshiftlmr

    rfpi = sqrt(16.0d0*atan(1.0d0))

    do ih = 1, natyp

      imt1 = ircut(1, ih)
      irc1 = irc(ih)
      irmin1 = irmin(ih)

      do ispin = 1, nspin
        ! shift potential spin dependent
        vshift = -dble(2*ispin-3)*hfield*inipol(ih)

        write (1337, *) 'SHIFTING OF THE POTENTIALS OF ATOM', ih, 'spin', ispin, ' BY', vshift, 'RY.'
        ipot = nspin*(ih-1) + ispin

        call rinit(irm*lmpot, pshiftlmr)
        call rinit(irm, pshiftr)
        do ir = 1, irc1
          pshiftlmr(ir, 1) = vshift
        end do

        if (kshape==0) then        ! ASA
          do ir = 1, irc1
            visp(ir, ipot) = visp(ir, ipot) + pshiftlmr(ir, 1)
          end do
        else                       ! Full-potential

          call convol(imt1, irc1, ntcell(ih), imaxsh(lmpot), ilm_map, ifunm, lmpot, gsh, thetas, thesme, 0.0_dp, rfpi, rmesh(1,ih), pshiftlmr, pshiftr, lmsp)

          do ir = 1, irc1
            visp(ir, ipot) = visp(ir, ipot) + pshiftlmr(ir, 1)
          end do

          do lm = 2, lmpot
            do ir = irmin1, irc1
              vins(ir, lm, ipot) = vins(ir, lm, ipot) + pshiftlmr(ir, lm)*rfpi
            end do
          end do
        end if                     ! (kshape.eq.0)
      end do
    end do

  end subroutine bshift_ns

  !-------------------------------------------------------------------------------
  !> Summary: 
  !> Author: 
  !> Category: KKRhost, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> 
  !-------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------
  ! SUBROUTINE: init_misc_variables
  !> @brief Set default values for misc variables for the calculation
  !> @details The idea behind this kind of routine is to separate the initialization
  !> of variables such that one can use this to modularize the code
  !> @author Jonathan Chico
  !> @date 22.12.2017
  ! -------------------------------------------------------------------------
  subroutine init_misc_variables()

    implicit none

    ipe = 0
    ipf = 0
    ipfe = 0
    kcor = 0
    kefg = 0
    khyp = 0
    kpre = 0
    nsra = 1
    iend = 1
    nvirt = 0
    kvmad = 0
    itscf = 0
    nsheld = 301                   ! Number of blocks of the GF matrix that need to be calculated (NATYPD + off-diagonals in case of impurity)
    invmod = 2                     ! Corner band matrix inversion scheme
    ielast = 0
    ishift = 0
    kfrozn = 0
    nsymat = 0
    nqcalc = 0
    kforce = 0                     ! Calculation of the forces
    para = .true.
    lrhosym = .false.

  end subroutine init_misc_variables

  !-------------------------------------------------------------------------------
  !> Summary: 
  !> Author: 
  !> Category: KKRhost, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> 
  !-------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------
  ! subroutine: init_relativistic_variables
  !> @brief set default values for the relativistic variables
  !> @details The idea behind this kind of routine is to separate the initialization
  !> of variables such that one can use this to modularize the code
  !> @author Jonathan Chico
  !> @date 26.12.2017
  ! -------------------------------------------------------------------------
  subroutine init_relativistic_variables()

    implicit none

    krel = 0                       ! Switch for non- (or scalar-) relativistic/relativistic (Dirac) program (0/1). Attention: several other parameters depend explicitly on KREL, they are set automatically Used for Dirac solver in ASA
    kvrel = 1                      ! Scalar-relativistic calculation
    korbit = 0                     ! No spin-orbit coupling
    lnc = .true.                   ! Coupled equations in two spins (switches true if KREL=1 or KORBIT=1 or KNOCO=1)

  end subroutine init_relativistic_variables

  !-------------------------------------------------------------------------------
  !> Summary: 
  !> Author: 
  !> Category: KKRhost, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> 
  !-------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------
  ! subroutine: init_cluster_variables
  !> @brief set default values for the cluster variables
  !> @details The idea behind this kind of routine is to separate the initialization
  !> of variables such that one can use this to modularize the code
  !> @author Jonathan Chico
  !> @date 26.12.2017
  ! -------------------------------------------------------------------------
  subroutine init_cluster_variables()

    implicit none

    nclsd = 2                      ! NAEZD + NEMBD maximum number of different TB-clusters
    naclsd = 500                   ! Maximum number of atoms in a TB-cluster
    nofgij = 2                     ! NATOMIMPD*NATOMIMPD+1 probably the same variable than NOFGIJD
    natomimp = 0                   ! Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
    natomimpd = 150                ! Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
    i25 = 'scoef                                   ' ! Default name of scoef file

  end subroutine init_cluster_variables

  !-------------------------------------------------------------------------------
  !> Summary: 
  !> Author: 
  !> Category: KKRhost, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> 
  !-------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------
  ! subroutine: init_io_variables
  !> @brief set default values for the I/O variables
  !> @details The idea behind this kind of routine is to separate the initialization
  !> of variables such that one can use this to modularize the code
  !> @author Jonathan Chico
  !> @date 26.12.2017
  ! -------------------------------------------------------------------------
  subroutine init_io_variables()

    implicit none

    igf = 0                        ! Not printing the Green functions
    icc = 0                        ! Not printing the Green functions
    wlength = 1                    ! Word length for direct access files, compiler dependent ifort/others (1/4)
#ifdef __GFORTRAN__
    wlength = 4
#endif
    i12 = '                                        '
    i40 = '                                        '

  end subroutine init_io_variables

  !-------------------------------------------------------------------------------
  !> Summary: 
  !> Author: 
  !> Category: KKRhost, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> 
  !-------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------
  ! subroutine: init_cell_variables
  !> @brief set default values for the unit cell variables
  !> @details The idea behind this kind of routine is to separate the initialization
  !> of variables such that one can use this to modularize the code
  !> @author Jonathan Chico
  !> @date 26.12.2017
  ! -------------------------------------------------------------------------
  subroutine init_cell_variables()

    implicit none

    naez = 1                       ! Number of atoms in the unit cell
    nemb = 1                       ! Number of embedded atoms
    ncls = 1                       ! Number of clusters
    natyp = 1                      ! Number of kinds of atoms in the unit cell
    nineq = 1
    alat = 1.0d0                   ! Lattice constant in a.u.
    abasis = 1.0d0                 ! Scaling factors for rbasis
    bbasis = 1.0d0                 ! Scaling factors for rbasis
    cbasis = 1.0d0                 ! Scaling factors for rbasis
    alatnew = 1.0d0
    volume0 = 1.0d0
    lcartesian = .false.           ! True: Basis in cartesian coords; false: in internal coords
    linterface = .false.           ! If True a matching with semi-inifinite surfaces must be performed

  end subroutine init_cell_variables

  !-------------------------------------------------------------------------------
  !> Summary: 
  !> Author: 
  !> Category: KKRhost, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> 
  !-------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------
  ! subroutine: init_slab_variables
  !> @brief set default values for the slab calculation variables
  !> @details The idea behind this kind of routine is to separate the initialization
  !> of variables such that one can use this to modularize the code
  !> @author Jonathan Chico
  !> @date 26.12.2017
  ! -------------------------------------------------------------------------
  subroutine init_slab_variables()

    implicit none

    nleft = 1
    nright = 1
    nlayer = 1                     ! Number of principal layer
    nlbasis = 0                    ! Number of basis layers of left host (repeated units)
    nrbasis = 0                    ! Number of basis layers of right host (repeated units)
    nprincd = 1                    ! Number of principle layers, set to a number >= NRPINC in output of main0
    nlayerd = 1                    ! (NAEZD/NPRINCD)

  end subroutine init_slab_variables

  !-------------------------------------------------------------------------------
  !> Summary: 
  !> Author: 
  !> Category: KKRhost, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> 
  !-------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------
  ! subroutine: init_energy_variables
  !> @brief set default values for the energy variables
  !> @details The idea behind this kind of routine is to separate the initialization
  !> of variables such that one can use this to modularize the code
  !> @author Jonathan Chico
  !> @date 26.12.2017
  ! -------------------------------------------------------------------------
  subroutine init_energy_variables()

    implicit none

    lly = 0                        ! No Lloyds formula
    kte = 1                        ! Calculate the total energy
    npol = 7                       ! Number of poles
    iemxd = 101                    ! Dimension for energy-dependent arrays
    npnt1 = 7                      ! Number of points in the energy contour
    npnt2 = 7                      ! Number of points in the energy contour
    npnt3 = 7                      ! Number of points in the energy contour
    n1semi = 0                     ! Number of energy points for the semicore contour
    n2semi = 0                     ! Number of energy points for the semicore contour
    n3semi = 0                     ! Number of energy points for the semicore contour
    ivshift = 0
    npolsemi = 0
    iesemicore = 0
    idosemicore = 0
    tk = 800.d0                    ! Temperature
    fcm = 20.0d0
    e2in = 0.0d0
    emin = -0.30d0                 ! Energies needed in EMESHT
    emax = 0.70d0                  ! Energies needed in EMESHT
    efermi = 0.d0                  ! Setting the Fermi energy to zero
    eshift = 0.d0
    tksemi = 800.d0                ! Temperature for the semi-core contour
    emusemi = 0.0d0
    ebotsemi = -0.5d0
    fsemicore = 0.0d0
    deltae = (1.d-5, 0.d0)

  end subroutine init_energy_variables

  !-------------------------------------------------------------------------------
  !> Summary: 
  !> Author: 
  !> Category: KKRhost, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> 
  !-------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------
  ! subroutine: init_convergence_variables
  !> @brief set default values for the convergence and solver variables
  !> @details The idea behind this kind of routine is to separate the initialization
  !> of variables such that one can use this to modularize the code
  !> @author Jonathan Chico
  !> @date 26.12.2017
  ! -------------------------------------------------------------------------
  subroutine init_convergence_variables()

    implicit none

    imix = 0                       ! Straight mixing
    ishld = 200000                 ! Paremeters for the Ewald summations
    nmaxd = 2000000                ! Paremeters for the Ewald summations
    ncheb = 10                     ! Number of Chebychev pannels for the new solver
    itdbry = 10
    ntrefd = 0                     ! Parameter in broyden subroutine MUST BE 0 for the host program
    ntperd = 1                     ! Parameter in broyden subroutines
    npan_eq = 30
    npan_log = 30
    scfsteps = 100                 ! Number of SCF steps
    rmax = 7.0d0                   ! Ewald summation cutoff parameter for real space summation
    gmax = 100.d0                  ! Ewald summation cutoff parameter for reciprocal space summation
    r_log = 0.1d0
    rcutz = 2.30d0                 ! Parameter for the screening cluster along the z-direction
    rcutxy = 2.30d0                ! Parameter for the screening cluster along the x-y plane
    qbound = 1.d-7                 ! Convergence parameter for the potential
    mixing = 0.001d0               ! Magnitude of the mixing parameter
    tolrdif = 1.0d0                ! Tolerance for r<tolrdif (a.u.) to handle vir. atoms
    solver = 'BS        '          ! Set the BS-solver as default

  end subroutine init_convergence_variables

  !-------------------------------------------------------------------------------
  !> Summary: 
  !> Author: 
  !> Category: KKRhost, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> 
  !-------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------
  ! subroutine: init_potential_variables
  !> @brief set default values for the potential variables
  !> @details The idea behind this kind of routine is to separate the initialization
  !> of variables such that one can use this to modularize the code
  !> @author Jonathan Chico
  !> @date 26.12.2017
  ! -------------------------------------------------------------------------
  subroutine init_potential_variables()

    implicit none

    kws = 1                        ! FP/ASA potential ()
    kxc = 2                        ! VWN potential
    ins = 1                        ! FP/ASA calculation ()
    icst = 2                       ! Number of Born approximation
    irid = 350                     ! Shape functions parameters in non-spherical part
    nref = 1                       ! Number of reference potentials
    lpot = 3                       ! Maximum l for expansion of the potential
    nfund = 289                    ! Shape functions parameters in non-spherical part
    ngshd = 60000                  ! Shape functions parameters in non-spherical part
    ipand = 50                     ! Number of panels in non-spherical part
    ntotd = 80                     ! IPAND+30
    ifile = 13                     ! File identifier of the potential file
    lmpot = 16                     ! (LPOT+1)**2
    insref = 0                     ! INS For reference potential
    knosph = 1                     ! Switch for spherical/non-spherical(0/1) program. Same obs. as for KREL applies.
    kshape = 2                     ! FP/ASA calculation (2/0)
    ncelld = 1                     ! Number of cells (shapes) in non-spherical part
    nspotd = 2                     ! Number of potentials for storing non-sph. potentials (2*KREL+(1-KREL)*NSPIND)*NSATYPD
    nsatypd = 1                    ! (NATYPD-1)*KNOSPH+1
    vconst = 0.d0                  ! Potential shift
    lambda_xc = 1.0d0              ! Scale magnetic moment (0 < Lambda_XC < 1, 0=zero moment, 1= full moment)
    i13 = 'potential                               ' ! Default name of potential file
    i19 = 'shapefun                                ' ! Default name of shapefunction file

  end subroutine init_potential_variables

  !-------------------------------------------------------------------------------
  !> Summary: 
  !> Author: 
  !> Category: KKRhost, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> 
  !-------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------
  ! subroutine: init_angular_momentum_variables
  !> @brief set default values for the angular momentum variables
  !> @details The idea behind this kind of routine is to separate the initialization
  !> of variables such that one can use this to modularize the code
  !> @author Jonathan Chico
  !> @date 26.12.2017
  ! -------------------------------------------------------------------------
  subroutine init_angular_momentum_variables()

    implicit none

    lmax = 3                       ! Maximum l for expansion
    ncleb = 784                    ! (LMAX*2+1)**2 * (LMAX+1)**2
    lmmax = 16                     ! (LMAX+1)**2
    lmmaxd = 32                    ! (KREL+KORBIT+1)*(LMAX+1)**2

  end subroutine init_angular_momentum_variables

  !-------------------------------------------------------------------------------
  !> Summary: 
  !> Author: 
  !> Category: KKRhost, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> 
  !-------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------
  ! subroutine: init_magnetization_variables
  !> @brief set default values for the magnetisation variables for the calculation
  !> @details the idea behind this kind of routine is to separate the initialization
  !> of variables such that one can use this to modularize the code
  !> @author Jonathan Chico
  !> @date 22.12.2017
  ! -------------------------------------------------------------------------
  subroutine init_magnetization_variables()

    implicit none

    kmrot = 0                      ! 0: no rotation of the magnetisation; 1: individual rotation of the magnetisation for every site
    knoco = 0                      ! Collinear calculation
    nspin = 2                      ! Spin polarised calculation
    nspind = 2                     ! KREL+(1-KREL)*(NSPIN+1)
    khfeld = 0                     ! No external magnetic field
    nspindd = 1                    ! NSPIND-KORBIT
    hfield = 0.d0                  ! External magnetic field, for initial potential shift in spin polarised case
    linipol = .false.              ! True: Initial spin polarization; false: no initial spin polarization

  end subroutine init_magnetization_variables

  !-------------------------------------------------------------------------------
  !> Summary: 
  !> Author: 
  !> Category: KKRhost, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> 
  !-------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------
  ! SUBROUTINE: init_CPA_variables
  !> @brief Set default values for the CPA variables for the calculation
  !> @details The idea behind this kind of routine is to separate the initialization
  !> of variables such that one can use this to modularize the code
  !> @author Jonathan Chico
  !> @date 22.12.2017
  ! -------------------------------------------------------------------------
  subroutine init_cpa_variables()

    implicit none

    ncpa = 0                       ! NCPA = 0/1 CPA flag
    itcpamax = 0                   ! Max. number of CPA iterations
    cpatol = 1d-4                  ! Convergency tolerance for CPA-cycle

  end subroutine init_cpa_variables

  !-------------------------------------------------------------------------------
  !> Summary: 
  !> Author: 
  !> Category: KKRhost, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> 
  !-------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------
  ! SUBROUTINE: init_CPA_variables
  !> @brief Set default values for the LDA+U variables for the calculation
  !> @details The idea behind this kind of routine is to separate the initialization
  !> of variables such that one can use this to modularize the code
  !> @author Jonathan Chico
  !> @date 22.12.2017
  ! -------------------------------------------------------------------------
  subroutine init_ldau_variables()

    implicit none

    ntldau = 0                     ! number of atoms on which LDA+U is applied
    idoldau = 0                    ! flag to perform LDA+U
    itrunldau = 0                  ! iteration index
    kreadldau = 0                  ! LDA+U arrays available

  end subroutine init_ldau_variables

  !-------------------------------------------------------------------------------
  !> Summary: 
  !> Author: 
  !> Category: KKRhost, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> 
  !-------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------
  ! subroutine: init_mesh_variables
  !> @brief set default values for the mesh variables
  !> @details The idea behind this kind of routine is to separate the initialization
  !> of variables such that one can use this to modularize the code
  !> @author Jonathan Chico
  !> @date 26.12.2017
  ! -------------------------------------------------------------------------
  subroutine init_mesh_variables()

    implicit none

    nrd = 20000                    ! Number of real space
    irmd = 900                     ! Number of radial mesh points in (0,...,RWS)
    irnsd = 890                    ! Number of radial mesh points in (RMT,...,RWS)
    kpoibz = 250000                ! Number of reciprocal space vectors
    intervx = 0                    ! Number of intervals in x-direction for k-net in IB of the BZ
    intervy = 0                    ! Number of intervals in y-direction for k-net in IB of the BZ
    intervz = 0                    ! Number of intervals in z-direction for k-net in IB of the BZ
    maxmesh = 1

  end subroutine init_mesh_variables


  !-------------------------------------------------------------------------------
  !> Summary: 
  !> Author: 
  !> Category: KKRhost, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> 
  !-------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------
  ! subroutine: init_all_wrapper
  !> @brief wrapper for initialization subroutines and allocation of test/opt arrays
  !> @author Philipp Ruessmann
  !> @date 22.08.2018
  ! -------------------------------------------------------------------------
  subroutine init_all_wrapper()
    use :: mod_wunfiles, only: t_params
    implicit none
    integer :: i_stat

    ! Set the default values for the I/O variables
    call init_io_variables()
    ! Set the defaults values for the CPA variables
    call init_cpa_variables()
    ! Set the default values for the Slab mode variables
    call init_slab_variables()
    ! Set the default values for the LDA+U variables
    call init_ldau_variables()
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
    allocate (t_params%optc(32), stat=i_stat) ! CHARACTER*8
    call memocc(i_stat, product(shape(t_params%optc))*kind(t_params%optc), 't_params%OPTC', 'main0')
    t_params%optc(1:32) = '        '
    allocate (t_params%testc(32), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%testc))*kind(t_params%testc), 't_params%TESTC', 'main0')
    t_params%testc(1:32) = '        '

  end subroutine init_all_wrapper


  !-------------------------------------------------------------------------------
  !> Summary: 
  !> Author: 
  !> Category: KKRhost, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> 
  !-------------------------------------------------------------------------------
  subroutine print_versionserial(iunit, version1, version2, version3, version4, serialnr)
    implicit none
    integer, intent (in) :: iunit
    character (len=*), intent (in) :: version1
    character (len=*), intent (in) :: version2
    character (len=*), intent (in) :: version3
    character (len=*), intent (in) :: version4
    character (len=*), intent (in) :: serialnr

    write (iunit, '(1A)') '     Screened Korringa-Kohn-Rostoker Electronic Structure Code'
    write (iunit, '(1A)') '                      for Bulk and Interfaces'
    write (iunit, '(1A)') '                    Juelich-Munich 2001 - 2018'
    write (iunit, '(1A)') ''
    write (iunit, '(2A)') '  Code version: ', trim(version1)
    write (iunit, '(6A)') '  Compile options: ', trim(version2), ' ', trim(version3), ' ', trim(version4)
    write (iunit, '(2A)') '  serial number for files: ', serialnr
  end subroutine print_versionserial

end module mod_main0
