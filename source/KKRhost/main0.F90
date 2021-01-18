!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Wrapper module for the reading and setup of the JM-KKR program
!> Author: Philipp Ruessmann, Bernd Zimmermann, Phivos Mavropoulos, R. Zeller,
!> and many others ...
!> Wrapper module for the reading and setup of the JM-KKR program, it also calculates
!> the electrostatic potential.
!------------------------------------------------------------------------------------
!> @todo JC: `NATOMIMP` and `NATOMIMPD` seem to be the same variable, however, right
!> now find no way to eliminate one of them.
!> @endtodo
!> @todo JC: There seem to be several repeated variables doing the same, e.g. `INS`,
!> `KNOSPH`, `KWS` and `KSHAPE`, all seem to dictate whether one has ASA or FP.
!> Maybe it would be good to consolidate and eliminate any unnecessary variables.
!> @endtodo
!> @todo JC: Several variables such as `IRMD` and `IRNSD` are actually determined in
!> the `startb1()` subroutine, maybe change the allocations such that they are done
!> there instead
!> @endtodo
!------------------------------------------------------------------------------------

! set CPP_OMPSTUFF if either HYBRID or OpenMP parallelization is chosen
#ifdef CPP_HYBRID
#define CPP_OMPSTUFF
#endif
#ifdef CPP_OMP
#define CPP_OMPSTUFF
#endif


!-------------------------------------------------------------------------------
!> Summary: Wrapper module for the reading and setup of the JM-KKR program
!> Author: 
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!>
!> @todo 
!> - JC: NATOMIMP and NATOMIMPD seem to be the same variable, however, right
!> now find no way to eliminate one of them.
!> - JC: There seem to be several repeated variables doing the same, e.g. INS,
!> KNOSPH, KWS and KSHAPE, all seem to dictate whether one has ASA or FP.
!> Maybe it would be good to consolidate and eliminate any unnecessary variables.
!> - JC: Several variables such as IRMD and IRNSD are actually determined in
!> the startb1 subroutine, maybe change the allocations such that they are done
!> there instead
!> @endtodo
!>
!> @note
!> - Jonathan Chico Jan. 2018: Removed inc.p dependencies and rewrote to Fortran90
!> @endnote
!-------------------------------------------------------------------------------
module mod_main0


  use :: mod_datatypes, only: dp
  use :: mod_constants, only: nsymaxd, pi

  implicit none

  private
  public :: main0, bshift_ns
  ! ------------- > scalars > ------------- 
  !integers
  public :: kte, kws, kxc, igf, icc, ins, irm, ipe, ipf, ipfe, kcor, kefg, khyp, kpre, nprinc, nsra, lpot, imix, iend, icst, &
    naez, nemb, lmax, ncls, nref, npol, npnt1, npnt2, npnt3, lmmax0d, nvirt, lmpot, kvmad, itscf, ncheb, nineq, natyp, ifile, &
    kvrel, nspin, nleft, nright, invmod, khfeld, itdbry, insref, kshape, ielast, ishift, ivshift, kfrozn, nsymat, nqcalc, kforce, n1semi, &
    n2semi, n3semi, nlayer, nlbasis, nrbasis, intervx, intervy, intervz, maxmesh, npan_eq, npan_log, npolsemi, scfsteps, natomimp, &
    iesemicore, idosemicore
  !real(kind=dp)
  public :: tk, fcm, e2in, emin, emax, alat, rmax, gmax, r_log, rcutz, rcutxy, qbound, vconst, hfield, mixing, abasis, bbasis, &
    cbasis, efermi, eshift, tksemi, tolrdif, alatnew, volume0, emusemi, ebotsemi, fsemicore, lambda_xc
  !character
  public :: solver, i12, i13, i19, i25, i40
  !logicals
  public :: lrhosym, linipol, lcartesian
  ! ------------- < scalars < ------------- 
  ! ------------- > arrays > ------------- 
  !integer
  public :: isymindex, cls, irc, imt, nfu, nsh1, nsh2, lmxc, ipan, irns, irws, kmesh, irmin, loflm, nacls, ncore, imaxsh, nshell, &
    inipol, ixipol, refpot, ntcell, iqcalc, iofgij, jofgij, atomimp, ijtabsh, ijtabsym, npan_tot, ijtabcalc, npan_eq_at, npan_log_at, &
    ijtabcalc_i, ish, jsh, ilm_map, kfg, atom, ezoa, lmsp, lcore, icleb, ircut, llmsp, lmsp1, kaoez, ifunm, ifunm1, ititle, icheck, &
    ipan_intervall, jend, kmrot, ncpa, itcpamax, noq, iqat, icpa, hostimp, zrel, jwsrel, irshift, nrrel, ntldau, idoldau, itrunldau, &
    kreadldau, lopt, itldau, lly, irrel
  !real
  public :: vbc, zperight, zperleft, recbv, bravais, rsymat, a, b, wg, gsh, zat, rmt, rws, vref, mtfac, rmtnew, rmtref, &
    rmtrefat, fpradius, socscale, rmesh, s, rr, drdi, dror, cleb, visp, cscl, rnew, ratom, ecore, tleft, tright, socscl, &
    rbasis, rclsimp, cmomhost, rpan_intervall, rs, yrg, vins, rcls, rrot, qmtet, qmphi, qmgam, qmgamtab, qmphitab, qmtettab, cpatol, &
    conc, fact, vtrel, btrel, rmrel, drdirel, r2drdirel, thesme, thetas, thetasnew, ueff, jeff, erefldau, wldau, uldau
  !complex
  public :: ez, dez, wez, dsymll, dsymll1, lefttinvll, righttinvll, rc, crel, rrel, srrel, drotq, phildau, deltae
  !character
  public :: txc
  !logical
  public :: vacflag, para, symunitary, emeshfile
  ! ------------- < arrays < ------------- 

  
  ! decalration of common variables

  integer :: kte = 1               !! Calculation of the total energy On/Off (1/0)
  integer :: kws = 1               !! 0 (MT), 1(ASA)
  integer :: kxc = 2               !! Type of xc-potential 0=MJW 1=vBH 2=VWN 3=PW91 4=PBE 5=PBEsol
  integer :: igf = 0               !! Do not print or print (0/1) the KKRFLEX_* files
  integer :: icc = 0               !! Enables the calculation of off-diagonal elements of the GF.(0=SCF/DOS; 1=cluster; -1=custom)
  integer :: ins = 2               !! 0 (MT), 1(ASA), 2(Full Potential)
  integer :: irm                   !! Maximum number of radial points
  integer :: ipe = 0               !! Not real used, IPFE should be 0
  integer :: ipf = 0               !! Not real used, IPFE should be 0
  integer :: ipfe = 0              !! Not real used, IPFE should be 0
  integer :: kcor = 0
  integer :: kefg = 0
  integer :: khyp = 0
  integer :: kpre = 0
  integer :: nprinc = 1            !! number of principal layers used for slab-inversion
  integer :: nsra = 1              !! scalar-relativistic (nsra==2) or non-relativistic (nsra==1)
  integer :: lpot = 2              !! Maximum l component in potential expansion
  integer :: imix = 0              !! Type of mixing scheme used (0=straight, 4=Broyden 2nd, 5=Anderson)
  integer :: iend = 1              !! Number of nonzero gaunt coefficients
  integer :: icst = 2              !! Number of Born approximation
  integer :: naez = 1              !! Number of atoms in unit cell
  integer :: nemb = 1              !! Number of 'embedding' positions
  integer :: lmax = 2              !! Maximum l component in wave function expansion
  integer :: ncls = 1              !! Number of reference clusters
  integer :: nref = 1              !! Number of diff. ref. potentials
  integer :: npol = 7              !! Number of Matsubara Poles (EMESHT)
  integer :: npnt1 = 3             !! number of E points (EMESHT) for the contour integration going up
  integer :: npnt2 = 10            !! number of E points (EMESHT) for the contour integration goind parallel to the real axis
  integer :: npnt3 = 4             !! number of E points (EMESHT) for the contour integration going down
  integer :: lmmax0d = 16          !! (LMAX+1)^2 wihtout spin doubling
  integer :: nvirt = 0
  integer :: lmpot = 16            !! (LPOT+1)**2
  integer :: kvmad = 0
  integer :: itscf = 0             !! counter scf iterations
  integer :: ncheb = 10            !! Number of Chebychev pannels for the new solver
  integer :: nineq = 1             !! Number of ineq. positions in unit cell
  integer :: natyp = 1             !! Number of kinds of atoms in unit cell
  integer :: ifile = 13            !! Unit specifier for potential card
  integer :: kvrel = 1             !! 0,1,2 : non / scalar relat. / full Dirac calculation
  integer :: nspin = 2             !! Counter for spin directions
  integer :: nleft = 1             !! Number of repeated basis for left host to get converged electrostatic potentials
  integer :: nright = 1            !! Number of repeated basis for right host to get converged electrostatic potentials
  integer :: invmod = -1           !! inversion mode, 0=full inversion, 1= banded matrix, 2= supercell, 3=godfrin; default signals that inversion mode was not set yet
  integer :: khfeld = 0            !! 0,1: no / yes external magnetic field
  integer :: itdbry = 10           !! Number of SCF steps to remember for the Broyden mixing
  integer :: insref = 0            !! INS for reference pot. (usual 0)
  integer :: kshape = 2            !! Exact treatment of WS cell
  integer :: ielast = 0            !! number of energy points in complex energy contour
  integer :: ishift = 0            !! Parameter controling the potential shift after mixing
  integer :: kfrozn = 0
  integer :: nsymat = 0
  integer :: nqcalc = 0
  integer :: kforce = 0            !! Calculation of the forces
  integer :: n1semi = 0            !! Number of energy points for the semicore contour
  integer :: n2semi = 0            !! Number of energy points for the semicore contour
  integer :: n3semi = 0            !! Number of energy points for the semicore contour
  integer :: nlayer = 1            !! Number of principal layer
  integer :: nlbasis = 0           !! Number of basis layers of left host (repeated units)
  integer :: nrbasis = 0           !! Number of basis layers of right host (repeated units)
  integer :: intervx = 10          !! Number of intervals in x-direction for k-net in IB of the BZ
  integer :: intervy = 10          !! Number of intervals in y-direction for k-net in IB of the BZ
  integer :: intervz = 10          !! Number of intervals in z-direction for k-net in IB of the BZ
  integer :: maxmesh = 1           !! Number of different k-meshes
  integer :: npan_eq = 30          !! Number of intervals from [R_LOG] to muffin-tin radius Used in conjunction with runopt NEWSOSOL
  integer :: npan_log = 30         !! Number of intervals from nucleus to [R_LOG] Used in conjunction with runopt NEWSOSOL
  integer :: npolsemi = 0          !! Number of poles for the semicore contour
  integer :: scfsteps = 1          !! number of scf iterations
  integer :: natomimp = 0          !! Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
  integer :: iesemicore = 0
  integer :: idosemicore = 0
  integer :: ivshift = 0           !! for selected potential shift: atom index of potentials to be shifted by VCONST
  integer :: verbosity = 0         !! verbosity level for timings and output: 0=old default, 1,2,3 = timing and ouput verbosity level the same (low,medium,high)
  integer :: MPI_scheme = 0        !! scheme for MPI parallelization: 0 = determine automatically (default), 1 = atoms, 2 = energies, 3 = after 2 runs with (1 and 2), select best option.
  integer :: special_straight_mixing = 0 !!id to specify modified straight mixing scheme: 0=normal, 1=alternating mixing factor (i.e. reduced mixing factor in every odd iteration), 2=charge-neurality based mixing factor (former: 'alt mix' and 'spec mix')
  real (kind=dp) :: tk = 800.0_dp       !! Temperature
  real (kind=dp) :: fcm = 20.0_dp       !! Factor for increased linear mixing of magnetic part of potential compared to non-magnetic part.
  real (kind=dp) :: e2in = 0.0_dp    
  real (kind=dp) :: emin = -0.30_dp     !! Lower value (in Ryd) for the energy contour
  real (kind=dp) :: emax = 0.70_dp      !! Maximum value (in Ryd) for the DOS calculation Controls also [NPT2] in some cases
  real (kind=dp) :: alat = 1.0_dp       !! Lattice constant in a.u.
  real (kind=dp) :: rmax = 10.0_dp      !! Ewald summation cutoff parameter for real space summation
  real (kind=dp) :: gmax = 100.0_dp     !! Ewald summation cutoff parameter for reciprocal space summation
  real (kind=dp) :: r_log = 0.5_dp      !! Radius up to which log-rule is used for interval width. Used in conjunction with runopt NEWSOSOL
  real (kind=dp) :: rcutz = 2.30_dp     !! Parameter for the screening cluster along the z-direction
  real (kind=dp) :: rcutxy = 2.30_dp    !! Parameter for the screening cluster along the x-y plane
  real (kind=dp) :: qbound = 1.0e-7_dp  !! Convergence parameter for the potential
  real (kind=dp) :: vconst = 0.0_dp     !! Value of potential shift in the first iteration in Ry
  real (kind=dp) :: hfield = 0.0_dp     !! Value of external magnetic field in Ry, for initial potential shift in spin polarised case
  real (kind=dp) :: mixing = 0.01_dp    !! Magnitude of the mixing parameter
  real (kind=dp) :: abasis = 1.0_dp     !! Scaling factors for rbasis
  real (kind=dp) :: bbasis = 1.0_dp     !! Scaling factors for rbasis
  real (kind=dp) :: cbasis = 1.0_dp     !! Scaling factors for rbasis
  real (kind=dp) :: efermi = 0.0_dp     !! Fermi energy
  real (kind=dp) :: eshift = 0.0_dp
  real (kind=dp) :: tksemi = 800.0_dp   !! Temperature of semi-core contour
  real (kind=dp) :: tolrdif = 1.0_dp    !! For distance between scattering-centers smaller than [<TOLRDIF>], free GF is set to zero. Units are Bohr radii.
  real (kind=dp) :: alatnew = 1.0_dp
  real (kind=dp) :: volume0 = 1.0_dp    !! Unit cell volume
  real (kind=dp) :: emusemi = 0.0_dp    !! Top of semicore contour in Ryd.
  real (kind=dp) :: ebotsemi = -0.50_dp !! Bottom of semicore contour in Ryd
  real (kind=dp) :: fsemicore = 0.00_dp !! Initial normalization factor for semicore states (approx. 1.)
  real (kind=dp) :: lambda_xc = 1.0_dp !! Scale magnetic moment (0 < Lambda_XC < 1, 0=zero moment, 1= full moment)
  character (len=10) :: solver = 'BS        ' !! Type of solver
  character (len=40) :: i12 = '                                        ' !! File identifiers
  character (len=40) :: i13 = 'potential                               ' !! Potential file name
  character (len=40) :: i19 = 'shapefun                                ' !! Shape function file name
  character (len=40) :: i25 = 'scoef                                   ' !! Default name of scoef file
  character (len=40) :: i40 = '                                        ' !! File identifiers
  logical :: lrhosym = .false.        !! 
  logical :: linipol = .false.        !! True: Initial spin polarization; false: no initial spin polarization
  logical :: lcartesian = .false.     !! True: Basis in cartesian coords; false: in internal coords
  ! ..
  ! .. Arrays ..
  integer, dimension (nsymaxd) :: isymindex
  integer, dimension (:), allocatable :: cls  !! Cluster around atomic sites
  integer, dimension (:), allocatable :: irc  !! R point for potential cutting
  integer, dimension (:), allocatable :: imt  !! R point at MT radius
  integer, dimension (:), allocatable :: nfu  !! number of shape function components in cell 'icell'
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
  integer, dimension (:, :), allocatable :: atom  !! Atom at site in cluster
  integer, dimension (:, :), allocatable :: ezoa  !! EZ of atom at site in cluster
  integer, dimension (:, :), allocatable :: lmsp  !! 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
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
  real (kind=dp), dimension (2) :: vbc              !! Potential constants
  real (kind=dp), dimension (3) :: zperight         !! Vector to define how to repeat the basis of the right host
  real (kind=dp), dimension (3) :: zperleft         !! Vector to define how to repeat the basis of the left host
  real (kind=dp), dimension (3, 3) :: recbv         !! Reciprocal basis vectors
  real (kind=dp), dimension (3, 3) :: bravais       !! Bravais lattice vectors
  real (kind=dp), dimension (64, 3, 3) :: rsymat
  real (kind=dp), dimension (:), allocatable :: a   !! Constants for exponential R mesh
  real (kind=dp), dimension (:), allocatable :: b   !! Constants for exponential R mesh
  real (kind=dp), dimension (:), allocatable :: wg  !! Integr. weights for Legendre polynomials
  real (kind=dp), dimension (:), allocatable :: gsh 
  real (kind=dp), dimension (:), allocatable :: zat !! Nuclear charge
  real (kind=dp), dimension (:), allocatable :: rmt !! Muffin-tin radius of true system
  real (kind=dp), dimension (:), allocatable :: rws !! Wigner Seitz radius
  real (kind=dp), dimension (:), allocatable :: vref
  real (kind=dp), dimension (:), allocatable :: mtfac       !! Scaling factor for radius MT
  real (kind=dp), dimension (:), allocatable :: rmtnew      !! Adapted muffin-tin radius
  real (kind=dp), dimension (:), allocatable :: rmtref      !! Muffin-tin radius of reference system
  real (kind=dp), dimension (:), allocatable :: rmtrefat
  real (kind=dp), dimension (:), allocatable :: fpradius !! R point at which full-potential treatment starts
  real (kind=dp), dimension (:), allocatable :: socscale !! Spin-orbit scaling
  real (kind=dp), dimension (:, :), allocatable :: rmesh !! Radial mesh ( in units a Bohr)
  real (kind=dp), dimension (:, :), allocatable :: s
  real (kind=dp), dimension (:, :), allocatable :: rr   !! Set of real space vectors (in a.u.)
  real (kind=dp), dimension (:, :), allocatable :: drdi !! Derivative dr/di
  real (kind=dp), dimension (:, :), allocatable :: dror
  real (kind=dp), dimension (:, :), allocatable :: cleb !! GAUNT coefficients (GAUNT)
  real (kind=dp), dimension (:, :), allocatable :: visp !! Spherical part of the potential
  real (kind=dp), dimension (:, :), allocatable :: cscl !! Speed of light scaling
  real (kind=dp), dimension (:, :), allocatable :: rnew
  real (kind=dp), dimension (:, :), allocatable :: ratom
  real (kind=dp), dimension (:, :), allocatable :: ecore  !! Core energies
  real (kind=dp), dimension (:, :), allocatable :: tleft  !! Vectors of the basis for the left host
  real (kind=dp), dimension (:, :), allocatable :: tright !! Vectors of the basis for the right host
  real (kind=dp), dimension (:, :), allocatable :: socscl
  real (kind=dp), dimension (:, :), allocatable :: rbasis !! Position of atoms in the unit cell in units of bravais vectors
  real (kind=dp), dimension (:, :), allocatable :: rclsimp
  real (kind=dp), dimension (:, :), allocatable :: cmomhost !! Charge moments of each atom of the (left/right) host
  real (kind=dp), dimension (:, :), allocatable :: rpan_intervall
  real (kind=dp), dimension (:, :, :), allocatable :: rs
  real (kind=dp), dimension (:, :, :), allocatable :: yrg  !! Spherical harmonics (GAUNT2)
  real (kind=dp), dimension (:, :, :), allocatable :: vins !! Non-spherical part of the potential
  real (kind=dp), dimension (:, :, :), allocatable :: rcls !! Real space position of atom in cluster
  real (kind=dp), dimension (:, :, :), allocatable :: rrot
  complex (kind=dp), dimension (:), allocatable :: ez  !! complex energy points
  complex (kind=dp), dimension (:), allocatable :: dez !! length of energy interval: (EF-EMIN)/NEPTS for uniform distribution of points
  complex (kind=dp), dimension (:), allocatable :: wez !! integration weights: wez(ie) = -2.0_dp/pi*dez(ie)
  complex (kind=dp), dimension (:, :, :), allocatable :: dsymll
  complex (kind=dp), dimension (:, :, :), allocatable :: dsymll1
  complex (kind=dp), dimension (:, :, :, :, :), allocatable :: lefttinvll
  complex (kind=dp), dimension (:, :, :, :, :), allocatable :: righttinvll
  character (len=124), dimension (6) :: txc
  logical, dimension (2) :: vacflag

  ! -------------------------------------------------------------------------
  ! Magnetisation angles -- description see RINPUT13
  ! -------------------------------------------------------------------------
  integer :: kmrot = 0                                    !! 0: no rotation of the magnetisation; 1: individual rotation of the magnetisation for every site
  real (kind=dp), dimension (:), allocatable :: qmtet !! \( \theta\) angle of the magnetization with respect to the z-axis
  real (kind=dp), dimension (:), allocatable :: qmphi !! \( \phi\) angle of the magnetization with respect to the z-axis
  ! -------------------------------------------------------------------------
  ! CPA variables
  ! -------------------------------------------------------------------------
  integer :: ncpa = 0                             !! NCPA = 0/1 CPA flag
  integer :: itcpamax = 0                         !! Max. number of CPA iterations
  integer, dimension (:), allocatable :: noq  !! Number of diff. atom types located
  integer, dimension (:), allocatable :: iqat !! The site on which an atom is located on a given site
  integer, dimension (:), allocatable :: icpa !! ICPA = 0/1 site-dependent CPA flag

  ! -------------------------------------------------------------------------
  !> @note ITERMDIR running option introduced Apr 2003 -- Munich
  !>              (H. Ebert + V. Popescu) allows a self-consistent
  !>              determination of the magnetic configuration in REL mode
  !> @endnote
  ! -------------------------------------------------------------------------
  real (kind=dp), dimension (:), allocatable :: qmgam
  real (kind=dp), dimension (:, :), allocatable :: qmgamtab
  real (kind=dp), dimension (:, :), allocatable :: qmphitab
  real (kind=dp), dimension (:, :), allocatable :: qmtettab
  ! -------------------------------------------------------------------------
  !> @note changes for impurity 20/02/2004 -- v.popescu according to
  !>                                          n.papanikolaou
  !> @endnote
  ! -------------------------------------------------------------------------
  integer, dimension (:), allocatable :: hostimp     !! 
  real (kind=dp) :: cpatol = 1e-4_dp                           !! Convergency tolerance for CPA-cycle
  real (kind=dp), dimension (:), allocatable :: conc !! Concentration of a given atom
  ! -------------------------------------------------------------------------------
  complex (kind=dp), dimension (:, :), allocatable :: rc       !! NREL REAL spher. harm. > CMPLX. spher. harm. NREL CMPLX. spher. harm. > REAL spher. harm.
  complex (kind=dp), dimension (:, :), allocatable :: crel     !! Non-relat. CMPLX. spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. CMPLX. spher. harm.
  complex (kind=dp), dimension (:, :), allocatable :: rrel     !! Non-relat. REAL spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. REAL spher. harm.
  complex (kind=dp), dimension (:, :, :), allocatable :: srrel
  complex (kind=dp), dimension (:, :, :), allocatable :: drotq !! Rotation matrices to change between LOCAL/GLOBAL frame of reference for magnetisation <> Oz or noncollinearity
  integer, dimension (:), allocatable :: zrel                  !! atomic number (cast integer)
  integer, dimension (:), allocatable :: jwsrel                !! index of the WS radius
  integer, dimension (:), allocatable :: irshift               !! shift of the REL radial mesh with respect no NREL
  integer, dimension (:, :), allocatable :: nrrel
  integer, dimension (:, :, :), allocatable :: irrel
  real (kind=dp), dimension (0:100) :: fact
  real (kind=dp), dimension (:, :), allocatable :: vtrel     !! potential (spherical part)
  real (kind=dp), dimension (:, :), allocatable :: btrel     !! magnetic field
  real (kind=dp), dimension (:, :), allocatable :: rmrel     !! radial mesh
  real (kind=dp), dimension (:, :), allocatable :: drdirel   !! derivative of radial mesh
  real (kind=dp), dimension (:, :), allocatable :: r2drdirel !! \( r^2 \frac{\partial}{\partial \mathbf{r}}\frac{\partial}{\partial i}\) (r**2 * drdi)
  real (kind=dp), dimension (:, :, :), allocatable :: thesme
  logical :: para = .true.
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
  !> @endnote
  ! -------------------------------------------------------------------------
  integer :: ntldau = 0    !! number of atoms on which LDA+U is applied
  integer :: idoldau = 0   !! flag to perform LDA+U
  integer :: itrunldau = 0 !! Iteration index for LDA+U
  integer :: kreadldau = 0 !! LDA+U arrays available
  integer, dimension (:), allocatable :: lopt            !! angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
  integer, dimension (:), allocatable :: itldau          !! integer pointer connecting the NTLDAU atoms to heir corresponding index in the unit cell
  real (kind=dp), dimension (:), allocatable :: ueff     !! input U parameter for each atom
  real (kind=dp), dimension (:), allocatable :: jeff     !! input J parameter for each atom
  real (kind=dp), dimension (:), allocatable :: erefldau !! the energies of the projector's wave functions (REAL)
  ! ..
  ! .. distinguish between spin-dependent and spin-independent
  ! .. quantities
  real (kind=dp), dimension (:, :, :, :), allocatable :: wldau    !! potential matrix
  real (kind=dp), dimension (:, :, :, :, :), allocatable :: uldau !! calculated Coulomb matrix elements (EREFLDAU)
  complex (kind=dp), dimension (:, :), allocatable :: phildau
  ! -------------------------------------------------------------------------
  ! LDA+U LDA+U LDA+U
  ! -------------------------------------------------------------------------
  ! Lloyds formula
  integer :: lly = 0              !! LLY <> 0 : apply Lloyds formula
  complex (kind=dp) :: deltae = (1.0e-5_dp, 0.0_dp) !! Energy difference for numerical derivative

  ! SUSC (BEGIN: modifications by Manuel and Benedikt)    ! susc
  ! LOGICAL THAT CHECKS WHETHER ENERGY MESH FILE EXISTS   ! susc
  logical :: emeshfile                                    ! susc
  ! SUSC (END:   modifications by Manuel and Benedikt)    ! susc

  ! allocations:
  real (kind=dp), dimension (:, :, :), allocatable :: thetas     !! shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
  real (kind=dp), dimension (:, :, :), allocatable :: thetasnew  !! shape function interpolated to Chebychev radial mesh


contains

  ! ----------------------------------------------------------------------------
  !> Summary: Main wrapper to handle input reading, allocation of arrays, and
  !> preparation of all the necessary data structures for a calculation.
  !> Author: Philipp Rüssmann, Bernd Zimmermann, Phivos Mavropoulos, R. Zeller,
  !> and many others ...
  !> Category: input-output, initialization, geometry, k-points, electrostatics, KKRhost
  !> Deprecated: False
  !> Main wrapper to handle input reading, allocation of arrays, and
  !> preparation of all the necessary data structures for a calculation. 
  ! ----------------------------------------------------------------------------
  subroutine main0() 

#ifdef CPP_OMPSTUFF
    use :: omp_lib ! necessary for omp functions
#endif
#ifdef CPP_MPI
    use :: mpi
#endif
    use :: mod_mympi, only: nranks
    use :: mod_runoptions, only: calc_DOS_Efermi, calc_GF_Efermi, relax_SpinAngle_Dirac, set_empty_system, use_Chebychev_solver, &
      use_decimation, use_ewald_2d, use_qdos, use_semicore, use_spherical_potential_only, use_virtual_atoms, write_deci_pot, &
      write_deci_tmat, write_energy_mesh, write_generalized_potential, write_green_host, write_green_imp, write_kkrimp_input, &
      write_kkrsusc_input, write_pkkr_input, write_pkkr_operators, write_potential_tests, write_rhoq_input, use_ldau, &
      disable_print_serialnumber, stop_1a, calc_wronskian
    use :: mod_version, only: version1, version2, version3, version4
    use :: mod_version_info, only: serialnr, version_print_header
    use :: mod_md5sums, only: get_md5sums, md5sum_potential, md5sum_shapefun
    use :: mod_wunfiles, only: wunfiles
    use :: mod_types, only: t_imp, t_inc, init_params_t_imp, init_t_imp
    use :: memoryhandling, only: memocc, allocate_cell, allocate_cpa, allocate_soc, allocate_ldau, allocate_magnetization, allocate_potential, &
      allocate_energies, allocate_relativistic, allocate_clusters, allocate_expansion, allocate_mesh, allocate_pannels, allocate_misc, &
      allocate_green, allocate_ldau_potential, allocate_rel_transformations, allocate_semi_inf_host
    use :: mod_create_newmesh, only: create_newmesh
    use :: mod_rhoqtools, only: rhoq_save_rmesh
    use :: rinput, only: rinput13
    use :: mod_addvirtual14, only: addviratoms14
    use :: mod_bzkint0, only: bzkint0
    use :: mod_calcrotmat, only: calcrotmat
    use :: mod_changerep, only: changerep
    use :: mod_cinit, only: cinit
    use :: mod_clsgen_tb, only: clsgen_tb
    use :: mod_deciopt, only: deciopt
    use :: mod_drvbastrans, only: drvbastrans
    use :: mod_epathtb, only: epathtb
    use :: mod_gaunt2, only: gaunt2
    use :: mod_gaunt, only: gaunt
    use :: mod_generalpot, only: generalpot
    use :: mod_getbr3, only: getbr3
    use :: mod_gfmask, only: gfmask
    use :: mod_lattix99, only: lattix99
    use :: mod_madelung2d, only: madelung2d
    use :: mod_madelung3d, only: madelung3d
    use :: mod_outpothost, only: outpothost
    use :: mod_outtmathost, only: outtmathost
    use :: mod_readimppot, only: readimppot
    use :: mod_relpotcvt, only: relpotcvt
    use :: mod_scalevec, only: scalevec
    use :: mod_setgijtab, only: setgijtab
    use :: mod_shape_corr, only: shape_corr
    use :: mod_startb1, only: startb1
    use :: mod_startldau, only: startldau
    use :: mod_testdim, only: testdim
    use :: mod_write_tbkkr_files, only: write_tbkkr_files
    use :: mod_writehoststructure, only: writehoststructure
    ! array dimensions
    use :: global_variables, only: krel, nspind, nrefd, irmd, ntotd, ipand, ncelld, nrmaxd, nchebd, natypd, naezd, lmaxd, alm, lmmaxd, &
      almgf0, lmgf0d, ndim_slabinv, nprincd, nembd, nembd1, nembd2, irmind, irnsd, nofgij, natomimpd, lpotd, lmpotd, npotd, nfund, &
      lmxspd, mmaxd, iemxd, ncleb, nclsd, nsheld, naclsd, lm2d, irid, lassld, nrd, nspind, nspindd, ngshd, linterface, nlayerd, knosph, &
      korbit, nmaxd, ishld, wlength, maxmshd, kpoibz, nspotd

    implicit none

    ! .. Local Scalars ..
    integer :: i
    integer :: j
    integer :: i1
    integer :: ie
    integer :: lm
    integer :: ns
    integer :: isvatom, nvatom
    integer :: i_stat, i_all
    integer :: irec
    integer :: lrecabmad
    integer :: i_commensurate !! counter to find closest divisor of naez for slab inversion (finds nprincd)
    integer :: ilayer !! loop counter for layer index, needed to find i_commensurate
    real (kind=dp) :: zattemp
    integer :: ierr
    real (kind=dp), dimension(:,:), allocatable :: tmp_rr
    ! for OPERATOR option
    logical :: lexist, operator_imp

#ifdef CPP_OMPSTUFF
    ! .. OMP ..
    integer :: nth, ith !! total number of threads, thread number
#endif
    ! ..
    ! .. External Functions ..



    ! -------------------------------------------------------------------------
    ! Write version info:
    ! -------------------------------------------------------------------------
    call print_versionserial(6, version1, version2, version3, version4, serialnr)
    call print_versionserial(1337, version1, version2, version3, version4, serialnr)

#ifdef CPP_OMPSTUFF
    !$omp parallel shared(nth) private(ith)
    ith = omp_get_thread_num()
    if (ith==0) then
      nth = omp_get_num_threads()
      write (*, '(/79("*")//1X,A,I5//79("*")/)') 'Number of OpenMP threads used:', nth
      write (1337, '(1X,A,I5)') 'Number of OpenMP threads used:', nth
    end if
    !$omp end parallel
#endif

#ifdef CPP_MPI
    write (*, '(1X,A,I5//79("*")/)') 'Number of MPI ranks used:', nranks
    write (1337, '(1X,A,I5//79("*")/)') 'Number of MPI ranks used:', nranks
#endif
    ! -------------------------------------------------------------------------
    ! End write version info
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! Reading of the inputcard, and allocation of several arrays
    !! @note JC: have added reading calls for the parameters that used to be in
    !! the `inc.p` and can now be modified via the `inputcard` directly
    !! @endnote
    ! -------------------------------------------------------------------------
    call rinput13(kte,igf,kxc,lly,icc,ins,kws,ipe,ipf,ipfe,icst,imix,lpot,naez,nemb,&
      nref,ncls,npol,lmax,kcor,kefg,khyp,kpre,kvmad,lmmax0d,lmpot,ncheb,nleft,ifile,  &
      kvrel,nspin,natyp,nineq,npnt1,npnt2,npnt3,kfrozn,ishift,n1semi,n2semi,n3semi, &
      scfsteps,insref,kshape,itdbry,nright,kforce,ivshift,khfeld,nlbasis,nrbasis,   &
      intervx,intervy,intervz,npan_eq,npan_log,npolsemi,tk,fcm,emin,emax,rmax,gmax, &
      alat,r_log,rcutz,rcutxy,eshift,qbound,hfield,mixing,abasis,bbasis,cbasis,     &
      vconst,tksemi,tolrdif,emusemi,ebotsemi,fsemicore,lambda_xc,deltae,lrhosym,    &
      linipol,lcartesian,imt,cls,lmxc,irns,irws,ntcell,refpot,inipol,ixipol,hostimp,&
      kfg,vbc,zperleft,zperight,bravais,rmt,zat,rws,mtfac,rmtref,rmtnew,rmtrefat,   &
      fpradius,tleft,tright,rbasis,socscale,cscl,socscl,solver,i12,i13,i19,i25,i40, &
      txc,drotq,ncpa,itcpamax,cpatol,noq,iqat,icpa,kaoez,conc,kmrot,qmtet,qmphi,    &
      kreadldau,lopt,ueff,jeff,erefldau,invmod,verbosity,MPI_scheme,                &
      special_straight_mixing)

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

    !--------------------------------------------------------------------------------
    ! Allocation calls
    !--------------------------------------------------------------------------------
    !! @note Jonathan Chico: The main idea here is to allocate all the needed arrays so that the `inc.p`
    !! file becomes irrelevant. In principle the philosophy would be to modularize
    !! the code such that each module has its own global variables and allocation routine
    !! e.g. a module called CPA_control could have defined all the needed CPA variables
    !! as well as the allocation calls, this module would be used in the needed routines
    !! and the arrays would only be allocated if a CPA calculation is actually performed
    !! in the current way ALL arrays are allocated which could cause an unnecessary memory
    !! consumption
    !! @endnote
    !--------------------------------------------------------------------------------
    ! Call to allocate the arrays associated with the potential
    call allocate_potential(1, irmd, natypd, npotd, ipand, nfund, lmxspd, lmpotd, irmind, nspotd, nfu, irc, ncore, irmin, lmsp, lmsp1, ircut, lcore, llmsp, ititle, visp, &
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

    !--------------------------------------------------------------------------------
    ! End of allocation calls
    !--------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------
    ! Deal with the lattice
    !--------------------------------------------------------------------------------

    call lattix99(linterface,alat,natyp,naez,conc,rws,bravais,recbv,volume0,rr,nrd, &
     natyp)

    ! rr has changed, fix allocation of array to new nrd size
    allocate(tmp_rr(3,0:nrd),stat=i_stat)
    call memocc(i_stat, product(shape(tmp_rr))*kind(tmp_rr), 'tmp_rr', 'main0')

    tmp_rr(:, :) = rr(1:3, 0:nrd)

    i_all = -product(shape(rr))*kind(rr)
    deallocate (rr, stat=i_stat)
    call memocc(i_stat, i_all, 'rr', 'main0')

    allocate(rr(3,0:nrd),stat=i_stat)
    call memocc(i_stat, product(shape(rr))*kind(rr), 'rr', 'main0')

    rr(:, :) = tmp_rr(:, :)
    
    i_all = -product(shape(tmp_rr))*kind(tmp_rr)
    deallocate (tmp_rr, stat=i_stat)
    call memocc(i_stat, i_all, 'tmp_rr', 'main0')

    call scalevec(lcartesian,rbasis,abasis,bbasis,cbasis,nlbasis,nrbasis,nleft,     &
      nright,zperleft,zperight,tleft,tright,linterface,naez,nemb,bravais,kaoez,noq, &
      naez,natyp,nemb)
    ! After SCALEVEC all basis positions are in cartesian coords.

    nvirt = 0
    if (use_virtual_atoms) then
      write (1337, *) 'Calling ADDVIRATOMS'
      call addviratoms14(linterface,nvirt,naez,naez,natyp,nemb,nemb,rbasis,.true.,  &
        bravais,ncls,nineq,refpot,kaoez,noq,nref,rmtrefat,i25)
    end if

    call clsgen_tb(naez,nemb,nvirt,rr,rbasis,kaoez,zat,cls,ncls,nacls,atom,ezoa,    &
      nlbasis,nrbasis,nleft,nright,zperleft,zperight,tleft,tright,rmtref,rmtrefat,  &
      vref,refpot,nref,rcls,rcutz,rcutxy,alat,natyp,nclsd,nrd,naclsd,nrefd,nembd,   &
      linterface,nprincd,nprinc)

    ! change nrefd to nref and reduce size of rmtre, vref accordingly
    ! do the same for ncls(d) with nacls and rcls arrays
    ! this is needed because in the 'allocate_clusters' call the maximal size was used
    call reduce_array_size(nref, nrefd, rmtref, vref, ncls, nclsd, nacls, rcls)

    nlayer = naez/nprinc
    ! overwrite nprincd if chosen too small (also updates array `icheck`)
    if (nprincd<nprinc) then
      ! find nprincd such that it is as big as it needs to be while being as
      ! small as commensurability with the number of layers etc. allows
      ! for this we loop over all layers and look for the divisors of naez
      i_commensurate = -1
      do ilayer = naez, 1, -1 ! go through loop backwards to find smallest divisor
        if (mod(naez, ilayer)==0) then
          if (naez/ilayer>=nprinc .and. i_commensurate==-1) then
            i_commensurate = naez/ilayer
          end if
        end if
      end do
      ! now we take the smallest divisor of naez that is >= nprinc
      if (i_commensurate>-1) nprinc = i_commensurate

      ! this is the fallback to reset it to the number of atoms
      if (nlayer*nprinc/=naez) then
        ! in this case we should actually do full inversion instead of slab inversion
        nprinc = naez
        ! this is enforced here automatically
        write (*, '(A)') 'WARNING: Found NPRINC==NAEZ!', 'Automatically overwriting inversion scheme with full inversion'
        write (1337, '(A)') 'WARNING: Found NPRINC==NAEZ!', 'Automatically overwriting inversion scheme with full inversion'
        invmod = 0
      end if

      ! now nprinc was found successfully, so nprincd can be set accordingly
      write (*, *) 'Automatically overwriting nprincd with ', nprinc
      write (1337, *) 'Automatically overwriting nprincd with ', nprinc
      nprincd = nprinc

      ! update parameter that depend on nprincd and change allocations of arrays that have nprincd
      ndim_slabinv = nprincd*lmmaxd

      i_all = -product(shape(icheck))*kind(icheck)
      deallocate (icheck, stat=i_stat)
      call memocc(i_stat, i_all, 'icheck', 'main0')
     
      allocate (icheck(naez/nprincd,naez/nprincd), stat=i_stat)
      call memocc(i_stat, product(shape(icheck))*kind(icheck), 'ICHECK', 'main0')
      icheck = 0
    end if ! nprincd<nprinc

    ! store nlayerd for later use
    nlayerd = nlayer

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

    call testdim(nspin,naez,nemb,natyp,ins,insref,nref,irns,nlayer,krel,nspind,     &
      nprincd,knosph,irnsd,korbit,invmod)

    if (ins>0) open (19, file=i19, status='old', form='formatted')
    if (ifile==13) open (ifile, file=i13, status='old', form='formatted')
    if (icc>0) open (25, file=i25, status='unknown', form='formatted')

    call startb1(ifile,1337,1337,ipe,krel,kws,lmax,1,natyp,alatnew,rmtnew,rmt,      &
      ititle,imt,irc,vconst,ins,irns,fpradius,nspin,vins,irmin,kshape,ntcell,ircut, &
      ipan,thetas,ifunm,nfu,llmsp,lmsp,e2in,vbc,dror,rs,s,visp,rws,ecore,lcore,     &
      ncore,drdi,rmesh,zat,a,b,irws,1,lmpot,irmind,irmd,lmxspd,ipand,irid,irnsd,    &
      natyp,ncelld,nfund,nspotd,ivshift, npotd)

    ! find md5sums for potential and shapefunction
    call get_md5sums(ins, i13, i19)
    write (1337, '(A,A)') 'Doing calculation with potential: ', md5sum_potential
    if (ins>0) then
      write (1337, '(A,A)') 'Doing calculation with shapefun: ', md5sum_shapefun
    end if

    if (write_rhoq_input) then
      call rhoq_save_rmesh(natyp,irm,ipand,irmin,irws,ipan,rmesh,ntcell,ircut,r_log,&
        npan_log,npan_eq)

      ! check consistency
      if (write_green_imp .or. write_green_host) then
        write (*, *) 'warning! rhoqtest cannot be used together with '
        write (*, *) '''GREENIMP'' or ''WRTGREEN'' options'
        stop
      end if

      ! ! enforce MPIatom here
      ! if (MPI_scheme!=1) stop 'rhoqtest assumes MPIatom'
      ! CALL ADDTEST('MPIatom')

      if (nranks>1) then
        write (*, *) 'at the moment rhoqtest does not work with MPI.'
        write (*, *) 'compile hybrid version and use OMP level only.'
        stop
      end if

    end if                         ! write_rhoq_input

    if (use_spherical_potential_only) then
      write (1337, *) 'TEST OPTION Vspher,', 'keeping only spherical component of potential.'
      vins(irmind:irmd, 2:lmpot, 1:nspotd) = 0.d0
    end if

    if (set_empty_system) then
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
    if (calc_GF_Efermi .or. calc_DOS_Efermi) then
      emin = e2in
      if (calc_GF_Efermi) then
        write (1337, fmt=130)
      else
        write (1337, fmt=140)
      end if
    end if

    if (abs(e2in-emax)>1d-10 .and. npol/=0) emax = e2in
    ! -------------------------------------------------------------------------
    if (write_generalized_potential) then
      rewind (3)
      call generalpot(3,1,natyp,nspin,zat,alat,rmt,rmtnew,rws,rmesh,drdi,visp,irws, &
        a,b,ins,irns,lpot,vins,qbound,irc,kshape,e2in,vbc,ecore,lcore,ncore,lmpot,  &
        irmd,irmind)
      close (3)
    end if
    ! -------------------------------------------------------------------------
    ! --> Apply external magnetic field
    ! -------------------------------------------------------------------------
    ! from startb1 moved here
    if (khfeld==1) then
      ! ---> maybe apply a magnetic field
      call bshift_ns(irmd,irid,ipand,lmpot,npotd,natyp,nspin,ngshd,nfund,ncelld,    &
        irmind,lmxspd,kshape,irc,irmin,inipol,ntcell,imaxsh,ilm_map,lmsp,ifunm,     &
        ircut,hfield,gsh,rmesh,thesme,thetas,visp,vins)
    end if
    if (write_potential_tests) then     ! ruess
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
    !if (krel+korbit==1) then
    if (krel==1 .or. use_Chebychev_solver) then
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
        call relpotcvt(1,visp,zat,rmesh,drdi,ircut,vtrel,btrel,zrel,rmrel,jwsrel,   &
          drdirel,r2drdirel,irshift,ipand,irmd,npotd,natyp)
      end if
    !end if ! KREL+KORBIT.EQ.1
    end if ! KREL==1 .or. opt('NEWSOSOL)
    ! -------------------------------------------------------------------------
    ! set up energy contour
    !--------------------------------------------------------------------------------
    idosemicore = 0
    if (use_semicore) idosemicore = 1

    call epathtb(ez,dez,e2in,ielast,iesemicore,idosemicore,emin,emax,tk,npol,npnt1, &
      npnt2,npnt3,ebotsemi,emusemi,tksemi,npolsemi,n1semi,n2semi,n3semi,iemxd)

    !--------------------------------------------------------------------------------
    ! SUSC (BEGIN: modifications by Manuel and Benedikt)                              
    !--------------------------------------------------------------------------------
    !                                                                                 ! susc
    if (write_energy_mesh) then                                                         ! susc
      ! write out the energy mesh and the corresponding                               ! susc
      ! weights to a file called 'emesh.scf'                                          ! susc
      write (*, '("main0: Runflag emesh is set.")')                                   ! susc
      write (*, '("       File emesh.scf will be written!")')                         ! susc
      write (*, *) 'writing emesh.scf file...'                                        ! susc
      open (file='emesh.scf', unit=12111984, status='replace')                        ! susc
      write (12111984, '(5x,i0)') ielast                                              ! susc
      do ie = 1, ielast                                                               ! susc
        write (12111984, '(4es16.8)') ez(ie), dez(ie)                                 ! susc
      end do                                                                          ! susc
      close (12111984)                                                                ! susc
      write (*, '("       Finished writing file emesh.scf.")')                        ! susc
    end if                                                                            ! susc
    !                                                                                 ! susc
    !                                                                                 ! susc
    if (write_kkrsusc_input) then                                                         ! susc
      ! read in 'emesh.dat' with new energy mesh-points                               ! susc
      inquire (file='emesh.dat', exist=emeshfile)                                     ! susc
      write (*, '("main0: Runflag KKRSUSC is set.")')                                 ! susc
      if (emeshfile) then                                                             ! susc
        write (*, '("main0: File emesh.dat exists and will ")', advance='no')         ! susc
        write (*, '("be read in.")')                                                  ! susc
        write (*, '("       Energy contour from inputcard ")', advance='no')          ! susc
        write (*, '("will be overwritten!")')                                         ! susc
        open (file='emesh.dat', unit=50)                                              ! susc
        read (50, *) ielast                                                           ! susc
        if (ielast>iemxd) stop 'ielast > iemxd!'                                      ! susc
        do ie = 1, ielast                                                             ! susc
          read (50, '(4es16.8)') ez(ie), dez(ie)                                      ! susc
          write (*, '(i8,4es16.8)') ie, ez(ie), dez(ie)                               ! susc
        end do                                                                        ! susc
        close (50)                                                                    ! susc
        write (*, '("       Finished reading in file emesh.dat.")')                   ! susc
      else                                                                            ! susc
        stop 'main0: Runflag KKRSUSC but cannot find file emesh.dat!'                 ! susc
      end if                                                                          ! susc
    end if                                                                            ! susc
    !                                                                                 ! susc
    !                                                                                 ! susc
    ! still missing: check here whether scfsteps is > 1                               ! susc
    ! if scfsteps>1 --> option a) stop program here                                   ! susc
    ! option b) set it to 1 and continue                                              ! susc
    if (write_kkrsusc_input .and. scfsteps>1) then                                        ! susc
      write (*, '("main0: Runflag KKRSUSC is set ")')                                 ! susc
      write (*, '("but scfsteps = ",i0)') scfsteps                                    ! susc
      write (*, '("       Here we enforce scfsteps = 1")')                            ! susc
      !------------------------------------------------------------------------------ ! susc
      scfsteps = 1                                                                    ! susc
      !------------------------------------------------------------------------------ ! susc
    end if                                                                            ! susc
    !--------------------------------------------------------------------------------
    ! SUSC (END:   modifications by Manuel and Benedikt)                              
    !--------------------------------------------------------------------------------

    do ie = 1, ielast
      wez(ie) = -2.0_dp/pi*dez(ie)
      if (ie<=iesemicore) wez(ie) = wez(ie)*fsemicore
    end do
    ! -------------------------------------------------------------------------
    ! update energy contour for Fermi-surface generation                       ! fswrt=fermi-surface write
    ! -------------------------------------------------------------------------
    if (write_pkkr_input) then      ! fswrt
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
    call gaunt(lmax,lpot,wg,yrg,cleb,loflm,icleb,iend,jend,ncleb,lmax,lmgf0d,lmpot)

    ! -------------------------------------------------------------------------
    ! set up of GAUNT coefficients C(l,m;l',m';l'',m'') for all
    ! nonvanishing (l'',m'')-components of the shape functions THETAS
    ! -------------------------------------------------------------------------
    if (kshape/=0) then
      call shape_corr(lpot,natyp,gsh,ilm_map,imaxsh,lmsp,ntcell,wg,yrg,lassld,lmpot,&
        natyp,ngshd)
    end if
    ! -------------------------------------------------------------------------
    ! calculate Madelung constants (needed only for SCF calculations)
    ! -------------------------------------------------------------------------
    ! fivos      IF ( SCFSTEPS.GT.1 .OR. ICC.GT.0 ) THEN
    if (npol/=0 .or. use_decimation) then ! No madelung calculation in case of DOS., needed for demination nevertheless
      ! OPEN(99,FILE='madelinfo.txt')

      !> @note Use option 'ewald2d' if the madelung summation is to be carried out in
      !> single-slab mode, otherwise it is carried out in repeated (periodic)
      !> slab mode.
      !> Reason: the 2d-mode gives wrong results sometimes [e.g. in diamond
      !> structure (110)].
      !> @endnote
      if (linterface .and. (use_ewald_2d .or. use_decimation)) then ! ewald2d
        write (*, *) 'Calling MADELUNG2D'
        ! -------------------------------------------------------------------
        ! 2D case
        ! -------------------------------------------------------------------
        call madelung2d(lpot,yrg,wg,naez,alat,volume0,bravais,recbv,rbasis,rmax,    &
          gmax,nlbasis,nleft,zperleft,tleft,nrbasis,nright,zperight,tright,lmxspd,  &
          lassld,lpot,lmpot,nmaxd,ishld,nembd1,wlength)
        write (*, *) 'Exited MADELUNG2D'
      else
        ! -------------------------------------------------------------------
        ! 3D case
        ! -------------------------------------------------------------------
        if (linterface) then
          call getbr3(nembd1,nlbasis,alat,tleft,nrbasis,tright,bravais,recbv,volume0)
        end if

        write (*, *) 'Calling MADELUNG3D'
        call madelung3d(lpot,yrg,wg,naez,alat,volume0,bravais,recbv,rbasis,rmax,    &
          gmax,naez,lmxspd,lassld,lpot,lmpot,nmaxd,ishld,nemb,wlength)
        write (*, *) 'Exited MADELUNG3D'
      end if

      ! CLOSE(99)
    else ! NPOL==0
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

    end if ! npol==0
    ! -------------------------------------------------------------------------
    ! fivos      END IF
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! Set up I,J pairs for ICC = -1
    ! -------------------------------------------------------------------------
    if (icc<0) then
      call setgijtab(linterface,icc,naez,iqat,rbasis,bravais,natomimp,atomimp,      &
        rclsimp,ijtabcalc,iofgij,jofgij,nqcalc,iqcalc,natomimpd,ijtabcalc_i)
    end if

    ! -------------------------------------------------------------------------

    dsymll = (0d0, 0d0)
    dsymll1 = (0d0, 0d0)

    call bzkint0(nshell,naez,natyp,noq,rbasis,kaoez,icc,bravais,recbv,atomimp,      &
      rsymat,isymindex,nsymat,i25,natomimp,nsh1,nsh2,rclsimp,ratom,ijtabsym,ijtabsh,&
      ijtabcalc,iofgij,jofgij,nofgij,ish,jsh,rrot,dsymll1,para,qmtet,qmphi,         &
      symunitary,hostimp,intervx,intervy,intervz,ielast,ez,kmesh,maxmesh,maxmshd,   &
      krel+korbit,lmax,lmmaxd,kpoibz,naez,natyp,natomimpd,nsheld,nemb,iemxd)

    ! -------------------------------------------------------------------------

    if (write_kkrimp_input) then

      call writehoststructure(bravais, natyp, rbasis, naez, nemb)

      open (58, file='kkrflex_atominfo', form='FORMATTED')
      call version_print_header(58, addition='; '//md5sum_potential//'; '//md5sum_shapefun, disable_print=disable_print_serialnumber)
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
    if (icc/=0 .and. .not. write_kkrimp_input) then
      open (58, file='shells.dat')
      write (1337, *) 'Writing out shells (also in shells.dat):'                 ! fivos
      write (1337, *) 'itype,jtype,iat,jat,r(iat),r(jat)'                        ! fivos
      write (1337, *) nshell(0), 'NSHELL(0)'                                     ! fivos
      write (58, *) nshell(0), 'NSHELL(0)'                                       ! fivos
      do i1 = 1, nshell(0)                                                       ! fivos
        write (1337, *) i1, nshell(i1), 'No. of shell, No. of atoms in shell'    ! fivos
        write (58, *) i1, nshell(i1), 'No. of shell, No. of atoms in shell'      ! fivos
        do lm = 1, nshell(i1)                                                    ! fivos
          write (1337, *) 'ish(i1,lm)', ish(i1, lm)
          if (ish(i1,lm)>0 .and. jsh(i1,lm)>0) then                              ! fix bernd
            write (1337, 100) nsh1(i1), nsh2(i1), ish(i1, lm), jsh(i1, lm) &     ! fivos
              , (rclsimp(i,ish(i1,lm)), i=1, 3), (rclsimp(i,jsh(i1,lm)), i=1, 3) ! fivos
            write (58, 100) nsh1(i1), nsh2(i1), ish(i1, lm), jsh(i1, lm), &      ! fivos
              (rclsimp(i,ish(i1,lm)), i=1, 3), (rclsimp(i,jsh(i1,lm)), i=1, 3)   ! fivos
          else                                                                   ! fix bernd
            write (1337, 110) nsh1(i1), nsh2(i1), ish(i1, lm), jsh(i1, lm)       ! fix bernd
            write (58, 110) nsh1(i1), nsh2(i1), ish(i1, lm), jsh(i1, lm)         ! fix bernd
          end if                                                                 ! fix bernd
100       format (4i5, 6f16.6)                                                   ! fivos
110       format (4i5)                                                           ! fix bernd
        end do                                                                   ! fivos
      end do                                                                     ! fivos
      write (1337, *) '###################'
      close (58)
    end if
    ! -------------------------------------------------------------------------
    ! end fivos
    ! -------------------------------------------------------------------------

    call gfmask(icheck,icc,invmod,nsh1,nsh2,naez,nshell(0),naez,nprincd)
    ! -------------------------------------------------------------------------
    ! set up transformation matrices between REL/NREL representations
    ! -------------------------------------------------------------------------
    if ((krel+korbit)==1) then
      call drvbastrans(rc,crel,rrel,srrel,nrrel,irrel,lmax+1,lmmaxd,2*(lmax+1),     &
        lmmaxd+2*(lmax+1),mmaxd,2*(lmax+1)*mmaxd)
    end if
    if (korbit==1) then
      do ns = 1, nsymat
        call changerep(dsymll1(1,1,ns),'REL>RLM',dsymll(1,1,ns),lmmaxd,lmmaxd,rc,   &
          crel,rrel,'DSYMLL',0)
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
        fact(i) = fact(i-1)*real(i, kind=dp)
      end do

      do i1 = 1, naez
        call calcrotmat(mmaxd,(krel+korbit)*3,qmphi(i1),qmtet(i1),0.0_dp,           &
          drotq(1,1,i1),fact,lmmaxd)
      end do
    end if
    ! -------------------------------------------------------------------------
    ! Treat decimation I/O cases
    ! -------------------------------------------------------------------------
    if (write_deci_pot) then
      call outpothost(alat,ins,krel+korbit,kmrot,nspin,naez,natyp,e2in,bravais,     &
        rbasis,qmtet,qmphi,noq,kaoez,iqat,zat,conc,ipan,ircut,solver,socscl,cscl,   &
        irws,rmtnew,rws,rmesh,drdi,visp,irshift,rmrel,drdirel,vtrel,btrel,lmax,     &
        natyp,naez,ipand,irmd)
    end if
    if (write_deci_tmat) then
      if (nranks>1) stop 'ERROR: deci-out does not work with MPI!'
      call outtmathost(alat,ins,krel+korbit,kmrot,nspin,naez,lmmax0d,bravais,rbasis,  &
        qmtet,qmphi,e2in,tk,npol,npnt1,npnt2,npnt3)
    end if
    if (use_decimation) then
      call deciopt(alat,ins,krel+korbit,kvrel,kmrot,nspin,naez,lmmax0d,bravais,tk,    &
        npol,npnt1,npnt2,npnt3,ez,ielast,kaoez,lefttinvll,righttinvll,vacflag,      &
        nlbasis,nrbasis,cmomhost,vref,rmtref,nref,refpot(naez),lmax,lmgf0d,lmmaxd,  &
        lm2d,nembd1,iemxd,nspindd,lmpot,natyp,irmd,ipand)
    end if
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! ITERMDIR  -- initialise
    ! -------------------------------------------------------------------------
    if (relax_SpinAngle_Dirac) then
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

    if (use_ldau) then
      call startldau(itrunldau,idoldau,kreadldau,lopt,ueff,jeff,erefldau,natyp,     &
        nspin,wldau,uldau,phildau,irws,ntldau,itldau,irmd,natyp,nspind,mmaxd)
    end if
    ! -------------------------------------------------------------------------
    ! LDA+U
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! write out information for the other program parts           =
    ! -------------------------------------------------------------------------

    ! new solver for full-potential, spin-orbit, initialise
    if (use_Chebychev_solver) then
      call create_newmesh(natyp,irmd,ipand,irid,ntotd,nfund,ncheb,ntotd*(ncheb+1),  &
        nspin,rmesh,irmin,ipan,ircut,r_log,npan_log,npan_eq,npan_log_at,npan_eq_at, &
        npan_tot,rnew,rpan_intervall,ipan_intervall,ncelld,ntcell,thetas,thetasnew)
    end if

    call wunfiles(npol, npnt1, npnt2, npnt3, ielast, tk, emin, emax, ez, wez, efermi, npolsemi, n1semi, n2semi, n3semi, iesemicore, tksemi, ebotsemi, emusemi, fsemicore, vins, &
      visp, vbc, vtrel, btrel, rmrel, drdirel, r2drdirel, zrel, jwsrel, irshift, itscf, scfsteps, cmomhost, ecore, lcore, ncore, qmtet, qmphi, qmphitab, qmtettab, qmgamtab, drotq, &
      nsra, ins, natypd, naezd, nineq, nref, nspin, ncls, icst, ipan, ircut, alat, zat, rmesh, drdi, refpot, rmtref, vref, iend, jend, cleb, icleb, atom, cls, rcls, nacls, loflm, &
      solver, socscl, cscl, icc, igf, nlbasis, nrbasis, ncpa, icpa, itcpamax, cpatol, rbasis, rr, ezoa, nshell, nsh1, nsh2, ijtabcalc, ijtabcalc_i, ish, jsh, ijtabsym, ijtabsh, &
      nofgij, nqcalc, iqcalc, kmrot, kaoez, iqat, noq, conc, kmesh, maxmesh, nsymat, symunitary, rrot, dsymll, invmod, icheck, natomimp, ratom, atomimp, rc, crel, rrel, srrel, &
      nrrel, irrel, lefttinvll, righttinvll, vacflag, a, b, ifunm, ifunm1, intervx, intervy, intervz, ititle, lmsp1, ntcell, thetas, lpotd, lmpotd, nright, nleft, linterface, imix, &
      mixing, qbound, fcm, itdbry, irns, kpre, kshape, kte, kvmad, kxc, lambda_xc, txc, ishift, ixipol, lrhosym, kforce, lmsp, llmsp, rmt, rmtnew, rws, imt, irc, irmin, irws, nfu, &
      hostimp, gsh, ilm_map, imaxsh, idoldau, itrunldau, ntldau, lopt, itldau, ueff, jeff, erefldau, uldau, wldau, phildau, iemxd, irmind, irmd, nspotd, npotd, nembd1, lmmaxd, &
      ipand, nembd2, lmax, ncleb, naclsd, nclsd, lm2d, lmax+1, mmaxd, nrd, nsheld, naez/nprincd, natomimpd, nspind, irid, nfund, ncelld, lmxspd, ngshd, krel, ntotd, ncheb, &
      npan_log, npan_eq, npan_log_at, npan_eq_at, r_log, npan_tot, rnew, rpan_intervall, ipan_intervall, nspindd, thetasnew, socscale, tolrdif, lly, deltae, rclsimp, verbosity, MPI_scheme, &
      special_straight_mixing )

    if (write_pkkr_input) then                                                          ! fswrt
      call write_tbkkr_files(lmax, nemb, ncls, natyp, naez, ielast, ins, alat, &        ! fswrt
        bravais, recbv, rbasis, cls, nacls, rcls, ezoa, atom, rr, nspin, nrd, korbit, & ! fswrt
        nclsd, naclsd)                                                                  ! fswrt
    end if                                                                              ! fswrt

    if (write_pkkr_operators) then
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

    if (write_green_imp .or. operator_imp) then                                       ! GREENIMP
      ! fill array dimensions and allocate arrays in t_imp                            ! GREENIMP
      call init_params_t_imp(t_imp,ipand,natyp,irmd,irid,nfund,nspin,irmind,lmpot)    ! GREENIMP
      call init_t_imp(t_inc, t_imp)                                                   ! GREENIMP
                                                                                      ! GREENIMP
      ! next read impurity potential and shapefunction                                ! GREENIMP
      call readimppot(natomimp,ins,1337,0,0,2,nspin,lpot,t_imp%ipanimp,             & ! GREENIMP
        t_imp%thetasimp,t_imp%ircutimp,t_imp%irwsimp,khfeld,hfield,t_imp%vinsimp,   & ! GREENIMP
        t_imp%vispimp,t_imp%irminimp,t_imp%rimp,t_imp%zimp,irmd,irnsd,irid,nfund,   & ! GREENIMP
        ipand)                                                                        ! GREENIMP
    end if                                                                            ! GREENIMP


    if (ishift==2) then                                                               ! fxf
      open (67, file='vmtzero', form='formatted')                                     ! fxf
      write (67, 120) vbc(1)                                                          ! fxf
      close (67)                                                                      ! fxf
120   format (d20.12)                                                                 ! fxf
    end if                                                                            ! fxf

    ! Check for inputcard consistency in case of qdos option
    if (use_qdos) then
      write (1337, *)
      write (1337, *) '     < QDOS > : consistency check '
      if ((npol/=0) .and. (npnt1==0) .and. (npnt3==0)) then
        stop 'For qdos calculation change enery contour to dos path'
      end if
      if (tk>50.d0) write (*, *) 'WARNING:  high energy smearing due to high value of TEMPR for energy contour integration could not be of advantage. Consider changeing ''TEMPR'' to lower value'
      if (tk>50.d0) write (1337, *) 'WARNING:  high energy smearing due to high value of TEMPR for energy contour integration could not be of advantage. Consider changeing ''TEMPR'' to lower value'
      write (1337, *) '       QDOS: consistecy check complete'
    end if

    ! Check consistency with Wronskian test calculation
    if (calc_wronskian) then
      write (1337, *)
      write (1337, *) '     < WRONSKIAN > : consistency check '
      write (1337, *) ' run wronskian calculation with single energy point only and then execute ''check_wronskian.py'' script.'
      if (.not. stop_1a) then
        stop_1a = .true.
        write (*, *) ' automatically adding ''stop_1a'' option.'
        write (1337, *) ' automatically adding ''stop_1a'' option.'
      end if
      if (npol/=0 .and. npnt1==0 .and. npnt3==0 .and. npnt2/=1) then
        write(*, *) 'Calculation of Wronskian only possible for a single energy point!'
        write(*, *) 'Otherwise files become too large.'
        stop 
      end if
      if (tk>10.0e-5_dp) then
        stop 'Calculation of Wronskian only works for real energies! Choose ''TEMPR=0'''
      end if
      if (nranks>1) then
        stop 'Calculation of Wronskian only works in serial!'
      end if
      write (1337, *) '       WRONSKIAN: consistecy check complete'
    end if

    ! -------------------------------------------------------------------------

    write (1337, '(79("="),/,31X,"< KKR0 finished >",/,79("="),/)')
130 format (5x, 'INFO:  Output of cluster Green function at E Fermi')
140 format (5x, 'INFO:  Determination of DOS at E Fermi')

  end subroutine main0

  !-------------------------------------------------------------------------------
  !> Summary: Adds a constant (=VSHIFT) to the potentials of atoms
  !> Author:
  !> Category: potential, KKRhost
  !> Deprecated: False 
  !> Adds a constant (=VSHIFT) to the potentials of atoms
  !-------------------------------------------------------------------------------
  subroutine bshift_ns(irm,irid,ipand,lmpot,npotd,natyp,nspin,ngshd,nfund,ncelld,   &
    irmind,lmxspd,kshape,irc,irmin,inipol,ntcell,imaxsh,ilm_map,lmsp, ifunm, ircut, &
    hfield, gsh, rmesh, thesme, thetas, visp, vins)

    use :: mod_datatypes, only: dp
    use :: global_variables, only: nspotd
    use :: mod_convol, only: convol
    use :: mod_rinit, only: rinit
    use :: mod_constants, only: pi
    implicit none

    ! .. Input variables
    integer, intent (in) :: irm     !! Maximum number of radial points
    integer, intent (in) :: irid    !! Shape functions parameters in non-spherical part
    integer, intent (in) :: ipand   !! Number of panels in non-spherical part
    integer, intent (in) :: lmpot   !! (LPOT+1)**2
    integer, intent (in) :: npotd   !! (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP)
    integer, intent (in) :: natyp   !! Number of kinds of atoms in unit cell
    integer, intent (in) :: nspin   !! Counter for spin directions
    integer, intent (in) :: ngshd   !! Shape functions parameters in non-spherical part
    integer, intent (in) :: nfund   !! Shape functions parameters in non-spherical part
    integer, intent (in) :: ncelld  !! Number of cells (shapes) in non-spherical part
    integer, intent (in) :: irmind  !! irmd - irnsd
    integer, intent (in) :: lmxspd
    integer, intent (in) :: kshape  !! exact treatment of WS cell
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

    ! .. In/Out variables
    real (kind=dp), dimension (irm, npotd), intent (inout) :: visp !! Spherical part of the potential
    real (kind=dp), dimension (irmind:irm, lmpot, nspotd), intent (inout) :: vins !! Non-spherical part of the potential

    ! .. Local variables
    integer :: ispin, ih, ipot, ir, lm, imt1, irc1, irmin1
    real (kind=dp) :: vshift
    real (kind=dp), dimension (irm) :: pshiftr
    real (kind=dp), dimension (irm, lmpot) :: pshiftlmr

    real (kind=dp), parameter :: rfpi = sqrt(4.0_dp*pi)

    do ih = 1, natyp

      imt1 = ircut(1, ih)
      irc1 = irc(ih)
      irmin1 = irmin(ih)

      do ispin = 1, nspin
        ! shift potential spin dependent
        vshift = -real(2*ispin-3, kind=dp)*hfield*inipol(ih)

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
          call convol(imt1,irc1,ntcell(ih),imaxsh(lmpot),ilm_map,ifunm,lmpot,gsh,   &
            thetas, thesme, 0.0_dp, rfpi, rmesh(1,ih), pshiftlmr, pshiftr, lmsp)
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
  !> Summary: Print the version info and header to the output file
  !> Author: Philipp Ruessmann 
  !> Category: input-output, KKRhost 
  !> Deprecated: False 
  !> Print the version info and header to the output file
  !-------------------------------------------------------------------------------  
  subroutine print_versionserial(iunit,version1,version2,version3,version4,serialnr)
    implicit none
    integer, intent (in) :: iunit !! Unit identifier for the output file
    character (len=*), intent (in) :: version1  !! Version of the code
    character (len=*), intent (in) :: version2  !! Compilation option
    character (len=*), intent (in) :: version3  !! Compilation option
    character (len=*), intent (in) :: version4  !! Compilation option
    character (len=*), intent (in) :: serialnr  !! File serial number

    write (iunit, '(1A)') '     Screened Korringa-Kohn-Rostoker Electronic Structure Code'
    write (iunit, '(1A)') '                      for Bulk and Interfaces'
    write (iunit, '(1A)') '                    Juelich-Munich 2001 - 2021'
    write (iunit, '(1A)') ''
    write (iunit, '(2A)') '  Code version: ', trim(version1)
    write (iunit, '(6A)') '  Compile options: ', trim(version2), ' ', trim(version3), ' ', trim(version4)
    write (iunit, '(2A)') '  serial number for files: ', serialnr
  end subroutine print_versionserial


  !-------------------------------------------------------------------------------  
  !> Summary: Reduce size of arrays depending on nrefd, ncsld 
  !> Author: Philipp Ruessmann
  !> Category: input-output, KKRhost 
  !> Deprecated: False 
  !> Overwrite nrefd and nclsd with actual values and change allocations accordingly
  !> Should be called after nrefd, nclsd have been determined in clsgen_tb
  !-------------------------------------------------------------------------------  
  subroutine reduce_array_size(nref, nrefd, rmtref, vref, ncls, nclsd, nacls, rcls)

    use mod_datatypes, only: dp

    implicit none

    ! interface
    integer, intent(in) :: nref      !! actual number of reference potentials
    integer, intent(in) :: ncls      !! actual number of clusters
    integer, intent(inout) :: nrefd  !! maximal number of reference potentials
    integer, intent(inout) :: nclsd  !! maximal number of clusters
    integer, dimension(:), allocatable, intent(inout) :: nacls !! number of atom in cluster
    real (kind=dp), dimension(:, :, :), allocatable, intent(inout) :: rcls !!real space position of atoms in cluster
    real (kind=dp), dimension (:), allocatable, intent(inout) :: rmtref !! Muffin-tin radius of reference system
    real (kind=dp), dimension (:), allocatable, intent(inout) :: vref !! reference potential
    ! local
    integer :: naclsd !! size of second dimension of rcls
    integer :: i_stat !! status of (de)allocations
    real (kind=dp), dimension (:), allocatable :: rmtref_temp    !! Muffin-tin radius of reference system
    real (kind=dp), dimension (:), allocatable :: vref_temp      !! reference potential
    integer, dimension(:), allocatable :: nacls_temp             !! number of atom in cluster
    real (kind=dp), dimension(:, :, :), allocatable :: rcls_temp !!real space position of atoms in cluster

    ! determine size of second rcls dimension
    naclsd = ubound(rcls, 2)

    ! Allocate temporary arrays
    allocate(rmtref_temp(nrefd),stat=i_stat)
    allocate(vref_temp(nrefd),stat=i_stat)
    allocate(nacls_temp(nclsd),stat=i_stat)
    allocate(rcls_temp(3,naclsd,nclsd),stat=i_stat)

    ! create temporary copy
    rmtref_temp(:) = rmtref(:)
    vref_temp(:) = vref(:)
    nacls_temp(:) = nacls(:)
    rcls_temp(:,:,:) = rcls(:,:,:)

    ! update array dimensions
    nrefd = nref
    nclsd = ncls

    ! Reallocate arrays
    deallocate (rmtref, vref, nacls, rcls, stat=i_stat)
    allocate(rmtref(nrefd),stat=i_stat)
    allocate(vref(nrefd),stat=i_stat)
    allocate(nacls(nclsd),stat=i_stat)
    allocate(rcls(3,naclsd,nclsd),stat=i_stat)

    ! copy value from temp arrays
    rmtref(:) = rmtref_temp(1:nref)
    vref(:) = vref_temp(1:nref)
    nacls(:) = nacls_temp(1:ncls)
    rcls(:,:,:) = rcls_temp(:,:,1:ncls)

    ! cleanup deallocations
    deallocate(rmtref_temp, vref_temp, nacls_temp, rcls_temp, stat=i_stat)

  end subroutine reduce_array_size

end module mod_main0
