!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Module responsible for storing the input variables and primary arrays so that they are distributed via MPI processes.
!> Author: Philipp Ruessmann and many others ...
!> Previously this routine wrote unformatted files to disk, so that they
!> would be used by the different executables. Since the advent of the single
!> executable mode, this routine creates a copy of most of the variables in the program
!> as special `type` parameters. This are then used in the MPI communication, and
!> in the rest of the variables used in the code.
!------------------------------------------------------------------------------------
!> @note Jonatan Chico: 02.01.2018 Modifications to ensure compatibility for the removal of
!> the inc.p file. Also added the memory profiling calls to the allocation/deallocation
!> of the arrays.
!> @endnote
!------------------------------------------------------------------------------------
module mod_wunfiles

  use :: mod_profiling
  use :: mod_datatypes
  use :: mod_constants, only: nsymaxd
  use :: mod_bfield, only: type_bfield

  implicit none

  !-------------------------------------------------------------------------------
  !> Summary: Type holding information of parameters for the communication of data
  !> Author:
  !> Category: communication, KKRhost
  !> Deprecated: False
  !> Type holding information of parameters for the communication of data
  !-------------------------------------------------------------------------------
  ! define type that replace wunfiles here, later define bcast routine
  type :: type_params

    integer :: nscalars = 126

    ! .. Scalars
    integer :: i1
    integer :: nr       !! Number of real space vectors rr
    integer :: irm      !! Maximum number of radial points
    integer :: lly      !! LLY!> 0 : apply Lloyds formula
    integer :: ins      !! 0 (MT), 1(ASA), 2(Full Potential)
    integer :: icc      !! Enables the calculation of off-diagonal elements of the GF.(0=SCF/DOS; 1=cluster; -1=custom)
    integer :: igf      !! Do not print or print (0/1) the KKRFLEX_* files
    integer :: kte      !! Calculation of the total energy On/Off (1/0)
    integer :: kxc      !! Type of xc-potential 0=vBH 1=MJW 2=VWN 3=PW91
    integer :: naez     !! Number of atoms in unit cell
    integer :: lmax     !! Maximum l component in wave function expansion
    integer :: nref     !! Number of diff. ref. potentials
    integer :: lm2d     !! (2*LMAX+1)**2
    integer :: irid     !! Shape functions parameters in non-spherical part
    integer :: krel     !! Switch for non-relativistic/relativistic (0/1) program. Attention: several other parameters depend explicitly on KREL, they are set automatically Used for Dirac solver in ASA
    integer :: kpre
    integer :: nsra
    integer :: nemb     !! Number of sites added to the slab in 2D calculations to extend the structure left and right (down and up)
    integer :: ncls     !! Number of reference clusters
    integer :: iend     !! Number of nonzero gaunt coefficients
    integer :: ncpa     !! NCPA = 0/1 CPA flag
    integer :: icst     !! Number of Born approximation
    integer :: imix     !! Type of mixing scheme used (0=straight, 4=Broyden 2nd, 5=Anderson)
    integer :: itab
    integer :: lpot     !! Maximum l component in potential expansion
    integer :: npol     !! Number of Matsubara Poles (EMESHT)
    integer :: npnt1    !! number of E points (EMESHT) for the contour integration
    integer :: npnt2    !! number of E points (EMESHT) for the contour integration
    integer :: npnt3    !! number of E points (EMESHT) for the contour integration
    integer :: itscf
    integer :: iemxd    !! Dimension for energy-dependent arrays
    integer :: npotd    !! (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP)
    integer :: natyp    !! Number of kinds of atoms in unit cell
    integer :: ipand    !! Number of panels in non-spherical part
    integer :: ncleb    !! Number of Clebsch-Gordon coefficients
    integer :: nclsd    !! Maximum number of different TB-clusters
    integer :: nfund    !! Shape functions parameters in non-spherical part
    integer :: ngshd    !! Shape functions parameters in non-spherical part
    integer :: mmaxd    !! 2*LMAX+1
    integer :: nineq    !! Number of ineq. positions in unit cell
    integer :: nspin    !! Counter for spin directions
    integer :: kmrot    !! 0: no rotation of the magnetisation; 1: individual rotation of the magnetisation for every site
    integer :: iltmp
    integer :: ncheb    !! Number of Chebychev pannels for the new solver
    integer :: ntotd
    integer :: kvmad
    integer :: irnsd    !! Number of radial mesh points in (RMT,...,RWS)
    integer :: knoco    !! (0/1) Collinear/Non-collinear magnetism (even in non-relativistic non-spin-orbit case)
    integer :: lmpot    !! (LPOT+1)**2
    integer :: nleft    !! Number of repeated basis for left host to get converged electrostatic potentials
    integer :: nright   !! Number of repeated basis for right host to get converged electrostatic potentials
    integer :: korbit   !! Spin-orbit/non-spin-orbit (1/0) added to the Schroedinger or SRA equations. Works with FP. KREL and KORBIT cannot be both non-zero.
    integer :: ntperd   !! Parameter in broyden subroutines
    integer :: ielast
    integer :: nrmaxd   !! NTOTD*(NCHEBD+1)
    integer :: ishift
    integer :: knosph   !! Switch for spherical/non-spherical (0/1) program. Same obs. as for KREL applies.
    integer :: kforce   !! Calculation of the forces
    integer :: itdbry   !! Number of SCF steps to remember for the Broyden mixing
    integer :: kshape   !! Exact treatment of WS cell
    integer :: nofgij   !! number of GF pairs IJ to be calculated as determined from IJTABCALC<>0
    integer :: nspind   !! KREL+(1-KREL)*(NSPIN+1)
    integer :: irmind   !! IRM-IRNSD
    integer :: nspotd   !! Number of potentials for storing non-sph. potentials
    integer :: nembd1   !! NEMB+1
    integer :: lmmaxd   !! (KREL+KORBIT+1)(LMAX+1)^2
    integer :: nembd2
    integer :: naclsd   !! Maximum number of atoms in a TB-cluster
    integer :: lmaxd1
    integer :: nsheld   !! Number of blocks of the GF matrix that need to be calculated (NATYP + off-diagonals in case of impurity)
    integer :: ncelld   !! = naez in main0, could/should be replaced to make the code more understandable
    integer :: lmxspd   !! (2*LPOT+1)**2
    integer :: nsymat
    integer :: nprinc   !! Number of atoms in one principal layer
    integer :: n1semi   !! Number of energy points for the semicore contour
    integer :: n2semi   !! Number of energy points for the semicore contour
    integer :: n3semi   !! Number of energy points for the semicore contour
    integer :: invmod   !! Inversion scheme
    integer :: nqcalc
    integer :: ntldau   !! number of atoms on which LDA+U is applied
    integer :: kpoibz   !! Number of reciprocal space vectors
    integer :: nsatypd  !! (NATYP-1)*NSPIN+1
    integer :: idoldau  !! flag to perform LDA+U
    integer :: nlayerd  !! Number of principal layers (NAEZD/NPRINCD) used in the inversion routines (independent on NATYPD)
    integer :: intervx  !! Number of intervals in x-direction for k-net in IB of the BZ
    integer :: intervy  !! Number of intervals in y-direction for k-net in IB of the BZ
    integer :: intervz  !! Number of intervals in z-direction for k-net in IB of the BZ
    integer :: nlbasis  !! Number of basis layers of left host (repeated units)
    integer :: nrbasis  !! Number of basis layers of right host (repeated units)
    integer :: wlength  !! Word length for direct access files, compiler dependent ifort/others (1/4)
    integer :: naezdpd
    integer :: maxmesh
    integer :: itmpdir
    integer :: nspindd  !! NSPIND-KORBIT
    integer :: npan_eq  !! Variables for the pannels for the new solver
    integer :: npan_log !! Variables for the pannels for the new solver
    integer :: scfsteps !! number of scf iterations
    integer :: itcpamax !! Max. number of CPA iterations
    integer :: natomimp !! Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
    integer :: nmvecmax
    integer :: npolsemi !! Number of poles for the semicore contour
    integer :: natomimpd !! Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
    integer :: itrunldau !! Iteration index for LDA+U
    integer :: iesemicore
    integer :: special_straight_mixing !!id to specify modified straight mixing scheme: 0=normal, 1=alternating mixing factor (i.e. reduced mixing factor in every odd iteration), 2=charge-neurality based mixing factor (former: 'alt mix' and 'spec mix')
    real (kind=dp) :: tk           !! Temperature
    real (kind=dp) :: fcm
    real (kind=dp) :: emin         !! Energies needed in EMESHT
    real (kind=dp) :: emax         !! Energies needed in EMESHT
    real (kind=dp) :: alat         !! Lattice constant in a.u.
    real (kind=dp) :: r_log
    real (kind=dp) :: efold
    real (kind=dp) :: denef
    real (kind=dp) :: efermi       !! Fermi energy
    real (kind=dp) :: cpatol       !! Convergency tolerance for CPA-cycle
    real (kind=dp) :: mixing       !! Magnitude of the mixing parameter
    real (kind=dp) :: qbound       !! Convergence parameter for the potential
    real (kind=dp) :: tksemi       !! Temperature of semi-core contour
    real (kind=dp) :: chrgold
    real (kind=dp) :: tolrdif      !! Tolerance for r<tolrdif (a.u.) to handle vir. atoms
    real (kind=dp) :: lasterr
    real (kind=dp) :: emusemi
    real (kind=dp) :: ebotsemi
    real (kind=dp) :: fsemicore
    real (kind=dp) :: lambda_xc    !! Scale magnetic moment (0! Lambda_XC! 1, 0=zero moment, 1= full moment)
    real (kind=dp) :: chrgsemicore
    complex (kind=dp) :: deltae    !! Energy difference for numerical derivative
    logical :: lnc                 !! Coupled equations in two spins (switches true if KREL=1 or KORBIT=1 or KNOCO=1)
    logical :: lrhosym
    logical :: linterface          !! If True a matching with semi-inifinite surfaces must be performed
    character (len=10) :: solver   !! Type of solver

    character (len=80) :: tmpdir

    ! .. Arrays
    complex (kind=dp), dimension (:), allocatable :: ez
    complex (kind=dp), dimension (:), allocatable :: wez
    complex (kind=dp), dimension (:, :), allocatable :: rc !! NREL REAL spher. harm. > CMPLX. spher. harm. NREL CMPLX. spher. harm. > REAL spher. harm.
    complex (kind=dp), dimension (:, :), allocatable :: crel !! Non-relat. CMPLX. spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. CMPLX. spher. harm.
    complex (kind=dp), dimension (:, :), allocatable :: rrel !! Non-relat. REAL spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. REAL spher. harm.
    complex (kind=dp), dimension (:, :), allocatable :: phildau
    complex (kind=dp), dimension (:, :, :), allocatable :: srrel
    complex (kind=dp), dimension (:, :, :), allocatable :: drotq !! Rotation matrices to change between LOCAL/GLOBAL frame of reference for magnetisation!> Oz or noncollinearity
    complex (kind=dp), dimension (:, :, :), allocatable :: dsymll
    complex (kind=dp), dimension (:, :, :, :, :), allocatable :: lefttinvll
    complex (kind=dp), dimension (:, :, :, :, :), allocatable :: righttinvll
    complex (kind=dp), dimension (:, :, :), allocatable :: mvevi
    complex (kind=dp), dimension (:, :, :), allocatable :: mvevief
    real (kind=dp), dimension (:), allocatable :: a !! Constants for exponential R mesh
    real (kind=dp), dimension (:), allocatable :: b !! Constants for exponential R mesh
    real (kind=dp), dimension (:), allocatable :: eu
    real (kind=dp), dimension (:), allocatable :: edc
    real (kind=dp), dimension (:), allocatable :: vbc !! Potential constants
    real (kind=dp), dimension (:), allocatable :: zat !! Nuclear charge
    real (kind=dp), dimension (:), allocatable :: rmt !! Muffin-tin radius of true system
    real (kind=dp), dimension (:), allocatable :: rws !! Wigner Seitz radius
    real (kind=dp), dimension (:), allocatable :: gsh
    real (kind=dp), dimension (:), allocatable :: phi !! Phi angle for a non-collinear calculation
    real (kind=dp), dimension (:), allocatable :: ueff !! input U parameter for each atom
    real (kind=dp), dimension (:), allocatable :: jeff !! input J parameter for each atom
    real (kind=dp), dimension (:), allocatable :: vref
    real (kind=dp), dimension (:), allocatable :: conc !! Concentration of a given atom
    real (kind=dp), dimension (:), allocatable :: theta !! Theta angle for a non-collinear calculation (don't confuse with thetas!)
    real (kind=dp), dimension (:), allocatable :: volbz
    real (kind=dp), dimension (:), allocatable :: qmtet !! \(\theta\) angle of the agnetization with respect to the z-axis
    real (kind=dp), dimension (:), allocatable :: qmphi !! \(\phi\) angle of the agnetization with respect to the z-axis
    real (kind=dp), dimension (:), allocatable :: rmtref !! Muffin-tin radius of reference system
    real (kind=dp), dimension (:), allocatable :: rmtnew !! Adapted muffin-tin radius
    real (kind=dp), dimension (:), allocatable :: denefat
    real (kind=dp), dimension (:), allocatable :: erefldau !! the energies of the projector's wave functions (REAL)
    real (kind=dp), dimension (:), allocatable :: socscale !! Spin-orbit scaling
    real (kind=dp), dimension (:, :, :), allocatable :: vins !! Non-spherical part of the potential
    real (kind=dp), dimension (:, :), allocatable :: rmesh !! Radial mesh ( in units a Bohr)
    real (kind=dp), dimension (:, :), allocatable :: rr
    real (kind=dp), dimension (:, :), allocatable :: drdi !! Derivative dr/di
    real (kind=dp), dimension (:, :), allocatable :: cscl !! Speed of light scaling
    real (kind=dp), dimension (:, :), allocatable :: cleb !! GAUNT coefficients (GAUNT)
    real (kind=dp), dimension (:, :), allocatable :: rhoc
    real (kind=dp), dimension (:, :), allocatable :: espv
    real (kind=dp), dimension (:, :), allocatable :: rnew
    real (kind=dp), dimension (:, :), allocatable :: visp !! Spherical part of the potential
    real (kind=dp), dimension (:, :), allocatable :: vtrel !! potential (spherical part)
    real (kind=dp), dimension (:, :), allocatable :: btrel !! magnetic field
    real (kind=dp), dimension (:, :), allocatable :: ecore !! Core energies
    real (kind=dp), dimension (:, :), allocatable :: rmrel !! radial mesh
    real (kind=dp), dimension (:, :), allocatable :: ratom
    real (kind=dp), dimension (:, :), allocatable :: rbasis !! Position of atoms in the unit cell in units of bravais vectors
    real (kind=dp), dimension (:, :), allocatable :: socscl
    real (kind=dp), dimension (:, :), allocatable :: volcub
    real (kind=dp), dimension (:, :), allocatable :: rhoorb
    real (kind=dp), dimension (:, :), allocatable :: rclsimp
    real (kind=dp), dimension (:, :), allocatable :: drdirel !! derivative of radial mesh
    real (kind=dp), dimension (:, :), allocatable :: ecorerel
    real (kind=dp), dimension (:, :), allocatable :: cmomhost !! Charge moments of each atom of the (left/right) host
    real (kind=dp), dimension (:, :), allocatable :: qmphitab
    real (kind=dp), dimension (:, :), allocatable :: qmtettab
    real (kind=dp), dimension (:, :), allocatable :: qmgamtab
    real (kind=dp), dimension (:, :), allocatable :: r2drdirel !! \( r^2 \frac{\partial}{\partial \mathbf{r}}\frac{\partial}{\partial i}\) (r**2 * drdi)
    real (kind=dp), dimension (:, :), allocatable :: rpan_intervall

    real (kind=dp), dimension (:, :, :), allocatable :: rcls !! Real space position of atom in cluster
    real (kind=dp), dimension (:, :, :), allocatable :: rrot
    real (kind=dp), dimension (:, :, :), allocatable :: bzkp
    real (kind=dp), dimension (:, :, :), allocatable :: thetas !! shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
    real (kind=dp), dimension (:, :, :), allocatable :: thetasnew
    real (kind=dp), dimension (:, :, :, :), allocatable :: r2nef
    real (kind=dp), dimension (:, :, :, :), allocatable :: wldau !! potential matrix
    real (kind=dp), dimension (:, :, :, :), allocatable :: rho2ns
    real (kind=dp), dimension (:, :, :, :, :), allocatable :: uldau !! calculated Coulomb matrix elements (EREFLDAU)
    integer, dimension (:), allocatable :: cls !! Cluster around atomic sites
    integer, dimension (:), allocatable :: noq !! Number of diff. atom types located
    integer, dimension (:), allocatable :: imt !! R point at MT radius
    integer, dimension (:), allocatable :: irc !! R point for potential cutting
    integer, dimension (:), allocatable :: nfu
    integer, dimension (:), allocatable :: zrel !! atomic number (cast integer)
    integer, dimension (:), allocatable :: lopt !! angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
    integer, dimension (:), allocatable :: ipan !! Number of panels in non-MT-region
    integer, dimension (:), allocatable :: iqat !! The site on which an atom is located on a given site
    integer, dimension (:), allocatable :: icpa !! ICPA = 0/1 site-dependent CPA flag
    integer, dimension (:), allocatable :: irns !! Position of atoms in the unit cell in units of bravais vectors
    integer, dimension (:), allocatable :: irws !! R point at WS radius
    integer, dimension (:), allocatable :: nsh1 !! Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
    integer, dimension (:), allocatable :: nsh2 !! Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
    integer, dimension (:), allocatable :: irmin !! Max R for spherical treatment
    integer, dimension (:), allocatable :: ncore !! Number of core states
    integer, dimension (:), allocatable :: nacls !! Number of atoms in cluster
    integer, dimension (:), allocatable :: nofks
    integer, dimension (:), allocatable :: loflm !! l of lm=(l,m) (GAUNT)
    integer, dimension (:), allocatable :: kmesh
    integer, dimension (:), allocatable :: itldau !! integer pointer connecting the NTLDAU atoms to heir corresponding index in the unit cell
    integer, dimension (:), allocatable :: nshell !! Index of atoms/pairs per shell (ij-pairs); nshell(0) = number of shells
    integer, dimension (:), allocatable :: iqcalc
    integer, dimension (:), allocatable :: refpot !! Ref. pot. card  at position
    integer, dimension (:), allocatable :: ntcell !! Index for WS cell
    integer, dimension (:), allocatable :: ixipol !! Constraint of spin pol.
    integer, dimension (:), allocatable :: jwsrel !! index of the WS radius
    integer, dimension (:), allocatable :: imaxsh
    integer, dimension (:), allocatable :: atomimp
    integer, dimension (:), allocatable :: hostimp
    integer, dimension (:), allocatable :: irshift !! shift of the REL radial mesh with respect no NREL
    integer, dimension (:), allocatable :: ijtabsh !! Linear pointer, assigns pair (i,j) to a shell in the array GS(*,*,*,NSHELD)
    integer, dimension (:), allocatable :: npan_tot
    integer, dimension (:), allocatable :: ijtabsym !! Linear pointer, assigns pair (i,j) to the rotation bringing GS into Gij
    integer, dimension (:), allocatable :: ijtabcalc !! Linear pointer, specifying whether the block (i,j) has to be calculated needs set up for ICC=-1, not used for ICC=1
    integer, dimension (:), allocatable :: npan_eq_at
    integer, dimension (:), allocatable :: npan_log_at
    integer, dimension (:), allocatable :: ijtabcalc_i
    integer, dimension (:), allocatable :: qdos_atomselect
    integer, dimension (:, :), allocatable :: ish
    integer, dimension (:, :), allocatable :: jsh
    integer, dimension (:, :), allocatable :: ilm_map
    integer, dimension (:, :), allocatable :: ezoa !! EZ of atom at site in cluster
    integer, dimension (:, :), allocatable :: atom !! Atom at site in cluster
    integer, dimension (:, :), allocatable :: lmsp !! 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
    integer, dimension (:, :), allocatable :: icleb !! Pointer array
    integer, dimension (:, :), allocatable :: lcore !! Angular momentum of core states
    integer, dimension (:, :), allocatable :: ircut !! R points of panel borders
    integer, dimension (:, :), allocatable :: kaoez !! Kind of atom at site in elem. cell
    integer, dimension (:, :), allocatable :: nrrel
    integer, dimension (:, :), allocatable :: lmsp1
    integer, dimension (:, :), allocatable :: ifunm
    integer, dimension (:, :), allocatable :: llmsp !! lm=(l,m) of 'nfund'th nonvanishing component of non-spherical pot.
    integer, dimension (:, :), allocatable :: icheck
    integer, dimension (:, :), allocatable :: ifunm1
    integer, dimension (:, :), allocatable :: ititle
    integer, dimension (:, :), allocatable :: nkcore
    integer, dimension (:, :), allocatable :: kapcore
    integer, dimension (:, :), allocatable :: ipan_intervall
    integer, dimension (:, :, :), allocatable :: jend !! Pointer array for icleb()
    integer, dimension (:, :, :), allocatable :: irrel
    logical, dimension (:), allocatable :: vacflag
    logical, dimension (:), allocatable :: symunitary !! unitary/antiunitary symmetry flag
    !character (len=8), dimension (:), allocatable :: optc
    !character (len=8), dimension (:), allocatable :: testc
    character (len=124), dimension (:), allocatable :: txc
    type (type_bfield) :: bfield

  end type type_params

  type (type_params), save :: t_params

contains

  !-------------------------------------------------------------------------------
  !> Summary: This routine takes the read parameters from the `inputcard` and stores them in the `t_params` type to be distributed via MPI
  !> Author: Philipp Rüssmann and many others ...
  !> Category: communication, input-output, KKRhost
  !> Deprecated: False
  !> This routine was oiginally meant to write unformated files to then
  !> be read by other executables, now it does the same job via storing types instead
  !> reducing I/O and allowing for MPI communication.
  !-------------------------------------------------------------------------------
  subroutine wunfiles(npol,npnt1,npnt2,npnt3,ielast,tk,emin,emax,ez,wez,efermi,     &
    npolsemi,n1semi,n2semi,n3semi,iesemicore,tksemi,ebotsemi,emusemi,fsemicore,vins,&
    visp,vbc,vtrel,btrel,rmrel,drdirel,r2drdirel,zrel,jwsrel,irshift,itscf,scfsteps,&
    cmomhost,ecore,lcore,ncore,qmtet,qmphi,qmphitab,qmtettab,qmgamtab,drotq,nsra,   &
    ins,natyp,naez,nineq,nref,nspin,ncls,icst,ipan,ircut,alat,zat,r,drdi,refpot,    &
    rmtref,vref,iend,jend,cleb,icleb,atom,cls,rcls,nacls,loflm,solver,socscl,cscl,  &
    icc,igf,nlbasis,nrbasis,ncpa,icpa,itcpamax,cpatol,rbasis,rr,ezoa,nshell,nsh1,   &
    nsh2,ijtabcalc,ijtabcalc_i,ish,jsh,ijtabsym,ijtabsh,nofgij,nqcalc,iqcalc,kmrot, &
    kaoez,iqat,noq,conc,kmesh,maxmesh,nsymat,symunitary,rrot,dsymll,invmod,icheck,  &
    natomimp,ratom,atomimp,rc,crel,rrel,srrel,nrrel,irrel,lefttinvll,righttinvll,   &
    vacflag,a,b,ifunm,ifunm1,intervx,intervy,intervz,ititle,lmsp1,ntcell,thetas,    &
    lpot,lmpot,nright,nleft,linterface,imix,mixing,qbound,fcm,itdbry,irns,kpre,     &
    kshape,kte,kvmad,kxc,lambda_xc,txc,ishift,ixipol,lrhosym,kforce,lmsp,llmsp,rmt, &
    rmtnew,rws,imt,irc,irmin,irws,nfu,hostimp,gsh,ilm_map,imaxsh,idoldau,itrunldau, &
    ntldau,lopt,itldau,ueff,jeff,erefldau,uldau,wldau,phildau,iemxd,irmind,irm,     &
    nspotd,npotd,nembd1,lmmaxd,ipand,nembd2,lmax,ncleb,naclsd,nclsd,lm2d,lmaxd1,    &
    mmaxd,nr,nsheld,naezdpd,natomimpd,nspind,irid,nfund,ncelld,lmxspd,ngshd,        &
    krel,ntotd,ncheb,npan_log,npan_eq,npan_log_at,npan_eq_at,r_log,npan_tot,rnew,   &
    rpan_intervall,ipan_intervall,nspindd,thetasnew,socscale,tolrdif,lly,deltae,    &
    rclsimp,verbosity,MPI_scheme,special_straight_mixing,bfield)
    ! **********************************************************************
    ! *                                                                    *
    ! *  This subroutine is part of the MAIN0 program in the tbkkr package *
    ! *  It writes out different unformatted files meant to provide the    *
    ! *  communication between the other parts (MAIN1a, 1b, 1c and 2)      *
    ! *  during an SCF cycle                                               *
    ! *  v.popescu, munich 2004                                            *
    ! *                                                                    *
    ! **********************************************************************

    use :: mod_types, only: t_inc, t_tgmat, t_lloyd, t_cpa
    use :: mod_mympi, only: nranks, mpiatom, mpiadapt
    use :: mod_runoptions, only: impurity_operator_only, relax_SpinAngle_Dirac, use_Chebychev_solver, use_qdos, &
      write_cpa_projection_file, write_deci_tmat, write_gmat_file, write_green_host, write_green_imp, write_gref_file, &
      write_kkrimp_input, write_lloyd_cdos_file, write_lloyd_dgref_file, write_lloyd_dtmat_file, write_lloyd_file, &
      write_lloyd_g0tr_file, write_lloyd_tralpha_file, write_pkkr_input, write_pkkr_operators, write_tmat_file, set_cheby_nosoc

    implicit none
    ! ..
    ! .. Scalar arguments
    integer, intent (in) :: nr     !! Number of real space vectors rr
    integer, intent (in) :: irm    !! Maximum number of radial points
    integer, intent (in) :: ins    !! 0 (MT), 1(ASA), 2(Full Potential)
    integer, intent (in) :: icc    !! Enables the calculation of off-diagonal elements of the GF.(0=SCF/DOS; 1=cluster; -1=custom)
    integer, intent (in) :: igf    !! Do not print or print (0/1) the KKRFLEX_* files
    integer, intent (in) :: lly    !! LLY!> 0 : apply Lloyds formula
    integer, intent (in) :: kte    !! Calculation of the total energy On/Off (1/0)
    integer, intent (in) :: kxc    !! Type of xc-potential 0=vBH 1=MJW 2=VWN 3=PW91
    integer, intent (in) :: kpre
    integer, intent (in) :: lpot   !! Maximum l component in potential expansion
    integer, intent (in) :: imix   !! Type of mixing scheme used (0=straight, 4=Broyden 2nd, 5=Anderson)
    integer, intent (in) :: ncls   !! Number of reference clusters
    integer, intent (in) :: icst   !! Number of Born approximation
    integer, intent (in) :: iend   !! Number of nonzero gaunt coefficients
    integer, intent (in) :: lmax   !! Maximum l component in wave function expansion
    integer, intent (in) :: nref   !! Number of diff. ref. potentials
    integer, intent (in) :: naez   !! Number of atoms in unit cell
    integer, intent (in) :: irid   !! Shape functions parameters in non-spherical part
    integer, intent (in) :: krel   !! Switch for non-relativistic/relativistic (0/1) program. Attention: several other parameters depend explicitly on KREL, they are set automatically Used for Dirac solver in ASA
    integer, intent (in) :: ncpa   !! NCPA = 0/1 CPA flag
    integer, intent (in) :: lm2d   !! (2*LMAX+1)**2
    integer, intent (in) :: nsra
    integer, intent (in) :: npol   !! Number of Matsubara Poles (EMESHT)
    integer, intent (in) :: npnt1  !! number of E points (EMESHT) for the contour integration
    integer, intent (in) :: npnt2  !! number of E points (EMESHT) for the contour integration
    integer, intent (in) :: npnt3  !! number of E points (EMESHT) for the contour integration
    integer, intent (in) :: itscf
    integer, intent (in) :: kvmad
    integer, intent (in) :: lmpot  !! (LPOT+1)**2
    integer, intent (in) :: ngshd  !! Shape functions parameters in non-spherical part
    integer, intent (in) :: mmaxd  !! 2*LMAX+1
    integer, intent (in) :: npotd  !! (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP)
    integer, intent (in) :: natyp  !! Number of kinds of atoms in unit cell
    integer, intent (in) :: ipand  !! Number of panels in non-spherical part
    integer, intent (in) :: nclsd  !! Maximum number of different TB-clusters
    integer, intent (in) :: ncleb  !! Number of Clebsch-Gordon coefficients
    integer, intent (in) :: nfund  !! Shape functions parameters in non-spherical part
    integer, intent (in) :: iemxd  !! Dimension for energy-dependent arrays
    integer, intent (in) :: nineq  !! Number of ineq. positions in unit cell
    integer, intent (in) :: nspin  !! Counter for spin directions
    integer, intent (in) :: kmrot  !! 0: no rotation of the magnetisation; 1: individual rotation of the magnetisation for every site
    integer, intent (in) :: ntotd
    integer, intent (in) :: ncheb  !! Number of Chebychev pannels for the new solver
    integer, intent (in) :: nleft  !! Number of repeated basis for left host to get converged electrostatic potentials
    integer, intent (in) :: nright !! Number of repeated basis for right host to get converged electrostatic potentials
    integer, intent (in) :: itdbry !! Number of SCF steps to remember for the Broyden mixing
    integer, intent (in) :: kshape !! Exact treatment of WS cell
    integer, intent (in) :: ishift
    integer, intent (in) :: ielast
    integer, intent (in) :: kforce !! Calculation of the forces
    integer, intent (in) :: ncelld !! Number of cells (shapes) in non-spherical part
    integer, intent (in) :: lmxspd !! (2*LPOT+1)**2
    integer, intent (in) :: irmind !! IRM-IRNSD
    integer, intent (in) :: nspotd !! Number of potentials for storing non-sph. potentials
    integer, intent (in) :: nembd1 !! NEMB+1
    integer, intent (in) :: lmmaxd !! (KREL+KORBIT+1)(LMAX+1)^2
    integer, intent (in) :: nembd2
    integer, intent (in) :: naclsd !! Maximum number of atoms in a TB-cluster
    integer, intent (in) :: lmaxd1
    integer, intent (in) :: nofgij !! number of GF pairs IJ to be calculated as determined from IJTABCALC<>0
    integer, intent (in) :: nspind !! KREL+(1-KREL)*(NSPIN+1)
    integer, intent (in) :: nsheld !! Number of blocks of the GF matrix that need to be calculated (NATYP + off-diagonals in case of impurity)
    integer, intent (in) :: nsymat
    integer, intent (in) :: invmod !! Inversion scheme
    integer, intent (in) :: ntldau !! number of atoms on which LDA+U is applied
    integer, intent (in) :: nqcalc
    integer, intent (in) :: n1semi !! Number of energy points for the semicore contour
    integer, intent (in) :: n2semi !! Number of energy points for the semicore contour
    integer, intent (in) :: n3semi !! Number of energy points for the semicore contour
    integer, intent (in) :: npan_eq !! Variables for the pannels for the new solver
    integer, intent (in) :: naezdpd
    integer, intent (in) :: nlbasis !! Number of basis layers of left host (repeated units)
    integer, intent (in) :: nrbasis !! Number of basis layers of right host (repeated units)
    integer, intent (in) :: nspindd !! NSPIND-KORBIT
    integer, intent (in) :: maxmesh
    integer, intent (in) :: intervx !! Number of intervals in x-direction for k-net in IB of the BZ
    integer, intent (in) :: intervy !! Number of intervals in y-direction for k-net in IB of the BZ
    integer, intent (in) :: intervz !! Number of intervals in z-direction for k-net in IB of the BZ
    integer, intent (in) :: idoldau !! flag to perform LDA+U
    integer, intent (in) :: npolsemi !! Number of poles for the semicore contour
    integer, intent (inout) :: scfsteps !! number of scf iterations
    integer, intent (in) :: itcpamax !! Max. number of CPA iterations
    integer, intent (in) :: natomimp !! Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
    integer, intent (in) :: npan_log !! Variables for the pannels for the new solver
    integer, intent (in) :: natomimpd !! Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
    integer, intent (in) :: itrunldau !! Iteration index for LDA+U
    integer, intent (in) :: iesemicore
    integer, intent (in) :: verbosity !! verbosity level for timings and output: 0=old default, 1,2,3 = timing and ouput verbosity level the same (low,medium,high)
    integer, intent (in) :: MPI_scheme !!!! scheme for MPI parallelization: 0 = automatic (default), 1 = atoms, 2 = energies, 3 = select best of (1,2)
    integer, intent (in) :: special_straight_mixing !!id to specify modified straight mixing scheme: 0=normal, 1=alternating mixing factor (i.e. reduced mixing factor in every odd iteration), 2=charge-neurality based mixing factor (former: 'alt mix' and 'spec mix')
    ! .. nembd2 = NAEZ+NEMB, lmaxd1=lmaxd+1, naezdpd=NAEZ/nprincd)
    real (kind=dp), intent (in) :: tk !! Temperature
    real (kind=dp), intent (in) :: fcm
    real (kind=dp), intent (in) :: alat !! Lattice constant in a.u.
    real (kind=dp), intent (inout) :: emin !! Energies needed in EMESHT
    real (kind=dp), intent (in) :: emax !! Energies needed in EMESHT
    real (kind=dp), intent (in) :: r_log
    real (kind=dp), intent (in) :: efermi !! Fermi energy
    real (kind=dp), intent (in) :: cpatol !! Convergency tolerance for CPA-cycle
    real (kind=dp), intent (in) :: mixing !! Magnitude of the mixing parameter
    real (kind=dp), intent (in) :: qbound !! Convergence parameter for the potential
    real (kind=dp), intent (in) :: tksemi !! Temperature of semi-core contour
    real (kind=dp), intent (in) :: emusemi
    real (kind=dp), intent (in) :: tolrdif !! Tolerance for r<tolrdif (a.u.) to handle vir. atoms
    real (kind=dp), intent (in) :: ebotsemi
    real (kind=dp), intent (in) :: fsemicore
    real (kind=dp), intent (in) :: lambda_xc !! Scale magnetic moment (0! Lambda_XC! 1, 0=zero moment, 1= full moment)
    logical, intent (in) :: lrhosym
    logical, intent (in) :: linterface !! If True a matching with semi-inifinite surfaces must be performed
    character (len=10), intent (in) :: solver                           !! Type of solver

    complex (kind=dp), intent (in) :: deltae !! Energy difference for numerical derivative
    ! ..
    ! .. Array arguments
    complex (kind=dp), dimension (iemxd), intent (in) :: ez
    complex (kind=dp), dimension (iemxd), intent (in) :: wez
    complex (kind=dp), dimension (lmmaxd, lmmaxd), intent (in) :: rc !! NREL REAL spher. harm. > CMPLX. spher. harm. NREL CMPLX. spher. harm. > REAL spher. harm.
    complex (kind=dp), dimension (lmmaxd, lmmaxd), intent (in) :: crel !! Non-relat. CMPLX. spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. CMPLX. spher. harm.
    complex (kind=dp), dimension (lmmaxd, lmmaxd), intent (in) :: rrel !! Non-relat. REAL spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. REAL spher. harm.
    complex (kind=dp), dimension (irm, natyp), intent (in) :: phildau

    complex (kind=dp), dimension (lmmaxd, lmmaxd, naez), intent (in) :: drotq !! Rotation matrices to change between LOCAL/GLOBAL frame of reference for magnetisation!> Oz or noncollinearity
    complex (kind=dp), dimension (lmmaxd, lmmaxd, nsymaxd), intent (in) :: dsymll
    complex (kind=dp), dimension (2, 2, lmmaxd), intent (in) :: srrel
    complex (kind=dp), dimension (lmmaxd, lmmaxd, nembd1, nspindd, iemxd), intent (in) :: lefttinvll
    complex (kind=dp), dimension (lmmaxd, lmmaxd, nembd1, nspindd, iemxd), intent (in) :: righttinvll
    real (kind=dp), dimension (natyp), intent (in) :: a !! Constants for exponential R mesh
    real (kind=dp), dimension (natyp), intent (in) :: b !! Constants for exponential R mesh
    real (kind=dp), dimension (2), intent (in) :: vbc !! Potential constants
    real (kind=dp), dimension (natyp), intent (in) :: zat !! Nuclear charge
    real (kind=dp), dimension (natyp), intent (in) :: rmt !! Muffin-tin radius of true system
    real (kind=dp), dimension (natyp), intent (in) :: rws !! Wigner Seitz radius
    real (kind=dp), dimension (ngshd), intent (in) :: gsh
    real (kind=dp), dimension (natyp), intent (in) :: conc !! Concentration of a given atom
    real (kind=dp), dimension (nref), intent (in) :: vref
    real (kind=dp), dimension (natyp), intent (in) :: ueff !! input U parameter for each atom
    real (kind=dp), dimension (natyp), intent (in) :: jeff !! input J parameter for each atom
    real (kind=dp), dimension (naez), intent (in) :: qmtet !! \( \theta\) angle of the agnetization with respect to the z-axis
    real (kind=dp), dimension (naez), intent (in) :: qmphi !! \( \phi\) angle of the agnetization with respect to the z-axis
    real (kind=dp), dimension (nref), intent (in) :: rmtref !! Muffin-tin radius of reference system
    real (kind=dp), dimension (natyp), intent (in) :: rmtnew !! Adapted muffin-tin radius
    real (kind=dp), dimension (natyp), intent (in) :: erefldau !! the energies of the projector's wave functions (REAL)
    real (kind=dp), dimension (natyp), intent (in) :: socscale !! Spin-orbit scaling
    real (kind=dp), dimension (irm, natyp), intent (in) :: r !! Radial mesh ( in units a Bohr)
    real (kind=dp), dimension (3, 0:nr), intent (in) :: rr
    real (kind=dp), dimension (irm, natyp), intent (in) :: drdi !! Derivative dr/di
    real (kind=dp), dimension (ncleb, 2), intent (in) :: cleb !! GAUNT coefficients (GAUNT)
    ! real (kind=dp), dimension(LMAXD1,NATYP), intent(in)          :: CSCL     !! Speed of light scaling
    real (kind=dp), dimension (krel*lmax+1, krel*natyp+(1-krel)), intent (inout) :: cscl !! Speed of light scaling
    real (kind=dp), dimension (ntotd*(ncheb+1), natyp), intent (in) :: rnew
    real (kind=dp), dimension (irm, npotd), intent (in) :: visp !! Spherical part of the potential
    ! real (kind=dp), dimension(IRM,NATYP), intent(in)             :: VTREL    !! potential (spherical part)
    real (kind=dp), dimension (irm*krel+(1-krel), natyp), intent (inout) :: vtrel !! potential (spherical part)
    ! real (kind=dp), dimension(IRM,NATYP), intent(in)             :: BTREL    !! magnetic field
    real (kind=dp), dimension (irm*krel+(1-krel), natyp), intent (inout) :: btrel !! magnetic field
    ! real (kind=dp), dimension(IRM,NATYP), intent(in)             :: RMREL    !! radial mesh
    real (kind=dp), dimension (irm*krel+(1-krel), natyp), intent (inout) :: rmrel !! radial mesh
    real (kind=dp), dimension (3, nsheld), intent (in) :: ratom
    real (kind=dp), dimension (20, npotd), intent (in) :: ecore !! Core energies
    real (kind=dp), dimension (3, nembd2), intent (in) :: rbasis !! Position of atoms in the unit cell in units of bravais vectors
    ! real (kind=dp), dimension(LMAXD1,NATYP), intent(in)          :: SOCSCL
    real (kind=dp), dimension (krel*lmax+1, krel*natyp+(1-krel)), intent (inout) :: socscl
    ! real (kind=dp), dimension(IRM,NATYP), intent(in)             :: DRDIREL  !! derivative of radial mesh
    real (kind=dp), dimension (irm*krel+(1-krel), natyp), intent (in) :: drdirel
    real (kind=dp), dimension (naez, 3), intent (in) :: qmphitab
    real (kind=dp), dimension (naez, 3), intent (in) :: qmtettab
    real (kind=dp), dimension (naez, 3), intent (in) :: qmgamtab
    real (kind=dp), dimension (3, natomimpd), intent (in) :: rclsimp
    real (kind=dp), dimension (lmpot, nembd1), intent (in) :: cmomhost !! Charge moments of each atom of the (left/right) host
    ! real (kind=dp), dimension(IRM,NATYP), intent(in)             :: R2DRDIREL   !! \( r^2 \frac{\partial}{\partial \mathbf{r}}\frac{\partial}{\partial i}\) (r**2 * drdi)
    real (kind=dp), dimension (irm*krel+(1-krel), natyp), intent (in) :: r2drdirel
    real (kind=dp), dimension (0:ntotd, natyp), intent (in) :: rpan_intervall
    real (kind=dp), dimension (48, 3, nsheld), intent (in) :: rrot
    real (kind=dp), dimension (3, naclsd, nclsd), intent (in) :: rcls !! Real space position of atom in cluster
    real (kind=dp), dimension (irmind:irm, lmpot, nspotd), intent (in) :: vins !! Non-spherical part of the potential
    real (kind=dp), dimension (irid, nfund, ncelld), intent (in) :: thetas !! shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
    real (kind=dp), dimension (ntotd*(ncheb+1), nfund, ncelld), intent (in) :: thetasnew
    real (kind=dp), dimension (mmaxd, mmaxd, nspind, natyp), intent (in) :: wldau !! potential matrix
    real (kind=dp), dimension (mmaxd, mmaxd, mmaxd, mmaxd, natyp), intent (in) :: uldau !! calculated Coulomb matrix elements (EREFLDAU)
    ! ..
    integer, dimension (naez), intent (in) :: noq !! Number of diff. atom types located
    integer, dimension (natyp), intent (in) :: imt !! R point at MT radius
    integer, dimension (natyp), intent (in) :: irc !! R point for potential cutting
    integer, dimension (natyp), intent (in) :: nfu
    integer, dimension (nembd2), intent (in) :: cls !! Cluster around atomic sites
    integer, dimension (naez), intent (in) :: icpa !! ICPA = 0/1 site-dependent CPA flag
    integer, dimension (natyp), intent (in) :: ipan !! Number of panels in non-MT-region
    integer, dimension (natyp), intent (in) :: zrel !! atomic number (cast integer)
    integer, dimension (natyp), intent (in) :: iqat !! The site on which an atom is located on a given site
    integer, dimension (natyp), intent (in) :: lopt !! angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
    integer, dimension (natyp), intent (in) :: irns !! Position of atoms in the unit cell in units of bravais vectors
    integer, dimension (nsheld), intent (in) :: nsh1 !! Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
    integer, dimension (nsheld), intent (in) :: nsh2 !! Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
    integer, dimension (natyp), intent (in) :: irws !! R point at WS radius
    integer, dimension (lm2d), intent (in) :: loflm !! l of lm=(l,m) (GAUNT)
    integer, dimension (nclsd), intent (in) :: nacls !! Number of atoms in cluster
    integer, dimension (iemxd), intent (in) :: kmesh
    integer, dimension (npotd), intent (in) :: ncore !! Number of core states
    integer, dimension (naez), intent (in) :: iqcalc
    integer, dimension (natyp), intent (in) :: irmin !! Max R for spherical treatment
    integer, dimension (natyp), intent (in) :: ntcell !! Index for WS cell
    integer, dimension (natyp), intent (in) :: ixipol !! Constraint of spin pol.
    integer, dimension (natyp), intent (in) :: jwsrel !! index of the WS radius
    integer, dimension (natyp), intent (in) :: itldau !! integer pointer connecting the NTLDAU atoms to heir corresponding index in the unit cel
    integer, dimension (nembd2), intent (in) :: refpot !! Ref. pot. card  at position
    integer, dimension (0:lmpot), intent (in) :: imaxsh
    integer, dimension (0:nsheld), intent (in) :: nshell !! Index of atoms/pairs per shell (ij-pairs); nshell(0) = number of shells
    integer, dimension (0:natyp), intent (in) :: hostimp
    integer, dimension (natyp), intent (in) :: irshift !! shift of the REL radial mesh with respect no NREL
    integer, dimension (nofgij), intent (in) :: ijtabsh !! Linear pointer, assigns pair (i,j) to a shell in the array GS(*,*,*,NSHELD)
    integer, dimension (natomimpd), intent (in) :: atomimp
    integer, dimension (natyp), intent (in) :: npan_tot
    integer, dimension (nofgij), intent (in) :: ijtabsym !! Linear pointer, assigns pair (i,j) to the rotation bringing GS into Gij
    integer, dimension (nofgij), intent (in) :: ijtabcalc !! Linear pointer, specifying whether the block (i,j) has to be calculated needs set up for ICC=-1, not used for ICC=1
    integer, dimension (natyp), intent (in) :: npan_eq_at
    integer, dimension (natyp), intent (in) :: npan_log_at
    integer, dimension (nofgij), intent (in) :: ijtabcalc_i
    integer, dimension (ngshd, 3), intent (in) :: ilm_map
    integer, dimension (nsheld, 2*nsymaxd), intent (in) :: ish
    integer, dimension (nsheld, 2*nsymaxd), intent (in) :: jsh
    integer, dimension (naclsd, nembd2), intent (in) :: atom !! Atom at site in cluster
    integer, dimension (naclsd, nembd2), intent (in) :: ezoa !! EZ of atom at site in cluster
    integer, dimension (natyp, lmxspd), intent (in) :: lmsp !! 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
    integer, dimension (natyp, nfund), intent (in) :: llmsp !! lm=(l,m) of 'nfund'th nonvanishing component of non-spherical pot.
    integer, dimension (lmxspd, natyp), intent (in) :: lmsp1
    integer, dimension (natyp, lmxspd), intent (in) :: ifunm
    integer, dimension (natyp, nembd2), intent (in) :: kaoez !! Kind of atom at site in elem. cell
    integer, dimension (0:ipand, natyp), intent (in) :: ircut !! R points of panel borders
    integer, dimension (ncleb, 4), intent (in) :: icleb !! Pointer array
    integer, dimension (20, npotd), intent (in) :: lcore !! Angular momentum of core states
    integer, dimension (2, lmmaxd), intent (in) :: nrrel
    integer, dimension (naezdpd, naezdpd), intent (in) :: icheck
    integer, dimension (lmxspd, natyp), intent (in) :: ifunm1
    integer, dimension (20, npotd), intent (in) :: ititle
    integer, dimension (0:ntotd, natyp), intent (in) :: ipan_intervall
    integer, dimension (lmpot, 0:lmax, 0:lmax), intent (in) :: jend !! Pointer array for icleb()
    integer, dimension (2, 2, lmmaxd), intent (in) :: irrel
    logical, dimension (2), intent (in) :: vacflag
    logical, dimension (nsymaxd), intent (in) :: symunitary !! unitary/antiunitary symmetry flag
    character (len=124), dimension (6), intent (in) :: txc
    ! Non-collinear Bfield options
    type (type_bfield), intent (in) :: bfield
    ! .. Local scalars
    integer :: i1
    integer :: ic, naclsmin, naclsmax, nqdos, irmdnew ! variables for t_inc filling
    integer :: itmpdir, iltmp
    real (kind=dp), dimension (natyp) :: phi
    real (kind=dp), dimension (natyp) :: theta
    character (len=80) :: tmpdir
    ! .. External Functions

    itmpdir = 0
    iltmp = 0

    ! put information about scf steps into t_inc
    t_inc%i_iteration = itscf
    t_inc%n_iteration = scfsteps
    ! put information for save_wavefun also in:
    t_inc%nsra = nsra
    irmdnew = 0
    do i1 = 1, natyp
      if (npan_tot(i1)*(ncheb+1)>irmdnew) then
        irmdnew = npan_tot(i1)*(ncheb+1)
      end if
    end do
    t_inc%irmdnew = irmdnew

    ! -------------------------------------------------------------------------
    ! itermdir
    ! -------------------------------------------------------------------------
    if (relax_SpinAngle_Dirac) then
      i1 = 0
      emin = 0.0_dp
      t_params%qmtet = qmtet
      t_params%qmphi = qmphi
      t_params%qmphitab = qmphitab
      t_params%qmtettab = qmtettab
      t_params%qmgamtab = qmgamtab
      t_params%i1 = i1
      t_params%emin = emin
      t_params%drotq = drotq
    end if
    ! -------------------------------------------------------------------------
    ! lda+u data in this file change
    ! -------------------------------------------------------------------------
    if (idoldau==1) then
      open (67, file='ldau.unformatted', form='unformatted')
      write (67) itrunldau, wldau, uldau, phildau
      close (67)
    end if
    ! -------------------------------------------------------------------------
    ! nonco_angle file
    ! -------------------------------------------------------------------------
    call read_angles(t_params, natyp, theta, phi)
    ! -------------------------------------------------------------------------
    ! fill t_inc type meant for simpler passing of basic parameter to other routines
    ! -------------------------------------------------------------------------

    ! find maximal cluster information (taken from main1b)
    naclsmax = 1
    do ic = 1, ncls
      if (nacls(ic)>naclsmax) naclsmax = nacls(ic)
    end do

    ! find NQDOS (taken from main1b)
    if (use_qdos) then
      open (67, file='qvec.dat')
      read (67, *) nqdos
      close (67)
    else
      nqdos = 1
    end if

    !--------------------------------------------------------------------------------
    ! t_inc t_inc t_inc t_inc t_inc t_inc t_inc t_inc t_inc t_inc
    !--------------------------------------------------------------------------------
    ! fill t_inc
    t_inc%ielast = ielast
    t_inc%nqdos = nqdos
    t_inc%naclsmax = naclsmax
    t_inc%nshell0 = nshell(0)
    if (use_Chebychev_solver) t_inc%newsosol = .true.
    if (set_cheby_nosoc) t_inc%nosoc = .true.
    if (write_deci_tmat) t_inc%deci_out = .true.
    !--------------------------------------------------------------------------------
    ! t_inc t_inc t_inc t_inc t_inc t_inc t_inc t_inc t_inc t_inc
    !--------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------
    ! writeout flags writeout flags writeout flags writeout flags writeout flags writeout flags writeout flags writeout flags
    !--------------------------------------------------------------------------------
    ! set logical switches in t_tgmat which control if tmat, gmat and gref are written to files or stored in memory
    if (write_tmat_file) t_tgmat%tmat_to_file = .true.
    if (write_gmat_file) t_tgmat%gmat_to_file = .true.
    if (write_gref_file) t_tgmat%gref_to_file = .true.
    if (write_cpa_projection_file) t_cpa%dmatproj_to_file = .true.


    !--------------------------------------------------------------------------------
    ! bug bug bug bug bug
    !--------------------------------------------------------------------------------
    ! in case of ASA DIRAC solver (KREL==1) then gmat file has to be written out otherwise something is going wrong.
    if (krel>0) t_tgmat%gmat_to_file = .true.
    ! bug bug bug bug bug

    ! some special run options:
    if (write_kkrimp_input) t_tgmat%tmat_to_file = .true. ! for KKRFLEX option tmat must be written to file
    if (use_qdos) t_tgmat%gmat_to_file = .true. ! for qdos write gmat to file since it scales with NQDOS and can become huge

    ! set logical switches in t_lloyd which control if files are written to files or stored in memory
    ! if(write_tmat_file.or.write_lloyd_file)
    if (write_lloyd_dtmat_file .or. write_lloyd_file) t_lloyd%dtmat_to_file = .true.
    if (write_lloyd_tralpha_file .or. write_lloyd_file) t_lloyd%tralpha_to_file = .true.
    if (write_lloyd_cdos_file .or. write_lloyd_file) t_lloyd%cdos_diff_lly_to_file = .true.
    if (write_lloyd_dgref_file .or. write_lloyd_file) t_lloyd%dgref_to_file = .true.
    if (write_lloyd_g0tr_file .or. write_lloyd_file) t_lloyd%g0tr_to_file = .true.

    ! set verbosity level in t_inc%i_write = 0,1,2 for default, verbose1, verbose2
    t_inc%i_write = 0                   ! default: write only output.000.txt and reset file after each iteration
    if (verbosity==2) t_inc%i_write = 1 ! write on all processors but only the latest iteration
    if (verbosity==3) t_inc%i_write = 2 ! write everything
    ! and t_inc_i_time for timing writeout
    t_inc%i_time = 1               ! default: only timings from master, all iterations
    if (verbosity==1) t_inc%i_time = 0 ! only timings from master, only the last iteration
    if (verbosity==3) t_inc%i_time = 2 ! all timing files, all iterations
    ! writeout flags writeout flags writeout flags writeout flags writeout flags writeout flags writeout flags writeout flags

    !--------------------------------------------------------------------------------
    ! MPI communication scheme
    !--------------------------------------------------------------------------------
    ! set switch for MPIatom test option (see mod_types and mod_mympi)
    ! default values for MPIadapt and MPIatom
    if (natyp<=ielast) then
      ! for more energy points than atoms we want energy parallelization
      mpiatom = .false.
    else
      ! for for a large number of atoms we want to use atom parallelization
      ! since this makes the k-loop and the atom loops slow
      mpiatom = .true.
    end if
    mpiadapt = 1                   ! 1 means check timings and then reshuffle ranks if necessary
    ! change default behaviour for the corresponding test flags are found
    if (MPI_scheme==1) then
      mpiatom = .true.
      mpiadapt = 0                 ! 0 means no change in communication scheme
    end if
    if (MPI_scheme==2) then
      mpiatom = .false.
      mpiadapt = 0
    end if
    if (MPI_scheme==3) then
      mpiadapt = 2                 ! 2 means force run with MPIatom then with MPIenerg and then compare to choose optimal
    end if
    ! so far changing does not work yet, so turn this off always:
    mpiadapt = 0

    if (write_pkkr_input .or. write_green_host .or. write_green_imp .or. write_pkkr_operators) then ! fswrt
      mpiatom = .true.             ! fswrt
      if (scfsteps>1) then         ! fswrt
        write (*, *) 'Warning: Setting SCFSTEPS=1 for FERMIOUT option' ! fswrt
        scfsteps = 1               ! fswrt
      end if                       ! fswrt
      if (nranks>natyp) then       ! fswrt
        write (*, *) 'FERMIOUT/WRTGREEN/GREENIMP option chosen' ! fswrt
        write (*, *) 'Nranks>NATYP', nranks, natyp ! fswrt
        stop 'Please choose Nranks<=NATYP' ! fswrt
      end if                       ! fswrt
      naclsmin = minval(nacls(1:ncls)) ! fswrt
      if (.not. write_green_imp .and. naclsmin<150) then ! fswrt
        write (*, *) ' !!!  WARNING  !!!' ! fswrt
        write (*, *) '   FERMIOUT/WRTGREEN option chosen' ! fswrt
        write (*, *) '   minimal cluster size smaller than 150 atoms!!!' ! fswrt
        write (*, *) '   should be increased to least 200-300 atoms' ! fswrt
      end if                       ! fswrt
      if (MPI_scheme==2) then   ! fswrt
        write (*, *) 'FERMIOUT/WRTGREEN/GREENIMP option chosen' ! fswrt
        write (*, *) 'found unsupported test option MPIenerg' ! fswrt
        stop 'Please choose MPIatom instead' ! fswrt
      end if                       ! fswrt
      if (write_pkkr_operators .and. write_pkkr_input) then ! fswrt
        write (*, *) 'OPERATOR and FERMIOUT cannot be used together' ! fswrt
        stop 'Please chose only one of the two' ! fswrt
      end if                       ! fswrt
      if (write_pkkr_operators .and. ielast/=1) then ! fswrt
        write (*, *) 'OPERATOR option chosen' ! fswrt
        write (*, *) 'energy contour should contain a single point' ! fswrt
        write (*, *) 'on the real axis only' ! fswrt
        stop 'Please correct energy contour' ! fswrt
      end if                       ! fswrt
    end if                         ! fswrt
    if (.not. write_pkkr_operators .and. impurity_operator_only) then ! fswrt
      write (*, *) 'test option "IMP_ONLY" can only be used' ! fswrt
      write (*, *) 'in combination with option "OPERATOR"' ! fswrt
      stop                         ! fswrt
    end if                         ! fswrt
    

    !--------------------------------------------------------------------------------
    ! Non-collinear magnetic field and constraining fields
    !--------------------------------------------------------------------------------
    ! store magnetic field information
    t_params%bfield = bfield

    !--------------------------------------------------------------------------------
    ! MPI communication scheme
    !--------------------------------------------------------------------------------

    ! all parameters are stored in t_params fomr mod_wunfiles
    ! first fill scalar values
    call fill_t_params_scalars(iemxd,irmind,irm,lmpot,nspotd,npotd,natyp,nembd1,    &
      lmmaxd,naez,ipand,nembd2,nref,lmax,ncleb,naclsd,nclsd,lm2d,lmaxd1,nr,nsheld,  &
      naezdpd,natomimpd,nofgij,nspind,nspindd,irid,nfund,ncelld,lmxspd,             &
      ngshd,krel,mmaxd,ielast,npol,npnt1,npnt2,npnt3,itscf,scfsteps,lly,nsra,ins,   &
      nineq,nspin,ncls,icst,iend,icc,igf,nlbasis,nrbasis,ncpa,itcpamax,kmrot,       &
      maxmesh,nsymat,natomimp,invmod,nqcalc,intervx,intervy,intervz,lpot,nright,    &
      nleft,imix,itdbry,kpre,kshape,kte,kvmad,kxc,ishift,kforce,idoldau,itrunldau,  &
      ntldau,npolsemi,n1semi,n2semi,n3semi,iesemicore,ebotsemi,emusemi,tksemi,      &
      fsemicore,r_log,emin,emax,tk,efermi,alat,cpatol,mixing,qbound,fcm,lambda_xc,  &
      tolrdif,linterface,lrhosym,solver,tmpdir,itmpdir,iltmp,ntotd,ncheb,deltae,    &
      special_straight_mixing,                                                      &
      t_params)

    ! initialize allocatable arrays
    call init_t_params(t_params)

    ! now fill arrays that have just been allocated
    call fill_t_params_arrays(t_params,iemxd,lmmaxd,naez,nembd1,nspindd,            &
      irmind,irm,lmpot,nspotd,npotd,natyp,nr,nembd2,nref,ncleb,nclsd,naclsd,nsheld, &
      ngshd,nfund,irid,ncelld,mmaxd,lm2d,lmxspd,lmaxd1,nspind,ntotd,ncheb,ipand,    &
      lmax,nofgij,naezdpd,natomimpd,ez,wez,drotq,dsymll,lefttinvll,righttinvll,crel,&
      rc,rrel,srrel,phildau,vins,visp,vbc,vtrel,btrel,socscale,drdirel,r2drdirel,   &
      rmrel,cmomhost,ecore,qmtet,qmphi,qmphitab,qmtettab,qmgamtab,zat,r,drdi,rmtref,&
      vref,cleb,rcls,socscl,cscl,rbasis,rr,conc,rrot,ratom,a,b,thetas,rmt,rmtnew,   &
      rws,gsh,erefldau,ueff,jeff,uldau,wldau,rpan_intervall,rnew,thetasnew,lopt,    &
      itldau,irshift,jwsrel,zrel,lcore,ncore,ipan,ircut,jend,icleb,atom,cls,nacls,  &
      loflm,ezoa,kaoez,iqat,icpa,noq,kmesh,nshell,nsh1,nsh2,ijtabcalc,ijtabcalc_i,  &
      ijtabsym,ijtabsh,ish,jsh,iqcalc,icheck,atomimp,refpot,irrel,nrrel,ifunm1,     &
      ititle,lmsp1,ntcell,ixipol,irns,ifunm,llmsp,lmsp,imt,irc,irmin,irws,nfu,      &
      hostimp,ilm_map,imaxsh,npan_log,npan_eq,npan_log_at,npan_eq_at,npan_tot,      &
      ipan_intervall,symunitary,vacflag,txc,rclsimp,krel)

    ! save information about the energy mesh
    call save_emesh(ielast,ez,wez,emin,emax,iesemicore,fsemicore,npol,tk,npnt1,     &
      npnt2,npnt3,ebotsemi,emusemi,tksemi,npolsemi,n1semi,n2semi,n3semi,iemxd,      &
      t_params)


  end subroutine wunfiles


  !-------------------------------------------------------------------------------
  !> Summary: Allocate initial parameters to be broadcasted via mpi
  !> Author: Philipp Ruessmann
  !> Category: memory-management, profiling, KKRhost
  !> Deprecated: False
  !> Allocate initial parameters to be broadcasted via mpi. allocate arrays, has to
  !> be done after `bcast t_params_scalars` for myrank<>master otherwise are the parameters not set
  !-------------------------------------------------------------------------------
  subroutine init_t_params(t_params)

    implicit none

    type (type_params), intent (inout) :: t_params

    integer :: i_stat

    ! -------------------------------------------------------------------------
    ! Allocate complex (kind=dp) arrays
    ! -------------------------------------------------------------------------
    allocate (t_params%ez(t_params%iemxd), stat=i_stat) ! complex (kind=dp)
    call memocc(i_stat, product(shape(t_params%ez))*kind(t_params%ez), 't_params%EZ', 'init_t_params')
    allocate (t_params%wez(t_params%iemxd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%wez))*kind(t_params%wez), 't_params%WEZ', 'init_t_params')
    allocate (t_params%drotq(t_params%lmmaxd,t_params%lmmaxd,t_params%naez), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%drotq))*kind(t_params%drotq), 't_params%DROTQ', 'init_t_params')
    allocate (t_params%dsymll(t_params%lmmaxd,t_params%lmmaxd,nsymaxd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%dsymll))*kind(t_params%dsymll), 't_params%DSYMLL', 'init_t_params')
    allocate (t_params%lefttinvll(t_params%lmmaxd,t_params%lmmaxd,t_params%nembd1,t_params%nspindd,t_params%iemxd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%lefttinvll))*kind(t_params%lefttinvll), 't_params%LEFTTINVLL', 'init_t_params')
    allocate (t_params%righttinvll(t_params%lmmaxd,t_params%lmmaxd,t_params%nembd1,t_params%nspindd,t_params%iemxd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%righttinvll))*kind(t_params%righttinvll), 't_params%RIGHTTINVLL', 'init_t_params')
    allocate (t_params%crel(t_params%lmmaxd,t_params%lmmaxd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%crel))*kind(t_params%crel), 't_params%CREL', 'init_t_params')
    allocate (t_params%rc(t_params%lmmaxd,t_params%lmmaxd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%rc))*kind(t_params%rc), 't_params%RC', 'init_t_params')
    allocate (t_params%rrel(t_params%lmmaxd,t_params%lmmaxd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%rrel))*kind(t_params%rrel), 't_params%RREL', 'init_t_params')
    allocate (t_params%srrel(2,2,t_params%lmmaxd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%srrel))*kind(t_params%srrel), 't_params%SRREL', 'init_t_params')
    allocate (t_params%phildau(t_params%irm,t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%phildau))*kind(t_params%phildau), 't_params%PHILDAU', 'init_t_params')
    ! -------------------------------------------------------------------------
    ! End of allocation of complex (kind=dp) arrays
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! Allocate real (kind=dp) arrays
    ! -------------------------------------------------------------------------
    allocate (t_params%vins(t_params%irmind:t_params%irm,t_params%lmpot,t_params%nspotd), stat=i_stat) ! real (kind=dp)
    call memocc(i_stat, product(shape(t_params%vins))*kind(t_params%vins), 't_params%VINS', 'init_t_params')
    allocate (t_params%visp(t_params%irm,t_params%npotd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%visp))*kind(t_params%visp), 't_params%VISP', 'init_t_params')
    allocate (t_params%vbc(2), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%vbc))*kind(t_params%vbc), 't_params%VBC', 'init_t_params')
    allocate (t_params%vtrel(t_params%irm,t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%vtrel))*kind(t_params%vtrel), 't_params%VTREL', 'init_t_params')
    allocate (t_params%btrel(t_params%irm,t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%btrel))*kind(t_params%btrel), 't_params%BTREL', 'init_t_params')
    allocate (t_params%socscale(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%socscale))*kind(t_params%socscale), 't_params%SOCSCALE', 'init_t_params')
    allocate (t_params%drdirel(t_params%irm,t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%drdirel))*kind(t_params%drdirel), 't_params%DRDIREL', 'init_t_params')
    allocate (t_params%r2drdirel(t_params%irm,t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%r2drdirel))*kind(t_params%r2drdirel), 't_params%R2DRDIREL', 'init_t_params')
    allocate (t_params%rmrel(t_params%irm,t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%rmrel))*kind(t_params%rmrel), 't_params%RMREL', 'init_t_params')
    allocate (t_params%cmomhost(t_params%lmpot,t_params%nembd1), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%cmomhost))*kind(t_params%cmomhost), 't_params%CMOMHOST', 'init_t_params')
    allocate (t_params%ecore(20,t_params%npotd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%ecore))*kind(t_params%ecore), 't_params%ECORE', 'init_t_params')
    allocate (t_params%qmtet(t_params%naez), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%qmtet))*kind(t_params%qmtet), 't_params%QMTET', 'init_t_params')
    allocate (t_params%qmphi(t_params%naez), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%qmphi))*kind(t_params%qmphi), 't_params%QMPHI', 'init_t_params')
    allocate (t_params%qmphitab(t_params%naez,3), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%qmphitab))*kind(t_params%qmphitab), 't_params%QMPHITAB', 'init_t_params')
    allocate (t_params%qmtettab(t_params%naez,3), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%qmtettab))*kind(t_params%qmtettab), 't_params%QMTETTAB', 'init_t_params')
    allocate (t_params%qmgamtab(t_params%naez,3), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%qmgamtab))*kind(t_params%qmgamtab), 't_params%QMGAMTAB', 'init_t_params')
    allocate (t_params%zat(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%zat))*kind(t_params%zat), 't_params%ZAT', 'init_t_params')
    allocate (t_params%rmesh(t_params%irm,t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%rmesh))*kind(t_params%rmesh), 't_params%RMESH', 'init_t_params')
    allocate (t_params%drdi(t_params%irm,t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%drdi))*kind(t_params%drdi), 't_params%DRDI', 'init_t_params')
    allocate (t_params%rmtref(t_params%nref), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%rmtref))*kind(t_params%rmtref), 't_params%RMTREF', 'init_t_params')
    allocate (t_params%vref(t_params%nref), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%vref))*kind(t_params%vref), 't_params%VREF', 'init_t_params')
    allocate (t_params%cleb(t_params%ncleb,2), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%cleb))*kind(t_params%cleb), 't_params%CLEB', 'init_t_params')
    allocate (t_params%rcls(3,t_params%naclsd,t_params%nclsd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%rcls))*kind(t_params%rcls), 't_params%RCLS', 'init_t_params')
    allocate (t_params%socscl(t_params%lmaxd1,t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%socscl))*kind(t_params%socscl), 't_params%SOCSCL', 'init_t_params')
    allocate (t_params%cscl(t_params%lmaxd1,t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%cscl))*kind(t_params%cscl), 't_params%CSCL', 'init_t_params')
    allocate (t_params%rbasis(3,t_params%nembd2), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%rbasis))*kind(t_params%rbasis), 't_params%RBASIS', 'init_t_params')
    allocate (t_params%rr(3,0:t_params%nr), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%rr))*kind(t_params%rr), 't_params%RR', 'init_t_params')
    allocate (t_params%conc(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%conc))*kind(t_params%conc), 't_params%CONC', 'init_t_params')
    allocate (t_params%rrot(48,3,t_params%nsheld), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%rrot))*kind(t_params%rrot), 't_params%RROT', 'init_t_params')
    allocate (t_params%ratom(3,t_params%nsheld), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%ratom))*kind(t_params%ratom), 't_params%RATOM', 'init_t_params')
    allocate (t_params%a(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%a))*kind(t_params%a), 't_params%A', 'init_t_params')
    allocate (t_params%b(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%b))*kind(t_params%b), 't_params%B', 'init_t_params')
    allocate (t_params%thetas(t_params%irid,t_params%nfund,t_params%ncelld), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%thetas))*kind(t_params%thetas), 't_params%THETAS', 'init_t_params')
    allocate (t_params%rmt(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%rmt))*kind(t_params%rmt), 't_params%RMT', 'init_t_params')
    allocate (t_params%rmtnew(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%rmtnew))*kind(t_params%rmtnew), 't_params%RMTNEW', 'init_t_params')
    allocate (t_params%rws(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%rws))*kind(t_params%rws), 't_params%RWS', 'init_t_params')
    allocate (t_params%gsh(t_params%ngshd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%gsh))*kind(t_params%gsh), 't_params%GSH', 'init_t_params')
    allocate (t_params%erefldau(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%erefldau))*kind(t_params%erefldau), 't_params%EREFLDAU', 'init_t_params')
    allocate (t_params%ueff(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%ueff))*kind(t_params%ueff), 't_params%UEFF', 'init_t_params')
    allocate (t_params%jeff(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%jeff))*kind(t_params%jeff), 't_params%JEFF', 'init_t_params')
    allocate (t_params%uldau(t_params%mmaxd,t_params%mmaxd,t_params%mmaxd,t_params%mmaxd,t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%uldau))*kind(t_params%uldau), 't_params%ULDAU', 'init_t_params')
    allocate (t_params%wldau(t_params%mmaxd,t_params%mmaxd,t_params%nspind,t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%wldau))*kind(t_params%wldau), 't_params%WLDAU', 'init_t_params')
    allocate (t_params%rpan_intervall(0:t_params%ntotd,t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%rpan_intervall))*kind(t_params%rpan_intervall), 't_params%RPAN_INTERVALL', 'init_t_params')
    allocate (t_params%rnew(t_params%ntotd*(t_params%ncheb+1),t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%rnew))*kind(t_params%rnew), 't_params%RNEW', 'init_t_params')
    allocate (t_params%mvevi(t_params%natyp,3,t_params%nmvecmax), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%mvevi))*kind(t_params%mvevi), 't_params%MVEVI', 'init_t_params')
    allocate (t_params%mvevief(t_params%natyp,3,t_params%nmvecmax), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%mvevief))*kind(t_params%mvevief), 't_params%MVEVIEF', 'init_t_params')
    allocate (t_params%thetasnew(t_params%ntotd*(t_params%ncheb+1),t_params%nfund,t_params%ncelld), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%thetasnew))*kind(t_params%thetasnew), 't_params%THETASNEW', 'init_t_params')
    allocate (t_params%rho2ns(t_params%irm,t_params%lmpot,t_params%natyp,2), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%rho2ns))*kind(t_params%rho2ns), 't_params%RHO2NS', 'init_t_params')
    allocate (t_params%r2nef(t_params%irm,t_params%lmpot,t_params%natyp,2), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%r2nef))*kind(t_params%r2nef), 't_params%R2NEF', 'init_t_params')
    allocate (t_params%rhoc(t_params%irm,t_params%npotd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%rhoc))*kind(t_params%rhoc), 't_params%RHOC', 'init_t_params')
    allocate (t_params%denefat(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%denefat))*kind(t_params%denefat), 't_params%DENEFAT', 'init_t_params')
    allocate (t_params%espv(0:t_params%lmaxd1,t_params%npotd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%espv))*kind(t_params%espv), 't_params%ESPV', 'init_t_params')
    allocate (t_params%edc(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%edc))*kind(t_params%edc), 't_params%EDC', 'init_t_params')
    allocate (t_params%eu(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%eu))*kind(t_params%eu), 't_params%EU', 'init_t_params')
    allocate (t_params%rhoorb(t_params%irm*t_params%krel+(1-t_params%krel),t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%rhoorb))*kind(t_params%rhoorb), 't_params%RHOORB', 'init_t_params')
    allocate (t_params%ecorerel(t_params%krel*20+(1-t_params%krel),t_params%npotd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%ecorerel))*kind(t_params%ecorerel), 't_params%ECOREREL', 'init_t_params')
    allocate (t_params%rclsimp(3,t_params%natomimpd), stat=i_stat) ! real (kind=dp)
    call memocc(i_stat, product(shape(t_params%rclsimp))*kind(t_params%rclsimp), 't_params%ECOREREL', 'init_t_params')
    if (.not. allocated(t_params%theta)) then
      allocate (t_params%theta(t_params%natyp), stat=i_stat) ! real (kind=dp)
      call memocc(i_stat, product(shape(t_params%theta))*kind(t_params%theta), 't_params%THETA', 'init_t_params')
      allocate (t_params%phi(t_params%natyp), stat=i_stat) ! real (kind=dp)
      call memocc(i_stat, product(shape(t_params%phi))*kind(t_params%phi), 't_params%PHI', 'init_t_params')
    end if
    ! -------------------------------------------------------------------------
    ! End of allocation of real (kind=dp) arrays
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! Allocate the integer arrays
    ! -------------------------------------------------------------------------
    allocate (t_params%lopt(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%lopt))*kind(t_params%lopt), 't_params%LOPT', 'init_t_params')
    allocate (t_params%itldau(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%itldau))*kind(t_params%itldau), 't_params%ITLDAU', 'init_t_params')
    allocate (t_params%irshift(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%irshift))*kind(t_params%irshift), 't_params%IRSHIFT', 'init_t_params')
    allocate (t_params%jwsrel(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%jwsrel))*kind(t_params%jwsrel), 't_params%JWSREL', 'init_t_params')
    allocate (t_params%zrel(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%zrel))*kind(t_params%zrel), 't_params%ZREL', 'init_t_params')
    allocate (t_params%lcore(20,t_params%npotd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%lcore))*kind(t_params%lcore), 't_params%LCORE', 'init_t_params')
    allocate (t_params%ncore(t_params%npotd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%npotd))*kind(t_params%npotd), 't_params%NPOTD', 'init_t_params')
    allocate (t_params%ipan(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%ipan))*kind(t_params%ipan), 't_params%IPAN', 'init_t_params')
    allocate (t_params%ircut(0:t_params%ipand,t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%ircut))*kind(t_params%ircut), 't_params%IRCUT', 'init_t_params')
    allocate (t_params%jend(t_params%lmpot,0:t_params%lmax,0:t_params%lmax), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%jend))*kind(t_params%jend), 't_params%JEND', 'init_t_params')
    allocate (t_params%icleb(t_params%ncleb,4), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%icleb))*kind(t_params%icleb), 't_params%ICLEB', 'init_t_params')
    allocate (t_params%atom(t_params%naclsd,t_params%nembd2), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%atom))*kind(t_params%atom), 't_params%ATOM', 'init_t_params')
    allocate (t_params%cls(t_params%nembd2), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%cls))*kind(t_params%cls), 't_params%CLS', 'init_t_params')
    allocate (t_params%nacls(t_params%nclsd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%nacls))*kind(t_params%nacls), 't_params%NACLS', 'init_t_params')
    allocate (t_params%loflm(t_params%lm2d), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%loflm))*kind(t_params%loflm), 't_params%LOFLM', 'init_t_params')
    allocate (t_params%ezoa(t_params%naclsd,t_params%nembd2), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%ezoa))*kind(t_params%ezoa), 't_params%EZOA', 'init_t_params')
    allocate (t_params%kaoez(t_params%natyp,t_params%nembd2), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%kaoez))*kind(t_params%kaoez), 't_params%KAOEZ', 'init_t_params')
    allocate (t_params%iqat(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%iqat))*kind(t_params%iqat), 't_params%IQAT', 'init_t_params')
    allocate (t_params%icpa(t_params%naez), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%icpa))*kind(t_params%icpa), 't_params%ICPA', 'init_t_params')
    allocate (t_params%noq(t_params%naez), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%noq))*kind(t_params%noq), 't_params%NOQ', 'init_t_params')
    allocate (t_params%kmesh(t_params%iemxd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%kmesh))*kind(t_params%kmesh), 't_params%KMESH', 'init_t_params')
    allocate (t_params%nshell(0:t_params%nsheld), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%nshell))*kind(t_params%nshell), 't_params%NSHELL', 'init_t_params')
    allocate (t_params%nsh1(t_params%nsheld), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%nsh1))*kind(t_params%nsh1), 't_params%NSH1', 'init_t_params')
    allocate (t_params%nsh2(t_params%nsheld), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%nsh2))*kind(t_params%nsh2), 't_params%NSH2', 'init_t_params')
    allocate (t_params%ijtabcalc(t_params%nofgij), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%ijtabcalc))*kind(t_params%ijtabcalc), 't_params%IJTABCALC', 'init_t_params')
    allocate (t_params%ijtabcalc_i(t_params%nofgij), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%ijtabcalc_i))*kind(t_params%ijtabcalc_i), 't_params%IJTABCALC_I', 'init_t_params')
    allocate (t_params%ijtabsym(t_params%nofgij), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%ijtabsym))*kind(t_params%ijtabsym), 't_params%IJTABSYM', 'init_t_params')
    allocate (t_params%ijtabsh(t_params%nofgij), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%ijtabsh))*kind(t_params%ijtabsh), 't_params%IJTABSH', 'init_t_params')
    allocate (t_params%ish(t_params%nsheld,2*nsymaxd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%ish))*kind(t_params%ish), 't_params%ISH', 'init_t_params')
    allocate (t_params%jsh(t_params%nsheld,2*nsymaxd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%jsh))*kind(t_params%jsh), 't_params%JSH', 'init_t_params')
    allocate (t_params%iqcalc(t_params%naez), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%iqcalc))*kind(t_params%iqcalc), 't_params%IQCALC', 'init_t_params')
    allocate (t_params%icheck(t_params%naezdpd,t_params%naezdpd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%icheck))*kind(t_params%icheck), 't_params%ICHECK', 'init_t_params')
    allocate (t_params%atomimp(t_params%natomimpd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%atomimp))*kind(t_params%atomimp), 't_params%ATOMIMP', 'init_t_params')
    allocate (t_params%refpot(t_params%nembd2), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%refpot))*kind(t_params%refpot), 't_params%REFPOT', 'init_t_params')
    allocate (t_params%irrel(2,2,t_params%lmmaxd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%irrel))*kind(t_params%irrel), 't_params%IRREL', 'init_t_params')
    allocate (t_params%nrrel(2,t_params%lmmaxd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%nrrel))*kind(t_params%nrrel), 't_params%NRREL', 'init_t_params')
    allocate (t_params%ifunm1(t_params%lmxspd,t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%ifunm1))*kind(t_params%ifunm1), 't_params%IFUNM1', 'init_t_params')
    allocate (t_params%ititle(20,t_params%npotd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%ititle))*kind(t_params%ititle), 't_params%ITITLE', 'init_t_params')
    allocate (t_params%lmsp1(t_params%lmxspd,t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%lmsp1))*kind(t_params%lmsp1), 't_params%LMSP1', 'init_t_params')
    allocate (t_params%ntcell(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%ntcell))*kind(t_params%ntcell), 't_params%NTCELL', 'init_t_params')
    allocate (t_params%ixipol(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%ixipol))*kind(t_params%ixipol), 't_params%IXIPOL', 'init_t_params')
    allocate (t_params%irns(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%irns))*kind(t_params%irns), 't_params%IRNS', 'init_t_params')
    allocate (t_params%ifunm(t_params%natyp,t_params%lmxspd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%ifunm))*kind(t_params%ifunm), 't_params%IFUNM', 'init_t_params')
    allocate (t_params%llmsp(t_params%natyp,t_params%nfund), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%llmsp))*kind(t_params%llmsp), 't_params%LLMSP', 'init_t_params')
    allocate (t_params%lmsp(t_params%natyp,t_params%lmxspd), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%lmsp))*kind(t_params%lmsp), 't_params%LMSP', 'init_t_params')
    allocate (t_params%imt(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%imt))*kind(t_params%imt), 't_params%IMT', 'init_t_params')
    allocate (t_params%irc(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%irc))*kind(t_params%irc), 't_params%IRC', 'init_t_params')
    allocate (t_params%irmin(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%irmin))*kind(t_params%irmin), 't_params%IRMIN', 'init_t_params')
    allocate (t_params%irws(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%irws))*kind(t_params%irws), 't_params%IRWS', 'init_t_params')
    allocate (t_params%nfu(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%nfu))*kind(t_params%nfu), 't_params%NFU', 'init_t_params')
    allocate (t_params%hostimp(0:t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%hostimp))*kind(t_params%hostimp), 't_params%HOSTIMP', 'init_t_params')
    allocate (t_params%ilm_map(t_params%ngshd,3), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%ilm_map))*kind(t_params%ilm_map), 't_params%ILM_MAP', 'init_t_params')
    allocate (t_params%imaxsh(0:t_params%lmpot), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%imaxsh))*kind(t_params%imaxsh), 't_params%IMAXSH', 'init_t_params')
    allocate (t_params%npan_log_at(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%npan_log_at))*kind(t_params%npan_log_at), 't_params%NPAN_LOG_AT', 'init_t_params')
    allocate (t_params%npan_eq_at(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%npan_eq_at))*kind(t_params%npan_eq_at), 't_params%NPAN_EQ_AT', 'init_t_params')
    allocate (t_params%npan_tot(t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%npan_tot))*kind(t_params%npan_tot), 't_params%NPAN_EQ_AT', 'init_t_params')
    allocate (t_params%ipan_intervall(0:t_params%ntotd,t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%ipan_intervall))*kind(t_params%ipan_intervall), 't_params%IPAN_INTERVALL', 'init_t_params')
    allocate (t_params%nkcore(20,t_params%natyp), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%nkcore))*kind(t_params%nkcore), 't_params%NKCORE', 'init_t_params')
    allocate (t_params%kapcore(20,t_params%npotd), stat=i_stat) ! INTEGER
    call memocc(i_stat, product(shape(t_params%kapcore))*kind(t_params%kapcore), 't_params%KAPCORE', 'init_t_params')
    if (.not. allocated(t_params%qdos_atomselect)) then
      allocate (t_params%qdos_atomselect(t_params%natyp), stat=i_stat) ! INTEGER
      call memocc(i_stat, product(shape(t_params%qdos_atomselect))*kind(t_params%qdos_atomselect), 't_params%qdos_atomselect', 'init_t_params')
    end if
    ! -------------------------------------------------------------------------
    ! End of allocation of integer arrays
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! Allocate the logical arrays
    ! -------------------------------------------------------------------------
    allocate (t_params%symunitary(nsymaxd), stat=i_stat) ! LOGICALS
    call memocc(i_stat, product(shape(t_params%symunitary))*kind(t_params%symunitary), 't_params%SYMUNITARY', 'init_t_params')
    allocate (t_params%vacflag(2), stat=i_stat)
    call memocc(i_stat, product(shape(t_params%vacflag))*kind(t_params%vacflag), 't_params%VACFLAG', 'init_t_params')
    ! -------------------------------------------------------------------------
    ! End allocation of logical arrays
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! Allocate the character arrays
    ! -------------------------------------------------------------------------
    allocate (t_params%txc(6), stat=i_stat) ! CHARACTER*124
    call memocc(i_stat, product(shape(t_params%txc))*kind(t_params%txc), 't_params%TXC', 'init_t_params')

!    if (.not. allocated(t_params%testc)) then
!      allocate (t_params%testc(32), stat=i_stat)
!      call memocc(i_stat, product(shape(t_params%testc))*kind(t_params%testc), 't_params%TESTC', 'init_t_params')
!      allocate (t_params%optc(32), stat=i_stat) ! CHARACTER*8
!      call memocc(i_stat, product(shape(t_params%optc))*kind(t_params%optc), 't_params%OPTC', 'init_t_params')
!    end if
    ! -------------------------------------------------------------------------
    ! End allocation of character arrays
    ! -------------------------------------------------------------------------

    if (.not. allocated(t_params%nofks)) then
      allocate (t_params%bzkp(3,t_params%kpoibz,t_params%maxmesh), stat=i_stat) ! real (kind=dp)
      call memocc(i_stat, product(shape(t_params%bzkp))*kind(t_params%bzkp), 't_params%BZKP', 'init_t_params')
      allocate (t_params%volcub(t_params%kpoibz,t_params%maxmesh), stat=i_stat) ! real (kind=dp)
      call memocc(i_stat, product(shape(t_params%volcub))*kind(t_params%volcub), 't_params%VOLCUB', 'init_t_params')
      allocate (t_params%volbz(t_params%maxmesh), stat=i_stat) ! real (kind=dp)
      call memocc(i_stat, product(shape(t_params%volbz))*kind(t_params%volbz), 't_params%VOLBZ', 'init_t_params')
      allocate (t_params%nofks(t_params%maxmesh), stat=i_stat) ! integer
      call memocc(i_stat, product(shape(t_params%nofks))*kind(t_params%nofks), 't_params%NOFKS', 'init_t_params')
    end if
      
    ! -------------------------------------------------------------------------------
    ! Non-collinear magnetic field 
    ! -------------------------------------------------------------------------------
    if (.not. allocated(t_params%bfield%theta)) then
      allocate (t_params%bfield%theta(t_params%natyp), stat=i_stat)
      call memocc(i_stat, product(shape(t_params%bfield%theta))*kind(t_params%bfield%theta), 't_params%bfield%theta', 'init_t_params')
      allocate (t_params%bfield%phi(t_params%natyp), stat=i_stat)
      call memocc(i_stat, product(shape(t_params%bfield%phi))*kind(t_params%bfield%phi), 't_params%bfield%phi', 'init_t_params')
      allocate (t_params%bfield%bfield(t_params%natyp,3), stat=i_stat)
      call memocc(i_stat, product(shape(t_params%bfield%bfield))*kind(t_params%bfield%bfield), 't_params%bfield%bfield', 'init_t_params')
      allocate (t_params%bfield%bfield_strength(t_params%natyp), stat=i_stat)
      call memocc(i_stat, product(shape(t_params%bfield%bfield_strength))*kind(t_params%bfield%bfield_strength), 't_params%bfield%bfield_strength', 'init_t_params')
      allocate (t_params%bfield%bfield_constr(t_params%natyp,3), stat=i_stat)
      call memocc(i_stat, product(shape(t_params%bfield%bfield_constr))*kind(t_params%bfield%bfield_constr), 't_params%bfield%bfield_constr', 'init_t_params')
    end if

  end subroutine init_t_params



#ifdef CPP_MPI
  !-------------------------------------------------------------------------------
  !> Summary: Broadcast scalar parameters via MPI
  !> Author: Philipp Ruessmann
  !> Category: communication, KKRhost
  !> Deprecated: False !
  !> Broadcast scalar parameters via MPI. Broadcast scalar parameters, deal with arrays later
  !-------------------------------------------------------------------------------
  subroutine bcast_t_params_scalars(t_params)

    use :: mpi
    use :: mod_mympi, only: master

    implicit none

    type (type_params), intent (inout) :: t_params
    integer :: ierr
    ! integer :: myMPItype1
    ! integer, dimension(t_params%Nscalars) :: blocklen1
    ! integer, dimension(t_params%Nscalars) :: etype1
    ! integer(kind=MPI_ADDRESS_KIND) :: disp1(t_params%Nscalars), base

    ! !INTEGER
    ! call MPI_Get_address(t_params%IEMXD,      disp1(1), ierr)
    ! call MPI_Get_address(t_params%IRMIND,     disp1(2), ierr)
    ! call MPI_Get_address(t_params%IRM,        disp1(3), ierr)
    ! call MPI_Get_address(t_params%LMPOT,      disp1(4), ierr)
    ! call MPI_Get_address(t_params%NSPOTD,     disp1(5), ierr)
    ! call MPI_Get_address(t_params%NPOTD,      disp1(6), ierr)
    ! call MPI_Get_address(t_params%NATYP,      disp1(7), ierr)
    ! call MPI_Get_address(t_params%NEMBD1,     disp1(8), ierr)
    ! call MPI_Get_address(t_params%LMMAXD,     disp1(9), ierr)
    ! call MPI_Get_address(t_params%NAEZ,       disp1(10), ierr)
    ! call MPI_Get_address(t_params%IPAND,      disp1(11), ierr)
    ! call MPI_Get_address(t_params%NEMBD2,     disp1(12), ierr)
    ! call MPI_Get_address(t_params%NREF,       disp1(13), ierr)
    ! call MPI_Get_address(t_params%LMAX,       disp1(14), ierr)
    ! call MPI_Get_address(t_params%NCLEB,      disp1(15), ierr)
    ! call MPI_Get_address(t_params%NACLSD,     disp1(16), ierr)
    ! call MPI_Get_address(t_params%NCLSD,      disp1(17), ierr)
    ! call MPI_Get_address(t_params%LM2D,       disp1(18), ierr)
    ! call MPI_Get_address(t_params%LMAXD1,     disp1(19), ierr)
    ! call MPI_Get_address(t_params%NR,         disp1(20), ierr)
    ! call MPI_Get_address(t_params%NSHELD,     disp1(21), ierr)
    ! call MPI_Get_address(t_params%NSYMAXD,    disp1(22), ierr)
    ! call MPI_Get_address(t_params%NAEZDPD,    disp1(23), ierr)
    ! call MPI_Get_address(t_params%NATOMIMPD,  disp1(24), ierr)
    ! call MPI_Get_address(t_params%NOFGIJ,     disp1(25), ierr)
    ! call MPI_Get_address(t_params%NSPIND,     disp1(26), ierr)
    ! call MPI_Get_address(t_params%NSPINDD,    disp1(27), ierr)
    ! call MPI_Get_address(t_params%IRID,       disp1(28), ierr)
    ! call MPI_Get_address(t_params%NFUND,      disp1(29), ierr)
    ! call MPI_Get_address(t_params%NCELLD,     disp1(30), ierr)
    ! call MPI_Get_address(t_params%LMXSPD,     disp1(31), ierr)
    ! call MPI_Get_address(t_params%NGSHD,      disp1(32), ierr)
    ! call MPI_Get_address(t_params%KREL,       disp1(33), ierr)
    ! call MPI_Get_address(t_params%MMAXD,      disp1(34), ierr)
    ! call MPI_Get_address(t_params%IELAST,     disp1(35), ierr)
    ! call MPI_Get_address(t_params%NPOL,       disp1(36), ierr)
    ! call MPI_Get_address(t_params%NPNT1,      disp1(37), ierr)
    ! call MPI_Get_address(t_params%NPNT2,      disp1(38), ierr)
    ! call MPI_Get_address(t_params%NPNT3,      disp1(39), ierr)
    ! call MPI_Get_address(t_params%ITSCF,      disp1(40), ierr)
    ! call MPI_Get_address(t_params%SCFSTEPS,   disp1(41), ierr)
    ! call MPI_Get_address(t_params%LLY,        disp1(42), ierr)
    ! call MPI_Get_address(t_params%NSRA,       disp1(43), ierr)
    ! call MPI_Get_address(t_params%INS,        disp1(44), ierr)
    ! call MPI_Get_address(t_params%KORBIT,     disp1(45), ierr)
    ! call MPI_Get_address(t_params%KNOCO,      disp1(46), ierr)
    ! call MPI_Get_address(t_params%NINEQ,      disp1(47), ierr)
    ! call MPI_Get_address(t_params%KNOSPH,     disp1(48), ierr)
    ! call MPI_Get_address(t_params%NSPIN,      disp1(49), ierr)
    ! call MPI_Get_address(t_params%IRNSD,      disp1(50), ierr)
    ! call MPI_Get_address(t_params%NPRINC,     disp1(51), ierr)
    ! call MPI_Get_address(t_params%NCLS,       disp1(52), ierr)
    ! call MPI_Get_address(t_params%ICST,       disp1(53), ierr)
    ! call MPI_Get_address(t_params%IEND,       disp1(54), ierr)
    ! call MPI_Get_address(t_params%ICC,        disp1(55), ierr)
    ! call MPI_Get_address(t_params%IGF,        disp1(56), ierr)
    ! call MPI_Get_address(t_params%NLBASIS,    disp1(57), ierr)
    ! call MPI_Get_address(t_params%NRBASIS,    disp1(58), ierr)
    ! call MPI_Get_address(t_params%NCPA,       disp1(59), ierr)
    ! call MPI_Get_address(t_params%ITCPAMAX,   disp1(60), ierr)
    ! call MPI_Get_address(t_params%KMROT,      disp1(61), ierr)
    ! call MPI_Get_address(t_params%MAXMESH,    disp1(62), ierr)
    ! call MPI_Get_address(t_params%NSYMAT,     disp1(63), ierr)
    ! call MPI_Get_address(t_params%NATOMIMP,   disp1(64), ierr)
    ! call MPI_Get_address(t_params%INVMOD,     disp1(65), ierr)
    ! call MPI_Get_address(t_params%NQCALC,     disp1(66), ierr)
    ! call MPI_Get_address(t_params%INTERVX,    disp1(67), ierr)
    ! call MPI_Get_address(t_params%INTERVY,    disp1(68), ierr)
    ! call MPI_Get_address(t_params%INTERVZ,    disp1(69), ierr)
    ! call MPI_Get_address(t_params%LPOT,       disp1(70), ierr)
    ! call MPI_Get_address(t_params%NRIGHT,     disp1(71), ierr)
    ! call MPI_Get_address(t_params%NLEFT,      disp1(72), ierr)
    ! call MPI_Get_address(t_params%IMIX,       disp1(73), ierr)
    ! call MPI_Get_address(t_params%ITDBRY,     disp1(74), ierr)
    ! call MPI_Get_address(t_params%KPRE,       disp1(75), ierr)
    ! call MPI_Get_address(t_params%KSHAPE,     disp1(76), ierr)
    ! call MPI_Get_address(t_params%KTE,        disp1(77), ierr)
    ! call MPI_Get_address(t_params%KVMAD,      disp1(78), ierr)
    ! call MPI_Get_address(t_params%KXC,        disp1(79), ierr)
    ! call MPI_Get_address(t_params%ISHIFT,     disp1(80), ierr)
    ! call MPI_Get_address(t_params%KFORCE,     disp1(81), ierr)
    ! call MPI_Get_address(t_params%IDOLDAU,    disp1(82), ierr)
    ! call MPI_Get_address(t_params%ITRUNLDAU,  disp1(83), ierr)
    ! call MPI_Get_address(t_params%NTLDAU,     disp1(84), ierr)
    ! call MPI_Get_address(t_params%NPOLSEMI,   disp1(85), ierr)
    ! call MPI_Get_address(t_params%N1SEMI,     disp1(86), ierr)
    ! call MPI_Get_address(t_params%N2SEMI,     disp1(87), ierr)
    ! call MPI_Get_address(t_params%N3SEMI,     disp1(88), ierr)
    ! call MPI_Get_address(t_params%IESEMICORE, disp1(89), ierr)
    ! call MPI_Get_address(t_params%ITMPDIR,    disp1(90), ierr)
    ! call MPI_Get_address(t_params%ILTMP,      disp1(91), ierr)
    ! call MPI_Get_address(t_params%NCHEB,      disp1(92), ierr)
    ! call MPI_Get_address(t_params%NTOTD,      disp1(93), ierr)
    ! call MPI_Get_address(t_params%WLENGTH,    disp1(94), ierr)
    ! call MPI_Get_address(t_params%NTPERD,     disp1(95), ierr)
    ! !DOUBPLE PRECISION
    ! call MPI_Get_address(t_params%EBOTSEMI,   disp1(96), ierr)
    ! call MPI_Get_address(t_params%EMUSEMI,    disp1(97), ierr)
    ! call MPI_Get_address(t_params%TKSEMI,     disp1(98), ierr)
    ! call MPI_Get_address(t_params%FSEMICORE,  disp1(99), ierr)
    ! call MPI_Get_address(t_params%R_LOG,      disp1(100), ierr)
    ! call MPI_Get_address(t_params%EMIN,       disp1(101), ierr)
    ! call MPI_Get_address(t_params%EMAX,       disp1(102), ierr)
    ! call MPI_Get_address(t_params%TK,         disp1(103), ierr)
    ! call MPI_Get_address(t_params%EFERMI,     disp1(104), ierr)
    ! call MPI_Get_address(t_params%ALAT,       disp1(105), ierr)
    ! call MPI_Get_address(t_params%CPATOL,     disp1(106), ierr)
    ! call MPI_Get_address(t_params%MIXING,     disp1(107), ierr)
    ! call MPI_Get_address(t_params%QBOUND,     disp1(108), ierr)
    ! call MPI_Get_address(t_params%FCM,        disp1(109), ierr)
    ! call MPI_Get_address(t_params%LAMBDA_XC,  disp1(110), ierr)
    ! call MPI_Get_address(t_params%TOLRDIF,    disp1(111), ierr)
    ! call MPI_Get_address(t_params%EFOLD,      disp1(112), ierr)
    ! call MPI_Get_address(t_params%CHRGOLD,    disp1(113), ierr)
    ! !!complex (kind=dp)
    ! call MPI_Get_address(t_params%DELTAE,     disp1(114), ierr)
    ! !LOGICAL
    ! call MPI_Get_address(t_params%LINTERFACE, disp1(115), ierr)
    ! call MPI_Get_address(t_params%LRHOSYM,    disp1(116), ierr)
    ! !CHARACTER*10
    ! call MPI_Get_address(t_params%SOLVER,     disp1(117), ierr)
    ! !CHARACTER*80
    ! call MPI_Get_address(t_params%TMPDIR,     disp1(118), ierr)
    ! !INTEGER
    ! call MPI_Get_address(t_params%I1,         disp1(119), ierr)
    ! call MPI_Get_address(t_params%NMVECMAX,   disp1(120), ierr)
    ! call MPI_Get_address(t_params%ITAB,       disp1(121), ierr)
    ! real (kind=dp)
    ! call MPI_Get_address(t_params%LASTERR,      disp1(122), ierr)
    ! call MPI_Get_address(t_params%DENEF,        disp1(123), ierr)
    ! call MPI_Get_address(t_params%CHRGSEMICORE, disp1(124), ierr)
    ! INTEGER
    ! call MPI_Get_address(t_params%NPAN_LOG    , disp1(125), ierr)
    ! call MPI_Get_address(t_params%NPAN_EQ     , disp1(126), ierr)

    ! base  = disp1(1)
    ! disp1 = disp1 - base

    ! blocklen1(1:116)=1
    ! blocklen1(117)=10
    ! blocklen1(118)=80
    ! blocklen1(119:126)=1

    ! etype1(1:95) = MPI_INTEGER
    ! etype1(96:113) = MPI_DOUBLE_PRECISION
    ! etype1(114) = MPI_DOUBLE_COMPLEX
    ! etype1(115:116) = MPI_LOGICAL
    ! etype1(117:118) = MPI_CHARACTER
    ! etype1(119:121) = MPI_INTEGER
    ! etype1(122:124) = MPI_DOUBLE_PRECISION
    ! etype1(125:126) = MPI_INTEGER

    ! call MPI_Type_create_struct(t_params%Nscalars, blocklen1, disp1,  &
    ! etype1, myMPItype1, ierr)
    ! if(ierr/=MPI_SUCCESS) stop 'Problem in create_mpimask_t_params'

    ! call MPI_Type_commit(myMPItype1, ierr)
    ! if(ierr/=MPI_SUCCESS) stop 'error comiting create_mpimsk_t_params'

    ! call MPI_Bcast(t_params%Nscalars, 1, myMPItype1, master,MPI_COMM_WORLD, ierr)
    ! if(ierr/=MPI_SUCCESS) stop 'error brodcasting t_params'

    ! call MPI_Type_free(myMPItype1, ierr)

    ! somehow this parameter gets overlooked in the communication, possibly a but somewhere, but for now this workaround does the job
    ! call MPI_Bcast(t_params%NCHEB, 1, MPI_INTEGER, master,MPI_COMM_WORLD, ierr)


    ! try bcast for every variables instead:

    ! integer
    call mpi_bcast(t_params%iemxd, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%irmind, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%irm, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%lmpot, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nspotd, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%npotd, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%natyp, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nembd1, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%lmmaxd, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%naez, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ipand, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nembd2, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nref, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%lmax, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ncleb, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%naclsd, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nclsd, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%lm2d, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%lmaxd1, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nr, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nsheld, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%naezdpd, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%natomimpd, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nofgij, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nspind, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nspindd, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%irid, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nfund, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ncelld, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%lmxspd, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ngshd, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%krel, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%mmaxd, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ielast, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%npol, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%npnt1, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%npnt2, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%npnt3, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%itscf, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%scfsteps, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%lly, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nsra, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ins, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%korbit, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%knoco, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nineq, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%knosph, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nspin, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%irnsd, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nprinc, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ncls, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%icst, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%iend, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%icc, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%igf, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nlbasis, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nrbasis, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ncpa, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%itcpamax, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%kmrot, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%maxmesh, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nsymat, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%natomimp, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%invmod, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nqcalc, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%intervx, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%intervy, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%intervz, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%lpot, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nright, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nleft, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%imix, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%itdbry, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%kpre, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%kshape, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%kte, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%kvmad, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%kxc, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ishift, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%kforce, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%idoldau, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%itrunldau, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ntldau, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%npolsemi, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%n1semi, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%n2semi, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%n3semi, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%iesemicore, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%itmpdir, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%iltmp, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ncheb, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ntotd, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%wlength, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ntperd, 1, mpi_integer, master, mpi_comm_world, ierr)
    ! forgot to bcast kpoibz?
    call mpi_bcast(t_params%kpoibz, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%special_straight_mixing, 1, mpi_integer, master, mpi_comm_world, ierr)

    ! double precision
    call mpi_bcast(t_params%ebotsemi, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%emusemi, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%tksemi, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%fsemicore, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%r_log, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%emin, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%emax, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%tk, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%efermi, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%alat, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%cpatol, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%mixing, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%qbound, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%fcm, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%lambda_xc, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%tolrdif, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%efold, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%chrgold, 1, mpi_double_precision, master, mpi_comm_world, ierr)

    ! double complex
    call mpi_bcast(t_params%deltae, 1, mpi_double_complex, master, mpi_comm_world, ierr)

    ! logical
    call mpi_bcast(t_params%linterface, 1, mpi_logical, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%lrhosym, 1, mpi_logical, master, mpi_comm_world, ierr)

    ! character
    call mpi_bcast(t_params%solver, 10, mpi_character, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%tmpdir, 80, mpi_character, master, mpi_comm_world, ierr)

    ! integer
    call mpi_bcast(t_params%i1, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nmvecmax, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%itab, 1, mpi_integer, master, mpi_comm_world, ierr)

    ! double precision
    call mpi_bcast(t_params%lasterr, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%denef, 1, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%chrgsemicore, 1, mpi_double_precision, master, mpi_comm_world, ierr)

    ! integer
    call mpi_bcast(t_params%npan_log, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%npan_eq, 1, mpi_integer, master, mpi_comm_world, ierr)


    ! -------------------------------------------------------------------------
    ! Non-collinear magnetic field
    ! -------------------------------------------------------------------------
    ! logicals
    call mpi_bcast(t_params%bfield%lbfield, 1, mpi_logical, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%bfield%lbfield_constr, 1, mpi_logical, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%bfield%lbfield_all, 1, mpi_logical, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%bfield%lbfield_trans, 1, mpi_logical, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%bfield%lbfield_mt, 1, mpi_logical, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%bfield%ltorque, 1, mpi_logical, master, mpi_comm_world, ierr)

    ! integer
    call mpi_bcast(t_params%bfield%ibfield, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%bfield%ibfield_constr, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%bfield%itscf0, 1, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%bfield%itscf1, 1, mpi_integer, master, mpi_comm_world, ierr)


  end subroutine bcast_t_params_scalars

  !-------------------------------------------------------------------------------
  !> Summary: Broadcast arrays via MPI
  !> Author: Philipp Ruessmann
  !> Category: communication, KKRhost
  !> Deprecated: False
  !> Broadcast arrays via MP. Broadcast arrays from t_params
  !-------------------------------------------------------------------------------
  subroutine bcast_t_params_arrays(t_params)

    use :: mpi
    use :: mod_mympi, only: master
    implicit none

    type (type_params), intent (inout) :: t_params
    integer :: ierr

    ! -------------------------------------------------------------------------
    ! complex (kind=dp) arrays
    ! -------------------------------------------------------------------------
    call mpi_bcast(t_params%ez, t_params%iemxd, mpi_double_complex, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%wez, t_params%iemxd, mpi_double_complex, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%drotq, t_params%lmmaxd*t_params%lmmaxd*t_params%naez, mpi_double_complex, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%dsymll, t_params%lmmaxd*t_params%lmmaxd*nsymaxd, mpi_double_complex, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%lefttinvll, t_params%lmmaxd*t_params%lmmaxd*t_params%nembd1*t_params%nspindd*t_params%iemxd, mpi_double_complex, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%righttinvll, t_params%lmmaxd*t_params%lmmaxd*t_params%nembd1*t_params%nspindd*t_params%iemxd, mpi_double_complex, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%crel, t_params%lmmaxd*t_params%lmmaxd, mpi_double_complex, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%rc, t_params%lmmaxd*t_params%lmmaxd, mpi_double_complex, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%rrel, t_params%lmmaxd*t_params%lmmaxd, mpi_double_complex, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%srrel, 2*2*t_params%lmmaxd, mpi_double_complex, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%phildau, t_params%irm*t_params%natyp, mpi_double_complex, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%mvevi, t_params%natyp*3*t_params%nmvecmax, mpi_double_complex, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%mvevief, t_params%natyp*3*t_params%nmvecmax, mpi_double_complex, master, mpi_comm_world, ierr)

    ! -------------------------------------------------------------------------
    ! real (kind=dp) arrays
    ! -------------------------------------------------------------------------
    call mpi_bcast(t_params%vins, ((t_params%irm-t_params%irmind+1)*t_params%lmpot*t_params%nspotd), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%visp, (t_params%irm*t_params%npotd), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%vbc, 2, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%vtrel, (t_params%irm*t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%btrel, (t_params%irm*t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%socscale, (t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%drdirel, (t_params%irm*t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%r2drdirel, (t_params%irm*t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%rmrel, (t_params%irm*t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%cmomhost, (t_params%lmpot*t_params%nembd1), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ecore, (20*t_params%npotd), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%qmtet, (t_params%naez), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%qmphi, (t_params%naez), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%qmphitab, (t_params%naez*3), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%qmtettab, (t_params%naez*3), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%qmgamtab, (t_params%naez*3), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%zat, (t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%rmesh, (t_params%irm*t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%drdi, (t_params%irm*t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%rmtref, (t_params%nref), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%vref, (t_params%nref), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%cleb, (t_params%ncleb*2), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%rcls, (3*t_params%naclsd*t_params%nclsd), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%socscl, (t_params%lmaxd1*t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%cscl, (t_params%lmaxd1*t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%rbasis, (3*t_params%nembd2), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%rr, (3*(t_params%nr+1)), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%conc, (t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%rrot, (48*3*t_params%nsheld), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ratom, (3*t_params%nsheld), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%a, t_params%natyp, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%b, t_params%natyp, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%thetas, (t_params%irid*t_params%nfund*t_params%ncelld), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%rmt, (t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%rmtnew, (t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%rws, (t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%gsh, (t_params%ngshd), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%erefldau, (t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ueff, (t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%jeff, (t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%uldau, (t_params%mmaxd*t_params%mmaxd*t_params%mmaxd*t_params%mmaxd*t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%wldau, (t_params%mmaxd*t_params%mmaxd*t_params%nspind*t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%rpan_intervall, ((t_params%ntotd+1)*t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%rnew, (t_params%ntotd*(t_params%ncheb+1)*t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%thetasnew, (t_params%ntotd*(t_params%ncheb+1)*t_params%nfund*t_params%ncelld), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%rho2ns, (t_params%irm*t_params%lmpot*t_params%natyp*2), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%r2nef, (t_params%irm*t_params%lmpot*t_params%natyp*2), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%rhoc, (t_params%irm*t_params%npotd), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%denefat, (t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%espv, (t_params%lmaxd1+1)*t_params%npotd, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%edc, (t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%eu, (t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%rhoorb, (t_params%irm*t_params%krel+(1-t_params%krel)*t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ecorerel, (t_params%krel*20+(1-t_params%krel)*t_params%npotd), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%theta, t_params%natyp, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%phi, t_params%natyp, mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%rclsimp, 3*t_params%natomimpd, mpi_double_precision, master, mpi_comm_world, ierr)

    ! -------------------------------------------------------------------------
    ! INTEGERS arrays
    ! -------------------------------------------------------------------------
    call mpi_bcast(t_params%lopt, (t_params%natyp), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%itldau, (t_params%natyp), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%irshift, (t_params%natyp), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%jwsrel, (t_params%natyp), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%zrel, (t_params%natyp), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%lcore, (20*t_params%npotd), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ncore, (t_params%npotd), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ipan, (t_params%natyp), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ircut, ((t_params%ipand+1)*t_params%natyp), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%jend, (t_params%lmpot*(t_params%lmax+1)*(t_params%lmax+1)), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%icleb, (t_params%ncleb*4), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%atom, (t_params%naclsd*t_params%nembd2), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%cls, (t_params%nembd2), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nacls, (t_params%nclsd), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%loflm, (t_params%lm2d), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ezoa, (t_params%naclsd*t_params%nembd2), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%kaoez, (t_params%natyp*t_params%nembd2), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%iqat, (t_params%natyp), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%icpa, (t_params%naez), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%noq, (t_params%naez), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%kmesh, (t_params%iemxd), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nshell, ((t_params%nsheld+1)), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nsh1, (t_params%nsheld), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nsh2, (t_params%nsheld), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ijtabcalc, (t_params%nofgij), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ijtabcalc_i, (t_params%nofgij), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ijtabsym, (t_params%nofgij), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ijtabsh, (t_params%nofgij), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ish, (t_params%nsheld*2*nsymaxd), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%jsh, (t_params%nsheld*2*nsymaxd), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%iqcalc, (t_params%naez), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%icheck, (t_params%naezdpd*t_params%naezdpd), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%atomimp, (t_params%natomimpd), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%refpot, (t_params%nembd2), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%irrel, (2*2*t_params%lmmaxd), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nrrel, (2*t_params%lmmaxd), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ifunm1, (t_params%lmxspd*t_params%natyp), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ititle, (20*t_params%npotd), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%lmsp1, (t_params%lmxspd*t_params%natyp), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ntcell, (t_params%natyp), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ixipol, (t_params%natyp), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%irns, (t_params%natyp), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ifunm, (t_params%natyp*t_params%lmxspd), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%llmsp, (t_params%natyp*t_params%nfund), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%lmsp, (t_params%natyp*t_params%lmxspd), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%imt, (t_params%natyp), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%irc, (t_params%natyp), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%irmin, (t_params%natyp), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%irws, (t_params%natyp), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nfu, (t_params%natyp), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%hostimp, ((t_params%natyp+1)), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ilm_map, (t_params%ngshd*3), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%imaxsh, ((t_params%lmpot+1)), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%npan_log_at, (t_params%natyp), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%npan_eq_at, (t_params%natyp), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%npan_tot, (t_params%natyp), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%ipan_intervall, ((t_params%ntotd+1)*t_params%natyp), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%nkcore, (20*t_params%natyp), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%kapcore, (20*t_params%npotd), mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%qdos_atomselect, (t_params%natyp), mpi_integer, master, mpi_comm_world, ierr)

    ! -------------------------------------------------------------------------
    ! LOGICAL arrays
    ! -------------------------------------------------------------------------
    call mpi_bcast(t_params%symunitary, nsymaxd, mpi_logical, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%vacflag, 2, mpi_logical, master, mpi_comm_world, ierr)

    ! -------------------------------------------------------------------------
    ! CHARACTER arrays
    ! -------------------------------------------------------------------------
    call mpi_bcast(t_params%txc, 6*124, & ! 6 entries of length 124
      mpi_character, master, mpi_comm_world, ierr) ! CHARACTER*124
    !call mpi_bcast(t_params%testc, 32*8, & ! 32 entries of length 8
    !  mpi_character, master, mpi_comm_world, ierr) ! CHARACTER*8
    !call mpi_bcast(t_params%optc, 32*8, & ! 32 entries of length 8
    !  mpi_character, master, mpi_comm_world, ierr) ! CHARACTER*8

    ! -------------------------------------------------------------------------
    ! K-points arrays
    ! -------------------------------------------------------------------------
    call mpi_bcast(t_params%bzkp, (3*t_params%kpoibz*t_params%maxmesh), mpi_double_precision, master, mpi_comm_world, ierr) ! real (kind=dp)
    call mpi_bcast(t_params%volcub, (t_params%kpoibz*t_params%maxmesh), mpi_double_precision, master, mpi_comm_world, ierr) ! real (kind=dp)
    call mpi_bcast(t_params%volbz, (t_params%maxmesh), mpi_double_precision, master, mpi_comm_world, ierr) ! real (kind=dp)
    call mpi_bcast(t_params%nofks, (t_params%maxmesh), mpi_integer, master, mpi_comm_world, ierr) ! integer
    

    ! -------------------------------------------------------------------------
    ! Non-collinear magnetic field
    ! -------------------------------------------------------------------------
    call mpi_bcast(t_params%bfield%bfield, (t_params%natyp*3), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%bfield%bfield_constr, (t_params%natyp*3), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%bfield%bfield_strength, (t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%bfield%theta, (t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    call mpi_bcast(t_params%bfield%phi, (t_params%natyp), mpi_double_precision, master, mpi_comm_world, ierr)
    
    call mpi_bcast(t_params%ntcell, (t_params%natyp), mpi_integer, master, mpi_comm_world, ierr)

  end subroutine bcast_t_params_arrays
#endif

  !-------------------------------------------------------------------------------
  !> Summary: Set the values of the `t_params` scalars with the input values
  !> Author: Philipp Ruessmann
  !> Category: initialization, communication, KKRhost
  !> Deprecated: False
  !> Set the values of the `t_params` scalars with the input values
  !-------------------------------------------------------------------------------
  subroutine fill_t_params_scalars(iemxd,irmind,irm,lmpot,nspotd,npotd,natyp,nembd1,&
    lmmaxd,naez,ipand,nembd2,nref,lmax,ncleb,naclsd,nclsd,lm2d,lmaxd1,nr,nsheld,    &
    naezdpd,natomimpd,nofgij,nspind,nspindd,irid,nfund,ncelld,lmxspd,ngshd,         &
    krel,mmaxd,ielast,npol,npnt1,npnt2,npnt3,itscf,scfsteps,lly,nsra,ins,nineq,     &
    nspin,ncls,icst,iend,icc,igf,nlbasis,nrbasis,ncpa,itcpamax,kmrot,maxmesh,nsymat,&
    natomimp,invmod,nqcalc,intervx,intervy,intervz,lpot,nright,nleft,imix,itdbry,   &
    kpre,kshape,kte,kvmad,kxc,ishift,kforce,idoldau,itrunldau,ntldau,npolsemi,      &
    n1semi,n2semi,n3semi,iesemicore,ebotsemi,emusemi,tksemi,fsemicore,r_log,emin,   &
    emax,tk,efermi,alat,cpatol,mixing,qbound,fcm,lambda_xc,tolrdif,linterface,      &
    lrhosym,solver,tmpdir,itmpdir,iltmp,ntotd,ncheb,deltae,special_straight_mixing, &
    t_params)
    ! fill scalars into t_params
    implicit none

    type (type_params), intent (inout) :: t_params
    ! ..
    ! .. Scalar arguments
    integer, intent (in) :: nr     !! Number of real space vectors rr
    integer, intent (in) :: irm    !! Maximum number of radial points
    integer, intent (in) :: ins    !! 0 (MT), 1(ASA), 2(Full Potential)
    integer, intent (in) :: lly    !! LLY!> 0 : apply Lloyds formula
    integer, intent (in) :: icc    !! Enables the calculation of off-diagonal elements of the GF.(0=SCF/DOS; 1=cluster; -1=custom)
    integer, intent (in) :: igf    !! Do not print or print (0/1) the KKRFLEX_* files
    integer, intent (in) :: kte    !! Calculation of the total energy On/Off (1/0)
    integer, intent (in) :: kxc    !! Type of xc-potential 0=vBH 1=MJW 2=VWN 3=PW91
    integer, intent (in) :: nref   !! Number of diff. ref. potentials
    integer, intent (in) :: lm2d   !! (2*LMAX+1)**2
    integer, intent (in) :: krel   !! Switch for non-relativistic/relativistic (0/1) program. Attention: several other parameters depend explicitly on KREL, they are set automatically Used for Dirac solver in ASA
    integer, intent (in) :: irid   !! Shape functions parameters in non-spherical part
    integer, intent (in) :: lmax   !! Maximum l component in wave function expansion
    integer, intent (in) :: ncls   !! Number of reference clusters
    integer, intent (in) :: icst   !! Number of Born approximation
    integer, intent (in) :: iend   !! Number of nonzero gaunt coefficients
    integer, intent (in) :: nsra
    integer, intent (in) :: lpot   !! Maximum l component in potential expansion
    integer, intent (in) :: kpre
    integer, intent (in) :: imix   !! Type of mixing scheme used (0=straight, 4=Broyden 2nd, 5=Anderson)
    integer, intent (in) :: ncpa   !! NCPA = 0/1 CPA flag
    integer, intent (in) :: naez   !! Number of atoms in unit cell
    integer, intent (in) :: npol   !! Number of Matsubara Poles (EMESHT)
    integer, intent (in) :: npnt1  !! number of E points (EMESHT) for the contour integration
    integer, intent (in) :: npnt2  !! number of E points (EMESHT) for the contour integration
    integer, intent (in) :: npnt3  !! number of E points (EMESHT) for the contour integration
    integer, intent (in) :: itscf
    integer, intent (in) :: iemxd  !! Dimension for energy-dependent arrays
    integer, intent (in) :: lmpot  !! (LPOT+1)**2
    integer, intent (in) :: npotd  !! (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP)
    integer, intent (in) :: natyp  !! Number of kinds of atoms in unit cell
    integer, intent (in) :: ipand  !! Number of panels in non-spherical part
    integer, intent (in) :: ncleb  !! Number of Clebsch-Gordon coefficients
    integer, intent (in) :: nclsd  !! Maximum number of different TB-clusters
    integer, intent (in) :: nfund  !! Shape functions parameters in non-spherical part
    integer, intent (in) :: ngshd  !! Shape functions parameters in non-spherical part
    integer, intent (in) :: mmaxd  !! 2*LMAX+1
    integer, intent (in) :: nineq  !! Number of ineq. positions in unit cell
    integer, intent (in) :: nspin  !! Counter for spin directions
    integer, intent (in) :: kmrot  !! 0: no rotation of the magnetisation; 1: individual rotation of the magnetisation for every site
    integer, intent (in) :: ntotd
    integer, intent (in) :: ncheb  !! Number of Chebychev pannels for the new solver
    integer, intent (in) :: kvmad
    integer, intent (in) :: iltmp
    integer, intent (in) :: nleft  !! Number of repeated basis for left host to get converged electrostatic potentials
    integer, intent (in) :: nright !! Number of repeated basis for right host to get converged electrostatic potentials
    integer, intent (in) :: itdbry !! Number of SCF steps to remember for the Broyden mixing
    integer, intent (in) :: kshape !! Exact treatment of WS cell
    integer, intent (in) :: ishift
    integer, intent (in) :: kforce !! Calculation of the forces
    integer, intent (in) :: irmind !! IRM-IRNSD
    integer, intent (in) :: nspotd !! Number of potentials for storing non-sph. potentials
    integer, intent (in) :: nembd1 !! NEMB+1
    integer, intent (in) :: lmmaxd !! (KREL+KORBIT+1)(LMAX+1)^2
    integer, intent (in) :: nembd2
    integer, intent (in) :: naclsd !! Maximum number of atoms in a TB-cluster
    integer, intent (in) :: lmaxd1
    integer, intent (in) :: nsheld !! Number of blocks of the GF matrix that need to be calculated (NATYPD + off-diagonals in case of impurity)
    integer, intent (in) :: nofgij !! number of GF pairs IJ to be calculated as determined from IJTABCALC<>0
    integer, intent (in) :: nspind !! KREL+(1-KREL)*(NSPIN+1)
    integer, intent (in) :: ncelld !! Number of cells (shapes) in non-spherical part
    integer, intent (in) :: lmxspd !! (2*LPOT+1)**2
    integer, intent (in) :: ielast
    integer, intent (in) :: nsymat
    integer, intent (in) :: invmod !! Inversion scheme
    integer, intent (in) :: nqcalc
    integer, intent (in) :: n1semi !! Number of energy points for the semicore contour
    integer, intent (in) :: n2semi !! Number of energy points for the semicore contour
    integer, intent (in) :: n3semi !! Number of energy points for the semicore contour
    integer, intent (in) :: ntldau !! number of atoms on which LDA+U is applied
    integer, intent (in) :: itmpdir
    integer, intent (in) :: maxmesh
    integer, intent (in) :: naezdpd
    integer, intent (in) :: nspindd !! NSPIND-KORBIT
    integer, intent (in) :: nlbasis !! Number of basis layers of left host (repeated units)
    integer, intent (in) :: nrbasis !! Number of basis layers of right host (repeated units)
    integer, intent (in) :: intervx !! Number of intervals in x-direction for k-net in IB of the BZ
    integer, intent (in) :: intervy !! Number of intervals in y-direction for k-net in IB of the BZ
    integer, intent (in) :: intervz !! Number of intervals in z-direction for k-net in IB of the BZ
    integer, intent (in) :: idoldau !! flag to perform LDA+U
    integer, intent (in) :: scfsteps !! number of scf iterations
    integer, intent (in) :: itcpamax !! Max. number of CPA iterations
    integer, intent (in) :: natomimp !! Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
    integer, intent (in) :: npolsemi !! Number of poles for the semicore contour
    integer, intent (in) :: natomimpd !! Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
    integer, intent (in) :: itrunldau !! Iteration index for LDA+U
    integer, intent (in) :: iesemicore
    integer, intent (in) :: special_straight_mixing
    real (kind=dp), intent (in) :: tk !! Temperature
    real (kind=dp), intent (in) :: fcm
    real (kind=dp), intent (in) :: emin !! Energies needed in EMESHT
    real (kind=dp), intent (in) :: emax !! Energies needed in EMESHT
    real (kind=dp), intent (in) :: alat !! Lattice constant in a.u.
    real (kind=dp), intent (in) :: r_log
    real (kind=dp), intent (in) :: tksemi !! Temperature of semi-core contour
    real (kind=dp), intent (in) :: efermi !! Fermi energy
    real (kind=dp), intent (in) :: cpatol !! Convergency tolerance for CPA-cycle
    real (kind=dp), intent (in) :: mixing !! Magnitude of the mixing parameter
    real (kind=dp), intent (in) :: qbound !! Convergence parameter for the potential
    real (kind=dp), intent (in) :: emusemi
    real (kind=dp), intent (in) :: tolrdif !! Tolerance for r<tolrdif (a.u.) to handle vir. atoms
    real (kind=dp), intent (in) :: ebotsemi
    real (kind=dp), intent (in) :: fsemicore
    real (kind=dp), intent (in) :: lambda_xc !! Scale magnetic moment (0! Lambda_XC! 1, 0=zero moment, 1= full moment)
    complex (kind=dp), intent (in) :: deltae !! Energy difference for numerical derivative
    logical, intent (in) :: lrhosym
    logical, intent (in) :: linterface !! If True a matching with semi-inifinite surfaces must be performed
    character (len=10), intent (in) :: solver                           !! Type of solver

    character (len=80), intent (in) :: tmpdir
    ! ..
    ! fill scalars:
    ! Integer
    t_params%nr = nr
    t_params%irm = irm
    t_params%lly = lly
    t_params%ins = ins
    t_params%icc = icc
    t_params%igf = igf
    t_params%kte = kte
    t_params%kxc = kxc
    t_params%lm2d = lm2d
    t_params%lpot = lpot
    t_params%imix = imix
    t_params%kpre = kpre
    t_params%nsra = nsra
    t_params%nref = nref
    t_params%lmax = lmax
    t_params%ncls = ncls
    t_params%icst = icst
    t_params%iend = iend
    t_params%ncpa = ncpa
    t_params%krel = krel
    t_params%irid = irid
    t_params%naez = naez
    t_params%npol = npol
    t_params%npnt1 = npnt1
    t_params%npnt2 = npnt2
    t_params%npnt3 = npnt3
    t_params%npotd = npotd
    t_params%natyp = natyp
    t_params%itscf = itscf
    t_params%nspin = nspin
    t_params%nineq = nineq
    t_params%iltmp = iltmp
    t_params%ncheb = ncheb
    t_params%ntotd = ntotd
    t_params%ncheb = ncheb
    t_params%lmpot = lmpot
    t_params%kmrot = kmrot
    t_params%kvmad = kvmad
    t_params%ngshd = ngshd
    t_params%mmaxd = mmaxd
    t_params%iemxd = iemxd
    t_params%lmpot = lmpot
    t_params%ipand = ipand
    t_params%ncleb = ncleb
    t_params%nclsd = nclsd
    t_params%nfund = nfund
    t_params%nleft = nleft
    t_params%nright = nright
    t_params%itdbry = itdbry
    t_params%kshape = kshape
    t_params%ishift = ishift
    t_params%kforce = kforce
    t_params%irmind = irmind
    t_params%nspotd = nspotd
    t_params%nembd1 = nembd1
    t_params%lmmaxd = lmmaxd
    t_params%nembd2 = nembd2
    t_params%naclsd = naclsd
    t_params%lmaxd1 = lmaxd1
    t_params%nsheld = nsheld
    t_params%nofgij = nofgij
    t_params%nspind = nspind
    t_params%ncelld = ncelld
    t_params%lmxspd = lmxspd
    t_params%ielast = ielast
    t_params%n1semi = n1semi
    t_params%n2semi = n2semi
    t_params%n3semi = n3semi
    t_params%invmod = invmod
    t_params%nqcalc = nqcalc
    t_params%nsymat = nsymat
    t_params%ntldau = ntldau
    t_params%nlbasis = nlbasis
    t_params%nrbasis = nrbasis
    t_params%maxmesh = maxmesh
    t_params%naezdpd = naezdpd
    t_params%nspindd = nspindd
    t_params%intervx = intervx
    t_params%intervy = intervy
    t_params%intervz = intervz
    t_params%idoldau = idoldau
    t_params%itmpdir = itmpdir
    t_params%scfsteps = scfsteps
    t_params%itcpamax = itcpamax
    t_params%natomimp = natomimp
    t_params%npolsemi = npolsemi
    t_params%natomimpd = natomimpd
    t_params%itrunldau = itrunldau
    t_params%iesemicore = iesemicore
    t_params%special_straight_mixing = special_straight_mixing
    ! Double precision
    t_params%tk = tk
    t_params%fcm = fcm
    t_params%emin = emin
    t_params%emax = emax
    t_params%alat = alat
    t_params%r_log = r_log
    t_params%efermi = efermi
    t_params%cpatol = cpatol
    t_params%mixing = mixing
    t_params%qbound = qbound
    t_params%tksemi = tksemi
    t_params%emusemi = emusemi
    t_params%tolrdif = tolrdif
    t_params%ebotsemi = ebotsemi
    t_params%fsemicore = fsemicore
    t_params%lambda_xc = lambda_xc
    ! Double complex
    t_params%deltae = deltae
    ! Logical
    t_params%lrhosym = lrhosym
    t_params%linterface = linterface
    ! Character
    t_params%solver = solver
    t_params%tmpdir = tmpdir

    t_params%nmvecmax = 4

  end subroutine fill_t_params_scalars

  !-------------------------------------------------------------------------------
  !> Summary: Set the values of the t_params arrays with the input values of the arrays
  !> Author: Who wrote this subroutine
  !> Category: initialization, communication, KKRhost
  !> Deprecated: False
  !> Set the values of the t_params arrays with the input values of the arrays.
  !> Fill arrays after they have been allocated in `init_t_params`
  !-------------------------------------------------------------------------------
  subroutine fill_t_params_arrays(t_params,iemxd,lmmaxd,naez,nembd1,nspindd,        &
    irmind,irm,lmpot,nspotd,npotd,natyp,nr,nembd2,nref,ncleb,nclsd,naclsd,nsheld,   &
    ngshd,nfund,irid,ncelld,mmaxd,lm2d,lmxspd,lmaxd1,nspind,ntotd,ncheb,ipand,lmax, &
    nofgij,naezdpd,natomimpd,ez,wez,drotq,dsymll,lefttinvll,righttinvll,crel,rc,    &
    rrel,srrel,phildau,vins,visp,vbc,vtrel,btrel,socscale,drdirel,r2drdirel,rmrel,  &
    cmomhost,ecore,qmtet,qmphi,qmphitab,qmtettab,qmgamtab,zat,r,drdi,rmtref,vref,   &
    cleb,rcls,socscl,cscl,rbasis,rr,conc,rrot,ratom,a,b,thetas,rmt,rmtnew,rws,gsh,  &
    erefldau,ueff,jeff,uldau,wldau,rpan_intervall,rnew,thetasnew,lopt,itldau,       &
    irshift,jwsrel,zrel,lcore,ncore,ipan,ircut,jend,icleb,atom,cls,nacls,loflm,ezoa,&
    kaoez,iqat,icpa,noq,kmesh,nshell,nsh1,nsh2,ijtabcalc,ijtabcalc_i,ijtabsym,      &
    ijtabsh,ish,jsh,iqcalc,icheck,atomimp,refpot,irrel,nrrel,ifunm1,ititle,lmsp1,   &
    ntcell,ixipol,irns,ifunm,llmsp,lmsp,imt,irc,irmin,irws,nfu,hostimp,ilm_map,     &
    imaxsh,npan_log,npan_eq,npan_log_at,npan_eq_at,npan_tot,ipan_intervall,         &
    symunitary,vacflag,txc,rclsimp,krel)
    ! ..
    implicit none

    type (type_params), intent (inout) :: t_params
    ! ..
    ! .. Scalars for array dimensions
    integer, intent (in) :: nr     !! Number of real space vectors rr
    integer, intent (in) :: irm    !! Maximum number of radial points
    integer, intent (in) :: lmax   !! Maximum l component in wave function expansion
    integer, intent (in) :: nref   !! Number of diff. ref. potentials
    integer, intent (in) :: naez   !! Number of atoms in unit cell
    integer, intent (in) :: irid   !! Shape functions parameters in non-spherical part
    integer, intent (in) :: lm2d   !! (2*LMAX+1)**2
    integer, intent (in) :: ncleb  !! Number of Clebsch-Gordon coefficients
    integer, intent (in) :: nclsd  !! Maximum number of different TB-clusters
    integer, intent (in) :: npotd  !! (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP)
    integer, intent (in) :: natyp  !! Number of kinds of atoms in unit cell
    integer, intent (in) :: lmpot  !! (LPOT+1)**2
    integer, intent (in) :: ngshd  !! Shape functions parameters in non-spherical part
    integer, intent (in) :: nfund  !! Shape functions parameters in non-spherical part
    integer, intent (in) :: mmaxd  !! 2*LMAX+1
    integer, intent (in) :: ntotd
    integer, intent (in) :: ncheb  !! Number of Chebychev pannels for the new solver
    integer, intent (in) :: ipand  !! Number of panels in non-spherical part
    integer, intent (in) :: iemxd  !! Dimension for energy-dependent arrays
    integer, intent (in) :: lmmaxd !! (KREL+KORBIT+1)(LMAX+1)^2
    integer, intent (in) :: nembd1 !! NEMB+1
    integer, intent (in) :: irmind !! IRM-IRNSD
    integer, intent (in) :: nspotd !! Number of potentials for storing non-sph. potentials
    integer, intent (in) :: nembd2
    integer, intent (in) :: naclsd !! Maximum number of atoms in a TB-cluster
    integer, intent (in) :: nsheld !! Number of blocks of the GF matrix that need to be calculated (NATYPD + off-diagonals in case of impurity)
    integer, intent (in) :: ncelld !! Number of cells (shapes) in non-spherical part
    integer, intent (in) :: lmxspd !! (2*LPOT+1)**2
    integer, intent (in) :: lmaxd1
    integer, intent (in) :: nofgij !! number of GF pairs IJ to be calculated as determined from IJTABCALC<>0
    integer, intent (in) :: nspind !! KREL+(1-KREL)*(NSPIN+1)
    integer, intent (in) :: nspindd !! NSPIND-KORBIT
    integer, intent (in) :: naezdpd
    integer, intent (in) :: npan_eq
    integer, intent (in) :: npan_log
    integer, intent (in) :: natomimpd !! Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
    integer, intent (in) :: krel
    ! .. Array arguments
    complex (kind=dp), dimension (iemxd), intent (in) :: ez
    complex (kind=dp), dimension (iemxd), intent (in) :: wez
    complex (kind=dp), dimension (lmmaxd, lmmaxd), intent (in) :: rc !! NREL REAL spher. harm. > CMPLX. spher. harm. NREL CMPLX. spher. harm. > REAL spher. harm.
    complex (kind=dp), dimension (lmmaxd, lmmaxd), intent (in) :: crel !! Non-relat. CMPLX. spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. CMPLX. spher. harm.
    complex (kind=dp), dimension (lmmaxd, lmmaxd), intent (in) :: rrel !! Non-relat. REAL spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. REAL spher. harm.
    complex (kind=dp), dimension (irm, natyp), intent (in) :: phildau
    complex (kind=dp), dimension (lmmaxd, lmmaxd, naez), intent (in) :: drotq !! Rotation matrices to change between LOCAL/GLOBAL frame of reference for magnetisation!> Oz or noncollinearity
    complex (kind=dp), dimension (2, 2, lmmaxd), intent (in) :: srrel
    complex (kind=dp), dimension (lmmaxd, lmmaxd, nsymaxd), intent (in) :: dsymll
    complex (kind=dp), dimension (lmmaxd, lmmaxd, nembd1, nspindd, iemxd), intent (in) :: lefttinvll
    complex (kind=dp), dimension (lmmaxd, lmmaxd, nembd1, nspindd, iemxd), intent (in) :: righttinvll

    real (kind=dp), dimension (natyp), intent (in) :: a !! Constants for exponential R mesh
    real (kind=dp), dimension (natyp), intent (in) :: b !! Constants for exponential R mesh
    real (kind=dp), dimension (2), intent (in) :: vbc !! Potential constants
    real (kind=dp), dimension (natyp), intent (in) :: rmt !! Muffin-tin radius of true system
    real (kind=dp), dimension (natyp), intent (in) :: rws !! Wigner Seitz radius
    real (kind=dp), dimension (ngshd), intent (in) :: gsh
    real (kind=dp), dimension (natyp), intent (in) :: zat !! Nuclear charge
    real (kind=dp), dimension (natyp), intent (in) :: ueff !! input U parameter for each atom
    real (kind=dp), dimension (natyp), intent (in) :: jeff !! input J parameter for each atom
    real (kind=dp), dimension (natyp), intent (in) :: conc !! Concentration of a given atom
    real (kind=dp), dimension (nref), intent (in) :: vref
    real (kind=dp), dimension (naez), intent (in) :: qmtet !! \( \theta\) angle of the agnetization with respect to the z-axis
    real (kind=dp), dimension (naez), intent (in) :: qmphi !! \( \phi\) angle of the agnetization with respect to the z-axis
    real (kind=dp), dimension (nref), intent (in) :: rmtref !! Muffin-tin radius of reference system
    real (kind=dp), dimension (natyp), intent (in) :: rmtnew !! Adapted muffin-tin radius
    real (kind=dp), dimension (natyp), intent (in) :: erefldau !! the energies of the projector's wave functions (REAL)
    real (kind=dp), dimension (natyp), intent (in) :: socscale !! Spin-orbit scaling

    real (kind=dp), dimension (irm, natyp), intent (in) :: r !! Radial mesh ( in units a Bohr)
    real (kind=dp), dimension (3, 0:nr), intent (in) :: rr !! Set of real space vectors (in a.u.)
    real (kind=dp), dimension (ncleb, 2), intent (in) :: cleb !! GAUNT coefficients (GAUNT)
    real (kind=dp), dimension (irm, natyp), intent (in) :: drdi !! Derivative dr/di
    real (kind=dp), dimension (irm, npotd), intent (in) :: visp !! Spherical part of the potential
    ! real (kind=dp), dimension(LMAXD1,NATYP), intent(in)          :: CSCL       !! Speed of light scaling
    real (kind=dp), dimension (krel*lmax+1, krel*natyp+(1-krel)), intent (inout) :: cscl !! Speed of light scaling
    real (kind=dp), dimension (ntotd*(ncheb+1), natyp), intent (in) :: rnew
    ! real (kind=dp), dimension(IRM,NATYP), intent(in)             :: VTREL      !! potential (spherical part)
    real (kind=dp), dimension (irm*krel+(1-krel), natyp), intent (inout) :: vtrel !! potential (spherical part)
    ! real (kind=dp), dimension(IRM,NATYP), intent(in)             :: BTREL      !! magnetic field
    real (kind=dp), dimension (irm*krel+(1-krel), natyp), intent (inout) :: btrel !! magnetic field
    ! real (kind=dp), dimension(IRM,NATYP), intent(in)             :: RMREL      !! radial mesh
    real (kind=dp), dimension (irm*krel+(1-krel), natyp), intent (inout) :: rmrel !! radial mesh
    real (kind=dp), dimension (20, npotd), intent (in) :: ecore !! Core energies
    real (kind=dp), dimension (3, nsheld), intent (in) :: ratom
    real (kind=dp), dimension (3, nembd2), intent (in) :: rbasis !! Position of atoms in the unit cell in units of bravais vectors
    ! real (kind=dp), dimension(LMAXD1,NATYP), intent(in)          :: SOCSCL
    real (kind=dp), dimension (krel*lmax+1, krel*natyp+(1-krel)), intent (inout) :: socscl
    real (kind=dp), dimension (3, natomimpd), intent (in) :: rclsimp
    ! real (kind=dp), dimension(IRM,NATYP), intent(in)             :: DRDIREL    !! derivative of radial mesh
    real (kind=dp), dimension (irm*krel+(1-krel), natyp), intent (in) :: drdirel
    real (kind=dp), dimension (naez, 3), intent (in) :: qmphitab
    real (kind=dp), dimension (naez, 3), intent (in) :: qmtettab
    real (kind=dp), dimension (naez, 3), intent (in) :: qmgamtab
    real (kind=dp), dimension (lmpot, nembd1), intent (in) :: cmomhost !! Charge moments of each atom of the (left/right) host
    ! real (kind=dp), dimension(IRM,NATYP), intent(in)             :: R2DRDIREL  !! \( r^2 \frac{\partial}{\partial \mathbf{r}}\frac{\partial}{\partial i}\) (r**2 * drdi)
    real (kind=dp), dimension (irm*krel+(1-krel), natyp), intent (in) :: r2drdirel
    real (kind=dp), dimension (0:ntotd, natyp), intent (in) :: rpan_intervall

    real (kind=dp), dimension (3, naclsd, nclsd), intent (in) :: rcls !! Real space position of atom in cluster
    real (kind=dp), dimension (48, 3, nsheld), intent (in) :: rrot
    real (kind=dp), dimension (irmind:irm, lmpot, nspotd), intent (in) :: vins !! Non-spherical part of the potential
    real (kind=dp), dimension (irid, nfund, ncelld), intent (in) :: thetas !! shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
    real (kind=dp), dimension (ntotd*(ncheb+1), nfund, ncelld), intent (in) :: thetasnew
    real (kind=dp), dimension (mmaxd, mmaxd, nspind, natyp), intent (in) :: wldau !! potential matrix
    real (kind=dp), dimension (mmaxd, mmaxd, mmaxd, mmaxd, natyp), intent (in) :: uldau !! calculated Coulomb matrix elements (EREFLDAU)
    ! ..
    integer, dimension (naez), intent (in) :: noq !! Number of diff. atom types located
    integer, dimension (natyp), intent (in) :: imt !! R point at MT radius
    integer, dimension (natyp), intent (in) :: irc !! R point for potential cutting
    integer, dimension (natyp), intent (in) :: nfu
    integer, dimension (nembd2), intent (in) :: cls !! Cluster around atomic sites
    integer, dimension (natyp), intent (in) :: irws !! R point at WS radius
    integer, dimension (natyp), intent (in) :: irns !! Position of atoms in the unit cell in units of bravais vectors
    integer, dimension (natyp), intent (in) :: zrel !! atomic number (cast integer)
    integer, dimension (natyp), intent (in) :: iqat !! The site on which an atom is located on a given site
    integer, dimension (naez), intent (in) :: icpa !! ICPA = 0/1 site-dependent CPA flag
    integer, dimension (natyp), intent (in) :: ipan !! Number of panels in non-MT-region
    integer, dimension (natyp), intent (in) :: lopt !! angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
    integer, dimension (nsheld), intent (in) :: nsh1 !! Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
    integer, dimension (nsheld), intent (in) :: nsh2 !! Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
    integer, dimension (nclsd), intent (in) :: nacls !! Number of atoms in cluster
    integer, dimension (iemxd), intent (in) :: kmesh
    integer, dimension (natyp), intent (in) :: irmin !! Max R for spherical treatment
    integer, dimension (lm2d), intent (in) :: loflm !! l of lm=(l,m) (GAUNT)
    integer, dimension (npotd), intent (in) :: ncore !! Number of core states
    integer, dimension (naez), intent (in) :: iqcalc
    integer, dimension (natyp), intent (in) :: itldau !! integer pointer connecting the NTLDAU atoms to heir corresponding index in the unit cell
    integer, dimension (natyp), intent (in) :: jwsrel !! index of the WS radius
    integer, dimension (natyp), intent (in) :: ntcell !! Index for WS cell
    integer, dimension (natyp), intent (in) :: ixipol !! Constraint of spin pol.
    integer, dimension (nembd2), intent (in) :: refpot !! Ref. pot. card  at position
    integer, dimension (0:nsheld), intent (in) :: nshell !! Index of atoms/pairs per shell (ij-pairs); nshell(0) = number of shells
    integer, dimension (0:lmpot), intent (in) :: imaxsh
    integer, dimension (natyp), intent (in) :: irshift
    integer, dimension (natomimpd), intent (in) :: atomimp
    integer, dimension (nofgij), intent (in) :: ijtabsh !! Linear pointer, assigns pair (i,j) to a shell in the array GS(*,*,*,NSHELD)
    integer, dimension (0:natyp), intent (in) :: hostimp
    integer, dimension (nofgij), intent (in) :: ijtabsym !! Linear pointer, assigns pair (i,j) to the rotation bringing GS into Gij
    integer, dimension (natyp), intent (in) :: npan_tot
    integer, dimension (nofgij), intent (in) :: ijtabcalc !! Linear pointer, specifying whether the block (i,j) has to be calculated needs set up for ICC=-1, not used for ICC=1
    integer, dimension (natyp), intent (in) :: npan_eq_at
    integer, dimension (natyp), intent (in) :: npan_log_at
    integer, dimension (nofgij), intent (in) :: ijtabcalc_i
    integer, dimension (ngshd, 3), intent (in) :: ilm_map
    integer, dimension (nsheld, 2*nsymaxd), intent (in) :: ish
    integer, dimension (nsheld, 2*nsymaxd), intent (in) :: jsh
    integer, dimension (natyp, lmxspd), intent (in) :: lmsp !! 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
    integer, dimension (naclsd, nembd2), intent (in) :: atom !! Atom at site in cluster
    integer, dimension (naclsd, nembd2), intent (in) :: ezoa !! EZ of atom at site in cluster
    integer, dimension (ncleb, 4), intent (in) :: icleb !! Pointer array
    integer, dimension (20, npotd), intent (in) :: lcore !! Angular momentum of core states
    integer, dimension (2, lmmaxd), intent (in) :: nrrel
    integer, dimension (natyp, nfund), intent (in) :: llmsp !! lm=(l,m) of 'nfund'th nonvanishing component of non-spherical pot.
    integer, dimension (lmxspd, natyp), intent (in) :: lmsp1
    integer, dimension (natyp, lmxspd), intent (in) :: ifunm
    integer, dimension (0:ipand, natyp), intent (in) :: ircut !! R points of panel borders
    integer, dimension (natyp, nembd2), intent (in) :: kaoez !! Kind of atom at site in elem. cell
    integer, dimension (20, npotd), intent (in) :: ititle
    integer, dimension (lmxspd, natyp), intent (in) :: ifunm1
    integer, dimension (naezdpd, naezdpd), intent (in) :: icheck
    integer, dimension (0:ntotd, natyp), intent (in) :: ipan_intervall
    integer, dimension (lmpot, 0:lmax, 0:lmax), intent (in) :: jend !! Pointer array for icleb()
    integer, dimension (2, 2, lmmaxd), intent (in) :: irrel

    logical, dimension (2), intent (in) :: vacflag
    logical, dimension (nsymaxd), intent (in) :: symunitary !! unitary/antiunitary symmetry flag
    character (len=124), dimension (6), intent (in) :: txc
    ! ..
    ! fill arrays:

    t_params%ez = ez
    t_params%wez = wez
    t_params%drotq = drotq
    t_params%dsymll = dsymll
    t_params%lefttinvll = lefttinvll
    t_params%righttinvll = righttinvll
    t_params%crel = crel
    t_params%rc = rc
    t_params%rrel = rrel
    t_params%srrel = srrel
    t_params%phildau = phildau
    t_params%vins = vins
    t_params%visp = visp
    t_params%vbc = vbc
    if (t_params%krel>0) t_params%vtrel = vtrel
    if (t_params%krel>0) t_params%btrel = btrel
    t_params%socscale = socscale
    if (t_params%krel>0) t_params%drdirel = drdirel
    if (t_params%krel>0) t_params%r2drdirel = r2drdirel
    if (t_params%krel>0) t_params%rmrel = rmrel
    t_params%cmomhost = cmomhost
    t_params%ecore = ecore
    t_params%qmtet = qmtet
    t_params%qmphi = qmphi
    t_params%qmphitab = qmphitab
    t_params%qmtettab = qmtettab
    t_params%qmgamtab = qmgamtab
    t_params%zat = zat
    t_params%rmesh = r
    t_params%drdi = drdi
    t_params%rmtref = rmtref
    t_params%vref = vref
    t_params%cleb = cleb
    t_params%rcls = rcls
    t_params%socscl = socscl
    t_params%cscl = cscl
    t_params%rbasis = rbasis
    t_params%rr = rr
    t_params%conc = conc
    t_params%rrot = rrot
    t_params%ratom = ratom
    t_params%a = a
    t_params%b = b
    t_params%thetas = thetas
    t_params%rmt = rmt
    t_params%rmtnew = rmtnew
    t_params%rws = rws
    t_params%gsh = gsh
    t_params%erefldau = erefldau
    t_params%ueff = ueff
    t_params%jeff = jeff
    if (t_params%idoldau==1) t_params%uldau = uldau
    t_params%wldau = wldau
    t_params%rclsimp = rclsimp
    t_params%rpan_intervall = rpan_intervall
    t_params%rnew = rnew
    t_params%thetasnew = thetasnew
    t_params%lopt = lopt
    t_params%itldau = itldau
    t_params%irshift = irshift
    t_params%jwsrel = jwsrel
    t_params%zrel = zrel
    t_params%lcore = lcore
    t_params%ncore = ncore
    t_params%ipan = ipan
    t_params%ircut = ircut
    t_params%jend = jend
    t_params%icleb = icleb
    t_params%atom = atom
    t_params%cls = cls
    t_params%nacls = nacls
    t_params%loflm = loflm
    t_params%ezoa = ezoa
    t_params%kaoez = kaoez
    t_params%iqat = iqat
    t_params%icpa = icpa
    t_params%noq = noq
    t_params%kmesh = kmesh
    t_params%nshell = nshell
    t_params%nsh1 = nsh1
    t_params%nsh2 = nsh2
    t_params%ijtabcalc = ijtabcalc
    t_params%ijtabcalc_i = ijtabcalc_i
    t_params%ijtabsym = ijtabsym
    t_params%ijtabsh = ijtabsh
    t_params%ish = ish
    t_params%jsh = jsh
    t_params%iqcalc = iqcalc
    t_params%icheck = icheck
    t_params%atomimp = atomimp
    t_params%refpot = refpot
    t_params%irrel = irrel
    t_params%nrrel = nrrel
    t_params%ifunm1 = ifunm1
    t_params%ititle = ititle
    t_params%lmsp1 = lmsp1
    t_params%ntcell = ntcell
    t_params%ixipol = ixipol
    t_params%irns = irns
    t_params%ifunm = ifunm
    t_params%llmsp = llmsp
    t_params%lmsp = lmsp
    t_params%imt = imt
    t_params%irc = irc
    t_params%irmin = irmin
    t_params%irws = irws
    t_params%nfu = nfu
    t_params%hostimp = hostimp
    t_params%ilm_map = ilm_map
    t_params%imaxsh = imaxsh
    t_params%npan_log = npan_log
    t_params%npan_eq = npan_eq
    t_params%npan_log_at = npan_log_at
    t_params%npan_eq_at = npan_eq_at
    t_params%npan_tot = npan_tot
    t_params%ipan_intervall = ipan_intervall
    t_params%symunitary = symunitary
    t_params%vacflag = vacflag
    t_params%txc = txc
    t_params%ntotd = ntotd

  end subroutine fill_t_params_arrays

  !-------------------------------------------------------------------------------
  !> Summary: Set the values of the local variables according to the stored `t_params` so that they can be passed between different control modules, specifically for `main1a`
  !> Author: Philipp Ruessmann
  !> Category: communication, KKRhost
  !> Deprecated: False
  !> Set the values of the local variables according to the stored `t_params`
  !> so that they can be passed between different control modules, specifically for `main1a`
  !-------------------------------------------------------------------------------
  subroutine get_params_1a(t_params,ipand,natypd,irmd,naclsd,ielast,nclsd,nrefd,    &
    ncleb,nembd,naezd,lm2d,nsra,ins,nspin,icst,ipan,ircut,lmax,ncls,nineq,idoldau,  &
    lly,krel,atom,cls,icleb,loflm,nacls,refpot,irws,iend,ez,vins,irmin,itmpdir,     &
    iltmp,alat,drdi,rmesh,zat,rcls,iemxd,visp,rmtref,vref,cleb,cscl,socscale,       &
    socscl,erefldau,ueff,jeff,solver,tmpdir,deltae,tolrdif,npan_log,npan_eq,ncheb,  &
    npan_tot,ipan_intervall,rpan_intervall,rnew,ntotd,nrmaxd,r_log,ntldau,itldau,   &
    lopt,vtrel,btrel,drdirel,r2drdirel,rmrel,irmind,lmpot,nspotd,npotd,jwsrel,zrel, &
    itscf,natomimpd,natomimp,atomimp,iqat,naez,natyp,nref)
    ! get relevant parameters from t_params
    ! ..
#ifdef CPP_MPI
    use :: mpi
#endif
    implicit none

    type (type_params), intent (in) :: t_params
    integer, intent (in) :: irmd   !! Maximum number of radial points
    integer, intent (in) :: krel   !! Switch for non-relativistic/relativistic (0/1) program. Attention: several other parameters depend explicitly on KREL, they are set automatically Used for Dirac solver in ASA
    integer, intent (in) :: nembd  !! Number of 'embedding' positions
    integer, intent (in) :: lm2d   !! (2*LMAX+1)**2
    integer, intent (in) :: nclsd  !! Maximum number of different TB-clusters
    integer, intent (in) :: ncleb  !! Number of Clebsch-Gordon coefficients
    integer, intent (in) :: ntotd
    integer, intent (in) :: ipand  !! Number of panels in non-spherical part
    integer, intent (in) :: iemxd  !! Dimension for energy-dependent arrays
    integer, intent (in) :: lmpot  !! (LPOT+1)**2
    integer, intent (in) :: npotd  !! (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP)
    integer, intent (in) :: naclsd !! Maximum number of atoms in a TB-cluster
    integer, intent (in) :: nrmaxd
    integer, intent (in) :: irmind !! IRM-IRNSD
    integer, intent (in) :: nspotd !! Number of potentials for storing non-sph. potentials
    integer, intent (in) :: natomimpd !! Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
    integer, intent (inout) :: lly !! LLY!> 0 : apply Lloyds formula
    integer, intent (inout) :: ins !! 0 (MT), 1(ASA), 2(Full Potential)
    integer, intent (inout) :: nsra
    integer, intent (inout) :: lmax !! Maximum l component in wave function expansion
    integer, intent (in) :: nrefd  !! Number of diff. ref. potentials
    integer, intent (out) :: nref  !! Number of diff. ref. potentials
    integer, intent (inout) :: icst !! Number of Born approximation
    integer, intent (inout) :: ncls !! Number of reference clusters
    integer, intent (inout) :: iend !! Number of nonzero gaunt coefficients
    integer, intent (in) :: naezd  !! Number of atoms in unit cell
    integer, intent (out) :: naez  !! Number of atoms in unit cell
    integer, intent (in) :: natypd !! Number of kinds of atoms in unit cell
    integer, intent (out) :: natyp !! Number of kinds of atoms in unit cell
    integer, intent (inout) :: nspin !! Counter for spin directions
    integer, intent (inout) :: nineq !! Number of ineq. positions in unit cell
    integer, intent (inout) :: iltmp
    integer, intent (inout) :: itscf
    integer, intent (inout) :: ncheb !! Number of Chebychev pannels for the new solver
    integer, intent (inout) :: ntldau !! number of atoms on which LDA+U is applied
    integer, intent (inout) :: ielast
    integer, intent (inout) :: itmpdir
    integer, intent (inout) :: idoldau !! flag to perform LDA+U
    integer, intent (inout) :: natomimp !! Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
    integer, dimension (naezd+nembd), intent (inout) :: cls !! Cluster around atomic sites
    integer, dimension (natypd), intent (inout) :: lopt !! angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
    integer, dimension (natypd), intent (inout) :: irws !! R point at WS radius
    integer, dimension (natypd), intent (inout) :: ipan !! Number of panels in non-MT-region
    integer, dimension (natypd), intent (inout) :: zrel !! atomic number (cast integer)
    integer, dimension (natypd), intent (inout) :: iqat !! The site on which an atom is located on a given site
    integer, dimension (lm2d), intent (inout) :: loflm !! l of lm=(l,m) (GAUNT)
    integer, dimension (natypd), intent (inout) :: irmin !! Max R for spherical treatment
    integer, dimension (nclsd), intent (inout) :: nacls !! Number of atoms in cluster
    integer, dimension (naezd+nembd), intent (inout) :: refpot !! Ref. pot. card  at position
    integer, dimension (natypd), intent (inout) :: itldau !! integer pointer connecting the NTLDAU atoms to heir corresponding index in the unit cell
    integer, dimension (natypd), intent (inout) :: jwsrel !! index of the WS radius
    integer, dimension (natypd), intent (inout) :: npan_eq
    integer, dimension (natypd), intent (inout) :: npan_log
    integer, dimension (natypd), intent (inout) :: npan_tot
    integer, dimension (natomimpd), intent (inout) :: atomimp
    integer, dimension (naclsd, naezd+nembd), intent (inout) :: atom !! Atom at site in cluster
    integer, dimension (ncleb, 4), intent (inout) :: icleb !! Pointer array
    integer, dimension (0:ipand, natypd), intent (inout) :: ircut !! R points of panel borders
    integer, dimension (0:ntotd, natypd), intent (inout) :: ipan_intervall
    real (kind=dp), intent (inout) :: alat !! Lattice constant in a.u.
    real (kind=dp), intent (inout) :: r_log
    real (kind=dp), intent (inout) :: tolrdif !! Tolerance for r<tolrdif (a.u.) to handle vir. atoms
    real (kind=dp), dimension (natypd), intent (inout) :: zat !! Nuclear charge
    real (kind=dp), dimension (nrefd), intent (inout) :: vref
    real (kind=dp), dimension (natypd), intent (inout) :: ueff !! input U parameter for each atom
    real (kind=dp), dimension (natypd), intent (inout) :: jeff !! input J parameter for each atom
    real (kind=dp), dimension (nrefd), intent (inout) :: rmtref !! Muffin-tin radius of reference system
    real (kind=dp), dimension (natypd), intent (inout) :: erefldau !! the energies of the projector's wave functions (REAL)
    real (kind=dp), dimension (natypd), intent (inout) :: socscale !! Spin-orbit scaling
    real (kind=dp), dimension (ncleb, 2), intent (inout) :: cleb !! GAUNT coefficients (GAUNT)
    real (kind=dp), dimension (irmd, npotd), intent (inout) :: visp !! Spherical part of the potential
    real (kind=dp), dimension (irmd, natypd), intent (inout) :: drdi !! Derivative dr/di
    real (kind=dp), dimension (nrmaxd, natypd), intent (inout) :: rnew
    real (kind=dp), dimension (krel*lmax+1, krel*natypd+(1-krel)), intent (inout) :: cscl !! Speed of light scaling
    real (kind=dp), dimension (irmd, natypd), intent (inout) :: rmesh
    real (kind=dp), dimension (irmd*krel+(1-krel), natypd), intent (inout) :: vtrel !! potential (spherical part)
    real (kind=dp), dimension (irmd*krel+(1-krel), natypd), intent (inout) :: btrel !! magnetic field
    real (kind=dp), dimension (irmd*krel+(1-krel), natypd), intent (inout) :: rmrel !! radial mesh
    real (kind=dp), dimension (krel*lmax+1, krel*natypd+(1-krel)), intent (inout) :: socscl
    real (kind=dp), dimension (irmd*krel+(1-krel), natypd), intent (inout) :: drdirel !! derivative of radial mesh
    real (kind=dp), dimension (irmd*krel+(1-krel), natypd), intent (inout) :: r2drdirel !! \( r^2 \frac{\partial}{\partial \mathbf{r}}\frac{\partial}{\partial i}\) (r**2 * drdi)
    real (kind=dp), dimension (0:ntotd, natypd), intent (inout) :: rpan_intervall
    real (kind=dp), dimension (3, naclsd, nclsd), intent (inout) :: rcls !! Real space position of atom in cluster
    real (kind=dp), dimension (irmind:irmd, lmpot, nspotd), intent (inout) :: vins !! Non-spherical part of the potential

    complex (kind=dp), intent (inout) :: deltae !! Energy difference for numerical derivative
    complex (kind=dp), dimension (iemxd), intent (inout) :: ez
    character (len=10), intent (inout) :: solver                           !! Type of solver

    character (len=80), intent (inout) :: tmpdir

    nsra = t_params%nsra
    ins = t_params%ins
    naez = t_params%naez
    natyp = t_params%natyp
    nspin = t_params%nspin
    icst = t_params%icst
    ipan = t_params%ipan
    ircut = t_params%ircut
    lmax = t_params%lmax
    ncls = t_params%ncls
    nineq = t_params%nineq
    nref = t_params%nref
    idoldau = t_params%idoldau
    lly = t_params%lly

    ! -------------------------------------------------------------------------
    ! Consistency check
    ! -------------------------------------------------------------------------
    if ((krel==1) .and. (ins/=0)) then
      write (6, *) ' FULL-POTENTIAL RELATIVISTIC mode not implemented '
      stop ' set INS = 0 in the input'
    end if

    if (nsra<=2) then
      if (krel==1) stop ' KVREL!= 1 in input, but relativistic program used'
    else
      if (krel==0) stop ' KVREL > 1 in input, but non-relativistic program used'
    end if
    ! -------------------------------------------------------------------------
    ! End of consistency check
    ! -------------------------------------------------------------------------

    alat = t_params%alat
    zat = t_params%zat
    drdi = t_params%drdi
    rmesh = t_params%rmesh
    rmtref = t_params%rmtref
    vref = t_params%vref
    iend = t_params%iend
    cleb = t_params%cleb
    rcls = t_params%rcls
    atom = t_params%atom
    cls = t_params%cls
    icleb = t_params%icleb
    loflm = t_params%loflm
    nacls = t_params%nacls
    refpot = t_params%refpot
    irws = t_params%irws
    irmin = t_params%irmin
    tolrdif = t_params%tolrdif
    deltae = t_params%deltae
    socscale = t_params%socscale
    tmpdir = t_params%tmpdir
    itmpdir = t_params%itmpdir
    iltmp = t_params%iltmp
    npan_log = t_params%npan_log_at
    npan_eq = t_params%npan_eq_at
    ncheb = t_params%ncheb
    r_log = t_params%r_log
    npan_tot = t_params%npan_tot
    rnew = t_params%rnew
    rpan_intervall = t_params%rpan_intervall
    ipan_intervall = t_params%ipan_intervall

    if (krel==1) then
      solver = t_params%solver
      socscl = t_params%socscl
      cscl = t_params%cscl
    end if

    if (idoldau==1) then
      ntldau = t_params%ntldau
      itldau = t_params%itldau
      lopt = t_params%lopt
      ueff = t_params%ueff
      jeff = t_params%jeff
      erefldau = t_params%erefldau
    end if
    ! -------------------------------------------------------------------------
    ! Energy_mesh
    ! -------------------------------------------------------------------------
    ielast = t_params%ielast
    ez = t_params%ez
    ! -------------------------------------------------------------------------
    ! Input_potential
    ! -------------------------------------------------------------------------
    vins = t_params%vins
    visp = t_params%visp
    if (krel==1) then
      rmrel = t_params%rmrel
      drdirel = t_params%drdirel
      r2drdirel = t_params%r2drdirel
      zrel = t_params%zrel
      jwsrel = t_params%jwsrel
      vtrel = t_params%vtrel
      btrel = t_params%btrel
    end if
    itscf = t_params%itscf
    ! -------------------------------------------------------------------------
    ! Cluster atoms
    ! -------------------------------------------------------------------------
    natomimp = t_params%natomimp
    atomimp = t_params%atomimp
    iqat = t_params%iqat

  end subroutine get_params_1a

  !-------------------------------------------------------------------------------
  !> Summary: Set the values of the local variables according to the stored `t_params` so that they can be passed between different control modules, specifically for `main1b`
  !> Author: Philipp Ruessmann
  !> Category: communication, KKRhost
  !> Deprecated: False
  !> Set the values of the local variables according to the stored `t_params`
  !> so that they can be passed between different control modules, specifically for `main1b`
  !-------------------------------------------------------------------------------
  subroutine get_params_1b(t_params,natypd,naezd,natyp,naclsd,ielast,npol,nclsd,    &
    nrefd,nref,nembd,naez,nsra,ins,nspin,lmax,ncls,lly,krel,atom,cls,nacls,refpot,  &
    ez,itmpdir,iltmp,alat,rcls,iemxd,rmtref,vref,tmpdir,nsheld,nprincd,kpoibz,      &
    atomimp,natomimpd,icc,igf,nlbasis,nrbasis,ncpa,icpa,itcpamax,cpatol,nrd,ideci,  &
    rbasis,rr,ezoa,nshell,kmrot,kaoez,ish,jsh,nsh1,nsh2,noq,iqat,nofgij,natomimp,   &
    conc,kmesh,maxmesh,nsymat,nqcalc,ratom,rrot,drotq,ijtabcalc,ijtabcalc_i,        &
    ijtabsym,ijtabsh,iqcalc,dsymll,invmod,icheck,symunitary,rc,crel,rrel,srrel,     &
    nrrel,irrel,lefttinvll,righttinvll,vacflag,nofks,volbz,bzkp,volcub,wez,nembd1,  &
    lmmaxd,nspindd,maxmshd,rclsimp)
    ! get relevant parameters from t_params
    ! ..

    use :: mod_runoptions, only: relax_SpinAngle_Dirac, use_decimation, write_kpts_file

    implicit none

    type (type_params), intent (in) :: t_params

    integer, intent (in) :: natypd
    integer, intent (in) :: naezd
    integer, intent (in) :: nembd
    integer, intent (in) :: nrd
    integer, intent (in) :: nrefd
    integer, intent (in) :: krel
    integer, intent (in) :: iemxd
    integer, intent (in) :: nclsd
    integer, intent (in) :: naclsd
    integer, intent (in) :: nembd1
    integer, intent (in) :: lmmaxd
    integer, intent (in) :: nsheld
    integer, intent (in) :: kpoibz
    integer, intent (in) :: nspindd
    integer, intent (in) :: maxmshd
    integer, intent (in) :: nprincd
    integer, intent (in) :: natomimpd
    integer, intent (inout) :: lly
    integer, intent (inout) :: icc
    integer, intent (inout) :: igf
    integer, intent (inout) :: ins
    integer, intent (inout) :: npol
    integer, intent (inout) :: nsra
    integer, intent (inout) :: lmax
    integer, intent (inout) :: ncls
    integer, intent (out) :: nref
    integer, intent (inout) :: ncpa
    integer, intent (out) :: naez
    integer, intent (out) :: natyp
    integer, intent (inout) :: nspin
    integer, intent (inout) :: iltmp
    integer, intent (inout) :: ideci
    integer, intent (inout) :: kmrot
    integer, intent (inout) :: invmod
    integer, intent (inout) :: nsymat
    integer, intent (inout) :: ielast
    integer, intent (inout) :: nqcalc
    integer, intent (inout) :: nofgij
    integer, intent (inout) :: maxmesh
    integer, intent (inout) :: itmpdir
    integer, intent (inout) :: nlbasis
    integer, intent (inout) :: nrbasis
    integer, intent (inout) :: natomimp
    integer, intent (inout) :: itcpamax
    integer, dimension (naezd+nembd), intent (inout) :: cls
    integer, dimension (naezd), intent (inout) :: noq
    integer, dimension (naezd), intent (inout) :: icpa
    integer, dimension (natypd), intent (inout) :: iqat
    integer, dimension (nsheld), intent (inout) :: nsh1
    integer, dimension (nsheld), intent (inout) :: nsh2
    integer, dimension (iemxd), intent (inout) :: kmesh
    integer, dimension (nclsd), intent (inout) :: nacls
    integer, dimension (maxmshd), intent (inout) :: nofks
    integer, dimension (0:nsheld), intent (inout) :: nshell
    integer, dimension (naezd), intent (inout) :: iqcalc
    integer, dimension (naezd+nembd), intent (inout) :: refpot
    integer, dimension (natomimpd), intent (inout) :: atomimp
    integer, dimension (nofgij), intent (inout) :: ijtabsh
    integer, dimension (nofgij), intent (inout) :: ijtabsym
    integer, dimension (nofgij), intent (inout) :: ijtabcalc
    integer, dimension (nofgij), intent (inout) :: ijtabcalc_i
    integer, dimension (nsheld, 2*nsymaxd), intent (inout) :: ish
    integer, dimension (nsheld, 2*nsymaxd), intent (inout) :: jsh
    integer, dimension (naclsd, naezd+nembd), intent (inout) :: ezoa
    integer, dimension (naclsd, naezd+nembd), intent (inout) :: atom
    integer, dimension (2, lmmaxd), intent (inout) :: nrrel
    integer, dimension (natypd, naezd+nembd), intent (inout) :: kaoez
    integer, dimension (naezd/nprincd, naezd/nprincd), intent (inout) :: icheck
    integer, dimension (2, 2, lmmaxd), intent (inout) :: irrel
    real (kind=dp), intent (inout) :: alat
    real (kind=dp), intent (inout) :: cpatol
    real (kind=dp), dimension (nrefd), intent (inout) :: vref
    real (kind=dp), dimension (natypd), intent (inout) :: conc
    real (kind=dp), dimension (nrefd), intent (inout) :: rmtref
    real (kind=dp), dimension (maxmshd), intent (inout) :: volbz
    real (kind=dp), dimension (3, 0:nrd), intent (inout) :: rr
    real (kind=dp), dimension (3, nsheld), intent (inout) :: ratom
    real (kind=dp), dimension (3, naezd+nembd), intent (inout) :: rbasis
    real (kind=dp), dimension (kpoibz, maxmshd), intent (inout) :: volcub
    real (kind=dp), dimension (3, natomimpd), intent (inout) :: rclsimp
    real (kind=dp), dimension (48, 3, nsheld), intent (inout) :: rrot
    real (kind=dp), dimension (3, naclsd, nclsd), intent (inout) :: rcls
    real (kind=dp), dimension (3, kpoibz, maxmshd), intent (inout) :: bzkp
    complex (kind=dp), dimension (iemxd), intent (inout) :: ez
    complex (kind=dp), dimension (iemxd), intent (inout) :: wez
    complex (kind=dp), dimension (lmmaxd, lmmaxd), intent (inout) :: rc
    complex (kind=dp), dimension (lmmaxd, lmmaxd), intent (inout) :: crel
    complex (kind=dp), dimension (lmmaxd, lmmaxd), intent (inout) :: rrel
    complex (kind=dp), dimension (lmmaxd, lmmaxd, naezd), intent (inout) :: drotq
    complex (kind=dp), dimension (lmmaxd, lmmaxd, nsymaxd), intent (inout) :: dsymll
    complex (kind=dp), dimension (2, 2, lmmaxd), intent (inout) :: srrel
    complex (kind=dp), dimension (lmmaxd, lmmaxd, nembd1, nspindd, iemxd), intent (inout) :: lefttinvll
    complex (kind=dp), dimension (lmmaxd, lmmaxd, nembd1, nspindd, iemxd), intent (inout) :: righttinvll
    character (len=80), intent (inout) :: tmpdir
    logical, dimension (2), intent (inout) :: vacflag
    logical, dimension (nsymaxd), intent (inout) :: symunitary

    integer :: i, l, id

    ! .. External Functions ..

    nsra = t_params%nsra
    ins = t_params%ins
    natyp = t_params%natyp
    naez = t_params%naez
    nspin = t_params%nspin
    lmax = t_params%lmax
    nref = t_params%nref
    icc = t_params%icc
    igf = t_params%igf
    nlbasis = t_params%nlbasis
    nrbasis = t_params%nrbasis
    ncpa = t_params%ncpa
    icpa = t_params%icpa
    itcpamax = t_params%itcpamax
    cpatol = t_params%cpatol
    ! -------------------------------------------------------------------------
    ! Consistency check
    ! -------------------------------------------------------------------------

    if ((krel==1) .and. (ins/=0)) then
      write (6, *) ' FULL-POTENTIAL RELATIVISTIC mode not implemented '
      stop ' set INS = 0 in the input'
    end if

    if (nsra<=2) then
      if (krel==1) stop ' KVREL!= 1 in input, but relativistic program used'
    else
      if (krel==0) stop ' KVREL > 1 in input, but non-relativistic program used'
    end if
    ideci = 0
    ! -------------------------------------------------------------------------
    alat = t_params%alat
    rbasis = t_params%rbasis
    refpot = t_params%refpot
    rmtref = t_params%rmtref
    vref = t_params%vref
    rcls = t_params%rcls
    rr = t_params%rr
    atom = t_params%atom
    cls = t_params%cls
    ncls = t_params%ncls
    ezoa = t_params%ezoa
    nacls = t_params%nacls
    nshell = t_params%nshell
    kmrot = t_params%kmrot
    kaoez = t_params%kaoez
    iqat = t_params%iqat
    noq = t_params%noq
    conc = t_params%conc
    kmesh = t_params%kmesh
    maxmesh = t_params%maxmesh
    lly = t_params%lly
    nsymat = t_params%nsymat
    natomimp = t_params%natomimp
    nqcalc = t_params%nqcalc
    ratom = t_params%ratom
    rrot = t_params%rrot
    nsh1 = t_params%nsh1
    nsh2 = t_params%nsh2
    drotq = t_params%drotq
    ijtabcalc = t_params%ijtabcalc
    ijtabcalc_i = t_params%ijtabcalc_i
    do i = 1, nshell(0)
      do l = 1, 2*nsymaxd
        ish(i, l) = t_params%ish(i, l)
        jsh(i, l) = t_params%jsh(i, l)
      end do
    end do
    ijtabsym = t_params%ijtabsym
    ijtabsh = t_params%ijtabsh
    iqcalc = t_params%iqcalc
    dsymll = t_params%dsymll
    invmod = t_params%invmod
    icheck = t_params%icheck
    atomimp = t_params%atomimp
    symunitary = t_params%symunitary
    tmpdir = t_params%tmpdir
    itmpdir = t_params%itmpdir
    iltmp = t_params%iltmp
    rclsimp = t_params%rclsimp
    if (krel==1) then
      rc = t_params%rc
      crel = t_params%crel
      rrel = t_params%rrel
      srrel = t_params%srrel
      nrrel = t_params%nrrel
      irrel = t_params%irrel
    end if
    if (use_decimation) then
      lefttinvll = t_params%lefttinvll
      righttinvll = t_params%righttinvll
      vacflag = t_params%vacflag
      ideci = 1
    end if
    ! -------------------------------------------------------------------------
    ! K-points
    ! -------------------------------------------------------------------------
    if (write_kpts_file) then
      open (52, file='kpoints', form='formatted')
      rewind (52)
      write (1337, *) 'kpoints read from kpoints file due to test option "kptsfile"'
    end if
    do l = 1, maxmesh
      if (write_kpts_file) then
        read (52, fmt='(I8,f15.10)') nofks(l), volbz(l)
        write (1337, *) 'kpts:', nofks(l), volbz(l), t_params%nofks(l), t_params%volbz(l)
      else
        nofks(l) = t_params%nofks(l)
        volbz(l) = t_params%volbz(l)
      end if
      do i = 1, nofks(l)
        if (write_kpts_file) then
          read (52, fmt=*)(bzkp(id,i,l), id=1, 3), volcub(i, l)
        else
          do id = 1, 3
            ! write(*,*)'bzkp', BZKP(ID,I,L), t_params%BZKP(ID,I,L)
            bzkp(id, i, l) = t_params%bzkp(id, i, l)
          end do
          ! write(*,*)'volcub', VOLCUB(I,L), t_params%VOLCUB(I,L)
          volcub(i, l) = t_params%volcub(i, l)
        end if
        ! write(*,'(A,4F21.17,I5,F21.17)') 'bzkmesh input:',
        ! &             (BZKP(ID,I,L),ID=1,3),VOLCUB(I,L),NOFKS(L),VOLBZ(L)
      end do
    end do
    if (write_kpts_file) close (52)
    ! -------------------------------------------------------------------------
    ! Energy_mesh
    ! -------------------------------------------------------------------------
    ielast = t_params%ielast
    ez = t_params%ez
    wez = t_params%wez
    npol = t_params%npol
    ! -------------------------------------------------------------------------
    ! Itermdir
    ! -------------------------------------------------------------------------
    if (relax_SpinAngle_Dirac) then
      drotq = t_params%drotq
      if (kmrot==0) kmrot = 1
    end if

  end subroutine get_params_1b

  !-------------------------------------------------------------------------------
  !> Summary: Set the values of the local variables according to the stored `t_params` so that they can be passed between different control modules, specifically for `main1c`
  !> Author: Philipp Ruessmann
  !> Category: communication, KKRhost
  !> Deprecated: False
  !> Set the values of the local variables according to the stored `t_params`
  !> so that they can be passed between different control modules, specifically for `main1c`
  !-------------------------------------------------------------------------------
  subroutine get_params_1c(t_params,krel,naezd,natypd,ncleb,lm2d,ncheb,ipand,lmpotd,&
    lmaxd,lmxspd,nfund,npotd,ntotd,mmaxd,iemxd,irmd,nsra,ins,nspin,nacls1,icst,     &
    kmrot,iqat,idoldau,irws,ipan,ircut,iend,icleb,loflm,jend,ifunm1,lmsp1,nfu,llmsp,&
    lcore,ncore,ntcell,irmin,ititle,intervx,intervy,intervz,lly,itmpdir,iltmp,      &
    npan_eq,ipan_intervall,npan_log,npan_tot,ntldau,lopt,itldau,ielast,iesemicore,  &
    npol,irshift,jwsrel,zrel,itrunldau,qmtet,qmphi,conc,alat,zat,drdi,rmesh,a,b,    &
    cleb,thetas,socscale,rpan_intervall,cscl,rnew,socscl,thetasnew,efermi,erefldau, &
    ueff,jeff,emin,emax,tk,vins,visp,ecore,drdirel,r2drdirel,rmrel,vtrel,btrel,     &
    wldau,uldau,ez,wez,phildau,tmpdir,solver,nspind,nspotd,irmind,lmaxd1,ncelld,    &
    irid,r_log,naez,natyp,lmax)
    ! get relevant parameters from t_params
    ! ..
    use :: mod_runoptions, only: relax_SpinAngle_Dirac, use_spherical_potential_only
    implicit none

    type (type_params), intent (in) :: t_params

    integer, intent (in) :: irmd
    integer, intent (in) :: krel
    integer, intent (in) :: irid
    integer, intent (in) :: lm2d
    integer, intent (in) :: ncleb
    integer, intent (in) :: ntotd
    integer, intent (in) :: ipand
    integer, intent (in) :: lmpotd
    integer, intent (in) :: nfund
    integer, intent (in) :: npotd
    integer, intent (in) :: mmaxd
    integer, intent (in) :: iemxd
    integer, intent (in) :: lmxspd
    integer, intent (in) :: nspind
    integer, intent (in) :: nspotd
    integer, intent (in) :: irmind
    integer, intent (in) :: lmaxd1
    integer, intent (in) :: ncelld
    integer, intent (inout) :: ins
    integer, intent (inout) :: lly
    integer, intent (inout) :: nsra
    integer, intent (inout) :: icst
    integer, intent (in) :: naezd
    integer, intent (out) :: naez
    integer, intent (in) :: lmaxd
    integer, intent (out) :: lmax
    integer, intent (inout) :: npol
    integer, intent (inout) :: iend
    integer, intent (inout) :: ncheb
    integer, intent (in) :: natypd
    integer, intent (out) :: natyp
    integer, intent (inout) :: nspin
    integer, intent (inout) :: kmrot
    integer, intent (inout) :: iltmp
    integer, intent (inout) :: ntldau
    integer, intent (inout) :: ielast
    integer, intent (inout) :: nacls1
    integer, intent (inout) :: idoldau
    integer, intent (inout) :: intervx
    integer, intent (inout) :: intervy
    integer, intent (inout) :: intervz
    integer, intent (inout) :: itmpdir
    integer, intent (inout) :: itrunldau
    integer, intent (inout) :: iesemicore
    integer, dimension (natypd), intent (inout) :: nfu
    integer, dimension (natypd), intent (inout) :: lopt
    integer, dimension (natypd), intent (inout) :: iqat
    integer, dimension (natypd), intent (inout) :: irws
    integer, dimension (natypd), intent (inout) :: ipan
    integer, dimension (natypd), intent (inout) :: zrel
    integer, dimension (lm2d), intent (inout) :: loflm
    integer, dimension (npotd), intent (inout) :: ncore
    integer, dimension (natypd), intent (inout) :: irmin
    integer, dimension (natypd), intent (inout) :: ntcell
    integer, dimension (natypd), intent (inout) :: itldau
    integer, dimension (natypd), intent (inout) :: jwsrel
    integer, dimension (natypd), intent (inout) :: npan_eq
    integer, dimension (natypd), intent (inout) :: irshift
    integer, dimension (natypd), intent (inout) :: npan_log
    integer, dimension (natypd), intent (inout) :: npan_tot
    integer, dimension (lmpotd, 0:lmaxd, 0:lmaxd), intent (inout) :: jend
    integer, dimension (0:ipand, natypd), intent (inout) :: ircut
    integer, dimension (ncleb, 4), intent (inout) :: icleb
    integer, dimension (lmxspd, natypd), intent (inout) :: lmsp1
    integer, dimension (natypd, nfund), intent (inout) :: llmsp
    integer, dimension (20, npotd), intent (inout) :: lcore
    integer, dimension (lmxspd, natypd), intent (inout) :: ifunm1
    integer, dimension (20, npotd), intent (inout) :: ititle
    integer, dimension (0:ntotd, natypd), intent (inout) :: ipan_intervall

    real (kind=dp), intent (inout) :: tk
    real (kind=dp), intent (inout) :: emin
    real (kind=dp), intent (inout) :: emax
    real (kind=dp), intent (inout) :: alat
    real (kind=dp), intent (inout) :: r_log
    real (kind=dp), intent (inout) :: efermi

    real (kind=dp), dimension (natypd), intent (inout) :: a
    real (kind=dp), dimension (natypd), intent (inout) :: b
    real (kind=dp), dimension (natypd), intent (inout) :: zat
    real (kind=dp), dimension (natypd), intent (inout) :: conc
    real (kind=dp), dimension (natypd), intent (inout) :: ueff
    real (kind=dp), dimension (natypd), intent (inout) :: jeff
    real (kind=dp), dimension (naezd), intent (inout) :: qmphi
    real (kind=dp), dimension (naezd), intent (inout) :: qmtet
    real (kind=dp), dimension (natypd), intent (inout) :: socscale
    real (kind=dp), dimension (natypd), intent (inout) :: erefldau
    real (kind=dp), dimension (irmd, natypd), intent (inout) :: drdi
    real (kind=dp), dimension (irmd, natypd), intent (inout) :: rmesh
    real (kind=dp), dimension (ncleb, 2), intent (inout) :: cleb
    real (kind=dp), dimension (krel*lmaxd+1, krel*natypd+(1-krel)), intent (inout) :: cscl !! Speed of light scaling
    real (kind=dp), dimension (irmd, npotd), intent (inout) :: visp
    real (kind=dp), dimension (ntotd*(ncheb+1), natypd), intent (inout) :: rnew
    real (kind=dp), dimension (irmd*krel+(1-krel), natypd), intent (inout) :: rmrel !! radial mesh
    real (kind=dp), dimension (irmd*krel+(1-krel), natypd), intent (inout) :: vtrel !! potential (spherical part)
    real (kind=dp), dimension (irmd*krel+(1-krel), natypd), intent (inout) :: btrel !! magnetic field
    real (kind=dp), dimension (20, npotd), intent (inout) :: ecore
    real (kind=dp), dimension (krel*lmaxd+1, krel*natypd+(1-krel)), intent (inout) :: socscl
    real (kind=dp), dimension (irmd*krel+(1-krel), natypd), intent (inout) :: drdirel
    real (kind=dp), dimension (irmd*krel+(1-krel), natypd), intent (inout) :: r2drdirel !! \( r^2 \frac{\partial}{\partial \mathbf{r}}\frac{\partial}{\partial i}\) (r**2 * drdi)
    real (kind=dp), dimension (0:ntotd, natypd), intent (inout) :: rpan_intervall
    real (kind=dp), dimension (irmind:irmd, lmpotd, nspotd), intent (inout) :: vins
    real (kind=dp), dimension (irid, nfund, ncelld), intent (inout) :: thetas
    real (kind=dp), dimension (ntotd*(ncheb+1), nfund, ncelld), intent (inout) :: thetasnew
    real (kind=dp), dimension (mmaxd, mmaxd, nspind, natypd), intent (inout) :: wldau
    real (kind=dp), dimension (mmaxd, mmaxd, mmaxd, mmaxd, natypd), intent (inout) :: uldau
    complex (kind=dp), dimension (iemxd), intent (inout) :: ez
    complex (kind=dp), dimension (iemxd), intent (inout) :: wez
    complex (kind=dp), dimension (irmd, natypd), intent (inout) :: phildau
    character (len=10), intent (inout) :: solver
    character (len=80), intent (inout) :: tmpdir
    ! .. External Functions ..

    nsra = t_params%nsra
    ins = t_params%ins
    natyp = t_params%natyp
    naez = t_params%naez
    nspin = t_params%nspin
    icst = t_params%icst
    ipan = t_params%ipan
    ircut = t_params%ircut
    kmrot = t_params%kmrot
    iqat = t_params%iqat
    conc = t_params%conc
    qmtet = t_params%qmtet
    qmphi = t_params%qmphi
    idoldau = t_params%idoldau
    lmax = t_params%lmax
    irws = t_params%irws
    ! -------------------------------------------------------------------------
    ! Consistency check
    ! -------------------------------------------------------------------------
    if ((krel==1) .and. (ins/=0)) then
      write (6, *) ' FULL-POTENTIAL RELATIVISTIC mode not implemented '
      stop ' set INS = 0 in the input'
    end if

    if (nsra<=2) then
      if (krel==1) stop ' KVREL!= 1 in input, but relativistic program used'
    else
      if (krel==0) stop ' KVREL > 1 in input, but non-relativistic program used'
    end if
    ! -------------------------------------------------------------------------
    alat = t_params%alat
    zat = t_params%zat
    drdi = t_params%drdi
    rmesh = t_params%rmesh
    a = t_params%a
    b = t_params%b
    iend = t_params%iend
    cleb = t_params%cleb
    icleb = t_params%icleb
    loflm = t_params%loflm
    jend = t_params%jend
    thetas = t_params%thetas
    ifunm1 = t_params%ifunm1
    lmsp1 = t_params%lmsp1
    nfu = t_params%nfu
    llmsp = t_params%llmsp
    lcore = t_params%lcore
    ncore = t_params%ncore
    ntcell = t_params%ntcell
    irmin = t_params%irmin
    ititle = t_params%ititle
    intervx = t_params%intervx
    intervy = t_params%intervy
    intervz = t_params%intervz
    nacls1 = t_params%nacls(1)
    lly = t_params%lly
    socscale = t_params%socscale
    tmpdir = t_params%tmpdir
    itmpdir = t_params%itmpdir
    iltmp = t_params%iltmp
    npan_log = t_params%npan_log_at
    npan_eq = t_params%npan_eq_at
    ncheb = t_params%ncheb
    r_log = t_params%r_log
    npan_tot = t_params%npan_tot
    rnew = t_params%rnew
    rpan_intervall = t_params%rpan_intervall
    ipan_intervall = t_params%ipan_intervall
    thetasnew = t_params%thetasnew
    if (krel==1) then
      solver = t_params%solver
      socscl = t_params%socscl
      cscl = t_params%cscl
    end if
    if (idoldau==1) then
      ntldau = t_params%ntldau
      itldau = t_params%itldau
      lopt = t_params%lopt
      ueff = t_params%ueff
      jeff = t_params%jeff
      erefldau = t_params%erefldau
    end if
    ! -------------------------------------------------------------------------
    ! Energy_mesh
    ! -------------------------------------------------------------------------
    ielast = t_params%ielast
    ez = t_params%ez
    wez = t_params%wez
    emin = t_params%emin
    emax = t_params%emax
    iesemicore = t_params%iesemicore
    npol = t_params%npol
    tk = t_params%tk
    if (npol==0) efermi = t_params%efermi
    ! -------------------------------------------------------------------------
    ! Input_potential
    ! -------------------------------------------------------------------------
    vins = t_params%vins
    visp = t_params%visp
    ecore = t_params%ecore
    if (krel==1) then
      rmrel = t_params%rmrel
      drdirel = t_params%drdirel
      r2drdirel = t_params%r2drdirel
      zrel = t_params%zrel
      jwsrel = t_params%jwsrel
      irshift = t_params%irshift
      vtrel = t_params%vtrel
      btrel = t_params%btrel
    end if
    if (use_spherical_potential_only) vins(irmind:irmd, 2:lmpotd, 1:nspotd) = 0.0_dp
    ! -------------------------------------------------------------------------
    ! Itermdir
    ! -------------------------------------------------------------------------
    if (relax_SpinAngle_Dirac) then
      qmtet = t_params%qmtet
      qmphi = t_params%qmphi
    end if
    ! -------------------------------------------------------------------------
    ! LDA+U
    ! -------------------------------------------------------------------------
    if (idoldau==1) then
      itrunldau = t_params%itrunldau
      wldau = t_params%wldau
      uldau = t_params%uldau
      phildau = t_params%phildau
    end if

  end subroutine get_params_1c

  !-------------------------------------------------------------------------------
  !> Summary: Set the values of the local variables according to the stored `t_params` so that they can be passed between different control modules, specifically for `main2`
  !> Author: Philipp Ruessmann
  !> Category: communication, KKRhost
  !> Deprecated: False
  !> Set the values of the local variables according to the stored `t_params`
  !> so that they can be passed between different control modules, specifically for `main2`
  !-------------------------------------------------------------------------------
  subroutine get_params_2(t_params,krel,natyp,ipand,npotd,natomimpd,lmxspd,nfund,   &
    lmpot,ncelld,irmd,nembd1,nembd,irmind,nsra,ins,nspin,ipan,ircut,lcore,ncore,    &
    lmax,ntcell,lpot,nlbasis,nrbasis,nright,nleft,natomimp,atomimp,imix,qbound,fcm, &
    itdbry,irns,kpre,kshape,kte,kvmad,kxc,icc,ishift,ixipol,kforce,ifunm,lmsp,imt,  &
    irc,irmin,irws,llmsp,ititle,nfu,hostimp,ilm_map,imaxsh,ielast,npol,npnt1,npnt2, &
    npnt3,itscf,scfsteps,iesemicore,kaoez,iqat,noq,lly,npolsemi,n1semi,n2semi,      &
    n3semi,zrel,jwsrel,irshift,mixing,lambda_xc,a,b,thetas,drdi,r,zat,rmt,rmtnew,   &
    rws,emin,emax,tk,alat,efold,chrgold,cmomhost,conc,gsh,ebotsemi,emusemi,tksemi,  &
    vins,visp,rmrel,drdirel,vbc,fsold,r2drdirel,ecore,ez,wez,txc,linterface,lrhosym,&
    ngshd,naez,irid,nspotd,iemxd,special_straight_mixing)
    ! get relevant parameters from t_params
    ! ..
    use :: mod_runoptions, only: use_spherical_potential_only
    implicit none
    type (type_params), intent (in) :: t_params

    integer, intent (in) :: irmd
    integer, intent (in) :: irid
    integer, intent (in) :: nembd
    integer, intent (in) :: krel
    integer, intent (in) :: ipand
    integer, intent (in) :: npotd
    integer, intent (in) :: nfund
    integer, intent (in) :: ngshd
    integer, intent (in) :: iemxd
    integer, intent (in) :: lmxspd
    integer, intent (in) :: ncelld
    integer, intent (in) :: nembd1
    integer, intent (in) :: irmind
    integer, intent (in) :: nspotd
    integer, intent (in) :: natomimpd
    integer, intent (inout) :: kte
    integer, intent (inout) :: kxc
    integer, intent (inout) :: icc
    integer, intent (inout) :: lly
    integer, intent (inout) :: ins
    integer, intent (inout) :: nsra
    integer, intent (inout) :: imix
    integer, intent (inout) :: lpot
    integer, intent (inout) :: kpre
    integer, intent (inout) :: lmax
    integer, intent (inout) :: naez
    integer, intent (inout) :: npol
    integer, intent (inout) :: npnt1
    integer, intent (inout) :: npnt2
    integer, intent (inout) :: npnt3
    integer, intent (inout) :: itscf
    integer, intent (inout) :: nspin
    integer, intent (inout) :: natyp
    integer, intent (inout) :: kvmad
    integer, intent (inout) :: lmpot
    integer, intent (inout) :: nleft
    integer, intent (inout) :: nright
    integer, intent (inout) :: ishift
    integer, intent (inout) :: kforce
    integer, intent (inout) :: n1semi
    integer, intent (inout) :: n2semi
    integer, intent (inout) :: n3semi
    integer, intent (inout) :: ielast
    integer, intent (inout) :: kshape
    integer, intent (inout) :: itdbry
    integer, intent (inout) :: nlbasis
    integer, intent (inout) :: nrbasis
    integer, intent (inout) :: natomimp
    integer, intent (inout) :: scfsteps
    integer, intent (inout) :: npolsemi
    integer, intent (inout) :: iesemicore
    integer, intent (inout) :: special_straight_mixing
    integer, dimension (naez), intent (inout) :: noq
    integer, dimension (natyp), intent (inout) :: imt
    integer, dimension (natyp), intent (inout) :: irc
    integer, dimension (natyp), intent (inout) :: nfu
    integer, dimension (natyp), intent (inout) :: iqat
    integer, dimension (natyp), intent (inout) :: zrel
    integer, dimension (natyp), intent (inout) :: irws
    integer, dimension (natyp), intent (inout) :: ipan
    integer, dimension (natyp), intent (inout) :: irns
    integer, dimension (npotd), intent (inout) :: ncore
    integer, dimension (natyp), intent (inout) :: irmin
    integer, dimension (0:lmpot), intent (inout) :: imaxsh
    integer, dimension (natyp), intent (inout) :: jwsrel
    integer, dimension (natyp), intent (inout) :: ixipol
    integer, dimension (natyp), intent (inout) :: ntcell
    integer, dimension (0:natyp), intent (inout) :: hostimp
    integer, dimension (natyp), intent (inout) :: irshift
    integer, dimension (natomimpd), intent (inout) :: atomimp
    integer, dimension (ngshd, 3), intent (inout) :: ilm_map
    integer, dimension (natyp, lmxspd), intent (inout) :: lmsp
    integer, dimension (20, npotd), intent (inout) :: lcore
    integer, dimension (natyp, nfund), intent (inout) :: llmsp
    integer, dimension (natyp, lmxspd), intent (inout) :: ifunm
    integer, dimension (natyp, naez+nembd), intent (inout) :: kaoez
    integer, dimension (0:ipand, natyp), intent (inout) :: ircut
    integer, dimension (20, npotd), intent (inout) :: ititle
    real (kind=dp), intent (inout) :: tk
    real (kind=dp), intent (inout) :: fcm
    real (kind=dp), intent (inout) :: emin
    real (kind=dp), intent (inout) :: emax
    real (kind=dp), intent (inout) :: alat
    real (kind=dp), intent (inout) :: efold
    real (kind=dp), intent (inout) :: fsold
    real (kind=dp), intent (inout) :: qbound
    real (kind=dp), intent (inout) :: tksemi
    real (kind=dp), intent (inout) :: mixing
    real (kind=dp), intent (inout) :: chrgold
    real (kind=dp), intent (inout) :: emusemi
    real (kind=dp), intent (inout) :: ebotsemi
    real (kind=dp), intent (inout) :: lambda_xc
    real (kind=dp), dimension (natyp), intent (inout) :: a
    real (kind=dp), dimension (natyp), intent (inout) :: b
    real (kind=dp), dimension (2), intent (inout) :: vbc
    real (kind=dp), dimension (natyp), intent (inout) :: rws
    real (kind=dp), dimension (ngshd), intent (inout) :: gsh
    real (kind=dp), dimension (natyp), intent (inout) :: zat
    real (kind=dp), dimension (natyp), intent (inout) :: rmt
    real (kind=dp), dimension (natyp), intent (inout) :: conc
    real (kind=dp), dimension (natyp), intent (inout) :: rmtnew
    real (kind=dp), dimension (irmd, natyp), intent (inout) :: r
    real (kind=dp), dimension (irmd, natyp), intent (inout) :: drdi
    real (kind=dp), dimension (irmd, npotd), intent (inout) :: visp
    real (kind=dp), dimension (20, npotd), intent (inout) :: ecore
    real (kind=dp), dimension (irmd*krel+(1-krel), natyp), intent (inout) :: rmrel
    real (kind=dp), dimension (irmd*krel+(1-krel), natyp), intent (inout) :: drdirel
    real (kind=dp), dimension (irmd*krel+(1-krel), natyp), intent (inout) :: r2drdirel
    real (kind=dp), dimension (lmpot, nembd1), intent (inout) :: cmomhost
    real (kind=dp), dimension (irmind:irmd, lmpot, nspotd), intent (inout) :: vins
    real (kind=dp), dimension (irid, nfund, ncelld), intent (inout) :: thetas
    complex (kind=dp), dimension (iemxd), intent (inout) :: ez
    complex (kind=dp), dimension (iemxd), intent (inout) :: wez
    character (len=124), dimension (6), intent (inout) :: txc
    logical, intent (inout) :: lrhosym
    logical, intent (inout) :: linterface
    ! .. External Functions ..

    nsra = t_params%nsra
    ins = t_params%ins
    natyp = t_params%natyp
    naez = t_params%naez
    nspin = t_params%nspin
    ipan = t_params%ipan
    ircut = t_params%ircut
    lcore = t_params%lcore
    ncore = t_params%ncore
    ntcell = t_params%ntcell
    lmax = t_params%lmax
    lpot = t_params%lpot
    lmpot = t_params%lmpot
    nlbasis = t_params%nlbasis
    nrbasis = t_params%nrbasis
    nright = t_params%nright
    nleft = t_params%nleft
    linterface = t_params%linterface
    atomimp = t_params%atomimp
    natomimp = t_params%natomimp

    ! -------------------------------------------------------------------------
    ! Consistency check
    ! -------------------------------------------------------------------------
    if ((krel==1) .and. (ins/=0)) then
      write (6, *) ' FULL-POTENTIAL RELATIVISTIC mode not implemented '
      stop ' set INS = 0 in the input'
    end if
    if (nsra<=2) then
      if (krel==1) stop ' KVREL!= 1 in input, but relativistic program used'
    else
      if (krel==0) stop ' KVREL > 1 in input, but non-relativistic program used'
    end if
    ! -------------------------------------------------------------------------
    imix = t_params%imix
    special_straight_mixing = t_params%special_straight_mixing
    mixing = t_params%mixing
    qbound = t_params%qbound
    fcm = t_params%fcm
    itdbry = t_params%itdbry
    irns = t_params%irns
    kpre = t_params%kpre
    kshape = t_params%kshape
    kte = t_params%kte
    kvmad = t_params%kvmad
    kxc = t_params%kxc
    lambda_xc = t_params%lambda_xc
    txc = t_params%txc
    icc = t_params%icc
    ishift = t_params%ishift
    ixipol = t_params%ixipol
    lrhosym = t_params%lrhosym
    kforce = t_params%kforce
    a = t_params%a
    b = t_params%b
    drdi = t_params%drdi
    r = t_params%rmesh
    thetas = t_params%thetas
    zat = t_params%zat
    ifunm = t_params%ifunm
    lmsp = t_params%lmsp
    rmt = t_params%rmt
    rmtnew = t_params%rmtnew
    rws = t_params%rws
    imt = t_params%imt
    irc = t_params%irc
    irmin = t_params%irmin
    irws = t_params%irws
    ititle = t_params%ititle
    llmsp = t_params%llmsp
    nfu = t_params%nfu
    hostimp = t_params%hostimp
    alat = t_params%alat
    kaoez = t_params%kaoez
    iqat = t_params%iqat
    noq = t_params%noq
    conc = t_params%conc
    gsh = t_params%gsh
    ilm_map = t_params%ilm_map
    imaxsh = t_params%imaxsh
    lly = t_params%lly
    ! -------------------------------------------------------------------------
    ! Energy_mesh
    ! -------------------------------------------------------------------------
    ielast = t_params%ielast
    ez = t_params%ez
    wez = t_params%wez
    emin = t_params%emin
    emax = t_params%emax
    iesemicore = t_params%iesemicore
    fsold = t_params%fsemicore
    npol = t_params%npol
    tk = t_params%tk
    npnt1 = t_params%npnt1
    npnt2 = t_params%npnt2
    npnt3 = t_params%npnt3
    ebotsemi = t_params%ebotsemi
    emusemi = t_params%emusemi
    tksemi = t_params%tksemi
    npolsemi = t_params%npolsemi
    n1semi = t_params%n1semi
    n2semi = t_params%n2semi
    n3semi = t_params%n3semi
    ! -------------------------------------------------------------------------
    ! Input_potential
    ! -------------------------------------------------------------------------
    vins = t_params%vins
    visp = t_params%visp
    ecore = t_params%ecore
    vbc = t_params%vbc
    if (krel==1) then
      rmrel = t_params%rmrel
      drdirel = t_params%drdirel
      r2drdirel = t_params%r2drdirel
      zrel = t_params%zrel
      jwsrel = t_params%jwsrel
      irshift = t_params%irshift
    end if
    itscf = t_params%itscf
    scfsteps = t_params%scfsteps
    efold = t_params%efold
    chrgold = t_params%chrgold
    cmomhost = t_params%cmomhost
    if (use_spherical_potential_only) vins(irmind:irmd, 2:lmpot, 1:nspotd) = 0.0_dp

  end subroutine get_params_2

  !-------------------------------------------------------------------------------
  !> Summary: Store the values of the local variables related to the energy mesh, in the `t_params` data types
  !> Author: Philipp Ruessmann
  !> Category: communication, KKRhost
  !> Deprecated: False
  !> Store the values of the local variables related to the energy mesh,
  !> in the `t_params` data types
  !-------------------------------------------------------------------------------
  subroutine save_emesh(ielast,ez,wez,emin,emax,iesemicore,fsemicore,npol,tk,npnt1, &
    npnt2,npnt3,ebotsemi,emusemi,tksemi,npolsemi,n1semi,n2semi,n3semi,iemxd,t_params)
    ! save information of energy mesh in t_params
    implicit none

    type (type_params), intent (inout) :: t_params

    integer, intent (in) :: npol
    integer, intent (in) :: npnt1
    integer, intent (in) :: npnt2
    integer, intent (in) :: npnt3
    integer, intent (in) :: iemxd
    integer, intent (in) :: ielast
    integer, intent (in) :: n1semi
    integer, intent (in) :: n2semi
    integer, intent (in) :: n3semi
    integer, intent (in) :: npolsemi
    integer, intent (in) :: iesemicore
    real (kind=dp), intent (in) :: tk
    real (kind=dp), intent (in) :: emin
    real (kind=dp), intent (in) :: emax
    real (kind=dp), intent (in) :: tksemi
    real (kind=dp), intent (in) :: emusemi
    real (kind=dp), intent (in) :: ebotsemi
    real (kind=dp), intent (in) :: fsemicore
    complex (kind=dp), dimension (iemxd), intent (in) :: ez
    complex (kind=dp), dimension (iemxd), intent (in) :: wez

    t_params%ielast = ielast
    t_params%ez = ez
    t_params%wez = wez
    t_params%emin = emin
    t_params%emax = emax
    t_params%iesemicore = iesemicore
    t_params%fsemicore = fsemicore
    t_params%npol = npol
    t_params%tk = tk
    t_params%npnt1 = npnt1
    t_params%npnt2 = npnt2
    t_params%npnt3 = npnt3
    t_params%ebotsemi = ebotsemi
    t_params%emusemi = emusemi
    t_params%tksemi = tksemi
    t_params%npolsemi = npolsemi
    t_params%n1semi = n1semi
    t_params%n2semi = n2semi
    t_params%n3semi = n3semi

  end subroutine save_emesh

  !-------------------------------------------------------------------------------
  !> Summary: Store the values of the local variables related to the SCF parameters in the `t_params` data types
  !> Author: Philipp Ruessmann
  !> Category: communication, KKRhost
  !> Deprecated: False
  !> Store the values of the local variables related to the SCF parameters
  !> in the `t_params` data types. Save information that is needed in next iteration
  !> and that is changeing, i.e. potential etc.
  !-------------------------------------------------------------------------------
  subroutine save_scfinfo(t_params,vins,visp,ecore,vbc,rmrel,drdirel,r2drdirel,zrel,&
    jwsrel,irshift,vtrel,btrel,itscf,scfsteps,efold,chrgold,cmomhost,krel,irmind,   &
    irm,lmpot,nspotd,natyp,npotd,nembd1)
    ! save information that is needed in next iteration and that is changeing, i.e. potential etc.
    implicit none

    type (type_params), intent (inout) :: t_params

    integer, intent (in) :: irm
    integer, intent (in) :: krel
    integer, intent (in) :: lmpot
    integer, intent (in) :: natyp
    integer, intent (in) :: npotd
    integer, intent (in) :: itscf
    integer, intent (in) :: nembd1
    integer, intent (in) :: nspotd
    integer, intent (in) :: irmind
    integer, intent (in) :: scfsteps
    integer, dimension (natyp), intent (in) :: zrel
    integer, dimension (natyp), intent (in) :: jwsrel
    integer, dimension (natyp), intent (in) :: irshift
    real (kind=dp), intent (in) :: efold
    real (kind=dp), intent (in) :: chrgold
    real (kind=dp), dimension (2), intent (in) :: vbc
    real (kind=dp), dimension (irm, npotd), intent (in) :: visp
    real (kind=dp), dimension (irm*krel+(1-krel), natyp), intent (in) :: vtrel
    real (kind=dp), dimension (irm*krel+(1-krel), natyp), intent (in) :: btrel
    real (kind=dp), dimension (20, npotd), intent (in) :: ecore
    real (kind=dp), dimension (irm*krel+(1-krel), natyp), intent (in) :: rmrel
    real (kind=dp), dimension (irm*krel+(1-krel), natyp), intent (in) :: drdirel
    real (kind=dp), dimension (lmpot, nembd1), intent (in) :: cmomhost
    real (kind=dp), dimension (irm*krel+(1-krel), natyp), intent (in) :: r2drdirel
    real (kind=dp), dimension (irmind:irm, lmpot, nspotd), intent (in) :: vins

    t_params%vins = vins
    t_params%visp = visp
    t_params%ecore = ecore
    t_params%vbc = vbc
    if (krel==1) then
      t_params%rmrel = rmrel
      t_params%drdirel = drdirel
      t_params%r2drdirel = r2drdirel
      t_params%zrel = zrel
      t_params%jwsrel = jwsrel
      t_params%irshift = irshift
      t_params%vtrel = vtrel
      t_params%btrel = btrel
    end if
    t_params%itscf = itscf
    t_params%scfsteps = scfsteps
    t_params%efold = efold
    t_params%chrgold = chrgold
    t_params%cmomhost = cmomhost
  end subroutine save_scfinfo

  !-------------------------------------------------------------------------------
  !> Summary: Store the values of the local variables related to the electronic density in the `t_params` data types
  !> Author: Philipp Ruessmann
  !> Category: communication, physical-observables, KKRhost
  !> Deprecated: False
  !> Store the values of the local variables related to the electronic density
  !> in the `t_params` data types. Save density after it has been calculated in
  !> `main1c`, is further processed in `main2`
  !-------------------------------------------------------------------------------
  subroutine save_density(t_params,rho2ns,r2nef,rhoc,denef,denefat,espv,ecore,      &
    idoldau,lopt,eu,edc,chrgsemicore,rhoorb,ecorerel,nkcore,kapcore,krel,natyp,     &
    npotd,irm,lmpot,lmaxd1)

    implicit none

    type (type_params), intent (inout) :: t_params

    integer, intent (in) :: irm
    integer, intent (in) :: krel
    integer, intent (in) :: natyp
    integer, intent (in) :: npotd
    integer, intent (in) :: lmpot
    integer, intent (in) :: lmaxd1
    integer, intent (in) :: idoldau
    integer, dimension (natyp), intent (in) :: lopt
    integer, dimension (20, natyp), intent (in) :: nkcore
    integer, dimension (20, npotd), intent (in) :: kapcore
    real (kind=dp), intent (in) :: denef
    real (kind=dp), intent (in) :: chrgsemicore
    real (kind=dp), intent (in) :: eu(natyp)
    real (kind=dp), intent (in) :: edc(natyp)
    real (kind=dp), intent (in) :: denefat(natyp)
    real (kind=dp), intent (in) :: rhoc(irm, npotd)
    real (kind=dp), intent (in) :: espv(0:lmaxd1, npotd)
    real (kind=dp), intent (in) :: ecore(20, npotd)
    real (kind=dp), intent (in) :: rhoorb(irm*krel+(1-krel), natyp)
    real (kind=dp), intent (in) :: ecorerel(krel*20+(1-krel), npotd)
    real (kind=dp), intent (in) :: r2nef(irm, lmpot, natyp, 2)
    real (kind=dp), intent (in) :: rho2ns(irm, lmpot, natyp, 2)

    t_params%rho2ns = rho2ns
    t_params%r2nef = r2nef
    t_params%rhoc = rhoc
    t_params%denef = denef
    t_params%denefat = denefat
    t_params%espv = espv
    t_params%ecore = ecore
    t_params%idoldau = idoldau
    t_params%lopt = lopt
    t_params%eu = eu
    t_params%edc = edc
    t_params%chrgsemicore = chrgsemicore
    if (krel==1) then
      t_params%rhoorb = rhoorb
      t_params%ecorerel = ecorerel
      t_params%nkcore = nkcore
      t_params%kapcore = kapcore
    end if

  end subroutine save_density

  !-------------------------------------------------------------------------------
  !> Summary: Set the values of the local variables related to the electronic density in the `t_params` data types
  !> Author: Philipp Ruessmann
  !> Category: communication, physical-observables, KKRhost
  !> Deprecated: False
  !> Set the values of the local variables related to the electronic density
  !> in the `t_params` data types. Store the values of the density in `main2`
  !-------------------------------------------------------------------------------
  subroutine read_density(t_params,rho2ns,r2nef,rhoc,denef,denefat,espv,ecore,      &
    idoldau,lopt,eu,edc,chrgsemicore,rhoorb,ecorerel,nkcore,kapcore,krel,natyp,     &
    npotd,irm,lmpot,lmaxd1)

    implicit none

    type (type_params), intent (inout) :: t_params

    integer, intent (in) :: irm
    integer, intent (in) :: krel
    integer, intent (in) :: natyp
    integer, intent (in) :: npotd
    integer, intent (in) :: lmpot
    integer, intent (in) :: lmaxd1
    integer, intent (out) :: idoldau
    integer, dimension (natyp), intent (out) :: lopt
    integer, dimension (20, natyp), intent (out) :: nkcore
    integer, dimension (20, npotd), intent (out) :: kapcore
    real (kind=dp), intent (out) :: rho2ns(irm, lmpot, natyp, 2)
    real (kind=dp), intent (out) :: denef
    real (kind=dp), intent (out) :: chrgsemicore
    real (kind=dp), dimension (natyp), intent (out) :: eu
    real (kind=dp), dimension (natyp), intent (out) :: edc
    real (kind=dp), dimension (natyp), intent (out) :: denefat
    real (kind=dp), dimension (irm, npotd), intent (out) :: rhoc
    real (kind=dp), dimension (0:lmaxd1, npotd), intent (out) :: espv
    real (kind=dp), dimension (20, npotd), intent (out) :: ecore
    real (kind=dp), dimension (irm*krel+(1-krel), natyp), intent (out) :: rhoorb
    real (kind=dp), dimension (krel*20+(1-krel), npotd), intent (out) :: ecorerel
    real (kind=dp), dimension (irm, lmpot, natyp, 2), intent (out) :: r2nef

    rho2ns = t_params%rho2ns
    r2nef = t_params%r2nef
    rhoc = t_params%rhoc
    denef = t_params%denef
    denefat = t_params%denefat
    espv = t_params%espv
    ecore = t_params%ecore
    idoldau = t_params%idoldau
    lopt = t_params%lopt
    eu = t_params%eu
    edc = t_params%edc
    chrgsemicore = t_params%chrgsemicore

    if (krel==1) then
      rhoorb = t_params%rhoorb
      ecorerel = t_params%ecorerel
      nkcore = t_params%nkcore
      kapcore = t_params%kapcore
    end if

  end subroutine read_density

  !-------------------------------------------------------------------------------
  !> Summary: Read the angles variables associated with the angles of magnetic moments in a non-collinear calculation
  !> Author: Philipp Ruessmann
  !> Category: input-output, dirac, KKRhost
  !> Deprecated: False
  !>  Read the angles variables associated with the angles of magnetic
  !> moments in a non-collinear calculation. Read `nonco_angles`.
  !-------------------------------------------------------------------------------
  subroutine read_angles(t_params, natyp, theta, phi)
    ! read nonco_angles
    use :: mod_types, only: t_inc
    use :: mod_mympi, only: myrank, master
    use :: mod_version_info
    use :: mod_constants, only: pi

    implicit none

    type (type_params), intent (inout) :: t_params

    integer, intent (in) :: natyp
    real (kind=dp), dimension (natyp), intent (out) :: theta
    real (kind=dp), dimension (natyp), intent (out) :: phi

    logical :: lread, lcheckangles
    integer :: i1, i_stat
    real (kind=dp) :: th1, ph1
    real (kind=dp), parameter :: eps = 1.0e-5_dp
    ! if executed first in wunfiles theta is not allocated, thus read angles from file
    if (.not. allocated(t_params%theta)) then

      theta(:) = 0.0_dp
      phi(:) = 0.0_dp
      lread = .false.
      lcheckangles = .false.
      inquire (file='nonco_angle.dat', exist=lread)
      if (lread) then
        open (unit=10, file='nonco_angle.dat', form='FORMATTED')
        call version_check_header(10)
        do i1 = 1, natyp
          read (10, *) th1, ph1
          if ((abs(th1)<(pi+eps) .and. abs(th1)>eps) .or. (abs(ph1)<(2*pi+eps) .and. abs(ph1)>eps)) then
            lcheckangles = .true.
          end if
          theta(i1) = th1*(pi/180.0_dp)
          phi(i1) = ph1*(pi/180.0_dp)
        end do
        close (10)
        if (lcheckangles .and. ((t_inc%i_write>0) .or. (myrank==master))) then
          write (1337, *) 'WARNING: Check if your nonco_angels file is correct! Found only values that are smaller than pi for theta and 2pi for &
            &phi, respectively. But angles are given in degree (0...360)'
        end if
        write (1337, '(A)') '      I1  THETA[deg]  PHI[deg]'
        do i1 = 1, natyp
          write (1337, '(I8,2F12.6)') i1, theta(i1)*180.0_dp/pi, phi(i1)*180.0_dp/pi
        end do                     ! i1
      end if                       ! LREAD

      ! now save this also to t_params
      allocate (t_params%theta(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(t_params%theta))*kind(t_params%theta), 't_params%THETA', 'read_angles')

      allocate (t_params%phi(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(t_params%phi))*kind(t_params%phi), 't_params%PHI', 'read_angles')

      t_params%theta = theta
      t_params%phi = phi

    else                           ! not first run: information saved in t_params

      theta = t_params%theta
      phi = t_params%phi

    end if

  end subroutine read_angles

end module mod_wunfiles
