module projection
! Module that handles tasks related to the projected GFs
! Most global variables are private, not visible outside the module
! Public are: which atoms will enter the susceptibility calculations
!             how many basis functions for each l for each atom
!             at what energies should they be computed
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!              NOTE: special treatment of the radial mesh is done
!       ir = 1 not used, remember to CHANGE this for other implementations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!                         transverse    longitudinal
! Cartesian: i=1,4  <=>      x,     y,     z,      n
! Spin:      i=1,4  <=>  up dn, dn up, up up, dn, dn
!
! The internal storage is done in spin density matrix form for each site
! Whenever a given quantity is required conversion to x,y,z,n form is done with the Pauli matrices
! (G^\dagger - G)/(i2pi) plays the role of the spectral spin density
!
!
! List of subroutines:
! --> susc_input	DONE	reads inpsusc.dat
! --> outsusc_out	DONE    opens outsusc.dat and writes initial data
! --> outsusc_in	DONE    opens outsusc.dat and reads all of it, as well as the pots and wfns files
! --> close_outsusc	DONE	closes outsusc.dat
! --> init_param	DONE	initializes global variables
! --> init_arrays       DONE    allocates global arrays
! --> ymy_gaunts        DONE    read angular mesh from file, compute spherical harmonics and Gaunt numbers
! --> save_rmesh	DONE	stores the radial mesh, radial integration weight and potentials on a file
! --> read_rmesh	DONE	reads the radial mesh, radial integration weight and potentials from a file
! --> ref_wfn		DONE	stores a wavefunction to use as basis function in memory
! --> new_basis1	DONE	calls find_basis for a set of input basis functions for fixed l and s
! --> new_basis2	DONE	calls find_basis for a set of input basis functions for fixed l, both spins
! --> new_basis3	DONE	calls find_basis for a set of input basis functions for all l and s
! --> find_basis	DONE	Gram-Schmidt orthonormalization on input basis functions
! --> find_basis2	DONE	diagonalization of overlap matrix of input basis functions
! --> reg_coeffs	DONE	projection coefficients of regular scattering solution
! --> irr_coeffs	DONE	projection coefficients of irregular scattering solution
! --> out_coeffs	DONE	output projection coefficients to outsusc.dat
! --> in_coeffs		DONE	get projection coefficients from outsusc.dat
! --> out_gmat		DONE	output structural GF to outsusc.dat
! --> in_gmat		DONE	get structural GF from outsusc.dat
! --> store_gscoll	DONE	save structural GF in RAM
! --> save_wfns		DONE	saves basis functions to file
! --> read_wfns		DONE	reads basis functions from file
! --> out_overlap	DONE	output overlaps between different basis functions
! --> radint		DONE	integration on radial mesh
! --> orbmoment         DONE    initializes the matrix elements of the angular momentum operator and Pauli matrices
! --> projected_gf      DONE    puts together a GF block for given energy
! --> full_susc         DONE    computes the static Kohn-Sham susceptibility
! --> groundstate       TODO    valence charge, spin and orbital moments, densities, static KS susc
! --> soc_correction    DONE    computes the SOC kernel, the change in the t-matrices and updates the GF
! --> build_vsoc        DONE    constructs the radial SOC potential
! --> build_vsocb       DONE    constructs the matrix elements of the SOC potential in the basis
! --> build_vscf        TODO    constructs the matrix elements of the SCF potential in the basis

  implicit none


!  Kind definitions
!  different compilers might handle kind values differently
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

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     Global parameters that are read in by suscmain
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! --> constant shift of Re GF; fudge factor
  integer(kind=i4b), public  :: iregf = 0
  real(kind=r8b),    public  :: fudge = 0.d0
! --> DOS calculation; energy parameters
  integer(kind=i4b), public  :: idos = 0, nedos = 0
  real(kind=r8b),    public  :: e0dos = 0.d0, e1dos = 0.d0, eimdos = 0.d0
! --> which part of the susceptibility to compute
  integer(kind=i4b), public  :: isusc = 0
  integer(kind=i4b), public  :: nomega = 0
  real(kind=r8b),    public  :: omegamin = 0.d0, omegamax = 0.d0, domega = 0.d0
! --> whether to correct for SOC
  integer(kind=i4b), public  :: isoc = 0
  logical,           private :: soc_applied = .false.
! --> whether to compute the Hartree kernel
  integer(kind=i4b), public  :: ikha = 0
! --> whether and which xc kernels should be computed
  integer(kind=i4b), public  :: ikxc = 0
! --> cartesian or spin representation
  logical,           public  :: cartesian = .false.
! --> compute analytic and nonanalytic integrals
  logical,           public  :: analytic  = .false., nonanalytic = .false.
! --> KS or full susc
  logical,           public  :: enhanced  = .false.
! --> iterations to correct the denominator
  integer(kind=i4b), public  :: itermax = 0
! --> mixing for lambda min
  real(kind=r8b),    public  :: lambdamix = 0.d0
! --> degrees of numerator and denominator in GF fit; shift of energy
  integer(kind=i4b), public  :: numd = 0, dend = 0
  complex(kind=c8b), public  :: eshift = (0.d0,0.d0)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     Global parameters that must be initialized with init_param
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! --> number of energy points for susc output
  integer(kind=i4b), public  :: nesusc = 0
! --> number of atoms for susc output
  integer(kind=i4b), public  :: nasusc = 0
! --> maximum number of radial basis functions per l, s and atom
  integer(kind=i4b), public  :: nbmax = 0
! --> maximum number of points in radial meshes
  integer(kind=i4b), private :: nrmax = 0
! --> maximum angular momentum components for basis (one per l-channel)
  integer(kind=i4b), public  :: nlmax = -1, nlmax2 = -1, nlmax4 = -1
! --> total number of angular momentum components
  integer(kind=i4b), public  :: lmmax = 0, lmmax2 = 0, lmmax4 = 0
! --> total number of angular momentum components for susceptibility
  integer(kind=i4b), public  :: nlmax0 = -1, lmmax0 = 0
! --> maximum number of spin components
  integer(kind=i4b), private :: nsmax = 0
! --> maximum angular momentum in Lebedev Laikov mesh
  integer(kind=i4b), private :: nlmaxll = 0, lmmaxll = 0
! --> number of angular mesh points
  integer(kind=i4b), private :: nll = 0
! --> angular momentum and spin
  integer(kind=i4b), public  :: nlms  = 0
! --> angular momentum, spin and basis functions
  integer(kind=i4b), public  :: nlmsb = 0
! --> dimensions for storage of structural GF
  integer(kind=i4b), private :: nalms = 0
! --> combined index for the susceptibility
  integer(kind=i4b), private :: nalmsb = 0
! --> whether to keep the small components of the SRA
  integer(kind=i4b), public  :: isra  = -1
! --> this checks if the parameters were initialized
! --> init_param sets this to true
  logical,           private :: noparameters = .false.
! --> unit number for I/O, change if there's a problem
  integer(kind=i4b), private, parameter :: io = 99
! --> same, for debugging
  integer(kind=i4b), private, parameter :: iodb = 666
! --> unit number for opening, reading/writing and closing a file in a single subroutine
  integer(kind=i4b), private, parameter :: iofile = 300
! --> useful arrays to read and write susc data
!     ngroup <= nasusc says how many groups of similar atoms there are
!     igroup(1:ngroup) says how many atoms there are in each group
  integer(kind=i4b), private :: ngroup, igroup(1000)


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! This checks if the Gram-Schmidt process has found a zero
  real(kind=r8b),    private, parameter :: gstol = 1.d-10
! This checks if the structural GF element is to be set to zero
  real(kind=r8b),    private, parameter :: gfilter = 1.d-8
! This checks if the KS susc element is to be set to zero
  real(kind=r8b),    private, parameter :: susctol = 1.d-6
! This decides about non-zero Gaunt numbers
  real(kind=r8b),    private, parameter :: ylmtol = 1.d-10
! This decides about non-zero values in the observables
  real(kind=r8b),    private, parameter :: atol = 1.d-8
! This drops small elements of the SOC wfn
  real(kind=r8b),    private, parameter :: soctol = 1.d-6
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     Global arrays that are allocated with init_param
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Numbering of atoms:
!   # atoms = # host atoms + # atoms in cluster
!   natypd  = natref       + natomd
! --> which atoms for susc, as read from potential file
  integer(kind=i4b), public,  allocatable :: iasusc(:)
! --> treat the atom as magnetic or non-magnetic
  integer(kind=i4b), private, allocatable :: issusc(:)
! --> which wfns to use for projection (s,p,d,f) for each atom
!     n = 0 means don't project
!     n > 0 means project using n regular solutions computed at ewsusc
  integer(kind=i4b), public,  allocatable :: iwsusc(:,:,:)   ! l only
! --> energy for wfn calculation for projection
  complex(kind=c8b), public,  allocatable :: ewsusc(:,:,:,:)
! --> storage for reference regular wfns
  real(kind=r8b),    private, allocatable :: phiref(:,:,:,:,:), overlap(:,:,:)
  logical,           private, allocatable :: nowfns(:,:,:), nobasis(:,:,:)
! --> angular mesh, spherical harmonics, gaunt numbers
  real(kind=r8b),    private, allocatable :: ull(:,:), wll(:), cthetall(:), phill(:), ylm(:,:), gaunt(:,:,:)
! --> radial mesh and weights for radial integration with reference wfns
! --> the weight must ensure that the initial wnfs are integrable
  integer(kind=i4b), private, allocatable :: nrpts(:), nrpts0(:), nrpts1(:)
  real(kind=r8b),    private, allocatable :: rmesh(:,:), rsmesh(:,:,:), drmesh(:,:), drproj(:,:)
  logical,           private, allocatable :: normesh(:)
! --> info needed about the atoms
! --> atomic number
  real(kind=r8b),    private, allocatable :: zat(:)
! --> direction of magnetization
  real(kind=r8b),    private, allocatable :: magdir(:,:)
! --> scalar and magnetic potentials (nuclear pot not included)
  real(kind=r8b),    private, allocatable :: vr(:,:), br(:,:)
! --> potential in the lmsb basis
  complex(kind=c8b), private, allocatable :: vlmsb(:,:,:)
! --> core charge and magnetisation densities
  real(kind=r8b),    private, allocatable :: nrc(:,:), mrc(:,:)
! --> valence charge and magnetisation densities
  real(kind=r8b),    private, allocatable :: nrv(:,:,:)
! --> energy mesh for susc calculation
  complex(kind=c8b), public,  allocatable :: ensusc(:)
! --> coefficients of the projected scattering solutions
  complex(kind=c8b), private, allocatable :: pzc(:,:,:,:), fzc(:,:,:,:)
! --> the coefficients of the products of regular and irregular wfns
! --> these are computed in the impurity program, in ASA mode
  complex(kind=c8b), private, allocatable :: pqc(:,:,:,:,:), psc(:,:,:,:,:), fqc(:,:,:,:,:), fsc(:,:,:,:,:)
  logical,           private, allocatable :: noregcoeffs(:,:), noirrcoeffs(:,:)
! --> storage for the structural GF
  complex(kind=c8b), private, allocatable :: gstruct(:,:,:)
! --> storage for fitted GF
  complex(kind=c8b), private, allocatable :: gffit(:,:,:,:,:)
! --> pointers from lm to storage and back (gaunt numbers)
  integer(kind=i4b), public,  allocatable :: lm2i(:,:), i2lm(:,:)
! --> pointers from lms to storage and back
  integer(kind=i4b), public,  allocatable :: lms2i(:,:), i2lms(:,:)
! --> pointers from lmsb to storage and back
  integer(kind=i4b), public,  allocatable :: lmsb2i(:,:,:,:), i2lmsb(:,:,:)
! --> sizes of blocks of projected GF
  integer(kind=i4b), public,  allocatable :: nlmsba(:)
! --> pointers from alms to storage and back
  integer(kind=i4b), public,  allocatable :: alms2i(:,:), i2alms(:,:)
! --> pointers from almsb to storage and back
  integer(kind=i4b), public,  allocatable :: almsb2i(:,:,:), i2almsb(:,:)
! --> coefficients of the projected scattering solutions
  complex(kind=c8b), private, allocatable :: pzl(:,:,:,:), pzr(:,:,:,:), fzl(:,:,:,:), fzr(:,:,:,:)
! --> storage for single-site Green functions
  complex(kind=c8b), private, allocatable :: gfpq(:,:,:,:), gfps(:,:,:,:), gffq(:,:,:,:), gffs(:,:,:,:)
! --> energy mesh
  complex(kind=c8b), private, allocatable :: esusc(:), eksusc(:), desusc(:), wsusc(:)
! --> to compute SCF-like quantities
  integer(kind=i4b), private :: nescf
  complex(kind=c8b), private, allocatable :: escf(:), ekscf(:), descf(:)
! --> Fermi energy
  real(kind=r8b),    private :: efermi = 0.d0
! --> Pauli matrices
  complex(kind=c8b), private :: pauli(2,2,4)
! --> transformation matrices from cartesian to spin labels and back
! potential from cartesian to spin
  complex(kind=c8b), private :: pc2s(2,2,4)
! potential from spin to cartesian
  complex(kind=c8b), private :: ps2c(4,2,2)
! density from cartesian to spin
  complex(kind=c8b), private :: dc2s(2,2,4)
! density from spin to cartesian
  complex(kind=c8b), private :: ds2c(4,2,2)
! --> angular momentum matrices
  complex(kind=c8b), private, allocatable :: lorb(:,:,:)
! --> matrix elements for L.S
  complex(kind=c8b), private, allocatable :: ldots(:,:,:)
! --> storage for kronecker form of KS susc
  complex(kind=c8b), private, allocatable :: kssusc(:,:)
! --> storage for kronecker form of xc kernel
  complex(kind=c8b), private, allocatable :: kxcsusc(:,:)
! --> storage for kronecker form of susc denominator
  complex(kind=c8b), private, allocatable :: denominator(:,:)


  contains


!----------------------------------------------------------------------
  subroutine susc_input(lmax,natyp,nspin,nrad)
! Handles the input file

  implicit none

! --> values given in the main program
  integer(kind=i4b), intent(in) :: lmax, natyp, nspin, nrad
! -----------------------------------------------------------------
  integer(kind=i4b) :: ne, na, nb, nl, sra
  integer(kind=i4b) :: istart, iend, ng
  integer(kind=i4b) :: i, is, ib, il, iw(0:lmax), ntot
  real(kind=r8b)    :: ere, eim

! head of outsusc.dat can be used to make new inpsusc.dat
  open(file='inpsusc.dat',unit=io,status='old')
! number of energy points, number of atoms, max number of basis functions
  read(io,*)  ! comment line
  read(io,*)  ! comment line
  read(io,*) ne, na, nb, nl, sra
!  if (ne /= nepts) stop 'nesusc /= nepts'   ! messy in fpimpu
  if (nl > lmax)   stop 'susc_input: nlmax   >  lmax'
  if (na > natyp)  stop 'susc_input: nasusc  > natyp'
! -----------------------------------------------------------------
! initialize global variables and allocate arrays needed to read input
  call init_param(ne,na,nb,nrad,nl,nspin,sra)
! now global variables are in use
! -----------------------------------------------------------------
! read atom labels in groups
! if a set of atoms is to be treated in the same way, only one read is needed
  read(io,*)  ! comment line
  read(io,*) ngroup
  read(io,*) igroup(1:ngroup)
  ntot = sum(igroup(1:ngroup))
  if (ntot /= nasusc) stop 'susc_input: atom groups and nasusc'
! loop over atom groups
  iend = 0
  do ng=1,ngroup
    istart = iend + 1
    iend   = iend + igroup(ng)
    read(io,*)  ! comment line
!   atomic labels
    read(io,*) iasusc(istart:iend)
    read(io,*)  ! comment line
    read(io,*) is, iw(0:nlmax)
!   which atoms are magnetic (is == 2) or not (is == 1)
    issusc(istart:iend) = is
    if (is > nsmax) stop 'susc_input: is > nsmax'
    do il=0,nlmax
!     info on how many basis functions for each l channel
      iwsusc(il,1:nsmax,istart:iend) = iw(il)
!     read the projection energies for each l channel
!     if different energies for each spin are wanted change here
      do ib=1,iw(il)
        read(io,*) i, ere, eim
        ewsusc(ib,il,1:is,istart:iend) = cmplx(ere,eim)
      end do
    end do
  end do
  close(io)
  if (any(iasusc(1:nasusc) > natyp)) stop 'susc_input: atom labels'
! -----------------------------------------------------------------
! allocate global arrays according to input
  call init_arrays
! now global arrays are available
! All done!
  end subroutine susc_input
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine outsusc_out()
! write initial info to output file

  implicit none

  integer(kind=i4b) :: istart, iend, ng, il, ib

  if (noparameters) stop 'open_outsusc: run init_param first!'
  open(file='outsusc.dat',unit=io,status='replace')
  write(io,'(" Data for dynamical susceptibility calculation")')
  write(io,'(" ne, na, nbmax, nlmax, sra, nrmax, nsmax=")')
  write(io,'(8i6)') nesusc, nasusc, nbmax, nlmax, isra, nrmax, nsmax
  write(io,'(" ngroup, igroup(1:ngroup)=")')
  write(io,'(i4,/,1000i4)') ngroup, igroup(1:ngroup)
  iend = 0
  do ng=1,ngroup
    istart = iend + 1
    iend   = iend + igroup(ng)
    write(io,'(" atom labels in group")')  
    write(io,'(1000i6)') iasusc(istart:iend)
    write(io,'(" nspin, iwsusc(0:lmax), ewsusc(1:nbmax,0:lmax)")')
    write(io,'(i2,20i6)') issusc(istart), iwsusc(:,:,istart)
    do il=0,nlmax
      do ib=1,iwsusc(il,1,istart)
        write(io,'(i4,2f10.3)') ib, ewsusc(ib,il,1,istart)
      end do
    end do
  end do
! All done!
  end subroutine outsusc_out
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine outsusc_in()
! read all info from outsusc.dat

  implicit none

  integer(kind=i4b) :: ne, na, nb, nl, sra, nrad, nspin, iw(0:100,2)
  integer(kind=i4b) :: istart, iend, ng, il, ib, ib1, is, ie, ia
  integer(kind=i4b) :: icount, im
  real(kind=r8b)    :: ram, er, ei, der, dei
  complex(kind=c8b) :: e, ek, de
  complex(kind=c8b), allocatable :: gmat(:,:,:,:)
  character*60 :: header

  write(*,'(/"Opening outsusc.dat"/)')
  open(file='outsusc.dat',unit=io,status='old')
  read(io,'(a)') header ! data
!  write(iodb,*) header
  read(io,'(a)') header ! vars
!  write(iodb,*) header
  read(io,*) ne, na, nb, nl, sra, nrad, nspin
!  write(iodb,'(10i6)') ne, na, nb, nl, sra, nrad, nspin
! Global variables initialized
  call init_param(ne,na,nb,nrad,nl,nspin,sra)
  write(*,*) header
  write(*,'(10i6)') nesusc, nasusc, nbmax, nlmax, isra, nrmax, nsmax
  read(io,'(a)') header ! ngroup, igroup=
!  write(iodb,*) header
  read(io,*) ngroup
  read(io,*) igroup(1:ngroup)
  write(*,'(" ngroup, igroup=",1000i4)') ngroup, igroup(1:ngroup)
  iend = 0
  do ng=1,ngroup
    istart = iend + 1
    iend   = iend + igroup(ng)
    read(io,'(a)') header ! atom labels in group
!    write(iodb,*) header
    read(io,*) iasusc(istart:iend)
    write(*,'("iasusc",1000i4)') iasusc(istart:iend)
    read(io,'(a)') header ! nspin, iwsusc, ewwusc
!    write(iodb,*) header
    read(io,*) is, iw(0:nlmax,1:nsmax)
    issusc(istart:iend) = is
    do is=1,issusc(istart)
      do il=0,nlmax
        iwsusc(il,is,istart:iend) = iw(il,is)
      end do
    end do
    write(*,'("is,iwl=",1000i4)') issusc(istart), iwsusc(:,:,istart)
    do il=0,nlmax
      do ib=1,iwsusc(il,1,istart)
        read(io,'(a)') header ! projection energies are not needed
!        write(iodb,*) header
      end do
    end do
  end do
! Allocate global arrays according to input
  call init_arrays
! Set up real spherical harmonics and Gaunt coefficients
  call ymy_gaunts
! Set up the angular momentum matrices
  call orbmoment(nlmax,lmmax)
! Read radial mesh, potentials, core densities and basis functions
  do ia=1,nasusc
    call read_rmesh(ia)
    do is=1,issusc(ia)
      do il=0,nlmax
        nb = iwsusc(il,is,ia)
        if (nb > 0) call read_wfns(ia,il,is)
      end do
    end do
  end do
! Set up overlaps for the basis
  call overlaps
! Put the SCF potential in the basis
  call build_vscfb
! Now read projection coefficients and structural GFs
! auxiliary storage for reads
  allocate(gmat(lmmax,lmmax,nasusc,nasusc))
! Big loop to read everything from outsusc.dat
  do is=1,nsmax
    read(io,'(a)') header ! spin separator
!    write(iodb,*) header
    write(*,'(/,"ispin=",i4)') is
    do ie=1,nesusc
      write(*,*) "ie =", ie
!     Read coefficients
      call in_coeffs(is,ie,e,ek,de)
!     Save energy mesh; de should be dummy
      esusc(ie) = e; eksusc(ie) = ek; desusc(ie) = de
      write(*,'("e,ek,de=",6es16.8)') e, ek, de
!     Save projection coefficients
      call save_coeffs(ie,is)
!     Read structural GF
      call in_gmat(gmat)
!     Put it somewhere
!     This is the collinear case
      write(iodb,'("Saving GF in RAM",2i4)') ie, is
      call store_gscoll(gmat,is,ie)
    end do
  end do
  write(*,'(/"Closing outsusc.dat"/)')
  close(io)
  deallocate(gmat)
! Read SCF mesh from file
  write(*,*) "Reading E-mesh for integration"
  open(file='emesh.scf',unit=io)
  read(io,*) nescf
  write(*,*) "nescf=", nescf
  do ie=1,nescf
    read(io,*) er, ei, der, dei
    escf(ie)  = cmplx(er,ei)
    ekscf(ie) = sqrt(escf(ie))
    descf(ie) = cmplx(der,dei)
!    descf(ie) = -4.d0*atan(1.d0)*cmplx(der,dei)
    write(*,'("e,ek,de=",6es16.8)') escf(ie), ekscf(ie), descf(ie)
  end do
  write(*,'("de sums to=",2f12.6,/)') sum(descf(1:nescf))
  close(io)
! All done
  end subroutine outsusc_in
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine outsusc_in2()
! read all info from outsusc.dat

  implicit none

  integer(kind=i4b) :: ne, na, nb, nl, sra, nrad, nspin, iw(0:100,2)
  integer(kind=i4b) :: istart, iend, ng, il, ib, ib1, is, ie, ia
  integer(kind=i4b) :: icount, im
  real(kind=r8b)    :: ram, er, ei, der, dei
  complex(kind=c8b) :: e, ek, de
  complex(kind=c8b), allocatable :: gmat(:,:,:,:)
  character*60 :: header

  write(*,'(/"Opening outsusc.dat"/)')
  open(file='outsusc.dat',unit=io,status='old')
  read(io,'(a)') header ! data
!  write(iodb,*) header
  read(io,'(a)') header ! vars
!  write(iodb,*) header
  read(io,*) ne, na, nb, nl, sra, nrad, nspin
!  write(iodb,'(10i6)') ne, na, nb, nl, sra, nrad, nspin
! Global variables initialized
  call init_param(ne,na,nb,nrad,nl,nspin,sra)
  write(*,*) header
  write(*,'(10i6)') nesusc, nasusc, nbmax, nlmax, isra, nrmax, nsmax
  read(io,'(a)') header ! ngroup, igroup=
!  write(iodb,*) header
  read(io,*) ngroup
  read(io,*) igroup(1:ngroup)
  write(*,'(" ngroup, igroup=",1000i4)') ngroup, igroup(1:ngroup)
  iend = 0
  do ng=1,ngroup
    istart = iend + 1
    iend   = iend + igroup(ng)
    read(io,'(a)') header ! atom labels in group
!    write(iodb,*) header
    read(io,*) iasusc(istart:iend)
    write(*,'("iasusc",1000i4)') iasusc(istart:iend)
    read(io,'(a)') header ! nspin, iwsusc, ewwusc
!    write(iodb,*) header
    read(io,*) is, iw(0:nlmax,1:nsmax)
    issusc(istart:iend) = is
    do is=1,issusc(istart)
      do il=0,nlmax
        iwsusc(il,is,istart:iend) = iw(il,is)
      end do
    end do
    write(*,'("is,iwl=",1000i4)') issusc(istart), iwsusc(:,:,istart)
    do il=0,nlmax
      do ib=1,iwsusc(il,1,istart)
        read(io,'(a)') header ! projection energies are not needed
!        write(iodb,*) header
      end do
    end do
  end do
! Allocate global arrays according to input
  call init_arrays
! Set up real spherical harmonics and Gaunt coefficients
  call ymy_gaunts
! Set up the angular momentum matrices
  call orbmoment(nlmax,lmmax)
! Read radial mesh, potentials, core densities and basis functions
  do ia=1,nasusc
    call read_rmesh(ia)
    do is=1,issusc(ia)
      do il=0,nlmax
        nb = iwsusc(il,is,ia)
        if (nb > 0) call read_wfns(ia,il,is)
      end do
    end do
  end do
! Set up overlaps for the basis
  call overlaps
! Put the SCF potential in the basis
  call build_vscfb
! Now read projection coefficients and structural GFs
! auxiliary storage for reads
  allocate(gmat(lmmax,lmmax,nasusc,nasusc))
! Big loop to read everything from outsusc.dat
  do ie=1,nesusc
    !write(*,*) "ie =", ie
    do is=1,nsmax
      read(io,'(a)') header ! spin separator
!      write(iodb,*) header
      !write(*,'(/,"ispin=",i4)') is
!     Read coefficients
      call in_coeffs(is,ie,e,ek,de)
!     Save energy mesh; de should be dummy
      esusc(ie) = e; eksusc(ie) = ek; desusc(ie) = de
      write(*,'("ie=",i5," is=",i2,": e,ek,de=",6es16.8)') ie, is, e, ek, de
!     Save projection coefficients
      call save_coeffs(ie,is)
!     Read structural GF
      call in_gmat(gmat)
!     Put it somewhere
!     This is the collinear case
      write(iodb,'("Saving GF in RAM",2i4)') ie, is
      call store_gscoll(gmat,is,ie)
    end do
  end do
  ensusc = esusc
  write(*,'(/"Closing outsusc.dat"/)')
  close(io)
  deallocate(gmat)
! Read SCF mesh from file
  write(*,*) "Reading E-mesh for integration"
  open(file='emesh.scf',unit=io)
  read(io,*) nescf
  write(*,*) "nescf=", nescf
  do ie=1,nescf
    read(io,*) er, ei, der, dei
    escf(ie)  = cmplx(er,ei)
    ekscf(ie) = sqrt(escf(ie))
    descf(ie) = cmplx(der,dei)
!    descf(ie) = -4.d0*atan(1.d0)*cmplx(der,dei)
    write(*,'("e,ek,de=",6es16.8)') escf(ie), ekscf(ie), descf(ie)
  end do
  write(*,'("de sums to=",2f12.6,/)') sum(descf(1:nescf))
  efermi = real(escf(nescf))
  close(io)
! All done
  end subroutine outsusc_in2
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine close_outsusc()
  close(io)
  end subroutine close_outsusc
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine init_param(ne,na,nb,nr,nl,ns,sra)
! These values should be specified somewhere in the main program
! or read in from an input file

  implicit none

  integer(kind=i4b), intent(in) :: ne, na, nb, nr, nl, ns, sra

! Global parameters
  nesusc = ne                            ! exact number of energy points to be computed
  nasusc = na                            ! exact number of atoms for susc
  nlmax  = nl                            ! maximum angular momentum
  lmmax  = (nl+1)**2                     ! length of angular momentum arrays
  nlmax2 = 2*nl
  lmmax2 = (2*nl+1)**2                   ! L1 and L2 in Gaunt numbers
  nlmax4 = 4*nl
  lmmax4 = (4*nl+1)**2                   ! L3 in Gaunt numbers
!  nlmax0 = 1                            ! read only by suscmain
  lmmax0 = (nlmax0 + 1)**2               ! comes from above
  nbmax  = nb                            ! maximum number of wfns for each l, s and atom
  isra   = sra                           ! SRA small components
  nrmax  = nr                            ! set this to nrad from inc.p or parameters.file
  nsmax  = ns                            ! spin components; set to 2 if ANYTHING is magnetic
  nlms   = lmmax*nsmax                   ! total size of scattering matrices with spin
  nalms  = nasusc*nlms                   ! dimensions for structural GF
  nlmsb  = nbmax*nlms                    ! maximum size of block of projected GF
!***************************************************************************************************
! Combined index for the susceptibility
! HUGE, THINK TWICE BEFORE DECLARING A MATRIX WITH THIS
!  nalmsb = nasusc*nlmsb*nlmsb  >>>>  redefined in init_arrays, to save memory
!***************************************************************************************************
! Allocate storage for input parameters and control information
  allocate(iasusc(nasusc))                      ! which atoms for susc
  allocate(iwsusc(0:nlmax,nsmax,nasusc))        ! which l-blocks to use
  allocate(issusc(nasusc))                      ! atom magnetic or not
  allocate(ewsusc(nbmax,0:nlmax,nsmax,nasusc))  ! projection energies
  allocate(nowfns(0:nlmax,nsmax,nasusc))        ! were all wfns stored?
  allocate(nobasis(0:nlmax,nsmax,nasusc))       ! was a basis constructed?
! Flag these things as initialized
  noparameters = .false.
! All done!
  end subroutine init_param
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine init_arrays()
! Allocate the storage according to the susc input

  implicit none

  integer(kind=i4b) :: is, il, im, ib, i, lm, ia
  integer(kind=i4b) :: ilmsb, jlmsb, ilms, jlms
  real(kind=r8b)    :: ram

  if (noparameters) stop 'init_arrays: run init_param first!'
! How much memory is being allocated here:
  ram = 8.d0*nrmax*nbmax*(nlmax+1)*nsmax*nasusc
  ram = ram + 8.d0*(nlmax+8)*nrmax*nasusc
  ram = ram + 16.d0*nbmax*(nbmax+2)*(nlmax+1)*nsmax*nasusc
  ram = ram/(1024.d0**2)
  write(*,'(/"init_arrays:    RAM=",f16.3," MB"/)') ram
! The basis is constructed per l-channel
  allocate(phiref(nrmax,nbmax,0:nlmax,nsmax,nasusc))  ! reference wfns
  allocate(overlap(nlmsb,nlmsb,nasusc))               ! overlaps
  phiref = 0.d0; overlap = 0.d0
  allocate(nrpts(nasusc))                             ! number of points in each radial mesh
  allocate(nrpts0(nasusc))                            ! start of each radial mesh
  allocate(nrpts1(nasusc))                            ! end of radial mesh
  nrpts = 0; nrpts0 = 0; nrpts1 = 0
  allocate(rmesh(nrmax,nasusc))                       ! radial mesh for evaluation
  allocate(rsmesh(nrmax,0:nlmax,nasusc))              ! powers of r
  allocate(drmesh(nrmax,nasusc))                      ! weights for radial integration
  allocate(drproj(nrmax,nasusc))                      ! weights for radial projection
  rmesh = 0.d0; drmesh = 0.d0; drproj = 0.d0
  allocate(zat(nasusc))                               ! atomic numbers
  allocate(magdir(3,nasusc))                          ! magnetization direction
  allocate(vr(nrmax,nasusc))                          ! charge potentials
  allocate(br(nrmax,nasusc))                          ! magnetic potentials
  allocate(nrc(nrmax,nasusc))                         ! core charge densities
  allocate(mrc(nrmax,nasusc))                         ! core magnetization densities
  allocate(nrv(nrmax,0:6,nasusc))                     ! valence densities
  zat = 0.d0; vr = 0.d0; br = 0.d0; nrc = 0.d0; mrc = 0.d0
  allocate(normesh(nasusc))                           ! was the radial mesh stored?
  allocate(ensusc(nesusc))                            ! energies for susc calculation
  allocate(pzc(nbmax,0:nlmax,nsmax,nasusc))           ! coefficients of pz in radial basis
  allocate(pqc(nbmax,nbmax,0:nlmax,nsmax,nasusc))     ! coefficients of pzqz in radial basis
  pzc = 0.d0; pqc = 0.d0
! ---------------------------------------------------------------------------------------------
  if (isra == 1) then
    allocate(fzc(nbmax,0:nlmax,nsmax,nasusc))         ! coefficients of fz in radial basis
    allocate(psc(nbmax,nbmax,0:nlmax,nsmax,nasusc))   ! coefficients of pzsz in radial basis
    allocate(fqc(nbmax,nbmax,0:nlmax,nsmax,nasusc))   ! coefficients of fzqz in radial basis
    allocate(fsc(nbmax,nbmax,0:nlmax,nsmax,nasusc))   ! coefficients of fzsz in radial basis
    fzc = 0.d0; psc = 0.d0; fqc = 0.d0; fsc = 0.d0
  end if
! ---------------------------------------------------------------------------------------------
  allocate(noregcoeffs(nsmax,nesusc))                 ! were reg coeffs computed?
  allocate(noirrcoeffs(nsmax,nesusc))                 ! were irr coeffs computed?
  allocate(esusc(nesusc))                             ! energy mesh for susc
  allocate(eksusc(nesusc))                            ! sqrt of e
  allocate(desusc(nesusc))                            ! energy integration weights
  allocate(wsusc(nesusc))                             ! interpolation weights
  esusc = 0.d0; eksusc = 0.d0; desusc = 0.d0
  open(file='emesh.scf',unit=20141118)
    read(20141118,*) nescf
  close(20141118)
  allocate(escf(nescf))                               ! same things
  allocate(ekscf(nescf))                              ! but to store the energy contour
  allocate(descf(nescf))                              ! used in the SCF calculation
 !allocate(escf(nesusc))                              ! same things
 !allocate(ekscf(nesusc))                             ! but to store the energy contour
 !allocate(descf(nesusc))                             ! used in the SCF calculation
  escf = 0.d0; ekscf = 0.d0; descf = 0.d0
! Flag these things as not initialized
  normesh      = .true.
  nowfns       = .true.
  nobasis      = .true.
  noregcoeffs  = .true.
  noirrcoeffs  = .true.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                    Indices for packing and unpacking arrays
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! --> lm (big to accommodate Gaunt numbers)
  allocate(lm2i(-lmmax2:lmmax2,0:lmmax2))
  allocate(i2lm(2,lmmax2))
  i = 0
  do il=0,nlmax2
    do im=-il,il
      i = i + 1
      lm2i(im,il) = i
      i2lm(:,i) = (/im,il/)
    end do
  end do
! --> lms
  allocate(lms2i(lmmax,nsmax))
  allocate(i2lms(2,nlms))
  i = 0
  do lm=1,lmmax
    do is=1,nsmax
      i  = i + 1
      lms2i(lm,is) = i
      i2lms(:,i) = (/lm,is/)
    end do
  end do
  write(*,'("nlms   should be ",i8)') i
! --> alms
  allocate(alms2i(nlms,nasusc))
  allocate(i2alms(2,nalms))
  i = 0
  do ia=1,nasusc
    do ilms=1,nlms
      i = i + 1
      alms2i(ilms,ia) = i
      i2alms(:,i) = (/ilms,ia/)
    end do
  end do
  write(*,'("nalms  should be ",i8)') i
! --> lmsb, new version: block sizes atom dependent
  allocate(lmsb2i(nbmax,lmmax,nsmax,nasusc))
  allocate(i2lmsb(3,nlmsb,nasusc))
  allocate(nlmsba(nasusc))
  do ia=1,nasusc
    i = 0
    do lm=1,lmmax
      do is=1,nsmax
        do ib=1,iwsusc(i2lm(2,lm),is,ia)
          i = i + 1
          lmsb2i(ib,lm,is,ia) = i
          i2lmsb(:,i,ia) = (/ib,lm,is/)
        end do
      end do
    end do
    nlmsba(ia) = i
    write(*,'("nlmsb for ia=",i4," is ",i4)') ia, i
  end do
! --> almsb
  nalmsb = sum(nlmsba*nlmsba)
  allocate(almsb2i(nlmsb,nlmsb,nasusc))
  allocate(i2almsb(3,nalmsb))
  i = 0
  do ia=1,nasusc
    do jlmsb=1,nlmsba(ia)
      do ilmsb=1,nlmsba(ia)
        i = i + 1
        almsb2i(ilmsb,jlmsb,ia) = i
        i2almsb(:,i) = (/ilmsb,jlmsb,ia/)
      end do
    end do
  end do
  write(*,'("nalmsb should be ",i8)') i
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            Storage for structural GF and coefficients in full projection basis
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! how much ram to be used in these arrays
  ram = 16.d0*2*nlms*nlmsb*nasusc*nesusc              ! size of reg coeffs arrays (lhs and rhs)
  ram = ram + 16.d0*nlmsb*nlmsb*nasusc*nesusc         ! size of single-site GF array
  ram = ram + 16.d0*nalms*nalms*nesusc                ! size of backscattering GF
  ram = ram + 16.d0*8*nlmsb*nlmsb*nasusc              ! size of new potential and xc kernel arrays
! extra storage for small components
  if (isra == 1) then                                 
    ram = ram + 16.d0*2*nlmsb*nlms*nasusc*nesusc      ! fz left and right
    ram = ram + 16.d0*3*nlmsb*nlmsb*nasusc*nesusc     ! ps, fq, fs
  end if
  ram = ram/(1024.d0**2)
! time to cry
  write(*,'(/"init_arrays: GF parts  RAM=",f16.3," MB"/)') ram
  ram = 16.d0*(1+numd+dend)*nlmsb*nlmsb*nasusc*nasusc
  ram = ram/(1024.d0**2)
! time to cry again
  write(*,'(/"init_arrays: GF fit    RAM=",f16.3," MB"/)') ram
! ---------------------------------------------------------------------------------------------
  allocate(gffit(1+numd+dend,nlmsb,nlmsb,nasusc,nasusc)) ! storage for fitted GF
  allocate(gstruct(nalms,nalms,nesusc))               ! storage for structural GF
  gstruct = 0.d0
  allocate(vlmsb(nlmsb,nlmsb,nasusc))                 ! potential in the projection basis
  allocate(pzl(nlmsb,nlms,nasusc,nesusc))             ! storage for lhs pz in lms with basis labels
  allocate(pzr(nlmsb,nlms,nasusc,nesusc))             ! storage for rhs pz in lms with basis labels
  allocate(gfpq(nlmsb,nlmsb,nasusc,nesusc))           ! storage for onsite GF pq in lms with basis labels
  pzl = 0.d0; pzr = 0.d0;  gfpq = 0.d0
! ---------------------------------------------------------------------------------------------
  if (isra == 1) then
    allocate(fzl(nlmsb,nlms,nasusc,nesusc))           ! storage for lhs fz
    allocate(fzr(nlmsb,nlms,nasusc,nesusc))           ! storage for rhs fz
    allocate(gfps(nlmsb,nlmsb,nasusc,nesusc))         ! storage for ps
    allocate(gffq(nlmsb,nlmsb,nasusc,nesusc))         ! storage for fq
    allocate(gffs(nlmsb,nlmsb,nasusc,nesusc))         ! storage for fs
    fzl = 0.d0; fzr = 0.d0; gfps = 0.d0; gffq = 0.d0; gffs = 0.d0
  end if
! ---------------------------------------------------------------------------------------------
  allocate(lorb(lmmax,lmmax,3))                       ! storage for angular momentum matrices
! ---------------------------------------------------------------------------------------------
! KS susc
  if (isusc > 0) then
    ram = 16.d0*nalmsb**2
    ram = ram/(1024.d0**3)
    write(*,'(/"init_arrays: KS susc   RAM=",f16.3," GB"/)') ram
    if (ram > 4.d0) stop 'KS susc RAM too large'
    allocate(kssusc(nalmsb,nalmsb))                   ! storage for kronecker KS
    kssusc = 0.d0
  end if
! xc kernel
  if (ikxc > 0) then
    ram = 16.d0*nalmsb**2
    ram = ram/(1024.d0**3)
    write(*,'(/"init_arrays: xc kernel RAM=",f16.3," GB"/)') ram
    if (ram > 4.d0) stop 'xc kernel RAM too large'
    allocate(kxcsusc(nalmsb,nalmsb))                  ! storage for kronecker KS
    kxcsusc = 0.d0
  end if
! xc kernel
  if (isusc > 0 .and. ikxc > 0) then
    ram = 16.d0*nalmsb**2
    ram = ram/(1024.d0**3)
    write(*,'(/"init_arrays: susc den  RAM=",f16.3," GB"/)') ram
    if (ram > 4.d0) stop 'susc den RAM too large'
    allocate(denominator(nalmsb,nalmsb))              ! storage for kronecker KS
    kxcsusc = 0.d0
  end if
! All done!
  end subroutine init_arrays
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine ymy_gaunts()
! Opens file with the LL mesh, generates spherical harmonics and Gaunts
! Makes use of global arrays ull, wll, cthetall, phill, ylm and gaunt

  implicit none

  integer(kind=i4b) :: ill, ngaunt(0:nlmax2)
  integer(kind=i4b) :: ngauntl(0:nlmax,0:nlmax,0:nlmax2)
  integer(kind=i4b) :: ilm, jlm, klm, im, il, jm, jl, km, kl, i2(2)
  real(kind=r8b)    :: norm
  real(kind=r8b), allocatable :: work(:)

! read Lebedev-Laikov mesh
  open(file='lebedev_ascii.gga',unit=iofile)
  read(iofile,*) nll, nlmaxll  ! I added the lmax in the file: nll=434 -> lmax=16
  write(*,'("Lebedev-Laikov mesh: nll, lmax=",2i4)') nll, nlmaxll
  lmmaxll = (nlmaxll+1)**2
  if (nlmaxll < 2*nlmax) stop 'LL mesh: nlmaxll < 4*nlmax!'
  allocate(ull(3,nll),wll(nll),cthetall(nll),phill(nll))
  allocate(ylm(nll,lmmax2),work(lmmax2))
  allocate(gaunt(lmmax,lmmax,lmmax2))
  do ill=1,nll
    read(iofile,*) ull(:,ill), wll(ill), cthetall(ill), phill(ill)
    write(iodb,'("ill=",i6,6es10.2)') ill, ull(:,ill), wll(ill), cthetall(ill), phill(ill)
  end do
  write(*,'("LL weights integrate to",f8.3)') sum(wll)
!  wll = wll*16.d0*atan(1.d0)
!  write(*,'("LL weights multiplied by 4pi")')
  close(iofile)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! compute real spherical harmonics
  do ill=1,nll
    call ymy(ull(:,ill),nlmax2,lmmax2,work)
!    call ymy2(cthetall(ill),phill(ill),nlmax2,lmmax2,work)
    ylm(ill,:) = work
  end do
! check various things
  do jlm=1,lmmax2
    norm = sum(wll*ylm(:,jlm))
    write(iodb,'("<yilm>",2i4,es24.12)') i2lm(:,jlm), norm
    do ilm=1,lmmax2
      norm = sum(wll*ylm(:,jlm)*ylm(:,ilm))
      if (abs(norm) > ylmtol) then
        write(iodb,'("<yilm*yjlm>",4i4,es24.12)') i2lm(:,ilm), i2lm(:,ilm), norm
      end if
    end do
  end do
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! now Gaunt numbers
  gaunt = 0.d0; ngaunt = 0; ngauntl = 0
  do klm=1,lmmax2
    i2 = i2lm(:,klm)
    km = i2(1); kl = i2(2)
    do jlm=1,lmmax 
      i2 = i2lm(:,jlm)
      jm = i2(1); jl = i2(2)
      do ilm=1,lmmax 
        i2 = i2lm(:,ilm)
        im = i2(1); il = i2(2)
!     -------------------------------------------------------------
!                          selection rules
!                 Homeier et al, Theochem 368, p31 (1996)
!     -------------------------------------------------------------
!     m selection rules
        if (abs(km) == abs(im+jm) .or. abs(km) == abs(im-jm)) then
!     parity selection rule
        if (mod(il+jl+kl,2) == 0) then
!     triangle inequality for l
        if (kl <= il + jl .and. kl >= max(abs(il-jl),min(abs(im+jm),abs(im-jm)))) then
          norm = sum(wll*ylm(:,ilm)*ylm(:,jlm)*ylm(:,klm))  ! this is where the work is
!         numerical selection rule
          if (abs(norm) > ylmtol) then
!            write(iodb,'("Gaunt: ",6i6,es8.1)') im, il, jm, jl, km, kl, norm
            gaunt(ilm,jlm,klm) = norm
            ngaunt(kl) = ngaunt(kl) + 1
            if (ngauntl(il,jl,kl) == 0) ngauntl(il,jl,kl) = 1
          end if
        end if
        end if
        end if
!     -------------------------------------------------------------
      end do
    end do
  end do
  write(*,'("Number of non-zero Gaunt numbers",i4,/)') sum(ngaunt)
  write(iodb,'("Number of Gaunt number per l")')
  do kl=0,nlmax2
    write(iodb,'(10i6)') kl, ngaunt(kl)
    do jl=0,nlmax
      write(iodb,'(10i2)') ngauntl(:,jl,kl)
    end do
  end do
! All done!
  end subroutine ymy_gaunts
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine ymy2(ctheta,phi,nl,nlm,ylm)
! this subroutine calculates real spherical harmonics with the
! normalization : <y|y> =1
!
! pilfered and abused in 2012 by MdSD

  implicit none

! spherical angles
  real(kind=r8b),    intent(in)  :: ctheta, phi
! dimensions
  integer(kind=i4b), intent(in)  :: nl, nlm
! real spherical harmonics
  real(kind=r8b),    intent(out) :: ylm(nlm)
! -----------------------------------------------------------------
! parameters
  real(kind=r8b), parameter :: invr4pi = 0.25d0/sqrt(atan(1.d0))
  real(kind=r8b), parameter :: rtwo = sqrt(2.d0)
  real(kind=r8b), parameter :: snull = 1.d-10
! -----------------------------------------------------------------
  real(kind=r8b)    :: xy, xyz, cth, sth, cph, sph, sthm
  integer(kind=i4b) :: i, l, m
  real(kind=r8b)    :: p(0:nl,0:nl), c(0:nl), s(0:nl)
  real(kind=r8b)    :: dfac, root, root1, root2


! I can't be bothered doing the extra if's
  if (nl < 2) stop 'check ymy2'


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            calculate sin and cos of theta and phi
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  cth = ctheta
  sth = sqrt(1.d0 - ctheta*ctheta)
! pathological thetas
  if (abs(ctheta) < snull) then
    cth =  0.d0; sth = 1.d0
  else if (abs(ctheta - 1.d0) < snull) then
    cth =  1.d0; sth = 0.d0
  else if (abs(ctheta + 1.d0) < snull) then
    cth = -1.d0; sth = 0.d0
  end if
  cph = cos(phi)
  sph = sin(phi)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     generate normalized associated legendre functions for m >= 0
!                see NR section 6.7 (only 2007)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! loop over m values
  dfac = 1; root = 1.d0; sthm = 1.d0; p(0,0) = 1.d0!invr4pi
  do m=0,nl
    if (m > 0) then
      sthm = sthm*sth
      dfac = dfac*sqrt((2*m-1)**2/(2.d0*m*(2*m-1)))
      root = sqrt((2.d0*m+1.d0))!*invr4pi
      p(m,m) = dfac*root*sthm  ! there should be a minus sign here for Condon-Shortley
    end if
    if (m < nl) p(m+1,m) = sqrt(2*m+3.d0)*cth*p(m,m)
!   upward recursion in l
    do l=m+2,nl
      root1 = sqrt((4*l*l-1.d0)/(l*l-m*m))
      root2 = sqrt(((l-1)*(l-1)-m*m)/(4*(l-1)*(l-1)-1.d0))
      p(l,m) = root1*(cth*p(l-1,m) - root2*p(l-2,m))
    end do
  end do
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            recursions for sin(m phi) and cos(m phi)
!                     see NR section 5.4
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  s(0) = 0.d0; s(1) = sph
  c(0) = 1.d0; c(1) = cph
  do m=2,nl
    s(m) = sin(m*phi)
    c(m) = cos(m*phi)
  end do
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         multiply legendre functions with cosines and sines
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  i = 0
  do l=0,nl
    i = i + l + 1
    ylm(i) = p(l,0)  ! m = 0
    do m=1,l
      ylm(i+m) = rtwo*p(l,m)*c(m)  ! m > 0 cosine type
      ylm(i-m) = rtwo*p(l,m)*s(m)  ! m < 0 sine type
    end do
    i = i + l
  end do
! All done!
  end subroutine ymy2
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine ymy(uvec,nl,nlm,ylm)
! this subroutine calculates real spherical harmonics with the
! normalization : <y|y> =1
!
! pilfered and abused in 2012 by MdSD

  implicit none

! components of unit vector from which ylm is to be computed
  real(kind=r8b),    intent(in)  :: uvec(3)
! dimensions
  integer(kind=i4b), intent(in)  :: nl, nlm
! real spherical harmonics
  real(kind=r8b),    intent(out) :: ylm(nlm)
! -----------------------------------------------------------------
! parameters
  real(kind=r8b), parameter :: invr4pi = 0.25d0/sqrt(atan(1.d0))
  real(kind=r8b), parameter :: rtwo = sqrt(2.d0)
  real(kind=r8b), parameter :: snull = 1.d-20
! -----------------------------------------------------------------
  real(kind=r8b)    :: xy, xyz, cth, sth, cph, sph, sthm
  integer(kind=i4b) :: i, l, m
  real(kind=r8b)    :: p(0:nl,0:nl), c(0:nl), s(0:nl)
  real(kind=r8b)    :: dfac, root, root1, root2


! I can't be bothered doing the extra if's
  if (nl < 2) stop 'check ymy'


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            calculate sin and cos of theta and phi
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  xy = uvec(1)**2 + uvec(2)**2
  xyz = xy + uvec(3)**2
  if (xy > snull*xyz) then  ! theta not 0 or pi
    xy  = sqrt(xy)
    xyz = sqrt(xyz)
! cos and sin theta
    cth = uvec(3)/xyz
    sth = xy/xyz
! cos and sin phi
    cph = uvec(1)/xy
    sph = uvec(2)/xy
  else  ! theta 0 or pi
    sth = 0.d0; cth = 1.d0
    if (uvec(3) < 0.d0) cth = -1.d0
    cph = 1.d0; sph = 0.d0
  end if
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     generate normalized associated legendre functions for m >= 0
!                see NR section 6.7 (only 2007)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! loop over m values
  dfac = 1.d0; root = 1.d0; sthm = 1.d0; p(0,0) = 1.d0!invr4pi
  do m=0,nl
    if (m > 0) then
      sthm = sthm*sth
      dfac = dfac*sqrt((2*m-1)**2/(2.d0*m*(2*m-1)))
      root = sqrt((2.d0*m+1.d0))!*invr4pi
      p(m,m) = dfac*root*sthm  ! there should be a minus sign here for Condon-Shortley
    end if
    if (m < nl) p(m+1,m) = sqrt(2*m+3.d0)*cth*p(m,m)
! upward recursion in l
    do l=m+2,nl
      root1 = sqrt((4*l**2 - 1.d0)/(l**2 -m**2))
      root2 = sqrt(((l-1)**2 - m**2)/(4*(l-1)**2 - 1.d0))
      p(l,m) = root1*(cth*p(l-1,m) - root2*p(l-2,m))
    end do
  end do
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            recursions for sin(m phi) and cos(m phi)
!                     see NR section 5.4
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  s(0) = 0.d0; s(1) = sph
  c(0) = 1.d0; c(1) = cph
  do m=2,nl
    s(m) = 2.d0*cph*s(m-1) - s(m-2)
    c(m) = 2.d0*cph*c(m-1) - c(m-2)
  end do
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         multiply legendre functions with cosines and sines
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  i = 0
  do l=0,nl
    i = i + l + 1
    ylm(i) = p(l,0)  ! m = 0
    do m=1,l
      ylm(i+m) = rtwo*p(l,m)*c(m)  ! m > 0 cosine type
      ylm(i-m) = rtwo*p(l,m)*s(m)  ! m < 0 sine type
    end do
    i = i + l
  end do
! All done!
  end subroutine ymy
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine save_rmesh(ia,nr0,nr1,rr,rs,wr,z,vup,vdn,ncup,ncdn)
! Save radial mesh and potentials

  implicit none

! Start and end of radial mesh; susceptibility atom
  integer(kind=i4b), intent(in) :: nr0, nr1, ia
! Radial mesh, powers of r and radial integration weights
  real(kind=r8b),    intent(in) :: rr(nrmax), wr(nrmax), rs(nrmax,0:nlmax)
! Atomic number, spin up and spin down potentials, core densities
  real(kind=r8b),    intent(in) :: z, vup(nrmax), vdn(nrmax), ncup(nrmax), ncdn(nrmax)
!     -----------------------------------------------------------------
  integer(kind=i4b) :: nr, ir
  real(kind=r8b)    :: magdummy(3)
  complex(kind=c8b) :: work(nrmax), norm
  character*12      :: filename

  nr = nr1 - nr0 + 1
  nrpts(ia)  = nr
  nrpts0(ia) = nr0
  nrpts1(ia) = nr1
  magdummy = (/0.d0,0.d0,1.d0/)
! Check how the radial mesh is defined, what is the first point
  rmesh(1:nr,ia) = rr(nr0:nr1)
  rsmesh(1:nr,0:nlmax,ia) = rs(nr0:nr1,0:nlmax)
! Original radial integration weights
  drmesh(1:nr,ia) = wr(nr0:nr1)
! Modified weights for projection
  drproj(1:nr,ia) = drmesh(1:nr,ia)
! Save atomic number and potentials
  zat(ia) = z
  vr(1:nr,ia) = 0.5d0*(vup(nr0:nr1) + vdn(nr0:nr1))
!  vr(1:nr,ia) = vr(1:nr,ia) - 2.d0*z/rr(nr0:nr1)
  br(1:nr,ia) = 0.5d0*(vup(nr0:nr1) - vdn(nr0:nr1))
! Core charge and magnetization densities
  nrc(1:nr,ia) = ncup(nr0:nr1) + ncdn(nr0:nr1)
  mrc(1:nr,ia) = ncup(nr0:nr1) - ncdn(nr0:nr1)
! Write potential file for susc
  write(filename,'("susc",i4.4,".pot")') ia
  open(file=filename,unit=iofile,status='replace')
  write(iofile,'("# Susc potential: z, nr, lmax; then magdummy; then rmesh, rs, dr, vr, br, nrc, mrc")')
  write(iofile,'("# ",f6.1,2i8)') z, nr, nlmax  ! pot can be plotted
  write(iofile,'("# ",3f12.8)') magdummy  ! magnetization direction
  do ir=1,nr
    write(iofile,'(100es16.8)') rmesh(ir,ia), rsmesh(ir,0:nlmax,ia), drmesh(ir,ia), vr(ir,ia), br(ir,ia), nrc(ir,ia), mrc(ir,ia)
  end do
  close(iofile)
! Output mesh and potentials
  normesh(ia) = .false.
! All done!
  end subroutine save_rmesh
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine read_rmesh(ia)
! Read radial mesh and potentials

  implicit none

! Susceptibility atom
  integer(kind=i4b), intent(in) :: ia
! -----------------------------------------------------------------
  integer(kind=i4b) :: nr, ir, nlmax1
  real(kind=r8b)    :: z
  complex(kind=c8b) :: work(nrmax), norm
  character*1       :: dummy
  character*12      :: filename
  logical           :: exists

  write(filename,'("susc",i4.4,".scf")') ia
! Where is my mind?
  inquire(file=filename,exist=exists)
  if (.not.exists) then
    write(*,*) "read_rmesh: file ",filename," not found!"
    stop
  end if
  write(iodb,*) "Reading ", filename
! Read potential file for susc
  open(file=filename,unit=iofile,status='old')
  read(iofile,*) ! header
  read(iofile,*) dummy, z, nr, nlmax1  ! dummy is the symbol #
! Not sure how useful this is, as rs is not used anywhere so far
  if (nlmax1 > nlmax) then
    write(*,*) "read_rmesh: for ia=",ia," nlmax1 > nlmax"
    stop
  end if
! Direction of magnetization
  read(iofile,*) dummy, magdir(:,ia)
! Number of radial points and atomic number
  nrpts(ia) = nr
  zat(ia) = z
  write(*,'("nr,z=",i8,f6.1)') nrpts(ia), zat(ia)
! r, r^s, dr for integration, scalar and magnetic potentials
  do ir=1,nr
    read(iofile,*) rmesh(ir,ia), rsmesh(ir,0:nlmax,ia), drmesh(ir,ia), vr(ir,ia), br(ir,ia), nrc(ir,ia), mrc(ir,ia)
!    write(iodb,'(100es16.8)') rmesh(ir,ia), rsmesh(ir,0:nlmax,ia), drmesh(ir,ia), vr(ir,ia), br(ir,ia), nrc(ir,ia), mrc(ir,ia)
  end do
  close(iofile)
! Mesh and potentials were read
  normesh(ia) = .false.
! All done!
  end subroutine read_rmesh
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine ref_wfn(ia,il,is,ib,phi)
! Save wavefunction for projection
! phi is computed at real energy without normalization at WS radius
! Should be inside a loop over the energies used for projection
! If using irregular solutions put them at the end of the list

  implicit none

! --> susc atom
  integer(kind=i4b), intent(in) :: ia
! --> angular momentum
  integer(kind=i4b), intent(in) :: il
! --> spin channel
  integer(kind=i4b), intent(in) :: is
! --> basis function
  integer(kind=i4b), intent(in) :: ib
! --> regular scattering solution
  complex(kind=c8b), intent(in) :: phi(nrmax)
! -----------------------------------------------------------------
  integer(kind=i4b) :: nr0, nr1, nr
  complex(kind=c8b) :: proj(nrmax), work(nrmax), norm

! Good morning
  if (noparameters) stop 'ref_wfn: run init_param!'
  if (normesh(ia))  stop 'ref_wfn: save rmesh first!'
  if (ib > nbmax)   stop 'ref_wfn: ib > nbmax!'
! Store trial basis functions
  nr0 = nrpts0(ia); nr1 = nrpts1(ia); nr = nr1 - nr0 + 1
  phiref(1:nr,ib,il,is,ia) = real(phi(nr0:nr1))
! Make overall sign positive
  work = phiref(1:nr,ib,il,is,ia)
  norm = radint(nr,work(1:nr),drproj(1:nr,ia))
  if (real(norm) < 0.d0) phiref(1:nr,ib,il,is,ia) = -phiref(1:nr,ib,il,is,ia)
! Normalize to 1
  work = phiref(1:nr,ib,il,is,ia)*phiref(1:nr,ib,il,is,ia)
  norm = radint(nr,work(1:nr),drproj(1:nr,ia))
  phiref(1:nr,ib,il,is,ia) = phiref(1:nr,ib,il,is,ia)/sqrt(real(norm))
! Check if the basis is completely stored
  if (ib == iwsusc(il,is,ia)) nowfns(il,is,ia) = .false.
! All done!
  end subroutine ref_wfn
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine new_basis1(nb,ia,il,is,tol)
! The list of trial wavefunctions is reduced to a minimal basis
! New basis functions with norm below tol are discarded

  implicit none

! --> reduced number of basis functions
  integer(kind=i4b), intent(out) :: nb
! --> susc atom
  integer(kind=i4b), intent(in)  :: ia
! --> angular momentum
  integer(kind=i4b), intent(in)  :: il
! --> spin channel
  integer(kind=i4b), intent(in)  :: is
! --> tolerance to chop off basis functions
  real(kind=r8b),    intent(in)  :: tol
! -----------------------------------------------------------------
  integer(kind=i4b) :: nr, nb1, ib
  real(kind=r8b)    :: basis(nrmax,nbmax)

! Good morning
  if (noparameters)     stop 'new_basis1: run init_param!'
  if (nowfns(il,is,ia)) stop 'new_basis1: save wavefunctions first!'
  nr = nrpts(ia)
  nb = iwsusc(il,is,ia)
  write(iodb,'("new_basis1:",5i4)') ia, il, is, nr, nb
! Store trial basis functions and do GS
  nb1 = nb
  basis(1:nr,1:nb1) = phiref(1:nr,1:nb1,il,is,ia)
  call find_basis2(nb1,basis,nr,drproj(:,ia),tol)
  phiref(1:nr,1:nb,il,is,ia) = basis(1:nr,1:nb)
! The size of the basis is updated
  nb = nb1
  iwsusc(il,is,ia) = nb
! Basis constructed
  nobasis(il,is,ia) = .false.
! All done!
  end subroutine new_basis1
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine new_basis2(nb,ia,il,tol)
! The list of trial wavefunctions is reduced to a minimal basis
! New basis functions with norm below tol are discarded

  implicit none

! --> reduced number of basis functions
  integer(kind=i4b), intent(out) :: nb
! --> susc atom
  integer(kind=i4b), intent(in)  :: ia
! --> angular momentum
  integer(kind=i4b), intent(in)  :: il
! --> tolerance to chop off basis functions
  real(kind=r8b),    intent(in)  :: tol
! -----------------------------------------------------------------
  integer(kind=i4b) :: nr, nb1, ib, ib0, ib1, is, ns
  real(kind=r8b)    :: basis(nrmax,nbmax)

! Good morning
  if (noparameters) stop 'new_basis2: run init_param!'
  ns = issusc(ia)
  nr = nrpts(ia)
! Store trial basis functions and do GS
  ib1 = 0
  do is=1,ns
    nb1 = iwsusc(il,is,ia)
    if (nowfns(il,is,ia) .and. nb1 > 1) stop 'new_basis2: save wfns'
    ib0 = ib1 + 1
    ib1 = ib1 + nb1
    if (ib1 > nbmax) stop 'new_basis2: increase nbmax!'
    basis(1:nr,ib0:ib1) = phiref(1:nr,1:nb1,il,is,ia)
  end do
  nb1 = ib1
  nb  = nb1
  write(iodb,'("new_basis2:",5i4)') ia, il, ns, nr, nb
  call find_basis2(nb1,basis,nr,drproj(:,ia),tol)
! A bit redundant, but probably easier to use
  do is=1,ns
    phiref(1:nr,1:nb,il,is,ia) = basis(1:nr,1:nb)
  end do
! The size of the basis is updated
  nb = nb1
  iwsusc(il,1:ns,ia) = nb
! Basis constructed
  nobasis(il,1:ns,ia) = .false.
! All done!
  end subroutine new_basis2
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine new_basis3(nb,ia,tol)
! The list of trial wavefunctions is reduced to a minimal basis
! New basis functions with norm below tol are discarded

  implicit none

! --> reduced number of basis functions
  integer(kind=i4b), intent(out) :: nb
! --> susc atom
  integer(kind=i4b), intent(in)  :: ia
! --> tolerance to chop off basis functions
  real(kind=r8b),    intent(in)  :: tol
! -----------------------------------------------------------------
  integer(kind=i4b) :: nr, nb1, nb2, ib, ib0, ib1, il, is, ns
  real(kind=r8b)    :: basis(nrmax,nbmax)

! Good morning
  if (noparameters) stop 'new_basis3: run init_param!'
  ns = issusc(ia)
  nr = nrpts(ia)
! Store trial basis functions and do GS
  ib1 = 0
  do il=0,nlmax
    do is=1,ns
      nb1 = iwsusc(il,is,ia)
      if (nowfns(il,is,ia) .and. nb1 > 1) stop 'new_basis3: save wfns'
      ib0 = ib1 + 1
      ib1 = ib1 + nb1
      if (ib1 > nbmax) stop 'new_basis3: increase nbmax!'
      basis(1:nr,ib0:ib1) = phiref(1:nr,1:nb1,il,is,ia)
      write(iodb,'(3i4)') nb1, ib0, ib1
    end do
  end do
  nb1 = ib1
  nb  = nb1
  write(iodb,'("new_basis3:",5i4)') ia, nlmax, ns, nr, nb
  call find_basis2(nb1,basis,nr,drproj(:,ia),tol)
! A bit redundant, but probably easier to use
  do il=0,nlmax
    do is=1,ns
      phiref(1:nr,1:nb,il,is,ia) = basis(1:nr,1:nb)
    end do
  end do
! The size of the basis is updated
  nb = nb1
  iwsusc(0:nlmax,1:ns,ia) = nb
! Basis constructed
  nobasis(0:nlmax,1:ns,ia) = .false.
! All done!
  end subroutine new_basis3
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine find_basis(nb,phi,nr,wr,tol)
! Gram-Schmidt with basis reduction

  implicit none

! --> input number of basis function / output reduced number of bfs
  integer(kind=i4b), intent(inout) :: nb
! --> input trial basis functions / output independent bfs only
  real(kind=r8b),    intent(inout) :: phi(nrmax,nbmax)
! --> actual number of points in radial mesh
  integer(kind=i4b), intent(in)    :: nr
! --> radial mesh weights
  real(kind=r8b),    intent(in)    :: wr(nrmax)
! --> tolerance for discarding a basis function
  real(kind=r8b),    intent(in)    :: tol
! -----------------------------------------------------------------
  integer(kind=i4b) :: ib, jb, nb1, ikeep(nbmax), ibasis(nbmax)
  real(kind=r8b)    :: rnorm(nbmax)
  complex(kind=c8b) :: work(nrmax), norm

  if (nb > nbmax) stop 'find_basis: nb > nbmax'
  if (nr > nrmax) stop 'find_basis: nr > nrmax'
  write(iodb,'("find_basis:",2i4,es12.1)') nr, nb, tol
! begin by normalizing the initial basis function
  do ib=1,nb
    work = phi(:,ib)*phi(:,ib)
    norm = radint(nr,work(1:nr),wr(1:nr))
    phi(:,ib) = phi(:,ib)/sqrt(real(norm))
    write(iodb,'("find_basis:",i4,2es12.1)') ib, norm
  end do
! whether to keep a basis function
  ikeep(1:nb) = 1
! number of non-degenerate basis functions
  nb1 = 0
  do ib=1,nb
!   Gram-Schmidt
    do jb=1,ib-1
      if (ikeep(jb) == 1) then
!     subtract contributions from previous basis function
        work = phi(:,jb)*phi(:,ib)
        norm = radint(nr,work(1:nr),wr(1:nr))
        phi(:,ib) = phi(:,ib) - real(norm)*phi(:,jb)
        write(iodb,'("GS   ",i4,2es12.3)') jb, norm
      end if
    end do
!   norm of current basis function
    work = phi(:,ib)*phi(:,ib)
    norm = radint(nr,work(1:nr),wr(1:nr))
    rnorm(ib) = real(norm)
    write(iodb,'("Norm ",i4,es12.3)') ib, sqrt(rnorm(ib))
    if (rnorm(ib) < max(tol**2,gstol**2)) then
!   degenerate, discard it
      ikeep(ib) = 0
    else
!   one more
      nb1 = nb1 + 1
!   which is number...
      ibasis(nb1) = ib
!   normalize to the WS cell
      phi(:,ib) = phi(:,ib)/sqrt(rnorm(ib))
    end if
  end do
  write(iodb,'("ibasis=",100i4)') ibasis(1:nb1)
! how many non-degenerate basis functions
  nb = nb1
! fetch them to the front of the array
  do ib=1,nb
    phi(:,ib) = phi(:,ibasis(ib))
  end do
! the rest is history
  phi(:,nb+1:nbmax) = 0.d0
! All done!
  end subroutine find_basis
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine find_basis2(nb,phi,nr,wr,tol)
! Diagonalization of the overlap matrix

  implicit none

! --> input number of basis function / output reduced number of bfs
  integer(kind=i4b), intent(inout) :: nb
! --> input trial basis functions / output independent bfs only
  real(kind=r8b),    intent(inout) :: phi(nrmax,nbmax)
! --> actual number of points in radial mesh
  integer(kind=i4b), intent(in)    :: nr
! --> radial mesh weights
  real(kind=r8b),    intent(in)    :: wr(nrmax)
! --> tolerance for discarding a basis function
  real(kind=r8b),    intent(in)    :: tol
! -----------------------------------------------------------------
  integer(kind=i4b) :: ib, jb, nb1, ikeep(nbmax), ibasis(nbmax)
  real(kind=r8b)    :: rnorm(nbmax), overlaps(nb,nb)
  complex(kind=c8b) :: work(nrmax), norm
  real(kind=r8b)    :: evals(nb), rwork(2*nb*(1 + 3*nb) + 1)
  integer(kind=i4b) :: iwork(5*nb+3), info, lrwork, liwork
  real(kind=r8b)    :: invsqrt(nb,nb), newphi(nrmax,nbmax)

  if (nb > nbmax) stop 'find_basis2: nb > nbmax'
  if (nr > nrmax) stop 'find_basis2: nr > nrmax'
  write(iodb,'("find_basis2:",2i4,es12.1)') nr, nb, tol
! begin by computing the overlaps
  do jb=1,nb
    do ib=1,nb
      work = phi(:,ib)*phi(:,jb)
      norm = radint(nr,work(1:nr),wr(1:nr))
      overlaps(ib,jb) = real(norm)
    end do
    write(iodb,'("find_basis2:",i4,100es12.1)') jb, overlaps(1:nb,jb)
  end do
! diagonalize; eigenvalues in ascending order
  lrwork  = 2*nb*(1 + 3*nb) + 1
  liwork = 5*nb + 3
  call dsyevd('V','U',nb,overlaps,nb,evals,rwork,lrwork,iwork,liwork,info)
  write(iodb,'("find_basis: eigen")')
  do jb=1,nb
    if (overlaps(1,jb) < 0.d0) overlaps(1:nb,jb) = -overlaps(1:nb,jb)
    write(iodb,'(es12.3," | ",100f12.6)') evals(jb), overlaps(1:nb,jb)
  end do
! new basis
  nb1 = 0; newphi = 0.d0
  do jb=1,nb
    if (abs(evals(nb-jb+1)) > max(tol,gstol)*evals(nb)) then
      nb1 = nb1 + 1
      do ib=1,nb
        newphi(:,jb) = newphi(:,jb) + phi(:,ib)*overlaps(ib,nb-jb+1)
      end do
      work = newphi(:,jb)*newphi(:,jb)
      norm = radint(nr,work(1:nr),wr(1:nr))
!     this should return the eigenvalues
      write(iodb,'("find_basis2: newphi",i4,2es12.3)') jb, norm
      newphi(:,jb) = newphi(:,jb)/sqrt(norm)
    end if
  end do
! new overlaps
  do jb=1,nb
    do ib=1,nb
      work = newphi(:,ib)*newphi(:,jb)
      norm = radint(nr,work(1:nr),wr(1:nr))
      overlaps(ib,jb) = real(norm)
    end do
    write(iodb,'("find_basis2:",i4,100es12.1)') jb, overlaps(1:nb,jb)
  end do
! new non-degenerate basis functions
  nb  = nb1
  phi = newphi
! All done!
  end subroutine find_basis2
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine reg_coeffs(natomd,lmaxd,nrmaxd,pz,fz,is,ie,ek)
! Energy dependent projection coefficients of regular scattering wfns
! Expects solutions of the ASA SRA problem, no m-dependence of pz or fz
! BEWARE OF THE SPECIAL TREATMENT OF THE RADIAL MESH, IR=1 NOT USED

  implicit none

! --> dimensions of the wavefunction arrays
  integer(kind=i4b), intent(in) :: natomd, lmaxd, nrmaxd
! --> regular scattering solutions
  complex(kind=c8b), intent(in) :: pz(nrmaxd,0:lmaxd,natomd), fz(nrmaxd,0:lmaxd,natomd)
! --> spin channel
  integer(kind=i4b), intent(in) :: is
! --> current energy point
  integer(kind=i4b), intent(in) :: ie
! --> defined value of the square-root of the energy
  complex(kind=c8b), intent(in) :: ek
! -----------------------------------------------------------------
  integer(kind=i4b) :: ia, ih, nr, nr0, nr1, il, nb, ib
  real(kind=r8b)    :: dr(nrmax)
  complex(kind=c8b) :: work1(nrmax), work2(nrmax), norm1, norm2

  do ia=1,nasusc
    ih  = iasusc(ia)
    nr0 = nrpts0(ia); nr1 = nrpts1(ia); nr = nr1 - nr0 + 1
    dr  = drproj(:,ia)
    pzc(:,:,is,ia) = 0.d0; if (isra == 1)  fzc(:,:,is,ia) = 0.d0
    do il=0,nlmax
      nb = iwsusc(il,is,ia)
      if (nobasis(il,is,ia) .and. nb > 0) then
        write(*,'("reg_coeffs: no basis",3i4)') il, is, ia
        stop
      end if
!     projection on each basis function
!     basis functions assumed normalized to 1 with weight dr
      do ib=1,nb
        work1(1:nr) = pz(nr0:nr1,il,ih)*phiref(1:nr,ib,il,is,ia)
        norm1 = radint(nr,work1(1:nr),dr(1:nr))
        pzc(ib,il,is,ia) = norm1
!     -------------------------------------------------------------
        if (isra == 1) then
          work2(1:nr) = fz(nr0:nr1,il,ih)*phiref(1:nr,ib,il,is,ia)
          norm2 = radint(nr,work2(1:nr),dr(1:nr))
          fzc(ib,il,is,ia) = norm2
        end if
!     -------------------------------------------------------------
      end do
!     check how much of the function is left
      if (nb > 0) then
        work1(1:nr) = pz(nr0:nr1,il,ih)
        do ib=1,nb
          norm1 = pzc(ib,il,is,ia)
          work1(1:nr) = work1(1:nr) - norm1*phiref(1:nr,ib,il,is,ia)
        end do
        work1 = work1*conjg(work1)
        norm1 = radint(nr,work1(1:nr),dr(1:nr))
        work1(1:nr) = conjg(pz(nr0:nr1,il,ih))*pz(nr0:nr1,il,ih)
        norm1 = norm1/radint(nr,work1(1:nr),dr(1:nr))
!     -------------------------------------------------------------
        norm2 = 0.d0
      if (isra == 1) then
        work2(1:nr) = fz(nr0:nr1,il,ih)
        do ib=1,nb
          norm2 = fzc(ib,il,is,ia)
          work2(1:nr) = work2(1:nr) - norm2*phiref(1:nr,ib,il,is,ia)
        end do
        work2 = work2*conjg(work2)
        norm2 = radint(nr,work2(1:nr),dr(1:nr))
        work2(1:nr) = conjg(fz(nr0:nr1,il,ih))*fz(nr0:nr1,il,ih)
        norm2 = norm2/radint(nr,work2(1:nr),dr(1:nr))
      end if
!     -------------------------------------------------------------
        write(iodb,'("reg_coeffs",4i4,2es16.3)') ie, ia, is, il, sqrt(abs(norm1)), sqrt(abs(norm2))
      end if
    end do
  end do
! coefficients computed
  noregcoeffs(is,ie) = .false.
! All done
  end subroutine reg_coeffs
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine irr_coeffs(natomd,lmaxd,nrmaxd,pz,qz,fz,sz,is,ie,ek)
! Energy dependent projection coefficients of the on-site part of the GF
! Expects solutions of the ASA SRA problem, no m-dependence of pz or fz
! BEWARE OF THE SPECIAL TREATMENT OF THE RADIAL MESH, IR=1 NOT USED

  implicit none

! --> dimensions of the wavefunction arrays
  integer(kind=i4b), intent(in) :: natomd, lmaxd, nrmaxd
! --> regular scattering solutions
  complex(kind=c8b), intent(in) :: pz(nrmaxd,0:lmaxd,natomd), fz(nrmaxd,0:lmaxd,natomd)
! --> irregular scattering solutions
  complex(kind=c8b), intent(in) :: qz(nrmaxd,0:lmaxd,natomd), sz(nrmaxd,0:lmaxd,natomd)
! --> spin channel
  integer(kind=i4b), intent(in) :: is
! --> current energy point
  integer(kind=i4b), intent(in) :: ie
! --> defined value of the square-root of the energy
  complex(kind=c8b), intent(in) :: ek
! -----------------------------------------------------------------
  logical, parameter :: separable = .false.
! -----------------------------------------------------------------
  integer(kind=i4b) :: ia, ih, nr, nr0, nr1, il, nb
  integer(kind=i4b) :: ib, jb, ir
  real(kind=r8b)    :: dr(nrmax)
  complex(kind=c8b) :: work1(nrmax), work2(nrmax), norm, norm2

  do ia=1,nasusc
    ih = iasusc(ia)
    nr0 = nrpts0(ia); nr1 = nrpts1(ia); nr = nr1 - nr0 + 1
    dr = drproj(:,ia)
    pqc(:,:,:,:,ia) = 0.d0
    if (isra == 1) then
      psc(:,:,:,:,ia) = 0.d0
      fqc(:,:,:,:,ia) = 0.d0
      fsc(:,:,:,:,ia) = 0.d0
    end if
    do il=0,nlmax
      nb = iwsusc(il,is,ia)
      if (nobasis(il,is,ia) .and. nb > 0) then
        write(*,'("irr_coeffs: no basis",3i4)') il, is, ia
        stop
      end if
!   ====================================================================
      if (separable) then
!   expansion of the irregular solutions
      do jb=1,nb
        work1(1:nr) = qz(nr0:nr1,il,ih)*phiref(1:nr,jb,il,is,ia)
        norm = radint(nr,work1(1:nr),dr(1:nr))
        if (isra == 1) then
          work2(1:nr) = sz(nr0:nr1,il,ih)*phiref(1:nr,jb,il,is,ia)
          norm2 = radint(nr,work2(1:nr),dr(1:nr))
        end if
        do ib=1,nb
          pqc(ib,jb,il,is,ia) = norm*ek*pzc(ib,il,is,ia)
          if (isra == 1) then
            psc(ib,jb,il,is,ia) = norm2*ek*pzc(ib,il,is,ia)
            fqc(ib,jb,il,is,ia) = norm*ek*fzc(ib,il,is,ia)
            fsc(ib,jb,il,is,ia) = norm2*ek*fzc(ib,il,is,ia)
          end if
        end do
      end do
!   symmetrize
!      pqc(1:nb,1:nb,il,is,ia) = 0.5d0*(pqc(1:nb,1:nb,il,is,ia) - transpose(pqc(1:nb,1:nb,il,is,ia)))
!   ====================================================================
      else
!   ====================================================================
!   double expansion of the on-site product
!   basis functions assumed normalized to 1 with weight dr
      do jb=1,nb
      do ib=1,nb
!     integral over r' with switching between P and Q
!     ir <-> r; work2 stores data for r'
        do ir=1,nr
!       r' <= r --> Q(r)P(r')
          work2(1:ir) = qz(nr0+ir-1,il,ih)*pz(nr0:nr0+ir-1,il,ih)
!       r' > r --> P(r)Q(r')
          work2(ir+1:nr) = pz(nr0+ir-1,il,ih)*qz(nr0+ir:nr1,il,ih)
          work2(1:nr) = work2(1:nr)*phiref(1:nr,ib,il,is,ia)
!       integrate over r'
          norm = radint(nr,work2(1:nr),dr(1:nr))
          work1(ir) = norm
        end do
!     multiply by basis function and integrate over r
        work1(1:nr) = work1(1:nr)*phiref(1:nr,jb,il,is,ia)
        norm = radint(nr,work1(1:nr),dr(1:nr))
!     include the sqrt of the energy
        pqc(ib,jb,il,is,ia) = norm*ek
!     -------------------------------------------------------------
      if (isra == 1) then
!     integral over r' with switching between P and S
        do ir=1,nr
!       r' <= r --> S(r)P(r')
          work2(1:ir) = sz(nr0+ir-1,il,ih)*pz(nr0:nr0+ir-1,il,ih)
!       r' > r --> P(r)S(r')
          work2(ir+1:nr) = pz(nr0+ir-1,il,ih)*sz(nr0+ir:nr,il,ih)
          work2(1:nr) = work2(1:nr)*phiref(1:nr,ib,il,is,ia)
          norm = radint(nr,work2(1:nr),dr(1:nr))
          work1(ir) = norm
        end do
!     multiply by basis function and integrate over r'
        work1(1:nr) = work1(1:nr)*phiref(1:nr,jb,il,is,ia)
        norm = radint(nr,work1(1:nr),dr(1:nr))
!     include the sqrt of the energy
        psc(ib,jb,il,is,ia) = norm*ek
!     -------------------------------------------------------------
!     integral over r' with switching between F and Q
        do ir=1,nr
!       r' <= r --> Q(r)F(r')
          work2(1:ir) = qz(nr0+ir-1,il,ih)*fz(nr0:nr0+ir-1,il,ih)
!       r' > r --> F(r)Q(r')
          work2(ir+1:nr) = fz(nr0+ir-1,il,ih)*qz(nr0+ir:nr,il,ih)
          work2(1:nr) = work2(1:nr)*phiref(1:nr,ib,il,is,ia)
          norm = radint(nr,work2(1:nr),dr(1:nr))
          work1(ir) = norm
        end do
!     multiply by basis function and integrate over r'
        work1(1:nr) = work1(1:nr)*phiref(1:nr,jb,il,is,ia)
        norm = radint(nr,work1(1:nr),dr(1:nr))
!     include the sqrt of the energy
        fqc(ib,jb,il,is,ia) = norm*ek
!     -------------------------------------------------------------
!     integral over r' with switching between F and S
        do ir=1,nr
!       r' <= r --> S(r)F(r')
          work2(1:ir) = sz(nr0+ir-1,il,ih)*fz(nr0:nr0+ir-1,il,ih)
!       r' > r --> F(r)S(r')
          work2(ir+1:nr) = fz(nr0+ir-1,il,ih)*sz(nr0+ir:nr,il,ih)
          work2(1:nr) = work2(1:nr)*phiref(1:nr,ib,il,is,ia)
          norm = radint(nr,work2(1:nr),dr(1:nr))
          work1(ir) = norm
        end do
!     multiply by basis function and integrate over r'
        work1(1:nr) = work1(1:nr)*phiref(1:nr,jb,il,is,ia)
        norm = radint(nr,work1(1:nr),dr(1:nr))
!     include the sqrt of the energy
        fsc(ib,jb,il,is,ia) = norm*ek
      end if
!     -------------------------------------------------------------
      end do
      end do
!     ==================================================================
      end if
!     ==================================================================
    end do
  end do
! coefficients computed
  noirrcoeffs = .false.
! All done
  end subroutine irr_coeffs
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine out_coeffs(is,ie,e,ek,de)
! Output projection coefficients to outsusc.dat file

  implicit none

! --> spin channel
  integer(kind=i4b), intent(in) :: is
! --> energy point label
  integer(kind=i4b), intent(in) :: ie
! --> energy point value, its square-root, integration weight
  complex(kind=c8b), intent(in) :: e, ek, de
! -----------------------------------------------------------------
  integer(kind=i4b) :: ia, ih, nb, ib, il

  write(io,'(" ie, e, ek, de=",i4)') ie
  write(io,'(6es16.8)') e, ek, de
  write(io,'(" Projection coefficients:")')
  do ia=1,nasusc
    ih = iasusc(ia)
    write(io,'(" ia=",i4,"  i1=",i4)') ia, ih
    do il=0,nlmax
      nb = iwsusc(il,is,ia)
      if (nb > 0) then
        write(io,'(4i4)') il, is, nb
        write(io,*) " pz coeff"
        do ib=1,nb
          write(io,'(100es16.8)') pzc(ib,il,is,ia)
        end do
!     ---------------------------------------------------------------
        if (isra == 1) then
          write(io,*) " fz coeff"
          do ib=1,nb
            write(io,'(100es16.8)') fzc(ib,il,is,ia)
          end do
        end if
!     ---------------------------------------------------------------
      end if
    end do
  end do
  write(io,'(" Onsite coefficients:")')
  do ia=1,nasusc
    ih = iasusc(ia)
    write(io,'(" ia=",i4,"  i1=",i4)') ia, ih
    do il=0,nlmax
      nb = iwsusc(il,is,ia)
      if (nb > 0) then
        write(io,'(4i4)') il, is, nb
        write(io,*) " pq coeff"
        do ib=1,nb
          write(io,'(100es16.8)') pqc(1:nb,ib,il,is,ia)
        end do
!     ---------------------------------------------------------------
        if (isra == 1) then
          write(io,*) " ps coeff"
          do ib=1,nb
            write(io,'(100es16.8)') psc(1:nb,ib,il,is,ia)
          end do
          write(io,*) " fq coeff"
          do ib=1,nb
            write(io,'(100es16.8)') fqc(1:nb,ib,il,is,ia)
          end do
          write(io,*) " fs coeff"
          do ib=1,nb
            write(io,'(100es16.8)') fsc(1:nb,ib,il,is,ia)
          end do
        end if
!     ---------------------------------------------------------------
      end if
    end do
  end do
! All done
  end subroutine out_coeffs
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine in_coeffs(is,ie,e,ek,de)
! Get projection coefficients from outsusc.dat file

  implicit none

! --> spin channel
  integer(kind=i4b), intent(in)  :: is
! --> energy point label
  integer(kind=i4b), intent(in)  :: ie
! --> energy point value, its square-root, integration weight
  complex(kind=c8b), intent(out) :: e, ek, de
! -----------------------------------------------------------------
  integer(kind=i4b) :: ia, ih, nb, nb1, ib, jb, il, il1, is1
  real(kind=r8b)    :: x(1000)
  character*60 :: header

  read(io,'(a)') header ! ie, e, ek, de
!  write(iodb,*) header
  read(io,*) x(1:6)
!  write(iodb,'(6es16.8)') x(1:6)
  e  = cmplx(x(1),x(2))
  ek = cmplx(x(3),x(4))
  de = cmplx(x(5),x(6))
  read(io,'(a)') header ! proj coeffs
!  write(iodb,*) header
  do ia=1,nasusc
    read(io,'(a)') header ! ia, i1
!    write(iodb,*) header
    do il=0,nlmax
!   this should cross-check with read_wfns
      nb = iwsusc(il,is,ia)
      if (nb > 0) then
        read(io,*) il1, is1, nb
!        write(iodb,'(3i4)') il1, is1, nb
        read(io,'(a)') header ! pz coeffs
!        write(iodb,*) header
        do ib=1,nb
          read(io,*) x(1:2)
!          write(iodb,'(100es16.8)') x(1:2)
          pzc(ib,il,is,ia) = cmplx(x(1),x(2))
        end do
!     ---------------------------------------------------------------
        if (isra == 1) then
          read(io,'(a)') header ! fz coeffs
!          write(iodb,*) header
          do ib=1,nb
            read(io,*) x(1:2)
!            write(iodb,'(100es16.8)') x(1:2)
            fzc(ib,il,is,ia) = cmplx(x(1),x(2))
          end do
        end if
!     ---------------------------------------------------------------
      end if
    end do
  end do
  read(io,'(a)') header ! onsite coefficients
!  write(iodb,*) header
  do ia=1,nasusc
    read(io,'(a)') header ! ia, i1
!    write(iodb,*) header
    do il=0,nlmax
!   this was read in before
      nb = iwsusc(il,is,ia)
!      write(*,*) "nb=", nb
      if (nb > 0) then
        read(io,*) il1, is1, nb1
!        write(iodb,'(3i4)') il1, is1, nb
        read(io,'(a)') header ! pq coeff
!        write(iodb,*) header
        do ib=1,nb
          read(io,*) x(1:2*nb)
!          write(iodb,'(100es16.8)') x(1:2)
          do jb=1,nb
            pqc(jb,ib,il,is,ia) = cmplx(x(2*jb-1),x(2*jb))
          end do
        end do
!     -------------------------------------------------------------
        if (isra == 1) then
          read(io,'(a)') header ! ps coeff
!          write(iodb,*) header
          do ib=1,nb
            read(io,*) x(1:2*nb)
!            write(iodb,'(100es16.8)') x(1:2)
            do jb=1,nb
              psc(jb,ib,il,is,ia) = cmplx(x(2*jb-1),x(2*jb))
            end do
          end do
          read(io,'(a)') header ! fq coeff
!          write(iodb,*) header
          do ib=1,nb
            read(io,*) x(1:2*nb)
!            write(iodb,'(100es16.8)') x(1:2)
            do jb=1,nb
              fqc(jb,ib,il,is,ia) = cmplx(x(2*jb-1),x(2*jb))
            end do
          end do
          read(io,'(a)') ! fs coeff
!          write(iodb,*) header
          do ib=1,nb
            read(io,*) x(1:2*nb)
!            write(iodb,'(100es16.8)') x(1:2)
            do jb=1,nb
              fsc(jb,ib,il,is,ia) = cmplx(x(2*jb-1),x(2*jb))
            end do
          end do
        end if
!     -------------------------------------------------------------
      end if
    end do
  end do
! coefficients read in
  noregcoeffs = .false.
  noirrcoeffs = .false.
! All done
  end subroutine in_coeffs
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine save_coeffs(ie,is)
! Put projection coefficients in memory
! The onsite term is added to the big GF array
! For the collinear SRA case

  implicit none

! Current energy point
  integer(kind=i4b), intent(in) :: ie
! Current spin channel
  integer(kind=i4b), intent(in) :: is
! -----------------------------------------------------------------
  integer(kind=i4b) :: i, i1, j, il, lm, im, ia, ib, jb, nb

  do ia=1,nasusc
    do il=0,nlmax
      nb = iwsusc(il,is,ia)
!   Check if the angular momentum block is to keep
      if (nb > 0) then
        do im=-il,il
          lm = lm2i(im,il)
          do ib=1,nb
!       Saving wfn coefficients
            i = lmsb2i(ib,lm,is,ia)
            i1 = lms2i(lm,is)
!            write(*,'(10i4)') ia, is, lm, ib, i
            pzl(i,i1,ia,ie) = pzc(ib,il,is,ia)
            pzr(i,i1,ia,ie) = pzc(ib,il,is,ia)
            if (isra == 1) then
              fzl(i,i1,ia,ie) = fzc(ib,il,is,ia)
              fzr(i,i1,ia,ie) = fzc(ib,il,is,ia)
            end if
!           -------------------------------------------------------
!         Diagonal blocks coming from the RH term in the GF
            do jb=1,nb
              j = lmsb2i(jb,lm,is,ia)
              gfpq(i,j,ia,ie) = pqc(ib,jb,il,is,ia)
!          TEST TEST TEST TEST TEST TEST TEST TEST TEST
!              if (ib /= jb) gfpq(i,j,ia,ie) = 0.d0
!          TEST TEST TEST TEST TEST TEST TEST TEST TEST
              if (isra == 1) then
                gffq(i,j,ia,ie) = fqc(ib,jb,il,is,ia)
                gfps(i,j,ia,ie) = psc(ib,jb,il,is,ia)
                gffs(i,j,ia,ie) = fsc(ib,jb,il,is,ia)
              end if
            end do
!           -------------------------------------------------------
          end do
        end do
      end if
    end do
  end do
! All done
  end subroutine save_coeffs
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine out_gmat(nhost,lmmaxp,nsec,gmat)
! Output structural GF to outsusc.dat

  implicit none

! --> number of atoms for host
  integer(kind=i4b), intent(in) :: nhost
! --> size of the angular momentum matrix
  integer(kind=i4b), intent(in) :: lmmaxp
! --> size of the GF matrix
  integer(kind=i4b), intent(in) :: nsec
! --> the structural GF matrix
  complex(kind=c8b), intent(in) :: gmat(nsec,nsec)
! -----------------------------------------------------------------
  complex(kind=c8b) :: gtrim(lmmax,lmmax)
  real(kind=r8b)    :: gmax
  integer(kind=i4b) :: ia, ja, i1, i2, lm

  write(io,'(" GF: ia,ja, iasusc(ia),iasusc(ja), i1,i2")')
  do ja=1,nasusc
    i2 = (iasusc(ja)-nhost-1)*lmmaxp + 1
    do ia=1,nasusc
      i1 = (iasusc(ia)-nhost-1)*lmmaxp + 1
      write(io,'(6i8)') ia,ja, iasusc(ia),iasusc(ja), i1,i2
      gtrim = gmat(i1:i1+lmmax,i2:i2+lmmax)
      gmax = maxval(abs(gtrim))
      where (abs(gtrim) < gfilter*gmax) gtrim = 0.d0
      do lm=1,lmmax
        write(io,'(1000es16.8)') gtrim(1:lmmax,lm)
      end do
    end do
  end do
! All done
  end subroutine out_gmat
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine in_gmat(gmat)
! Get structural GF from outsusc.dat
! Reads the whole thing for a given spin channel and energy point

  implicit none

! --> the structural GF matrix for a given spin channel
  complex(kind=c8b), intent(out) :: gmat(lmmax,lmmax,nasusc,nasusc)
! -----------------------------------------------------------------
  integer(kind=i4b) :: ia, ja, lmi, lmj, i6(6)
  real(kind=r8b)    :: x(1000)
  character*60      :: header

!  write(iodb,*) "lmmax=", lmmax
  read(io,'(a)') header
!  write(iodb,*) header
  do ja=1,nasusc
    do ia=1,nasusc
      read(io,*) i6
!      write(iodb,'(6i4)') i6
      do lmj=1,lmmax
        read(io,*) x(1:2*lmmax)
        do lmi=1,lmmax
          gmat(lmi,lmj,ia,ja) = cmplx(x(2*lmi-1),x(2*lmi))
!          write(iodb,'(2i4,6es12.4)') lmi, lmj, gmat(lmi,lmj,ia,ja)
        end do
      end do
    end do
  end do
! All done
  end subroutine in_gmat
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine store_gscoll(gmat,is,ie)
! Puts the blocks of the structural GF in RAM
! Collinear spins here

  implicit none

  complex(kind=c8b), intent(in) :: gmat(lmmax,lmmax,nasusc,nasusc)
  integer(kind=i4b), intent(in) :: is, ie
! -----------------------------------------------------------------
  integer(kind=i4b) :: ja, ia, ilm, jlm, i, j

  do ja=1,nasusc
  do jlm=1,lmmax
    j = lms2i(jlm,is)
    j = alms2i(j,ja)
    do ia=1,nasusc
    do ilm=1,lmmax
      i = lms2i(ilm,is)
      i = alms2i(i,ia)
      gstruct(i,j,ie) = gmat(ilm,jlm,ia,ja)
!      write(iodb,'(4i4,6es12.4)') i2alms(:,i), i2alms(:,j), gmat(ilm,jlm,ia,ja)
    end do
    end do
  end do
  end do
! Make it structurally symmetric: Gij(L1,s1,L2,s2) = Gji(L2,s2,L1,s1)
  if (is == nsmax) then
  do j=1,nalms
    do i=j+1,nalms
      gstruct(i,j,ie) = 0.5d0*(gstruct(i,j,ie) + gstruct(j,i,ie))
      gstruct(j,i,ie) = gstruct(i,j,ie)
    end do
  end do
  end if
! All done
  end subroutine store_gscoll
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine save_wfns(ia,il,is)
! Output wavefunction for projection to file

  implicit none

! --> susc atom
  integer(kind=i4b), intent(in) :: ia
! --> angular momentum
  integer(kind=i4b), intent(in) :: il
! --> spin channel
  integer(kind=i4b), intent(in) :: is
! -----------------------------------------------------------------
  integer(kind=i4b) :: i, nb
  character*18      :: filename

  write(filename,'("ia",i4.4,"is",i2.2,"il",i2.2,".wfn")') ia,is,il
  open(file=filename,unit=iofile)
  nb = iwsusc(il,is,ia)
  write(iofile,'("# Susc basis set: ia, iasusc, il, is, nb, nr; then rmesh, dr, dr, phiref")')
  write(iofile,'("# ",6i8)') ia, iasusc(ia), il, is, nb, nrpts(ia)
  do i=1,nrpts(ia)
    write(iofile,'(1000es16.8)') rmesh(i,ia), drmesh(i,ia), drproj(i,ia), phiref(i,1:nb,il,is,ia)
  end do
  close(iofile)
! All done!
  end subroutine save_wfns
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine read_wfns(ia,il,is)
! Read basis functions from file

  implicit none

! --> susc atom
  integer(kind=i4b), intent(in) :: ia
! --> angular momentum
  integer(kind=i4b), intent(in) :: il
! --> spin channel
  integer(kind=i4b), intent(in) :: is
! -----------------------------------------------------------------
  integer(kind=i4b) :: i, ia1, ia2, il1, is1, nb
  character*1       :: dummy
  character*18      :: filename
  logical           :: exists

  write(filename,'("ia",i4.4,"is",i2.2,"il",i2.2,".wfn")') ia,is,il
! Where is my mind?
  inquire(file=filename,exist=exists)
  if (.not.exists) then
    write(*,*) "in_wfns: file ",filename," not found!"
    stop
  end if
! Was the potential file read first?
  if (normesh(ia)) stop 'in_wfns: read pot first!'
  write(iodb,*) "Reading ", filename
! Now read the basis set
  open(file=filename,unit=iofile,status='old')
  read(iofile,*) ! header
  read(iofile,*) dummy, ia1, ia2, il1, is1, nb
! The number of basis functions is updated here
  iwsusc(il,is,ia) = nb
  do i=1,nrpts(ia)
    read(iofile,*) rmesh(i,ia), drmesh(i,ia), drproj(i,ia), phiref(i,1:nb,il,is,ia)
  end do
  close(iofile)
! wfn was read in
  nowfns(il,is,ia) = .false.
! All done!
  end subroutine read_wfns
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine overlaps()
! Compute and store overlaps between basis functions

  implicit none

  integer(kind=i4b) :: ia, i3(3), nr
  integer(kind=i4b) :: i, j, is, js, ilm, jlm, il, jl, im, jm, ib, jb
  complex(kind=c8b) :: work(nrmax)

  overlap = 0.d0
  do ia=1,nasusc
    nr = nrpts(ia)
    write(iodb,*) "Overlaps for ia=", ia
    write(iodb,'("  ib  im  il  is  jb  jm  jl  js  overlap")')
    do j=1,nlmsba(ia)
      i3 = i2lmsb(:,j,ia)
      jb = i3(1); jlm = i3(2); js = i3(3)
      do i=1,nlmsba(ia)
        i3 = i2lmsb(:,i,ia)
        ib = i3(1); ilm = i3(2); is = i3(3)
        work = phiref(:,jb,i2lm(2,jlm),js,ia)*phiref(:,ib,i2lm(2,ilm),is,ia)
        overlap(i,j,ia) = radint(nr,work(1:nr),drmesh(1:nr,ia))
        write(iodb,'(8i4,es16.6)') ib, i2lm(:,ilm), is, jb, i2lm(:,jlm), js, overlap(i,j,ia)
      end do
    end do
  end do
! All done!
  end subroutine overlaps
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine overlaps_susc(ia,gfsum)
! Analysis of the product basis for the susceptibility

  implicit none

  integer(kind=i4b), intent(in) :: ia
  complex(kind=c8b), intent(in) :: gfsum(nlmsb,nlmsb)
! ---------------------------------------------------------------------
  real(kind=r8b), parameter :: tol = 1.d-3
  integer(kind=i4b) :: i3(3), nb, nr, i, j, lm, ib, ievals
  integer(kind=i4b) :: q1, b1, lm1, s1
  integer(kind=i4b) :: q2, b2, lm2, s2
  integer(kind=i4b) :: q3, b3, lm3, s3
  integer(kind=i4b) :: q4, b4, lm4, s4
  integer(kind=i4b) :: info, lrwork, liwork
  integer(kind=i4b), allocatable :: iwork(:), lmsb2i2(:,:)
  real(kind=r8b),    allocatable :: big_overlaps(:,:), rwork(:), evals(:)
  complex(kind=c8b) :: work(nrmax), norm, rho2, rho2proj, rho, rhoc
  real(kind=r8b)    :: dr(nrmax), doublegaunt, doublegaunt0


  nb = nlmsba(ia)**2
  allocate(big_overlaps(nb,nb),lmsb2i2(2,nb))
  big_overlaps = 0.d0; rho2 = 0.d0; rho2proj = 0.d0
  nr = nrpts(ia); dr = drmesh(:,ia)
  j = 0
  do q4=1,nlmsba(ia)
    i3 = i2lmsb(:,q4,ia)
    b4 = i3(1); lm4 = i3(2); s4 = i3(3)
    do q3=1,nlmsba(ia)
      i3 = i2lmsb(:,q3,ia)
      b3 = i3(1); lm3 = i3(2); s3 = i3(3)
      j = j + 1
      lmsb2i2(:,j) = (/q3,q4/)
      i = 0
      do q2=1,nlmsba(ia)
        i3 = i2lmsb(:,q2,ia)
        b2 = i3(1); lm2 = i3(2); s2 = i3(3)
        do q1=1,nlmsba(ia)
          i3 = i2lmsb(:,q1,ia)
          b1 = i3(1); lm1 = i3(2); s1 = i3(3)
          i = i + 1
!         -----------------------------------------------------------------------------
          if (s1 == s3 .and. s2 == s4) then
            doublegaunt0 = 0.d0
            do lm=1,lmmax0  ! truncated Ylm expansion
              doublegaunt0 = doublegaunt0 + gaunt(lm1,lm2,lm)*gaunt(lm3,lm4,lm)
            end do
            doublegaunt = doublegaunt0
            do lm=lmmax0+1,lmmax2  ! full Ylm expansion
              doublegaunt = doublegaunt + gaunt(lm1,lm2,lm)*gaunt(lm3,lm4,lm)
            end do
            if (abs(doublegaunt) > ylmtol) then
              write(iodb,'("big_overlaps: ia=",i4," nonzero q=",4i8,es16.8)') ia, q1, q2, q3, q4, doublegaunt
              work = phiref(:,b1,i2lm(2,lm1),s1,ia)*phiref(:,b2,i2lm(2,lm2),s2,ia)
              work = work*phiref(:,b3,i2lm(2,lm3),s3,ia)*phiref(:,b4,i2lm(2,lm4),s4,ia)
              work(1:nr) = work(1:nr)/rmesh(1:nr,ia)**2
              norm = radint(nr,work(1:nr),dr(1:nr))
              rho2 = rho2 + real(norm)*doublegaunt*gfsum(q1,q2)*conjg(gfsum(q4,q3))
              if (abs(doublegaunt0) > ylmtol) big_overlaps(i,j) = real(norm)*doublegaunt0
            end if
          end if
!         -----------------------------------------------------------------------------
        end do
      end do
    end do
  end do
! diagonalize; eigenvalues in ascending order
  lrwork = 2*nb*(1 + 3*nb) + 1
  liwork = 5*nb + 3
  allocate(evals(nb),rwork(lrwork),iwork(liwork))
  call dsyevd('V','U',nb,big_overlaps,nb,evals,rwork,lrwork,iwork,liwork,info)
  ievals = 0
  do ib=1,nb
    write(iodb,'("big_overlaps: ia,ib=",i4,i8," evals(ib)=",es16.8)') ia, ib, evals(ib)
    if (abs(evals(ib)) > tol*evals(nb)) then
      write(iodb,'("evec(ib)  ib ilm  is  jb jlm  js")')
      if (sum(big_overlaps(:,ib)) < 0.d0) big_overlaps(:,ib) = -big_overlaps(:,ib)
      do i=1,nb
        if (abs(big_overlaps(i,ib)) > 1.d-4) then
          write(iodb,'(f8.4,1000i4)') big_overlaps(i,ib), i2lmsb(:,lmsb2i2(1,i),ia), i2lmsb(:,lmsb2i2(2,i),ia)
        end if
      end do
      ievals = ievals + 1 
    end if
  end do
  write(iodb,'("big_overlaps: ia=",i4," evals > tol=",es10.1,i8)') ia, tol, ievals
! check how much of the density matrix is lost
  rho2proj = 0.d0
  do ib=nb-ievals+1,nb
    rho = 0.d0; rhoc = 0.d0
    i = 0
    do q2=1,nlmsba(ia)
      do q1=1,nlmsba(ia)
        i = i + 1
        rho  = rho  + big_overlaps(i,ib)*gfsum(q1,q2)
        rhoc = rhoc + big_overlaps(i,ib)*conjg(gfsum(q2,q1))
      end do
    end do
    rho2proj = rho2proj + evals(ib)*rho*rhoc
  end do
  write(iodb,'("big_overlaps: ia=",i4," rho2=",2es16.8,"  rho2proj=",2es16.8,"  diff=",es16.8)') ia, rho2, rho2proj, abs(rho2 - rho2proj)
  deallocate(big_overlaps,evals,rwork,iwork)
! All done!
  end subroutine overlaps_susc
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine out_overlap(ia,ja)
! Mutual overlaps between basis functions for ia and ja
! Only gives right answer if the same radial mesh applies to everyone

  implicit none

! --> susc atom
  integer(kind=i4b), intent(in) :: ia, ja
! -----------------------------------------------------------------
  integer(kind=i4b) :: ib, jb, il, jl, is, js, nr, nb, mb
  complex(kind=c8b) :: work(nrmax)
  complex(kind=c8b) :: pp(nbmax)

  nr = nrpts(ja)
  write(iodb,'("overlaps: ia, ja, is, js, il, jl")')
  do js=1,issusc(ja)
  do is=1,js
    do jl=0,nlmax
    do il=0,jl
      nb = iwsusc(jl,js,ja)
      mb = iwsusc(il,is,ia)
      if (nb > 0 .and. mb > 0) then
        write(iodb,'("block=",6i4)') ia, ja, is, js, il, jl
        do jb=1,nb
          pp = 0.d0
          do ib=1,mb
            work = phiref(:,jb,jl,js,ja)*phiref(:,ib,il,is,ia)
            pp(ib) = radint(nr,work(1:nr),drmesh(1:nr,ia))
          end do
          where (abs(pp(1:mb)) < gstol) pp = 0.d0
          write(iodb,'(100es10.3)') real(pp(1:mb))
        end do
      end if
    end do
    end do
  end do
  end do
! All done!
  end subroutine out_overlap
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  pure function radint(n,f,dx)
! Integrates array x with weights dx coming from radial mesh
! If x(n1:n2) and n = n2 - n1 + 1, can be used to integrate over an interval
! dx has the effect of the change of variable to an equidistant mesh

  implicit none

  integer(kind=i4b), intent(in) :: n
  complex(kind=c8b), intent(in) :: f(n)
  real(kind=r8b),    intent(in) :: dx(n)
  complex(kind=c8b) :: radint
! -----------------------------------------------------------------
  complex(kind=c8b) :: work(n)

  work = f*dx

  if (n < 2) then
! nothing
    radint = 0.d0
  else if (n == 2) then
! trapezoidal rule
    radint = 0.5d0*(work(1) + work(2))
  else if (mod(n,2) == 1) then
! normal extended simpson rule
    radint = 2.d0*sum(work(1:n:2))/3.d0
    radint = radint + 4.d0*sum(work(2:n:2))/3.d0
    radint = radint - (work(1) + work(n))/3.d0
  else
! trapezoidal rule for first two points, might be inaccurate
    radint = 0.5d0*(work(1) + work(2))
! extended simpson rule on n-1 points 
    radint = radint + 2.d0*sum(work(2:n:2))/3.d0
    radint = radint + 4.d0*sum(work(3:n:2))/3.d0
    radint = radint - (work(2) + work(n))/3.d0
  end if
! All done!
  end function radint
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine orbmoment(lmax,lmmax)
! Construction of the orbital and spin angular momentum matrices
! Also matrix elements in the lmsb basis

  implicit none

  integer(kind=i4b), intent(in) :: lmax, lmmax
! -----------------------------------------------------------------
  real(kind=r8b),    parameter :: invr2 = 1.d0/sqrt(2.d0)
  complex(kind=c8b), parameter :: iu = (0.d0,1.d0)
  integer(kind=i4b) :: l1, m1, lm1, l2, m2, lm2, i, j, i2(2)
  integer(kind=i4b) :: ib, ilm, is, jb, jlm, js, k, ia, i3(3)
  complex(kind=c8b) :: u(lmmax,lmmax), fac
  complex(kind=c8b) :: lp(lmmax,lmmax), lm(lmmax,lmmax)

! storage for matrix elements in lmsb basis
  allocate(ldots(nlmsb,nlmsb,nasusc))
! initialize pauli matrices
! x
  pauli(1,1,1) = ( 0.0d0, 0.0d0); pauli(1,2,1) = ( 1.0d0, 0.0d0)
  pauli(2,1,1) = ( 1.0d0, 0.0d0); pauli(2,2,1) = ( 0.0d0, 0.0d0)
! y
  pauli(1,1,2) = ( 0.0d0, 0.0d0); pauli(1,2,2) = ( 0.0d0, 1.0d0)
  pauli(2,1,2) = ( 0.0d0,-1.0d0); pauli(2,2,2) = ( 0.0d0, 0.0d0)
! z
  pauli(1,1,3) = (-1.0d0, 0.0d0); pauli(1,2,3) = ( 0.0d0, 0.0d0)
  pauli(2,1,3) = ( 0.0d0, 0.0d0); pauli(2,2,3) = ( 1.0d0, 0.0d0)
! n
  pauli(1,1,4) = ( 1.0d0, 0.0d0); pauli(1,2,4) = ( 0.0d0, 0.0d0)
  pauli(2,1,4) = ( 0.0d0, 0.0d0); pauli(2,2,4) = ( 1.0d0, 0.0d0)
! transformations from cartesian to spin labels and back
! potential cartesian to spin
  pc2s(1,1,1) = ( 0.0d0, 0.0d0); pc2s(1,1,2) = ( 0.0d0, 0.0d0); pc2s(1,1,3) = (-1.0d0, 0.0d0); pc2s(1,1,4) = ( 1.0d0, 0.0d0) 
  pc2s(2,1,1) = ( 1.0d0, 0.0d0); pc2s(2,1,2) = ( 0.0d0, 1.0d0); pc2s(2,1,3) = ( 0.0d0, 0.0d0); pc2s(2,1,4) = ( 0.0d0, 0.0d0) 
  pc2s(1,2,1) = ( 1.0d0, 0.0d0); pc2s(1,2,2) = ( 0.0d0,-1.0d0); pc2s(1,2,3) = ( 0.0d0, 0.0d0); pc2s(1,2,4) = ( 0.0d0, 0.0d0) 
  pc2s(2,2,1) = ( 0.0d0, 0.0d0); pc2s(2,2,2) = ( 0.0d0, 0.0d0); pc2s(2,2,3) = ( 1.0d0, 0.0d0); pc2s(2,2,4) = ( 1.0d0, 0.0d0) 
! potential spin to cartesian
  ps2c(1,1,1) = ( 0.0d0, 0.0d0); ps2c(2,1,1) = ( 0.0d0, 0.0d0); ps2c(3,1,1) = (-0.5d0, 0.0d0); ps2c(4,1,1) = ( 0.5d0, 0.0d0) 
  ps2c(1,2,1) = ( 0.5d0, 0.0d0); ps2c(2,2,1) = ( 0.0d0,-0.5d0); ps2c(3,2,1) = ( 0.0d0, 0.0d0); ps2c(4,2,1) = ( 0.0d0, 0.0d0) 
  ps2c(1,1,2) = ( 0.5d0, 0.0d0); ps2c(2,1,2) = ( 0.0d0, 0.5d0); ps2c(3,1,2) = ( 0.0d0, 0.0d0); ps2c(4,1,2) = ( 0.0d0, 0.0d0) 
  ps2c(1,2,2) = ( 0.0d0, 0.0d0); ps2c(2,2,2) = ( 0.0d0, 0.0d0); ps2c(3,2,2) = ( 0.5d0, 0.0d0); ps2c(4,2,2) = ( 0.5d0, 0.0d0) 
! density cartesian to spin
  dc2s(1,1,1) = ( 0.0d0, 0.0d0); dc2s(1,1,2) = ( 0.0d0, 0.0d0); dc2s(1,1,3) = (-0.5d0, 0.0d0); dc2s(1,1,4) = ( 0.5d0, 0.0d0) 
  dc2s(2,1,1) = ( 0.5d0, 0.0d0); dc2s(2,1,2) = ( 0.0d0,-0.5d0); dc2s(2,1,3) = ( 0.0d0, 0.0d0); dc2s(2,1,4) = ( 0.0d0, 0.0d0) 
  dc2s(1,2,1) = ( 0.5d0, 0.0d0); dc2s(1,2,2) = ( 0.0d0, 0.5d0); dc2s(1,2,3) = ( 0.0d0, 0.0d0); dc2s(1,2,4) = ( 0.0d0, 0.0d0) 
  dc2s(2,2,1) = ( 0.0d0, 0.0d0); dc2s(2,2,2) = ( 0.0d0, 0.0d0); dc2s(2,2,3) = ( 0.5d0, 0.0d0); dc2s(2,2,4) = ( 0.5d0, 0.0d0) 
! density spin to cartesian
  ds2c(1,1,1) = ( 0.0d0, 0.0d0); ds2c(2,1,1) = ( 0.0d0, 0.0d0); ds2c(3,1,1) = (-1.0d0, 0.0d0); ds2c(4,1,1) = ( 1.0d0, 0.0d0) 
  ds2c(1,2,1) = ( 1.0d0, 0.0d0); ds2c(2,2,1) = ( 0.0d0, 1.0d0); ds2c(3,2,1) = ( 0.0d0, 0.0d0); ds2c(4,2,1) = ( 0.0d0, 0.0d0) 
  ds2c(1,1,2) = ( 1.0d0, 0.0d0); ds2c(2,1,2) = ( 0.0d0,-1.0d0); ds2c(3,1,2) = ( 0.0d0, 0.0d0); ds2c(4,1,2) = ( 0.0d0, 0.0d0) 
  ds2c(1,2,2) = ( 0.0d0, 0.0d0); ds2c(2,2,2) = ( 0.0d0, 0.0d0); ds2c(3,2,2) = ( 1.0d0, 0.0d0); ds2c(4,2,2) = ( 1.0d0, 0.0d0) 

! orbital angular momentum
! 1 -> columns
  do lm1=1,lmmax
    i2 = i2lm(:,lm1)
    m1 = i2(1); l1 = i2(2)
! 2 -> rows
    do lm2=1,lmmax
      i2 = i2lm(:,lm2)
      m2 = i2(1); l2 = i2(2)
!   -----------------------------------------------------------
!   L+ in complex Ylm basis
!   row index larger than column index
      if (l2 == l1 .and. m2 == m1+1) then
        fac = l1*(l1+1) - m1*m2
        lp(lm2,lm1) = sqrt(fac)
!      write(iodb,'(6i4,2f6.1)') l2, m2, lm2, l1, m1, lm1, lp(lm2,lm1)
      else
        lp(lm2,lm1) = 0.d0
      end if
!   -----------------------------------------------------------
!   L- in complex Ylm basis
!   row index smaller than column index
      if (l2 == l1 .and. m2 == m1-1) then
        fac = l1*(l1+1) - m1*m2
        lm(lm2,lm1) = sqrt(fac)
!      write(iodb,'(6i4,2f6.1)') l2, m2, lm2, l1, m1, lm1, lm(lm2,lm1)
      else
        lm(lm2,lm1) = 0.d0
      end if
!   -----------------------------------------------------------
!   Lz in complex Ylm basis
!   column index equal to row index
      if (l2 == l1 .and. m2 == m1) then
        lorb(lm2,lm1,3) = m1
      else
        lorb(lm2,lm1,3) = 0.d0
      end if
!   -----------------------------------------------------------
!   Transformation matrix from complex to real Ylm
!   -----------------------------------------------------------
      if (l2 == l1 .and. m2 == 0 .and. m1 == 0) then
        u(lm2,lm1) = 1.d0
!  write(iodb,'(6i4,2f6.1)') l2, m2, lm2, l1, m1, lm1, u(lm2,lm1)
      else if (l2 == l1 .and. m2 ==  m1 .and. m1 < 0) then
        u(lm2,lm1) = iu*invr2
!  write(iodb,'(6i4,2f6.1)') l2, m2, lm2, l1, m1, lm1, u(lm2,lm1)
      else if (l2 == l1 .and. m2 == -m1 .and. m1 < 0) then
        u(lm2,lm1) = invr2
!  write(iodb,'(6i4,2f6.1)') l2, m2, lm2, l1, m1, lm1, u(lm2,lm1)
      else if (l2 == l1 .and. m2 ==  m1 .and. m1 > 0) then
        u(lm2,lm1) = invr2*(-1)**m1
!  write(iodb,'(6i4,2f6.1)') l2, m2, lm2, l1, m1, lm1, u(lm2,lm1)
      else if (l2 == l1 .and. m2 == -m1 .and. m1 > 0) then
        u(lm2,lm1) = -iu*invr2*(-1)**m1
!  write(iodb,'(6i4,2f6.1)') l2, m2, lm2, l1, m1, lm1, u(lm2,lm1)
      else
        u(lm2,lm1) = 0.d0
      end if
!   -----------------------------------------------------------
    end do
  end do
! -----------------------------------------------------------------
! Lx
  lorb(:,:,1) = 0.5d0*(lp + lm)
! Ly
  lorb(:,:,2) = 0.5d0*(lp - lm)/iu
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  write(iodb,*) "orbmom", lmax, lmmax
  write(iodb,*) "Angular momentum matrices in complex Ylm basis"
  write(iodb,*) "L+"
  do lm1=1,lmmax
    write(iodb,'(100f4.1)') lp(lm1,1:lmmax)
  end do
  write(iodb,*) "L-"
  do lm1=1,lmmax
    write(iodb,'(100f4.1)') lm(lm1,1:lmmax)
  end do
  write(iodb,*) "Lx"
  do lm1=1,lmmax
    write(iodb,'(100f4.1)') lorb(lm1,1:lmmax,1)
  end do
  write(iodb,*) "Ly"
  do lm1=1,lmmax
    write(iodb,'(100f4.1)') lorb(lm1,1:lmmax,2)
  end do
  write(iodb,*) "Lz"
  do lm1=1,lmmax
    write(iodb,'(100f4.1)') lorb(lm1,1:lmmax,3)
  end do
  lp = 0.d0
  do k=1,3
    lp = lp + matmul(lorb(:,:,k),lorb(:,:,k))
  end do
  write(iodb,*) "L^2"
  do lm1=1,lmmax
    write(iodb,'(100f4.1)') lp(1:lmmax,lm1)
  end do
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  write(iodb,*) "Angular momentum matrices in real Ylm basis"
  write(iodb,*) "U"
  do lm1=1,lmmax
    write(iodb,'(100f4.1)') u(lm1,1:lmmax)
  end do
  write(iodb,*) "U^dagger"
  do lm1=1,lmmax
    write(iodb,'(100f4.1)') conjg(u(1:lmmax,lm1))
  end do
  write(iodb,*) "U^dagger U"
  lp = matmul(conjg(transpose(u)),u)
  do lm1=1,lmmax
    write(iodb,'(100f4.1)') lp(1:lmmax,lm1)
  end do
  do k=1,3
    lorb(:,:,k) = matmul(conjg(u),matmul(lorb(:,:,k),transpose(u)))
  end do
  write(iodb,*) "Lx"
  do lm1=1,lmmax
    write(iodb,'(100f4.1)') lorb(lm1,1:lmmax,1)
  end do
  write(iodb,*) "Ly"
  do lm1=1,lmmax
    write(iodb,'(100f4.1)') lorb(lm1,1:lmmax,2)
  end do
  write(iodb,*) "Lz"
  do lm1=1,lmmax
    write(iodb,'(100f4.1)') lorb(lm1,1:lmmax,3)
  end do
  lp = 0.d0
  do k=1,3
    lp = lp + matmul(lorb(:,:,k),lorb(:,:,k))
  end do
  write(iodb,*) "L^2"
  do lm1=1,lmmax
    write(iodb,'(100f4.1)') lp(1:lmmax,lm1)
  end do
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Now put together the matrix elements of L.S
  ldots = 0.d0
  do ia=1,nasusc
    do j = 1,nlmsba(ia)
      i3 = i2lmsb(:,j,ia)
      jb = i3(1); jlm = i3(2); js = i3(3)
      do i = 1,nlmsba(ia)
        i3 = i2lmsb(:,i,ia)
        ib = i3(1); ilm = i3(2); is = i3(3)
!        if (ib == jb) then
        do k=1,3
          ldots(i,j,ia) = ldots(i,j,ia) + lorb(ilm,jlm,k)*pauli(is,js,k)
        end do
!        end if
      end do
    end do
  end do
! All done
  end subroutine orbmoment
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine projected_gf(ie,ia,ja,gf,onsite,struct)
! Block of the projected GF
! This version does not interpolate yet

  implicit none

! which energy
  integer(kind=i4b), intent(in)  :: ie
! which atoms
  integer(kind=i4b), intent(in)  :: ia, ja
! the block of the GF in the projection basis
  complex(kind=c8b), intent(out) :: gf(nlmsb,nlmsb)
! whether to include the onsite and the structural parts
  logical,           intent(in)  :: onsite, struct
! -------------------------------------------------------------------
  complex(kind=c8b) :: left(nlms), right(nlms)
  complex(kind=c8b) :: gs(nlms,nlms), block(4,4)
! decoding
  integer(kind=i4b) :: j, jlms, l
  integer(kind=i4b) :: i, ilms, k

  gf = 0.d0
! fetch the structural GF block if needed
  if (struct) then
    do jlms=1,nlms
      j = alms2i(jlms,ja)
      do ilms=1,nlms
        i = alms2i(ilms,ia)
        gs(ilms,jlms) = gstruct(i,j,ie)
      end do
    end do
  end if
! put together the projected GF
  do j=1,nlmsba(ja)
    do i=1,nlmsba(ia)
!   ****    add the onsite part    ****
      if (ia == ja .and. onsite) then
        gf(i,j) = gf(i,j) + gfpq(i,j,ia,ie)
      end if
!   ****  add the structural part  ****
      if (struct) then
        left  = pzl(j,:,ja,ie)
        right = pzr(i,:,ia,ie)
        do k=1,nlms
          do l=1,nlms
            gf(i,j) = gf(i,j) + right(k)*gs(k,l)*left(l)
          end do
        end do
      end if
!   ***********************************
    end do
  end do
! All done!
  end subroutine projected_gf 
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine symmetrize(n,a,tol)
! Compare all elements of a with each other

  implicit none

  integer(kind=i4b), intent(in)    :: n
  complex(kind=c8b), intent(inout) :: a(n,n)
  real(kind=r8b),    intent(in)    :: tol
! ---------------------------------------
  integer(kind=i4b) :: i1, j1, i2, j2
  real(kind=r8b)    :: re1, im1, re2, im2, avg

  do j2=1,n
  do i2=1,n
    do j1=1,n
    do i1=1,n
      re2 = real(a(i2,j2)); im2 = aimag(a(i2,j2))
      re1 = real(a(i1,j1)); im1 = aimag(a(i1,j1))
      if (abs(re1) < tol) re1 = 0.d0
      if (abs(im1) < tol) im1 = 0.d0
      if (abs(re2) < tol) re2 = 0.d0
      if (abs(im2) < tol) im2 = 0.d0
!      if (abs(abs(re1) - abs(re2)) < tol) then
!        avg = 0.5d0*(abs(re1) + abs(re2))
!        re1 = sign(avg,re1)
!        re2 = sign(avg,re2)
!      end if
!      if (abs(abs(im1) - abs(im2)) < tol) then
!        avg = 0.5d0*(abs(im1) + abs(im2))
!        im1 = sign(avg,im1)
!        im2 = sign(avg,im2)
!      end if
!      if (abs(abs(re1) - abs(im2)) < tol) then
!        avg = 0.5d0*(abs(re1) + abs(im2))
!        re1 = sign(avg,re1)
!        im2 = sign(avg,im2)
!      end if
!      if (abs(abs(im1) - abs(re2)) < tol) then
!        avg = 0.5d0*(abs(im1) + abs(re2))
!        im1 = sign(avg,im1)
!        re2 = sign(avg,re2)
!      end if
      a(i1,j1) = cmplx(re1,im1)
      a(i2,j2) = cmplx(re2,im2)
    end do
    end do
  end do
  end do
! All done!
  end subroutine symmetrize
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine groundstate()
! Ground state properties
! The matrix elements of the operators are filled in subroutine observables

  implicit none

! temporary parameters
!  integer(kind=i4b), parameter :: ne = 1001, nesusc = 300
! i 2pi
  real(kind=r8b),    parameter :: twopi = 8.d0*atan(1.d0)
  complex(kind=c8b), parameter :: i2pi = (0.d0,twopi)
!  real(kind=r8b),    parameter :: omegamax = 4.d-1
! auxiliary
  complex(kind=c8b) :: de, work(nrmax), ze, za, zb, tmp
  complex(kind=c8b) :: gf(nlmsb,nlmsb), gfsum(nlmsb,nlmsb)
  integer(kind=i4b) :: nr, imap(nlmsb,nlmsb)
  real(kind=r8b)    :: r(nrmax), zre, zim, norm, start, finish
  real(kind=r8b)    :: omega
! moments of charge density, spin and orbital moments
  real(kind=r8b)    :: phi1(nrmax), phi2(nrmax)
  complex(kind=c8b) :: rho(nrmax,0:6,lmmax,lmmax)
  complex(kind=c8b) :: rhoint(0:6,lmmax,lmmax)
  complex(kind=c8b) :: suscylm(4,4,lmmax0,lmmax0,nasusc,nasusc)
  complex(kind=c8b) :: suscy00(4,4,lmmax,lmmax,nasusc,nasusc)
  complex(kind=c8b) :: suscden(4,lmmax,nasusc)
  complex(kind=c8b) :: kxcylm(4,4,lmmax0,nasusc)
  complex(kind=c8b) :: kxcy00(4,4,lmmax,lmmax,nasusc)
  complex(kind=c8b) :: zdos(lmmax,lmmax,nsmax,nsmax,nasusc,nedos)
  complex(kind=c8b) :: zdosl(0:nlmax,nsmax,nasusc,nedos)
  real(kind=r8b)    :: rhotot(0:6,0:nlmax)
  real(kind=r8b)    :: rhoylm(0:6,lmmax2)
  real(kind=r8b)    :: kxc(nrmax,nasusc)
! looping
  integer(kind=i4b) :: ia, ja, ie, k, l, ir, klm
! decoding
  integer(kind=i4b) :: i2(2), i3(3), dsusc
  integer(kind=i4b) :: j, jlms, js, jlm, jl, jm, jb, jalmsb
  integer(kind=i4b) :: i, ilms, is, ilm, il, im, ib, ialmsb
  ! new character to save ldos file for different ia-values:
  character*15         :: filename
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Get fitting coefficients
  call cpu_time(start)
  do ja=1,nasusc
    do ia=1,nasusc
      call ratfit_gf(ia,ja,.true.,.true.)
    end do
  end do
  call cpu_time(finish)
  write(*,'(/,"GF rational fit time=",f10.3," s",/)') finish - start
! ----------------------------------------------------------------------
  call cpu_time(start)
  nrv = 0.d0; zdosl = 0.d0; zdos = 0.d0
  do ia=1,nasusc
!    write(*,*) "before GF"
    nr = nrpts(ia)
    gfsum = 0.d0
! use the energy mesh from SCF calculation here
! construct the energy integrated site-diagonal GF
    do ie=1,nescf
      de = descf(ie)
!     ****  site diagonal part for local quantities  ****
! JULEN one can choose between commenting one of the two lines below, 
!   projected_gf : does not interpolate, it uses the projection scheme and recalculate green function
!   ratval_gf: interpolates green functions, 
!      call projected_gf(ie,ia,ia,gf,.true.,.true.)
      call ratval_gf(ia,ia,escf(ie),gf)
!       (G^dagger - G)/(i 2pi)
      gf = (conjg(transpose(de*gf)) - de*gf)/i2pi
      gfsum = gfsum + gf
!     ***************************************************
    end do
! ================================================================
! Compute DOS
    if (idos == 1) then
!   ===================
    za = cmplx(e0dos,eimdos); zb = cmplx(e1dos,eimdos)
    do ie=1,nedos
      ze = za + (zb-za)*(ie-1.d0)/(nedos-1.d0)
!      call projected_gf(ie,ia,ia,gf,.true.,.true.)
      call ratval_gf(ia,ia,ze,gf)
!     DOS
      do j=1,nlmsba(ia)
        i3 = i2lmsb(:,j,ia)
        jb = i3(1); jlm = i3(2); js = i3(3)
        i2 = i2lm(:,jlm)
        jm = i2(1); jl = i2(2)
        phi1(1:nr) = phiref(1:nr,jb,jl,js,ia)
        do i=1,nlmsba(ia)
          i3 = i2lmsb(:,i,ia)
          ib = i3(1); ilm = i3(2); is = i3(3)
          i2 = i2lm(:,ilm)
          im = i2(1); il = i2(2)
          zdos(ilm,jlm,is,js,ia,ie) = zdos(ilm,jlm,is,js,ia,ie) + gf(i,j)*overlap(i,j,ia)
          if (ilm == jlm .and. is == js) zdosl(il,is,ia,ie) = zdosl(il,is,ia,ie) + gf(i,j)*overlap(i,j,ia)
        end do
      end do
    end do
! DOS to file
    write(filename,'("ldos_ia",i4.4,".dat")') ia
    open(file=filename,unit=iofile)
      write(iofile,'("# Re(energy)  totalDOS_down totalDOS_up  sDOS_down pDOS_down dDOS_down (fDOS_down) sDOS_up pDOS_up dDOS_up (fDOS_up)")')
      write(iofile,'("# imaginary part of energy is ",es16.8)') eimdos! == aimag(za) == aimag(zb)
      do ie=1,nedos
        ze = za + (zb-za)*(ie-1.d0)/(nedos-1.d0)
        write(iofile,'(100es16.8)') dreal(ze), sum(zdosl(:,1,ia,ie)), sum(zdosl(:,2,ia,ie)), zdosl(:,:,ia,ie)
      end do
    close(iofile)
!   ======
    end if
! ================================================================
! filter small elements of the GF away
!    imap = 0
    do j=1,nlmsba(ia)
      do i=1,nlmsba(ia)
        zre = real(gfsum(i,j)); zim = aimag(gfsum(i,j))
        if (abs(zre) < gstol) zre = 0.d0
        if (abs(zim) < gstol) zim = 0.d0
        gfsum(i,j) = cmplx(zre,zim)
      end do
    end do
!    write(*,*) "before density matrix"
! ----------------------------------------------------------------
! xc-potential and energy
    kxc = 0.d0
    call overlaps_susc(ia,gfsum)
!    call get_xc(ia,lmmax2,gfsum,kxc)
! ----------------------------------------------------------------
! density matrix and its radial integral (COMPLEX)
    rho = 0.d0; rhoint = 0.d0
    do j=1,nlmsba(ia)
      i3 = i2lmsb(:,j,ia)
      jb = i3(1); jlm = i3(2); js = i3(3)
      phi1(1:nr) = phiref(1:nr,jb,i2lm(2,jlm),js,ia)
      do i=1,nlmsba(ia)
        i3 = i2lmsb(:,i,ia)
        ib = i3(1); ilm = i3(2); is = i3(3)
        phi2(1:nr) = phiref(1:nr,ib,i2lm(2,ilm),is,ia)
        do k=1,4
          rho(1:nr,mod(k,4),ilm,jlm) = rho(1:nr,mod(k,4),ilm,jlm)   &
                                       + gfsum(i,j)*pauli(is,js,k)*phi1(1:nr)*phi2(1:nr)
          rhoint(mod(k,4),ilm,jlm)   = rhoint(mod(k,4),ilm,jlm)     &
                                       + gfsum(i,j)*pauli(is,js,k)*overlap(i,j,ia)
        end do
      end do
    end do
! extra: orbital densities and moments - only charge part (spin trace)
! ----------------------------------------------------------------
    do jlm=1,lmmax
    do ilm=1,lmmax
      do k=4,6
      do klm=1,lmmax
        rho(1:nr,k,ilm,jlm) = rho(1:nr,k,ilm,jlm) + lorb(ilm,klm,k-3)*rho(1:nr,0,klm,jlm)
        rhoint(k,ilm,jlm)   = rhoint(k,ilm,jlm)   + lorb(ilm,klm,k-3)*rhoint(0,klm,jlm)
      end do
      end do
    end do
    end do
!   ----------------------------------------------------------------
!    do i=1,nlmsb
!      write(*,'(1000i2)') imap(i,1:nlmsb)
!    end do
!    do i=1,nlmsb
!      write(*,'("gfsum ",1000es10.1)') gfsum(i,1:nlmsb)
!    end do
! this is the lm decomposition of the spherical part
    rhotot = 0.d0
    do ilm=1,lmmax
      i2 = i2lm(:,ilm)
      im = i2(1); il = i2(2)
      do k=0,6
        rhotot(k,il) = rhotot(k,il) + rhoint(k,ilm,ilm)
      end do
    end do
! this is the spherical harmonics resummation
    rhoylm = 0.d0
    do klm=1,lmmax2
      do jlm=1,lmmax
      do ilm=1,lmmax
        norm = gaunt(ilm,jlm,klm)
        if (abs(norm) > ylmtol) then
          rhoylm(:,klm) = rhoylm(:,klm) + rhoint(:,ilm,jlm)*norm
        end if
      end do
      end do
    end do
! and the winner is...
    write(*,'("GS quantities for ia =",i4,": sum, spdf")') ia
    write(*,'("  qe =",100f16.8)') sum(rhotot(0,:)), rhotot(0,:)
    write(*,'("  msx=",100f16.8)') sum(rhotot(1,:)), rhotot(1,:)
    write(*,'("  msy=",100f16.8)') sum(rhotot(2,:)), rhotot(2,:)
    write(*,'("  msz=",100f16.8)') sum(rhotot(3,:)), rhotot(3,:)
    write(*,'("  mox=",100f16.8)') sum(rhotot(4,:)), rhotot(4,:)
    write(*,'("  moy=",100f16.8)') sum(rhotot(5,:)), rhotot(5,:)
    write(*,'("  moz=",100f16.8)') sum(rhotot(6,:)), rhotot(6,:)
! output decomposition
    write(*,'("decomposition into spherical harmonics")')
    write(*,'("  qeylm=")')
    do ilm=1,lmmax2
      norm = rhoylm(0,ilm)
      if (abs(norm) > atol) then
        write(*,'(2i4,2f16.8)') i2lm(:,ilm), norm
      end if
    end do
    write(*,'("  msylm=")')
    do ilm=1,lmmax2
      norm = sqrt(dot_product(rhoylm(1:3,ilm),rhoylm(1:3,ilm)))
      if (norm > atol) then
        write(*,'(2i4,8f16.8)') i2lm(:,ilm), norm, rhoylm(1:3,ilm)
      end if
    end do
    write(*,'("  moylm=")')
    do ilm=1,lmmax2
      norm = sqrt(dot_product(rhoylm(4:6,ilm),rhoylm(4:6,ilm)))
      if (norm > atol) then
        write(*,'(2i4,8f16.8)') i2lm(:,ilm), norm, rhoylm(4:6,ilm)
      end if
    end do
! save spherical part of valence densities
    open(100+ia)
      nr = nrpts(ia)
      r = rmesh(1:nr,ia)
      do ilm=1,lmmax
        nrv(1:nr,:,ia) = nrv(1:nr,:,ia) + real(rho(1:nr,:,ilm,ilm))
      end do
      do ir=1,nr
        write(100+ia,'(8es16.8)') r(ir), nrv(ir,:,ia)!/r(ir)**2
      end do
    close(100+ia)
  end do
  call cpu_time(finish)
  write(*,'(/," GS quantities  time=",f10.3," s",/)') finish - start
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Now for the static kernel
  if (ikxc > 0) then
  call cpu_time(start)
  call build_kxcalda(kxc,kxcylm,kxcy00)
  write(*,'(/,"xc kernel   ia, m1, l1, m2, l2, spherical")')
  do ia=1,nasusc
    do jlm=1,lmmax
    do ilm=1,lmmax
      norm = sum(abs(kxcy00(:,:,ilm,jlm,ia)))
      if (norm > atol) then
        do j=1,4
          write(*,'(5i4,8f16.8)') ia, i2lm(:,ilm), i2lm(:,jlm), (kxcy00(i,j,ilm,jlm,ia),i=1,4)
        end do
      end if
    end do
    end do
  end do
  write(*,'("xc kernel   ia, im, il, jm, jl, kxc  ylm")')
  do ia=1,nasusc
    do jlm=1,lmmax0
      norm = sum(abs(kxcylm(:,:,jlm,ia)))
      if (norm > atol) then
      do j=1,4
        write(*,'(3i4,8f16.8)') ia, i2lm(:,jlm), (kxcylm(i,j,jlm,ia),i=1,4)
      end do
      end if
    end do
  end do
!  write(*,'("xc kernel   ia, ib, im, il, is, jb, jm, jl, js, kxcsusc(i,j,j,i)")')
!  do ia=1,nasusc
!    do j=1,nlmsb
!      i3 = i2lmsb(:,j)
!      jb = i3(1); jlm = i3(2); js = i3(3)
!      i2 = i2lm(:,jlm)
!      jm = i2(1); jl = i2(2)
!      if (jb <= iwsusc(jl,js,ia)) then
!      do i=1,nlmsb
!        i3 = i2lmsb(:,i)
!        ib = i3(1); ilm = i3(2); is = i3(3)
!        i2 = i2lm(:,ilm)
!        im = i2(1); il = i2(2)
!        if (ib <= iwsusc(il,is,ia)) then
!          if (abs(kxcsusc(i,j,j,i)) > atol) write(*,'(9i4,2f12.6)') ia, ib, im, il, is, jb, jm, jl, js, kxcsusc(i,j,j,i)
!        end if
!      end do
!      end if
!    end do
!  end do
  call cpu_time(finish)
  write(*,'(/," xc ALDA kernel time=",f10.3," s",/)') finish - start
  end if
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Now for the static KS susc
  if (isusc > 0) then  
! ++++++++++++++++++++++++++++++++++++
  call cpu_time(start)
! static susc
  call dyn_susc_real(0.d0,suscylm,suscy00,suscden,.true.,.false.,.false.)
  write(*,'("KS susc ia, ja, m1, l1, m2, l2, spherical")')
  do ja=1,nasusc
  do ia=1,nasusc
    do jlm=1,lmmax
    do ilm=1,lmmax
      norm = sum(abs(suscy00(:,:,ilm,jlm,ia,ja)))
      if (norm > atol) then
      do j=1,4
        write(*,'(6i4,8f16.8)') ia, ja, i2lm(:,ilm), i2lm(:,jlm), (suscy00(i,j,ilm,jlm,ia,ja),i=1,4)
      end do
      end if
    end do
    end do
  end do
  end do
  write(*,'(/,"den = KS susc * eff pot")')
  do ia=1,nasusc
  do ilm=1,lmmax
    write(*,'(3i4,8f16.8)') ia, i2lm(:,ilm), (suscden(i,ilm,ia),i=1,4)
  end do
  end do
  write(*,'(/,"KS susc ia, ja, im, il, jm, jl, susc ylm")')
  do ja=1,nasusc
  do ia=1,nasusc
    do jlm=1,lmmax0
    do ilm=1,lmmax0
      norm = sum(abs(suscylm(:,:,ilm,jlm,ia,ja)))
      if (norm > atol) then
      do j=1,4
        write(*,'(6i4,8f16.8)') ia, ja, i2lm(:,ilm), i2lm(:,jlm), (suscylm(i,j,ilm,jlm,ia,ja),i=1,4)
      end do
      end if
    end do
    end do
  end do
  end do
  call cpu_time(finish)
  write(*,'(/," Static  KS susc time=",f10.3," s",/)') finish - start
! ++++++
  end if
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Analyze the denominator
  if (isusc == 1 .and. ikxc == 1) then
  call cpu_time(start)
  call susc_denominator
  call cpu_time(finish)
  write(*,'(/," susc denominator time=",f10.3," s",/)') finish - start
  end if
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Now for the dynamic KS susc
! don't compute the susceptibility twice
  if (isoc == 0 .or. soc_applied) then  
! ++++++++++++++++++++++++++++++++++++
  call cpu_time(start)
! static susc
! Go over omegas
  open(file='dynsusc.dat',unit=iofile)
  do ie=1,nomega
    omega = omegamin + (omegamax - omegamin)*(ie - 1.d0)/(nomega - 1.d0)
    call dyn_susc_real(omega,suscylm,suscy00,suscden,analytic,nonanalytic,enhanced)
    tmp = 0.d0
    do k=1,nalmsb
      tmp = tmp + kssusc(k,k)
    end do
    write(iofile,'(1000es16.8)') omega, tmp, ((suscylm(i,j,1,1,1,1),i=1,4),j=1,4)
  end do
  close(iofile)
  call cpu_time(finish)
  write(*,'(/," Dynamic KS susc time=",f10.3," s",/)') finish - start
! ++++++
  end if
! ++++++
! All done
  end subroutine groundstate
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine soc_correction()
! SOC correction to the GF

  implicit none

! speed of light
  real(kind=r8b),    parameter :: c = 274.0720442d0
! complex parameters
  complex(kind=c8b), parameter :: cplus  = ( 1.d0, 0.d0)
  complex(kind=c8b), parameter :: cminus = (-1.d0, 0.d0)
  complex(kind=c8b), parameter :: czero  = ( 0.d0, 0.d0)
! change in t-matrix
  complex(kind=c8b) :: tmat(nlms,nlms,nasusc)
  complex(kind=c8b) :: tmatl(nlms,nlms), tmatr(nlms,nlms)
! SOC potential
  complex(kind=c8b) :: vsoc(nrmax), vsocb(nlmsb,nlmsb)
! auxiliary
  complex(kind=c8b) :: work1(nlmsb,nlmsb), work2(nalms,nalms)
  complex(kind=c8b) :: gfpqsoc(nlmsb,nlmsb)
  complex(kind=c8b) :: pzlsoc(nlmsb,nlms), pzrsoc(nlmsb,nlms)
  complex(kind=c8b) :: edummy
  real(kind=r8b)    :: maxelem, start, finish
! -----------------------------------------------------------------
  integer(kind=i4b) :: i, j, ia, is, js, ja, k
  integer(kind=i4b) :: jlm, ilm, ilms, jlms, klms
  integer(kind=i4b) :: ib, jb, ie, i3(3), j3(3)
  integer(kind=i4b) :: ipiv(nlmsb), jpiv(nalms), info

  call cpu_time(start)
! loop over energies
  do ie=1,nesusc
    write(iodb,*) "SOC correction ie=", ie
! loop over atoms
    do ia=1,nasusc
!   construct SOC potential
      edummy = escf(ie)
!      edummy = 0.d0
      call build_vsoc(ia,c,edummy,vsoc)
!   compute the matrix elements in the basis
      call build_vsocb(ia,vsoc,vsocb)
!   ********************************************************************
!              solve Lippmann-Schwinger and Dyson equations
!   ********************************************************************
!     --> construct 1 - G0.Vsoc
      work1 = 0.d0
      do i=1,nlmsba(ia)
        work1(i,i) = 1.d0
      end do
      call zgemm('N','N',nlmsba(ia),nlmsba(ia),nlmsba(ia),cminus,gfpq(:,:,ia,ie),nlmsb,vsocb,nlmsb,cplus,work1,nlmsb)
!     --> LU decomposition of 1 - G0.Vsoc
      call zgetrf(nlmsba(ia),nlmsba(ia),work1,nlmsb,ipiv,info)
      if (info /= 0) then
        write(*,*) "LU fail 1 - G0.Vsoc at ie,ia=", ie, ia, info
        stop
      end if
!     --> new single-site GF: (1 - G0.Vsoc).G = G0
      gfpqsoc = gfpq(:,:,ia,ie)
      call zgetrs('N',nlmsba(ia),nlmsba(ia),work1,nlmsb,ipiv,gfpqsoc,nlmsb,info)
      if (info /= 0) then
        write(*,*) "SOC single-site GF failed at ie,ia=", ie, ia
        stop
      end if
!     --> solve Lipp-Schw for rhs wfn: (1 - G0.Vsoc).R^r = R0^r
!       nlms are the boundary conditions
      pzrsoc = pzr(:,:,ia,ie)  ! the initial rhs solution is needed in t^l
      call zgetrs('N',nlmsba(ia),nlms,work1,nlmsb,ipiv,pzrsoc,nlmsb,info)
      if (info /= 0) then
        write(*,*) "SOC rhs Lipp Schw failed at ie,ia=", ie, ia
        stop
      end if
!      write(iodb,*) "non-zero elements of rhs wfn"
!      write(iodb,'("  ib ilm  is jlm  js      pzr       pzrsoc")')
      maxelem = maxval(abs(pzr(:,:,ia,ie)))
!      write(iodb,'("maxelem  ",es14.6)') maxelem
      do j=1,nlms
        do i=1,nlmsba(ia)
          if (abs(pzrsoc(i,j)) > soctol*maxelem) then
!        write(iodb,'(5i4,4es14.6)') i2lmsb(:,i), i2lms(:,j), pzr(i,j,ia,ie), pzrsoc(i,j)
          else
            pzrsoc(i,j) = 0.d0
          end if
        end do
      end do
!     --> change in t-matrix: t^r = R0^l.Vsoc.R^r
      tmatr = 0.d0; work1 = 0.d0
!     Vsoc.R^r -> work
      call zgemm('N','N',nlmsba(ia),nlms,nlmsba(ia),cplus,vsocb,nlmsb,pzrsoc,nlmsb,czero,work1,nlmsb)
!     R0^l.work -> t^r
      call zgemm('T','N',nlms,nlms,nlmsba(ia),cplus,pzl(:,:,ia,ie),nlmsb,work1,nlmsb,czero,tmatr,nlms)
!      write(iodb,*) "tmatr for ia=", ia
!      do j=1,nlms
!        write(iodb,'(100es10.2)') tmatr(:,j)
!      end do
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     --> construct 1 - G0^T.Vsoc^T
      work1 = 0.d0
      do i=1,nlmsba(ia)
        work1(i,i) = 1.d0
      end do
      call zgemm('T','T',nlmsba(ia),nlmsba(ia),nlmsba(ia),cminus,gfpq(:,:,ia,ie),nlmsb,vsocb,nlmsb,cplus,work1,nlmsb)
!     --> LU decomposition of 1 - G0^T.Vsoc^T
      call zgetrf(nlmsba(ia),nlmsba(ia),work1,nlmsb,ipiv,info)
      if (info /= 0) then
        write(*,*) "LU fail 1 - G0^T.Vsoc^T at ie,ia=", ie, ia, info
        stop
      end if
!     --> solve Lipp-Schw for left hand side wfn: (1 - G0^T.Vsoc^T).R^l = R0^l
!       nlms are the boundary conditions
      pzlsoc = pzl(:,:,ia,ie)
      call zgetrs('N',nlmsba(ia),nlms,work1,nlmsb,ipiv,pzlsoc,nlmsb,info)
      if (info /= 0) then
        write(*,*) "SOC lhs Lipp Schw failed at ie,ia=", ie, ia
        stop
      end if
!      write(iodb,*) "non-zero elements of lhs wfn"
!      write(iodb,'("  ib ilm  is jlm  js      pzl       pzlsoc")')
      maxelem = maxval(abs(pzl(:,:,ia,ie)))
!      write(iodb,'("maxelem  ",es14.6)') maxelem
      do j=1,nlms
        do i=1,nlmsba(ia)
          if (abs(pzlsoc(i,j)) > soctol*maxelem) then
!            write(iodb,'(5i4,4es14.6)') i2lmsb(:,i), i2lms(:,j), pzl(i,j,ia,ie), pzlsoc(i,j)
          else
            pzlsoc(i,j) = 0.d0
          end if
        end do
      end do
!     --> change in t-matrix: t^l = R^l.Vsoc.R0^r
      tmatl = 0.d0; work1 = 0.d0
!     R^l.Vsoc -> work
      call zgemm('T','N',nlms,nlmsba(ia),nlmsba(ia),cplus,pzlsoc,nlmsb,vsocb,nlmsb,czero,work1,nlmsb)
!     work.R0^r -> t^l
      call zgemm('N','N',nlms,nlms,nlmsba(ia),cplus,work1,nlmsb,pzr(:,:,ia,ie),nlmsb,czero,tmatl,nlms)
!      write(iodb,*) "tmatl for ia=", ia
!      do j=1,nlms
!        write(iodb,'(100es10.2)') tmatl(:,j)
!      end do
      maxelem = maxval(abs(tmatl - tmatr))
      write(iodb,'("|tmatl - tmatr| for ia=",i4," is ",es14.6)') ia, maxelem
!      do j=1,nlms
!        write(iodb,'(100es10.2)') abs(tmatl(:,j) - tmatr(:,j))
!      end do
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     --> update wfns
      pzl(:,:,ia,ie) = pzlsoc
      pzr(:,:,ia,ie) = pzrsoc
!     --> update t-matrix
!      if (ia == 1) tmat(:,:,ia) = tmatr
      tmat(:,:,ia) = tmatr
!     --> update single-site GF
      gfpq(:,:,ia,ie) = gfpqsoc
    end do
!   new structural GF: (1 - G0^s.t).G^s = G0^s
!   --> construct 1 - G0^s.t and G0^s
    do ja=1,nasusc
    do jlms=1,nlms
      j = alms2i(jlms,ja)
      do ia=1,nasusc
      do ilms=1,nlms
        i = alms2i(ilms,ia)
        work2(i,j) = 0.d0
        if (i == j) work2(i,j) = 1.d0
!     matrix multiplication
        do klms=1,nlms
          k = alms2i(klms,ja)
          work2(i,j) = work2(i,j) - gstruct(i,k,ie)*tmat(klms,jlms,ja)
        end do
      end do
      end do
    end do
    end do
!     --> solve Dyson equation
    call zgesv(nalms,nalms,work2,nalms,jpiv,gstruct(:,:,ie),nalms,info)
    if (info /= 0) then
      write(*,*) "SOC struct Dyson failed at ie=", ie
      stop
    end if
  end do
! Flag SOC as applied
  soc_applied = .true.
  call cpu_time(finish)
  write(*,'(/," SOC correction time=",f10.3," s",/)') finish - start
! All done!
  end subroutine soc_correction
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine build_vsoc(ia,c,e,vsoc)
! Assembles the radial part of the SOC potential
! vsoc = 1/(2Mc)^2 1/r dV/dr
! M = m + (E - V(r))/(2c^2)
! m = 1/2, e^2 = 2 in Ry units
! e-e potential     is repulsive  -> positive
! nuclear potential is attractive -> negative

  implicit none

! Which atom
  integer(kind=i4b), intent(in)  :: ia
! Speed of light
  real(kind=r8b),    intent(in)  :: c
! Complex energy
  complex(kind=c8b), intent(in)  :: e
! Complex SOC potential
  complex(kind=c8b), intent(out) :: vsoc(nrmax)
! -------------------------------------------------------------------
  complex(kind=c8b) :: mass(nrmax)
  real(kind=r8b)    :: dvdr(nrmax), r(nrmax), v(nrmax)
  integer(kind=i4b) :: ir, nr

  nr = nrpts(ia)
  r(1:nr) = rmesh(1:nr,ia)
  v(1:nr) = vr(1:nr,ia)
! Relativistic mass multiplied by c
  mass(1:nr) = c + (e - v(1:nr) + 2.d0*zat(ia)/r(1:nr))/c
! dV/dr for e-e potential
! forward difference
  dvdr(1) = (v(2)-v(1))/(r(2)-r(1))
! backward difference
  dvdr(nr) = (v(nr)-v(nr-1))/(r(nr)-r(nr-1))
! centered differences
  do ir=2,nr-1
    dvdr(ir) = 0.5d0*(v(ir+1)-v(ir))/(r(ir+1)-r(ir))
    dvdr(ir) = dvdr(ir) + 0.5d0*(v(ir)-v(ir-1))/(r(ir)-r(ir-1))
  end do
! 1/r dV/dr
  vsoc(1:nr) = 2.d0*zat(ia)/(r(1:nr)**3) + dvdr/r(1:nr)
! multiply by the inverse relativistic mass
  vsoc(1:nr) = vsoc(1:nr)/(mass(1:nr))**2
! All done!
  end subroutine build_vsoc
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine build_vsocb(ia,vsoc,vsocb)
! assembles the full SOC kernel in the projection basis

  implicit none

  integer(kind=i4b), intent(in)  :: ia
  complex(kind=c8b), intent(in)  :: vsoc(nrmax)
  complex(kind=c8b), intent(out) :: vsocb(nlmsb,nlmsb)
! -----------------------------------------------------------------
  integer(kind=i4b) :: is, js, ib, jb, nbi, nbj, i3(3), i2(2)
  integer(kind=i4b) :: jm, jl, im, il, jlm, ilm, i, j, nr
  real(kind=r8b)    :: dr(nrmax), maxnorm
  complex(kind=c8b) :: work(nrmax), norm
  integer(kind=i4b) :: imax, jmax

  vsocb = 0.d0
! use the structure of L.S
! also use the fact that the basis was constructed per l channel
  nr = nrpts(ia)
  dr(1:nr) = drproj(1:nr,ia)
  maxnorm = 0.d0
  do j=1,nlmsba(ia)
    i3 = i2lmsb(:,j,ia)
    jb = i3(1); jlm = i3(2); js = i3(3)
    i2 = i2lm(:,jlm)
    jm = i2(1); jl = i2(2)
    do i=1,nlmsba(ia)
      i3 = i2lmsb(:,i,ia)
      ib = i3(1); ilm = i3(2); is = i3(3)
      i2 = i2lm(:,ilm)
      im = i2(1); il = i2(2)
!   ----------------------------------------------------------------
!   selection rules
      if (abs(ldots(i,j,ia)) > atol) then
!     revise the basis
        work(1:nr) = phiref(1:nr,ib,il,is,ia)*vsoc(1:nr)
        work(1:nr) = work(1:nr)*phiref(1:nr,jb,il,js,ia)
        norm = radint(nr,work,dr)
        if (abs(norm) > maxnorm) then
          maxnorm = abs(norm)
          imax = i; jmax = j
        end if
        vsocb(i,j) = vsocb(i,j) + norm*ldots(i,j,ia)
      end if
!   ----------------------------------------------------------------
    end do
  end do
  write(iodb,'(" SOC kernel, ia=",i4)') ia
  write(iodb,'("vsocb norm:",6i4,2es16.8)') i2lmsb(:,imax,ia), i2lmsb(:,jmax,ia), maxnorm
!  write(iodb,'("  ib ilm  is, jb jlm  js, matrix element")')
  do j=1,nlmsba(ia)
    do i=1,nlmsba(ia)
      if (abs(vsocb(i,j)) > 1.d-8) then
!        write(iodb,'(6i4,2es10.2)') i2lmsb(:,i,ia), i2lmsb(:,j,ia), vsocb(i,j)
      end if
    end do
  end do
! All done
  end subroutine build_vsocb
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine build_vscfb()
! assembles the SCF potential in the projection basis
! converts from cartesian to spin representation

  implicit none

! -----------------------------------------------------------------
  integer(kind=i4b) :: ia, is, js, ib, jb, nbi, nbj, i3(3), i2(2)
  integer(kind=i4b) :: jm, jl, im, il, jlm, ilm, i, j, nr, k
  real(kind=r8b)    :: dr(nrmax), maxnorm, uvec(3)
  complex(kind=c8b) :: work(nrmax), norm
  real(kind=r8b)    :: potcart(nrmax,4)
  complex(kind=c8b) :: potspin(nrmax,nsmax,nsmax)
  integer(kind=i4b) :: imax, jmax

  vlmsb = 0.d0
! use the ASA
! use the fact that the basis was constructed per l channel
! ******************
  do ia=1,nasusc
! ******************
  nr = nrpts(ia)
  dr(1:nr) = drproj(1:nr,ia)
  maxnorm = 0.d0
! ----------------------------------------------------------------------
! potential in cartesian labels
! this is the real space orientation of the potential
!  uvec = (/ 0.d0, 1.d0, 0.d0/)
  potcart = 0.d0
! magnetic
  do k=1,3
    potcart(1:nr,k) = br(1:nr,ia)*magdir(k,ia)
  end do
! charge
!  potcart(1:nr,4) = vr(1:nr,ia) - 2.d0*zat(ia)*rmesh(1:nr,ia)
! potential from cartesian to spin labels
  potspin = 0.d0
  do js=1,2
    do is=1,2
      do k=1,4
        potspin(1:nr,is,js) = potspin(1:nr,is,js) + pc2s(is,js,k)*potcart(1:nr,k)
      end do
    end do
  end do
! ----------------------------------------------------------------------
  do j=1,nlmsba(ia)
    i3 = i2lmsb(:,j,ia)
    jb = i3(1); jlm = i3(2); js = i3(3)
    i2 = i2lm(:,jlm)
    jm = i2(1); jl = i2(2)
    do i=1,nlmsba(ia)
      i3 = i2lmsb(:,i,ia)
      ib = i3(1); ilm = i3(2); is = i3(3)
      i2 = i2lm(:,ilm)
      im = i2(1); il = i2(2)
!   ----------------------------------------------------------------
!   selection rules
      if (il == jl .and. im == jm) then
!     revise the basis
        work(1:nr) = phiref(1:nr,ib,il,is,ia)*potspin(1:nr,is,js)*phiref(1:nr,jb,jl,js,ia)
        norm = radint(nr,work,dr)
        vlmsb(i,j,ia) = norm
        if (abs(norm) > maxnorm) then
          maxnorm = abs(norm)
          imax = i; jmax = j
        end if
      end if
!   ----------------------------------------------------------------
    end do
  end do
  write(iodb,'(" SCF potential, ia=",i4)') ia
  write(iodb,'("vscfb norm:",6i4,2es16.8)') i2lmsb(:,imax,ia), i2lmsb(:,jmax,ia), maxnorm
!  write(iodb,'("  ib ilm  is, jb jlm  js, matrix element")')
!  do k=1,4
!  do jlm=1,lmmax
!  do ilm=1,lmmax
!    do js=1,nsmax
!    do is=1,nsmax
!      do jb=1,iwsusc(i2lm(2,jlm),js,ia)
!      do ib=1,iwsusc(i2lm(2,ilm),is,ia)
!        j = lmsb2i(jb,jlm,js,ia)
!        i = lmsb2i(ib,ilm,is,ia)
!        if (abs(vscfb(i,j,ia)) > 1.d-8) then
!          write(iodb,'(6i4,2es10.2)') i2lmsb(:,i,ia), i2lmsb(:,j,ia), vscfb(i,j,ia)
!        end if
!      end do
!      end do
!    end do
!    end do
!  end do
!  end do
!  end do
! ******************
  end do
! ******************
! All done
  end subroutine build_vscfb
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine regf_const_shift()
! shifts the real part of the GF by a constant

  implicit none

  integer(kind=i4b) :: ne, i, ie, ia, lmax
  integer(kind=i4b) :: i3(3), i2(2), ib, ilm, is, im, il
  real(kind=r8b), allocatable :: dummy(:,:,:)
  complex(kind=c8b) :: ze, gf(nlmsb,nlmsb), zdos(0:nlmax,nsmax), avgdev(0:nlmax,nsmax)
  real(kind=r8b)    :: zere, zeim

! read file with complex dos from impurity program
  open(file='imp.zdos',unit=iofile,status='old')
  read(iofile,*) ne, lmax
  allocate(dummy(2,-1:lmax,nsmax))
  do ia=1,nasusc
    avgdev = 0.d0
    call ratfit_gf(ia,ia,.true.,.true.)
    do ie=1,ne
      zdos = 0.d0
!     total zdos put into il = -1 component
      read(iofile,*) i, zere, zeim, dummy
      do is=1,2
        do il=0,nlmax
          zdos(il,is) = cmplx(dummy(1,il,is),dummy(2,il,is))
        end do
      end do
!      write(*,'(i4,100f8.3)') i, zere, zeim, zdos
      ze = cmplx(zere,zeim)
      call ratval_gf(ia,ia,ze,gf)
!     **** zdos from projected GF ****
      do i=1,nlmsba(ia)
        i3 = i2lmsb(:,i,ia)
        ib = i3(1); ilm = i3(2); is = i3(3)
        i2 = i2lm(:,ilm)
        im = i2(1); il = i2(2)
        zdos(il,is) = zdos(il,is) - gf(i,i)*overlap(i,i,ia)
      end do
!     ********************************
!     running averaged difference
      avgdev = avgdev + zdos(0:nlmax,:)
    end do
    avgdev = avgdev/ne
!   done with the differences; print them
    write(*,'("ia, il, avgdev(is=1:2)")')
    do il=0,nlmax
      write(*,'(2i4,4f16.8)') ia, il, avgdev(il,:)
    end do
!   now implement correction
    do ie=1,nesusc
      do i=1,nlmsba(ia)
        i3 = i2lmsb(:,i,ia)
        ib = i3(1); ilm = i3(2); is = i3(3)
        i2 = i2lm(:,ilm)
        im = i2(1); il = i2(2)
!     the location of the shift is a bit arbitrary
        if (ib == 1) then
!        if (ib <= iwsusc(il,is,ia)) then
!          gfpq(i,i,ia,ie) = gfpq(i,i,ia,ie) + fudge*real(avgdev(il,is))/((2*il+1)*iwsusc(il,is,ia))
          gfpq(i,i,ia,ie) = gfpq(i,i,ia,ie) + fudge*real(avgdev(il,is))/((2*il+1))
        end if
      end do
    end do
  end do
  close(iofile)
  deallocate(dummy)
! All done!
  end subroutine regf_const_shift
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine dyn_susc_real(omega,suscylm,suscy00,suscden,analytic,nonanalytic,enhanced)
! CHANGED: reordered indices in susceptibility
! Dynamical KS susceptibility with real frequency
! Radial integral and spherical harmonic resummation

  implicit none

  real(kind=r8b),    intent(in)  :: omega
  complex(kind=c8b), intent(out) :: suscylm(4,4,lmmax0,lmmax0,nasusc,nasusc)
  complex(kind=c8b), intent(out) :: suscy00(4,4,lmmax,lmmax,nasusc,nasusc)
  complex(kind=c8b), intent(out) :: suscden(4,lmmax,nasusc)
  logical,           intent(in)  :: analytic, nonanalytic, enhanced
! -----------------------------------------------------------------
!   i 2pi
  real(kind=r8b),    parameter :: twopi = 8.d0*atan(1.d0)
  complex(kind=c8b), parameter :: i2pi = (0.d0,twopi)
  complex(kind=c8b), parameter :: czero = (0.d0,0.d0), cminus = (-1.d0,0.d0)
! -----------------------------------------------------------------
  complex(kind=c8b) :: gfijw(nlmsb,nlmsb), gfjiw(nlmsb,nlmsb)
  complex(kind=c8b) :: gfij0(nlmsb,nlmsb), gfji0(nlmsb,nlmsb)
  complex(kind=c8b) :: norm, de, e
  integer(kind=i4b) :: q1, lm1, l1, m1, s1, p1
  integer(kind=i4b) :: q2, lm2, l2, m2, s2, p2
  integer(kind=i4b) :: q3, lm3, l3, m3, s3, p3
  integer(kind=i4b) :: q4, lm4, l4, m4, s4, p4
  integer(kind=i4b) :: i2(2), i3(3), ie, ia, ja, jlm, ilm, j, i, is2i(2,2), iq, jq
  integer(kind=i4b) :: ipiv(nalmsb), info
  real(kind=r8b)    :: gaunti, gauntj, doublegaunt, maxelem, efermi
  integer(kind=i4b) :: ne
!  real(kind=r8b)    :: r(nalmsb), c(nalmsb), rwork(2*nalmsb), rcond, ferr(nalmsb), berr(nalmsb)
!  complex(kind=c8b) :: work(2*nalmsb), temp(nalmsb,nalmsb)
!  complex(kind=c8b) :: x(nalmsb,nalmsb)
!  character*1       :: equed

  suscylm = 0.d0; suscy00 = 0.d0; maxelem = 0.d0; suscden = 0.d0
  is2i(1,1) = 4; is2i(2,1) = 1; is2i(1,2) = 2; is2i(2,2) = 3
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  do ja=1,nasusc    ! atom j
  do ia=1,nasusc    ! atom i
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    kssusc = 0.d0
!   ----------------------------------------------------------------
    if (analytic) then
!   Integration on usual box contour
    do ie=1,nescf     ! energy
!    write(*,*) "ie", ie
!   ----------------------------------------------------------------
      e  = escf(ie)
      de = descf(ie)
!     get the projected GFs
      call ratval_gf(ia,ja,e,gfij0)
!      gfij0 = 0.5d0*(gfij0 + transpose(gfij0))
      gfij0 = conjg(transpose(gfij0))
      call ratval_gf(ia,ja,e,gfji0)
!      gfji0 = 0.5d0*(gfji0 + transpose(gfji0))
      call ratval_gf(ia,ja,e+omega,gfijw)
!      gfijw = 0.5d0*(gfijw + transpose(gfijw))
      call ratval_gf(ia,ja,e-omega,gfjiw)
!      gfjiw = 0.5d0*(gfjiw + transpose(gfjiw))
      gfjiw = conjg(transpose(gfjiw))
      do q4=1,nlmsba(ja)
      do q3=1,nlmsba(ja)
        jq = almsb2i(q3,q4,ja)
        do q2=1,nlmsba(ia)
        do q1=1,nlmsba(ia)
          iq = almsb2i(q1,q2,ia)
          kssusc(iq,jq) = kssusc(iq,jq) + (conjg(de)*gfij0(q2,q3)*gfjiw(q4,q1) - de*gfijw(q2,q3)*gfji0(q4,q1))/i2pi
        end do
        end do
      end do
      end do
!   ----------------------------------------------------------------
    end do            ! energy
    end if
!   ----------------------------------------------------------------
!   ----------------------------------------------------------------
    if (nonanalytic) then
!   Integration along real axis
    efermi = real(escf(nescf))
    ne = max(int(abs(omega)/domega),11)
    do ie=1,ne        ! energy
!    write(*,*) "ie", ie
!     --------------------------------------------------------------
      e  = efermi + omega*((ie - 1)/(ne - 1.d0) - 1.d0)
      de = omega/(ne - 1) 
      if (ie == 1 .or. ie == ne) de = 0.5d0*de
!     get the projected GFs
      call ratval_gf(ia,ja,e,gfji0)
      gfji0 = conjg(transpose(gfji0))
!      call ratval_gf(ja,ia,e,gfji0)
      call ratval_gf(ia,ja,e+omega,gfijw)
!      call ratval_gf(ja,ia,e+omega,gfjiw)
!      gfjiw = conjg(transpose(gfjiw))
      do q4=1,nlmsba(ja)
      do q3=1,nlmsba(ja)
        jq = almsb2i(q3,q4,ja)
        do q2=1,nlmsba(ia)
        do q1=1,nlmsba(ia)
          iq = almsb2i(q1,q2,ia)
          kssusc(iq,jq) = kssusc(iq,jq) + de*gfijw(q2,q3)*gfji0(q4,q1)/i2pi
        end do
        end do
      end do
      end do
!   ----------------------------------------------------------------
    end do            ! energy
    end if
!   ----------------------------------------------------------------
    maxelem = maxval(abs(kssusc))
!    where (abs(kssusc) < susctol*maxelem) kssusc = 0.d0
    do q4=1,nlmsba(ja)
    i3 = i2lmsb(:,q4,ja)
    p4 = i3(1); lm4 = i3(2); s4 = i3(3)
    do q3=1,nlmsba(ja)
    i3 = i2lmsb(:,q3,ja)
    p3 = i3(1); lm3 = i3(2); s3 = i3(3)
    jq = almsb2i(q3,q4,ja)
      do q2=1,nlmsba(ia)
      i3 = i2lmsb(:,q2,ia)
      p2 = i3(1); lm2 = i3(2); s2 = i3(3)
      do q1=1,nlmsba(ia)
      i3 = i2lmsb(:,q1,ia)
      p1 = i3(1); lm1 = i3(2); s1 = i3(3)
      iq = almsb2i(q1,q2,ia)
      if (lm1 == lm2 .and. lm3 == lm4) then  ! density = susceptibility*potential
        if (cartesian) then
          do i=1,4
            suscden(i,lm1,ia) = suscden(i,lm1,ia) + ds2c(i,s1,s2)*overlap(q1,q2,ia)*kssusc(iq,jq)*vlmsb(q3,q4,ja)
          end do
        else
          i = is2i(s1,s2)
            suscden(i,lm1,ia) = suscden(i,lm1,ia) + overlap(q1,q2,ia)*kssusc(iq,jq)*vlmsb(q3,q4,ja)
        end if
      end if  ! density = susceptibility*potential
!      norm = 1.d0!overlap(q1,q2,ia)*overlap(q3,q4,ja)
!      if (abs(norm) > gstol) then
!        kssusc(iq,jq) = kssusc(iq,jq)*norm
!        norm = kssusc(iq,jq)
!        if (abs(norm) > maxelem) maxelem = abs(norm)
!      else
!        kssusc(iq,jq) = 0.d0
!      end if
      end do
      end do
    end do
    end do
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  end do            ! atom i
  end do            ! atom j
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ------------------------------------------------------------------
  if (enhanced) then
    call zgemm('N','N',nalmsb,nalmsb,nalmsb,cminus,kssusc,nalmsb,kxcsusc,nalmsb,czero,denominator,nalmsb)
    do iq=1,nalmsb
      denominator(iq,iq) = denominator(iq,iq) + 1.d0
    end do
    call zgesv(nalmsb,nalmsb,denominator,nalmsb,ipiv,kssusc,nalmsb,info)
    if (info /= 0) stop 'dyn_susc_real2: failure in zgesv'
!    call zgesvx('N','N',nalmsb,nalmsb,denominator,nalmsb,temp,nalmsb,ipiv,equed,r,c,kssusc,nalmsb,x,nalmsb,rcond,ferr,berr,work,rwork,info)
!    write(iodb,'("condition number of enhancement factor=",es16.3)') rcond
!    write(iodb,'("fwd & bkwd error in solution=",2es16.3)') maxval(abs(ferr)), maxval(abs(berr))
!    if (info /= 0) stop 'dyn_susc_real2: failure in zgesvx'
!    kssusc = x
  end if
! ------------------------------------------------------------------
! Symmetrize
!  call symmetrize(nalmsb,kssusc,susctol)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  do ja=1,nasusc    ! atom j
  do ia=1,nasusc    ! atom i
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   4444444444444444444444444444444444444444444444444444444444444444
    do q4=1,nlmsba(ja)  ! first set of labels of Gji
!   4444444444444444444444444444444444444444444444444444444444444444
    i3 = i2lmsb(:,q4,ja)
    p4 = i3(1); lm4 = i3(2); s4 = i3(3)
!   3333333333333333333333333333333333333333333333333333333333333333
    do q3=1,nlmsba(ja)  ! second set of labels of Gij
!   3333333333333333333333333333333333333333333333333333333333333333
    i3 = i2lmsb(:,q3,ja)
    p3 = i3(1); lm3 = i3(2); s3 = i3(3)
    jq = almsb2i(q3,q4,ja)
!     22222222222222222222222222222222222222222222222222222222222222
      do q2=1,nlmsba(ia)  ! first set of labels of Gij
!     22222222222222222222222222222222222222222222222222222222222222
      i3 = i2lmsb(:,q2,ia)
      p2 = i3(1); lm2 = i3(2); s2 = i3(3)
!     11111111111111111111111111111111111111111111111111111111111111
      do q1=1,nlmsba(ia)  ! second set of labels of Gji
!     11111111111111111111111111111111111111111111111111111111111111
      i3 = i2lmsb(:,q1,ia)
      p1 = i3(1); lm1 = i3(2); s1 = i3(3)
      iq = almsb2i(q1,q2,ia)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      norm  = kssusc(iq,jq)!*overlap(q1,q2,ia)*overlap(q3,q4,ja)
      doublegaunt = 0.d0
      if (abs(norm) > susctol*maxelem) then
!       ============================================================
        do jlm=1,lmmax0   ! resum for j
        gauntj = gaunt(lm3,lm4,jlm)
        if (abs(gauntj) > ylmtol) then
!       ============================================================
          do ilm=1,lmmax0   ! resum for i
          gaunti = gaunt(lm1,lm2,ilm)
          if (abs(gaunti) > ylmtol) then
!         ==========================================================
!            write(*,*) "jlm, ilm", jlm, ilm
            doublegaunt = doublegaunt + abs(gaunti)*abs(gauntj)
            if (cartesian) then
              do j=1,4
                do i=1,4
                  suscylm(i,j,ilm,jlm,ia,ja) = suscylm(i,j,ilm,jlm,ia,ja) + norm*gaunti*gauntj*ds2c(i,s1,s2)*pc2s(s3,s4,j)
                end do
              end do
            else
              j = is2i(s3,s4)
                i = is2i(s1,s2)
                  suscylm(i,j,ilm,jlm,ia,ja) = suscylm(i,j,ilm,jlm,ia,ja) + norm*gaunti*gauntj
            end if
!         ==========================================================
          end if
          end do            ! resum for i
!       ============================================================
        end if
        end do            ! resum for j
!       ============================================================
        if (abs(doublegaunt) < ylmtol) kssusc(iq,jq) = 0.d0
!       **** Spherical part ****
        if (lm1 == lm2 .and. lm3 == lm4) then
          if (cartesian) then
            do j=1,4
              do i=1,4
                suscy00(i,j,lm1,lm3,ia,ja) = suscy00(i,j,lm1,lm3,ia,ja) + norm*ds2c(i,s1,s2)*pc2s(s3,s4,j)
              end do
            end do
          else
            j = is2i(s3,s4)
              i = is2i(s1,s2)
                suscy00(i,j,lm1,lm3,ia,ja) = suscy00(i,j,lm1,lm3,ia,ja) + norm
          end if
        else
!       TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST
!          kssusc(iq,jq) = 0.d0
!       TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST
        end if
!       ************************
      end if
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     11111111111111111111111111111111111111111111111111111111111111
      end do            ! second set of labels of Gji
!     11111111111111111111111111111111111111111111111111111111111111
!     22222222222222222222222222222222222222222222222222222222222222
      end do            ! first set of labels of Gij
!     22222222222222222222222222222222222222222222222222222222222222
!   3333333333333333333333333333333333333333333333333333333333333333
    end do            ! second set of labels of Gij
!   3333333333333333333333333333333333333333333333333333333333333333
!   4444444444444444444444444444444444444444444444444444444444444444
    end do            ! first set of labels of Gji
!   4444444444444444444444444444444444444444444444444444444444444444
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  end do            ! atom i
  end do            ! atom j
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Final filtering (global)
!  maxelem = maxval(abs(suscylm))
!  where (abs(suscylm) < susctol*maxelem) suscylm = 0.d0
!  maxelem = maxval(abs(suscy00))
!  where (abs(suscy00) < susctol*maxelem) suscy00 = 0.d0
! All done!
  end subroutine dyn_susc_real
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine build_kxcalda(kxc_in,kxcylm,kxcy00)
! assembles the xc ALDA kernel in the projection basis

  implicit none

  real(kind=r8b),    intent(in)  :: kxc_in(nrmax,nasusc)
  complex(kind=c8b), intent(out) :: kxcylm(4,4,lmmax0,nasusc)
  complex(kind=c8b), intent(out) :: kxcy00(4,4,lmmax,lmmax,nasusc)
! -----------------------------------------------------------------
  real(kind=r8b),    parameter :: fourpi = 16.d0*atan(1.d0)
  integer(kind=i4b) :: q1, lm1, l1, m1, s1, p1
  integer(kind=i4b) :: q2, lm2, l2, m2, s2, p2
  integer(kind=i4b) :: q3, lm3, l3, m3, s3, p3
  integer(kind=i4b) :: q4, lm4, l4, m4, s4, p4
  integer(kind=i4b) :: i2(2), i3(3), ie, ia, ja, jlm, ilm, j, i, nr, is2i(2,2), iq, jq
  real(kind=r8b)    :: dr(nrmax), maxnorm, uveci(4), uvecj(4), cart(4,4), dummy(5)
  real(kind=r8b)    :: kxc(nrmax), gaunti, gauntj, maxelem, doublegaunt
  complex(kind=c8b) :: work(nrmax), norm, spin(2,2,2,2)
  real(kind=r8b)    :: kercart(nrmax,4,4)
  complex(kind=c8b) :: kerspin(nrmax,2,2,2,2)
  integer(kind=i4b) :: i1max, i2max, i3max, i4max

  kxcylm = 0.d0; kxcy00 = 0.d0
  is2i(1,1) = 4; is2i(2,1) = 1; is2i(1,2) = 2; is2i(2,2) = 3
! use the ASA
! use the fact that the basis was constructed per l channel
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  do ia=1,nasusc    ! atom i
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    nr = nrpts(ia)
    dr(1:nr) = drproj(1:nr,ia)
    maxnorm = 0.d0
!    uveci(1:3) = magdir(:,ia); uveci(4) = 0.d0
!    uvecj(1:3) = magdir(:,ia); uvecj(4) = 0.d0
    uveci = (/0.d0,0.d0,1.d0,0.d0/); uvecj = uveci
! for testing: cartesian structure of the transverse kernel
    cart = 0.d0
    do j=1,3
      cart(j,j) = 1.d0
    end do
    do j=1,4
      do i=1,4
        cart(i,j) = cart(i,j) - uveci(i)*uvecj(j)
      end do
    end do
    write(iodb,'("cart")')
    write(iodb,'(4f8.4)') cart
!   now convert to the spin representation
    spin = 0.d0
    do s4=1,nsmax
    do s3=1,nsmax
      do s2=1,nsmax
      do s1=1,nsmax
        do j=1,4
        do i=1,4
          spin(s1,s2,s3,s4) = spin(s1,s2,s3,s4) + pc2s(s1,s2,i)*cart(i,j)*ds2c(j,s3,s4)
        end do
        end do
      end do
      end do
    end do
    end do
    write(iodb,'("spin")')
    write(iodb,'(8f8.4)') ((((spin(s1,s2,s3,s4),s3=1,nsmax),s4=1,nsmax),s1=1,nsmax),s2=1,nsmax)
!   transverse kernel
    if (sum(abs(kxc_in(:,ia))) < 1.d-8) then
      kxc(1:nr) = br(1:nr,ia)/nrv(1:nr,3,ia)
    else
      kxc(1:nr) = kxc_in(1:nr,ia)
    end if
!    open(file="kerxc.scf",unit=iofile)
!    read(iofile,*)
!    do i=1,nr
!      read(iofile,*) dummy(1:5), kxc(i)
!    end do
!    kxc(1:nr) = kxc(1:nr)/sqrt(fourpi)
!    close(iofile)
!    where(rmesh(:,ia) < 0.1d0) kxc = 0.d0
!   kernel in cartesian labels
    kercart = 0.d0
    do j=1,4
      do i=1,4
        kercart(1:nr,i,j) = cart(i,j)*kxc(1:nr)!/fourpi
      end do
    end do
!   kernel from cartesian to spin labels
    do s4=1,nsmax
    do s3=1,nsmax
      do s2=1,nsmax
      do s1=1,nsmax
        kerspin(:,s1,s2,s3,s4) = 0.d0
        do j=1,4
        do i=1,4
          kerspin(1:nr,s1,s2,s3,s4) = kerspin(1:nr,s1,s2,s3,s4) + pc2s(s1,s2,i)*kercart(1:nr,i,j)*ds2c(j,s3,s4)
        end do
        end do
        if (sum(abs(kerspin(1:nr,s1,s2,s3,s4))) > atol) then
          do p4=1,iwsusc(2,s4,ia)
          do p3=1,iwsusc(2,s3,ia)
          do p2=1,iwsusc(2,s2,ia)
          do p1=1,iwsusc(2,s1,ia)
          work(1:nr) = phiref(1:nr,p1,2,s1,ia)*phiref(1:nr,p2,2,s2,ia)*phiref(1:nr,p3,2,s3,ia)*phiref(1:nr,p4,2,s4,ia)
          work(1:nr) = (work(1:nr)*br(1:nr,ia))/nrv(1:nr,3,ia)!*kerspin(1:nr,s1,s2,s3,s4)
          norm = radint(nr,work,dr)
          write(iodb,'(4es16.8)') br(nr,ia), nrv(nr,3,ia)
!          write(iodb,'(4es16.8)') phiref(nr,1,2,s1,ia), phiref(nr,1,2,s2,ia), phiref(nr,1,2,s3,ia), phiref(nr,1,2,s4,ia)
          write(iodb,'("U in eV:", 8i8,2f16.8)') s1, s2, s3, s4, p1, p2, p3, p4, norm*13.61d0*2.d0
          end do
          end do
          end do
          end do
          open(file="kerxc.dat",unit=400)
          do i=1,nr
            write(400,'(6es16.8)') rmesh(i,ia), phiref(i,1,2,1,ia), phiref(i,1,2,2,ia), br(i,ia), nrv(i,3,ia), kxc(i)
          end do
          close(400)
        end if
      end do
      end do
    end do
    end do
!   for my amusement: non-zero sums of gaunt products
    do lm4=5,9
    do lm3=5,9
      do lm2=5,9
      do lm1=5,9
        maxelem = 0.d0
        do ilm=1,lmmax2
          maxelem = maxelem + gaunt(lm1,lm2,ilm)*gaunt(lm3,lm4,ilm)
        end do
        if (abs(maxelem) > ylmtol) write(iodb,'(4i4,f10.3)') lm1, lm2, lm3, lm4, maxelem
      end do
      end do
    end do
    end do
!   444444444444444444444444444444444444444444444444444444444444444444
    do q4=1,nlmsba(ia)  ! first set of labels of Gji
!   444444444444444444444444444444444444444444444444444444444444444444
    i3 = i2lmsb(:,q4,ia)
    p4 = i3(1); lm4 = i3(2); s4 = i3(3)
!   333333333333333333333333333333333333333333333333333333333333333333
    do q3=1,nlmsba(ia)  ! second set of labels of Gij
!   333333333333333333333333333333333333333333333333333333333333333333
    i3 = i2lmsb(:,q3,ia)
    p3 = i3(1); lm3 = i3(2); s3 = i3(3)
    jq = almsb2i(q3,q4,ia)
!     2222222222222222222222222222222222222222222222222222222222222222
      do q2=1,nlmsba(ia)  ! first set of labels of Gij
!     2222222222222222222222222222222222222222222222222222222222222222
      i3 = i2lmsb(:,q2,ia)
      p2 = i3(1); lm2 = i3(2); s2 = i3(3)
!     1111111111111111111111111111111111111111111111111111111111111111
      do q1=1,nlmsba(ia)  ! second set of labels of Gji
!     1111111111111111111111111111111111111111111111111111111111111111
      i3 = i2lmsb(:,q1,ia)
      p1 = i3(1); lm1 = i3(2); s1 = i3(3)
      iq = almsb2i(q1,q2,ia)
      if (abs(spin(s1,s2,s3,s4)) > ylmtol) then
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       matrix elements of four basis functions
        work(1:nr) = phiref(1:nr,p1,i2lm(2,lm1),s1,ia)*phiref(1:nr,p2,i2lm(2,lm2),s2,ia)
        work(1:nr) = work(1:nr)*phiref(1:nr,p3,i2lm(2,lm3),s3,ia)*phiref(1:nr,p4,i2lm(2,lm4),s4,ia)
        work(1:nr) = work(1:nr)*kxc(1:nr)
!       multiply by spin factor
        norm = radint(nr,work,dr)*spin(s1,s2,s3,s4)
!        write(iodb,'(4i8,2es16.8)') q1, q2, q3, q4, norm
        if (abs(norm) > maxnorm) then
          maxnorm = abs(norm)
          i1max = q1; i2max = q2; i3max = q3; i4max = q4
        end if
!        if (abs(norm) > susctol*maxelem) then
!       ==============================================================
        doublegaunt = 0.d0
        do jlm=1,lmmax2   ! resum for j
        gauntj = gaunt(lm3,lm4,jlm)
        if (abs(gauntj) > ylmtol) then
          gaunti = gaunt(lm1,lm2,jlm)
          if (abs(gaunti) > ylmtol) then
!         ============================================================
!            write(*,*) "jlm, ilm", jlm, ilm
            doublegaunt = doublegaunt + gaunti*gauntj
            if (jlm <= lmmax0) then
              if (cartesian) then
                do j=1,4
                  do i=1,4
                    kxcylm(i,j,jlm,ia) = kxcylm(i,j,jlm,ia) + norm*gaunti*gauntj*ps2c(i,s1,s2)*dc2s(s3,s4,j)
                  end do
                end do
              else
                j = is2i(s3,s4)
                  i = is2i(s1,s2)
                    kxcylm(i,j,jlm,ia) = kxcylm(i,j,jlm,ia) + norm*gaunti*gauntj
              end if
            end if
!         ============================================================
          end if
        end if
        end do            ! resum for j
!       ==============================================================
        kxcsusc(iq,jq) = norm*doublegaunt
        if (abs(doublegaunt) < ylmtol) kxcsusc(iq,jq) = 0.d0
!        end if
!       **** Spherical part ****
        if (lm1 == lm2 .and. lm3 == lm4) then
          if (cartesian) then
            do j=1,4
              do i=1,4
                kxcy00(i,j,lm1,lm3,ia) = kxcy00(i,j,lm1,lm3,ia) + norm*ps2c(i,s1,s2)*dc2s(s3,s4,j)
              end do
            end do
          else
            j = is2i(s3,s4)
              i = is2i(s1,s2)
                kxcy00(i,j,lm1,lm3,ia) = kxcy00(i,j,lm1,lm3,ia) + norm
          end if
        else
!       TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST
!          kxcsusc(iq,jq) = 0.d0
!       TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST
        end if
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end if
!     11111111111111111111111111111111111111111111111111111111111111
      end do            ! second set of labels of Gji
!     11111111111111111111111111111111111111111111111111111111111111
!     22222222222222222222222222222222222222222222222222222222222222
      end do            ! first set of labels of Gij
!     22222222222222222222222222222222222222222222222222222222222222
!   3333333333333333333333333333333333333333333333333333333333333333
    end do            ! second set of labels of Gij
!   3333333333333333333333333333333333333333333333333333333333333333
!   4444444444444444444444444444444444444444444444444444444444444444
    end do            ! first set of labels of Gji
!   4444444444444444444444444444444444444444444444444444444444444444
    write(iodb,'(" xc kernel, ia=",i4)') ia
    write(iodb,'("kxc norm:",12i4,2es16.8)') i2lmsb(:,i1max,ia), i2lmsb(:,i2max,ia), i2lmsb(:,i3max,ia), i2lmsb(:,i4max,ia), maxnorm
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  end do            ! atom i
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  kxcsusc = real(kxcsusc)
! All done
  end subroutine build_kxcalda
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine susc_denominator
! Analysis of 1 - KS susc * xc kernel

  implicit none

!  integer(kind=i4b), parameter :: itermax = 1
  complex(kind=c8b), parameter :: czero = (0.d0,0.d0), cminus = (-1.d0,0.d0)
  integer(kind=i4b) :: iq, jq, info, iter
  complex(kind=c8b) :: w(nalmsb), vl, vr, work(2*nalmsb)
  real(kind=r8b)    :: start, finish, rwork(2*nalmsb), rw(nalmsb), minlambda

! KS susc eigendecomposition
  call cpu_time(start)
!  kssusc = 0.5d0*(kssusc + transpose(kssusc))
  denominator = kssusc
  call zgeev('N','N',nalmsb,denominator,nalmsb,w,vl,1,vr,1,work,2*nalmsb,rwork,info)
  if (info /= 0) stop 'susc_denominator: failure in zgeev'
  call cpu_time(finish)
  write(iodb,'(" denominator eigen time=",f10.3," s")') finish - start
!  rw = real(w)
! sort them into descending order (MKL utility)
!  call dlasrt('D',nlmsb*nlmsb,rw,info)
  do iq=1,nalmsb
    if (abs(w(iq)) > atol) write(iodb,'("kssusc eval=",i4,2es16.8)') iq, w(iq)
  end do
! xc kernel eigendecomposition
  call cpu_time(start)
!  kxcsusc = real(kxcsusc)
  denominator = kxcsusc
  call zgeev('N','N',nalmsb,denominator,nalmsb,w,vl,1,vr,1,work,2*nalmsb,rwork,info)
  if (info /= 0) stop 'susc_denominator: failure in zgeev'
  call cpu_time(finish)
  write(*,'(" denominator eigen time=",f10.3," s")') finish - start
!  rw = real(w)
! sort them into descending order (MKL utility)
!  call dlasrt('D',nlmsb*nlmsb,rw,info)
  do iq=1,nalmsb
    if (abs(w(iq)) > atol) write(iodb,'("kxcsusc eval=",i4,2es16.8)') iq, w(iq)
  end do
! enhancement factor
  call zgemm('N','N',nalmsb,nalmsb,nalmsb,cminus,kssusc,nalmsb,kxcsusc,nalmsb,czero,denominator,nalmsb)
  do iq=1,nalmsb
    denominator(iq,iq) = denominator(iq,iq) + 1.d0
  end do
! Eigendecomposition
  call cpu_time(start)
  call zgeev('N','N',nalmsb,denominator,nalmsb,w,vl,1,vr,1,work,2*nalmsb,rwork,info)
  if (info /= 0) stop 'susc_denominator: failure in zgeev'
  call cpu_time(finish)
  write(iodb,'(" denominator eigen time=",f10.3," s")') finish - start
!  rw = real(w)
! sort them into descending order (MKL utility)
!  call dlasrt('D',nlmsb*nlmsb,rw,info)
  minlambda = 1.d0
  do iq=1,nalmsb
    if (abs(w(iq)) < abs(minlambda)) minlambda = w(iq)
    if (abs(w(iq)) < 0.999d0 .or. abs(w(iq)) > 1.001d0) write(iodb,'("denominator eval=",i4,2es16.8)') iq, w(iq)
  end do
  write(*,'("susc_denominator: minlambda=",es8.1)') minlambda
! *******************
  if (isoc == 0) then
! *******************
  do iter=1,itermax
! Correction
  if (minlambda < 0.d0) kxcsusc = kxcsusc*(1.d0 + (1.d0 + lambdamix)*minlambda)
  if (minlambda > 0.d0) kxcsusc = kxcsusc*(1.d0 + (1.d0 - lambdamix)*minlambda)
! enhancement factor
  call zgemm('N','N',nalmsb,nalmsb,nalmsb,cminus,kssusc,nalmsb,kxcsusc,nalmsb,czero,denominator,nalmsb)
  do iq=1,nalmsb
    denominator(iq,iq) = denominator(iq,iq) + 1.d0
  end do
! Eigendecomposition
  call cpu_time(start)
  call zgeev('N','N',nalmsb,denominator,nalmsb,w,vl,1,vr,1,work,2*nalmsb,rwork,info)
  if (info /= 0) stop 'susc_denominator: failure in zgeev'
  call cpu_time(finish)
  write(iodb,'(" denominator eigen time=",f10.3," s")') finish - start
!  rw = real(w)
! sort them into descending order (MKL utility)
!  call dlasrt('D',nlmsb*nlmsb,rw,info)
  minlambda = 1.d0
  do iq=1,nalmsb
    if (abs(w(iq)) < abs(minlambda)) minlambda = w(iq)
    if (abs(w(iq)) < 0.999d0 .or. abs(w(iq)) > 1.001d0) write(iodb,'("denominator eval=",i4,2es16.8)') iq, w(iq)
  end do
  write(*,'("susc_denominator: minlambda=",es8.1)') minlambda
  end do
! ******
  end if
! ******
  end subroutine susc_denominator
!----------------------------------------------------------------------


! ----------------------------------------------------------------------
  subroutine baryweights(n,d,x,w)
! Weights for barycentric interpolation
! Rescaled by the factorial of the interpolation order

  implicit none

! Number of points, order of interpolation
  integer(kind=i4b), intent(in)  :: n, d
! Abcissas
  complex(kind=c8b), intent(in)  :: x(n)
! Weights
  complex(kind=c8b), intent(out) :: w(n)
! -----------------------------------------
  integer(kind=i4b) :: i, j, k, imin, imax, iset(d+1)
  real(kind=r8b)    :: minus, fact, hmin, h
  complex(kind=c8b) :: prod

  if (d > n) stop 'baryweights: d > n'
! factorial of the interpolation order
  fact = 1.d0
  do i=2,d
    fact = fact*i
  end do
  hmin = 1.d99
! smallest step size
  do k=1,n
    do i=k+1,n
      h = abs(x(k) - x(i))
      if (h < hmin) hmin = h
    end do
  end do
  hmin = hmin**d
! loop over interpolation points
  do k=1,n
!   form weights from neighbours
    imin = max(k-d,1)
    imax = min(k,n-d)
!   sum over products
    w(k) = 0.d0
    write(iodb,'("k imin imax:",3i6)') k, imin, imax
    do i=imin,imax
      minus = (-1.d0)**(i-1)
!     product
      prod = 1.d0
      do j=i,i+d
        if (j == k) cycle  ! skip 1/0
        write(iodb,'("k i j:",3i6)') k, i, j
        prod = prod/(x(k) - x(j))
      end do
      w(k) = w(k) + minus*prod
    end do
    w(k) = w(k)*fact*hmin
!    w(k) = (-1.d0)**(k-1)
  end do
!  w(1) = 0.5d0
!  w(n) = 0.5d0*(-1.d0)**(n-1)
! All done!
  end subroutine baryweights
! ----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine baryint(n,x,y,w,x0,y0)
! Rational barycentric interpolation
! Works also with equidistant points, as long as order is kept low

  implicit none

! Number of points
  integer(kind=i4b), intent(in)  :: n
! List of abcissas, function values and interpolation weights
  complex(kind=c8b), intent(in)  :: x(n), y(n), w(n)
! Desired point
  complex(kind=c8b), intent(in)  :: x0
! Result of interpolation
  complex(kind=c8b), intent(out) :: y0
! -----------------------------------------
  real(kind=r8b), parameter :: tol = 1.d-8
! -----------------------------------------
  integer(kind=i4b) :: k
  complex(kind=c8b) :: hstep, num, den

  num = 0.d0; den = 0.d0
  do k=1,n
    hstep = x0 - x(k)
    if (abs(hstep) < tol) then
      y0 = y(k)
      return
    end if
    hstep = w(k)/hstep
    num = num + hstep*y(k)
    den = den + hstep
  end do
  y0 = num/den
! All done!
  end subroutine baryint
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine baryint_gf(n,x,w,ia,ja,x0,gf0,onsite,struct)
! Rational barycentric interpolation
! Works also with equidistant points, as long as order is kept low

  implicit none

! Number of points
  integer(kind=i4b), intent(in)  :: n
! List of abcissas and interpolation weights
  complex(kind=c8b), intent(in)  :: x(n), w(n)
! Block of GF
  integer(kind=i4b), intent(in)  :: ia, ja
! Desired point
  complex(kind=c8b), intent(in)  :: x0
! Result of interpolation
  complex(kind=c8b), intent(out) :: gf0(nlmsb,nlmsb)
! To add onsite or structural parts of GF
  logical,           intent(in)  :: onsite, struct
! -----------------------------------------
  real(kind=r8b), parameter :: tol = 1.d-8
! -----------------------------------------
  integer(kind=i4b) :: k
  complex(kind=c8b) :: hstep, gf(nlmsb,nlmsb), den

  gf0 = 0.d0; den = 0.d0
  do k=1,n
    call projected_gf(k,ia,ja,gf,onsite,struct)
    hstep = x0 - x(k)
    if (abs(hstep) < tol) then
      gf0 = gf
      return
    end if
    hstep = w(k)/hstep
    gf0 = gf0 + hstep*gf
    den = den + hstep
  end do
  gf0 = gf0/den
! All done!
  end subroutine baryint_gf
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine ratfit_gf(ia,ja,onsite,struct)
! Rational function fit
! TEST VERSION

  implicit none

! Block of GF
  integer(kind=i4b), intent(in)  :: ia, ja
! To add onsite or structural parts of GF
  logical,           intent(in)  :: onsite, struct
! -----------------------------------------
  real(kind=r8b),    parameter :: tol = 1.d-8
  integer(kind=i4b), parameter :: itermax = 1
! -----------------------------------------
  integer(kind=i4b) :: ie, i, j, p, iter, bestiter, bestnum, bestden, maxie, maxi, maxj, inum, iden
  complex(kind=c8b) :: gf(nlmsb,nlmsb), gfdata(nesusc,nlmsb,nlmsb)
  complex(kind=c8b) :: phase, y, efac
  complex(kind=c8b) :: a(nesusc,1+numd+dend), b(nesusc), numroots(numd), denroots(dend)
  complex(kind=c8b) :: dev(nesusc), fit(nesusc), bestfit(nesusc)
  real(kind=r8b)    :: weights(nesusc), avgdev, maxdev, mindev, maxij, avgij, rms

  gffit(:,:,:,ia,ja) = 0.d0 ! bug solved,  line added 21.August2012
! Fill in array of data
  do ie=1,nesusc
    call projected_gf(ie,ia,ja,gf,onsite,struct)
!    call symmetrize(nlmsb,gf,gfilter)
    gfdata(ie,:,:) = gf
!    write(*,*) "esusc=", esusc(ie)
  end do
! Get coefficients of fit
  avgij = 0.d0; maxij = -1.d0
  do j=1,nlmsba(ja)
  do i=1,nlmsba(ia)
    avgdev = sum(abs(gfdata(1:nesusc,i,j)))/nesusc
!   **************************
    if (avgdev > gfilter) then
!   **************************
    rms = 1.d99
!    inum=numd
!    iden=dend
!    do inum=numd-5,numd
!    do iden=dend-5,dend
    do inum=numd/2,numd
    do iden=dend/2,dend
!    iden = inum+1
!   --------------------------------------------------------------------
!     first pass: all points weighted in the same way
      weights = 1.d0; avgdev = 0.d0; dev = 1.d0
      do iter=1,itermax
!     fill in design matrix; each equation is weighted
        do ie=1,nesusc
!          weights(ie) = 1.d0  ! unbiased fit
!          weights(ie) = 1.d0/(aimag(esusc(ie)) + 1.d-4)  ! small Im E is more important
          weights(ie) = 1.d0/(abs(esusc(ie) - efermi) + 1.d-4)  ! around Efermi is more important
!          weights(ie) = 1.d0/(minval(abs(esusc(ie) - escf(:))) + 1.d-2)  ! close to contour integration is more important
!          weights(ie) = aimag(esusc(ie))        ! large Im E is more important
!       rhs: weighted and shifted by average deviation
          phase = dev(ie)/(abs(dev(ie)) + tol)
          b(ie)  = weights(ie)*(gfdata(ie,i,j) + avgdev*phase)
!       numerator loop
          efac   = 1.d0
          a(ie,1) = weights(ie)*efac
          do p=1,inum
            efac = efac*(esusc(ie) - eshift)
            a(ie,1+p) = weights(ie)*efac
          end do
!       denominator loop
          efac   = 1.d0
          do p=1,iden
            efac = efac*(esusc(ie) - eshift)
            a(ie,1+inum+p) = -efac*b(ie)  ! weight is already included in b
          end do
        end do
!     fit
        call least_squares(nesusc,1+inum+iden,a,b,maxdev)
        fit = b
!     compute deviations on fitted points
        avgdev = 0.d0!; maxdev = -1.d0
        do ie=1,nesusc
          call ratval(nesusc,eshift,inum,iden,fit,esusc(ie),y)
          dev(ie) = y - gfdata(ie,i,j)
          weights(ie) = weights(ie)*abs(dev(ie))**2
!          if (weights(ie) > maxdev) then
!            maxdev = weights(ie); maxie = ie
!          end if
          avgdev = avgdev + weights(ie)
        end do
        avgdev  = sqrt(avgdev)/nesusc
!        call rat_roots(eshift,inum,iden,fit,numroots,denroots)
!        if (any(aimag(denroots(1:iden)) > 0.d0)) avgdev = 10.d0*avgdev
        if (avgdev < rms) then
          rms = avgdev; bestiter = iter; bestnum = inum; bestden = iden
          bestfit = 0.d0; bestfit(1:1+inum) = fit(1:1+inum); bestfit(2+numd:2+numd+iden) = fit(2+inum:2+inum+iden)
        end if
!        write(iodb,'("ratfit_gf: inum, iden=",2i4," maxdev, avgdev=",2es16.8)') inum, iden, maxdev, avgdev
        avgdev  = 0.d0
        weights = 1.d0
      end do
!   best fit
!      if (mindev > tol) write(iodb,'("ratfit_gf: ia,ja, i,j=",8i4," maxdev=",3i4,es16.8)') ia, ja, i2lmsb(:,i), i2lmsb(:,j), bestiter, bestnum, bestden, mindev
!   --------------------------------------------------------------------
    end do
    end do
    gffit(:,i,j,ia,ja) = bestfit(1:1+numd+dend)
    write(iodb,'("ratfit_gf: ia,ja, i,j=",8i4," num, den, rms=",3i4,2es16.8)') ia, ja, i2lmsb(:,i,ia), i2lmsb(:,j,ja), bestiter, bestnum, bestden, rms
!    if (mindev > maxij) then
!      maxi = i; maxj = j; maxij = mindev
!    end if
!    avgij = avgij + mindev
!   ******
    end if
!   ******
  end do
  end do
!  write(iodb,'("ratfit_gf: ia,ja, i,j=",8i4," maxij=",2es16.8)') ia, ja, i2lmsb(:,maxi,ia), i2lmsb(:,maxj,ja), maxij, avgij
! All done!
  end subroutine ratfit_gf
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine ratval(n,zero,numd,dend,fit,x0,y)
! Rational function evaluation

  implicit none

! Number of points
  integer(kind=i4b), intent(in)  :: n
! Shift of origin
  complex(kind=c8b), intent(in)  :: zero
! Degree of numerator and denominator of rational function in fit
  integer(kind=i4b), intent(in)  :: numd, dend
! Fitting coefficients
  complex(kind=c8b), intent(in)  :: fit(n)
! Desired point
  complex(kind=c8b), intent(in)  :: x0
! Result of fit
  complex(kind=c8b), intent(out) :: y
! -----------------------------------------
  real(kind=r8b), parameter :: tol = 1.d-12
! -----------------------------------------
  integer(kind=i4b) :: i, j, p
  complex(kind=c8b) :: xfac, num, den

! numerator loop
  num  = fit(1)
  xfac = 1.d0
  do p=1,numd
    xfac = xfac*(x0 - zero)
    num  = num + xfac*fit(1+p)
  end do
! denominator loop
  den  = 1.d0
  xfac = 1.d0
  do p=1,dend
    xfac = xfac*(x0 - zero)
    den  = den + xfac*fit(1+numd+p)
  end do
!  if (abs(den) < tol) stop 'ratval: denominator is zero!'
  y = num/(den + tol)
! All done!
  end subroutine ratval
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine rat_roots(zero,numd,dend,fit,numroots,denroots)
! Roots of a rational function's numerator and denominator

  implicit none

! Shift of origin
  complex(kind=c8b), intent(in)  :: zero
! Degree of numerator and denominator of rational function in fit
  integer(kind=i4b), intent(in)  :: numd, dend
! Fitting coefficients
  complex(kind=c8b), intent(in)  :: fit(1+numd+dend)
! Roots of numerator and denominator
  complex(kind=c8b), intent(out) :: numroots(numd), denroots(dend)
! -----------------------------------------
  real(kind=r8b), parameter :: tol = 1.d-12
! -----------------------------------------
  integer(kind=i4b) :: i, info
  complex(kind=c8b) :: w(numd+dend), vl, vr, work(2*(numd+dend)), a(numd+dend,numd+dend)
  real(kind=r8b)    :: rwork(2*(numd+dend))

! numerator loop
  numroots = 0.d0
  if (numd > 0) then
!   fill in upper Hessenber matrix
    a = 0.d0
    a(1,1) = -fit(numd)/fit(1+numd)
    do i=2,numd
      a(i,i-1) = 1.d0
      a(1,i) = -fit(1+numd-i)/fit(1+numd)
    end do
    call zgeev('N','N',numd,a,numd+dend,w,vl,1,vr,1,work,2*(numd+dend),rwork,info)
!    write(iodb,'("num roots: ",2es16.8)') w(1:numd)
    numroots = w(1:numd)
  end if
! denominator loop
  denroots = 0.d0
  if (dend > 0) then
! fill in upper Hessenber matrix
    a = 0.d0
    a(1,1) = -fit(numd+dend)/fit(1+numd+dend)
    do i=2,dend
      a(i,i-1) = 1.d0
      a(1,i) = -fit(1+numd+dend-i)/fit(1+numd+dend)
    end do
    a(1,dend) = -1.d0/fit(1+numd+dend)
    call zgeev('N','N',dend,a,numd+dend,w,vl,1,vr,1,work,2*(numd+dend),rwork,info)
!    write(iodb,'("den roots: ",2es16.8)') w(1:dend)
    denroots = w(1:dend)
  end if
! All done!
  end subroutine rat_roots
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine ratval_gf(ia,ja,e0,gf)
! Rational function evaluation
! TEST VERSION

  implicit none

! Which block of the GF
  integer(kind=i4b), intent(in)  :: ia, ja
! Desired point
  complex(kind=c8b), intent(in)  :: e0
! Result of fit
  complex(kind=c8b), intent(out) :: gf(nlmsb,nlmsb)
! -----------------------------------------
  real(kind=r8b), parameter :: tol = 1.d-16
! -----------------------------------------
  integer(kind=i4b) :: i, j, p
  complex(kind=c8b) :: efac, num, den

! Are we in the upper complex plane?
  if (aimag(e0) < 0.d0) stop 'ratval_gf: Im E < 0'
! Get coefficients
  do j=1,nlmsba(ja)
    do i=1,nlmsba(ia)
!   numerator loop
      num  = gffit(1,i,j,ia,ja)
      efac = 1.d0
      do p=1,numd
        efac = efac*(e0 - eshift)
        num  = num + efac*gffit(1+p,i,j,ia,ja)
      end do
!   denominator loop
      den  = 1.d0
      efac = 1.d0
      do p=1,dend
        efac = efac*(e0 - eshift)
        den  = den + efac*gffit(1+numd+p,i,j,ia,ja)
      end do
!      if (abs(den) < tol) stop 'ratval_gf: denominator is zero!'
      gf(i,j) = num/(den + tol)
    end do
  end do
! All done!
  end subroutine ratval_gf
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine least_squares(m,n,a,b,chi2)
! Solution of Ax = b by SVD: A = U sigma V^dagger, U and V unitary
! A is destroyed, x returned in b
! m >= n

  implicit none

! m data points, n unknowns
  integer(kind=i4b), intent(in)    :: m, n
! design matrix, rhs then x
  complex(kind=c8b), intent(inout) :: a(m,n), b(m)
! approximate chi-square
  real(kind=r8b),    intent(out)   :: chi2
! ------------------------------------------
! tolerance for small singular values
  real(kind=r8b),    parameter :: reltol = 1.d-8, abstol = 1.d-6
  complex(kind=c8b), parameter :: one = (1.d0,0.d0), zero = (0.d0,0.d0), minus = (-1.d0,0.d0)
  complex(kind=c8b) :: v(n,n), udummy, work(2*m+n), bsave(m)
  real(kind=r8b)    :: sigma(n), maxsigma, minsigma, rwork(5*n)
  integer(kind=i4b) :: info, i

  if (m < n) stop 'least_squares: m < n'
! computing the SVD of matrix A
!  write(iodb,'("least_squares: SVD")')
!  do i=1,m
!    write(*,'(100es16.8)') a(i,1:n), b(i)
!  end do
  call zgesvd('O','S',m,n,a,m,sigma,udummy,m,v,n,work,2*m+n,rwork,info)
!  write(*,*) "SVD info=", info
  if (info /= 0) stop 'least_squares: fail in SVD'
!  write(iodb,'("singular values:",100es16.8)') sigma
! backsubstitution
!  write(iodb,'("least_squares: backsubstitution")')
! save rhs and singular values
  bsave = b; rwork(1:n) = sigma
! set small singular values to zero
  maxsigma = maxval(sigma)
  minsigma = minval(sigma)
  where (rwork(1:n) < reltol*maxsigma) rwork = 0.d0
!  write(iodb,'("least_squares: inv cond number=",es10.1)') minsigma/maxsigma
  where (sigma > reltol*maxsigma)
    sigma = 1.d0/sigma
  elsewhere
    sigma = 0.d0
  end where
!  write(iodb,'("singular values:",100es16.8)') sigma
! multiply b by U^dagger, put result in work
  call zgemv('C',m,n,one,a,m,b,1,zero,work,1)
! now rescale with singular values
  work(1:n) = sigma*work(1:n)
! multiply by V, put result in b
  call zgemv('C',n,n,one,v,n,work,1,zero,b,1)
  b(n+1:m) = zero
! Now chi-square residual
! multiply x by V^dagger, put result in work
  call zgemv('N',n,n,one,v,n,b,1,zero,work,1)
! now rescale with singular values
  work(1:n) = rwork(1:n)*work(1:n)
! multiply by U, subtract from b
  call zgemv('N',m,n,one,a,m,work,1,minus,bsave,1)
! final result
  chi2 = dot_product(bsave,bsave)
!  write(iodb,'("least_squares: chi2=",es10.1)') chi2
! All done!
  end subroutine least_squares
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine get_xc(ia,lmmax,gfsum,kxc)
! xc potential and energy for atom ia
! The density matrix and the core densities are multiplied by 4pi
! This was factored out in the vxc_vwn subroutine

  implicit none

! which atom, angular momentum cutoff
  integer(kind=i4b), intent(in)  :: ia, lmmax
! density matrix
  complex(kind=c8b), intent(in)  :: gfsum(nlmsb,nlmsb)
! transverse kernel, spherical average
  real(kind=r8b),    intent(out) :: kxc(nrmax)
! ---------------------------------------------------
  real(kind=r8b), parameter :: fourpi = 16.d0*atan(1.d0)
  integer(kind=i4b) :: nr, ir, i2(2), i3(3), k, ill
  integer(kind=i4b) :: i, ilm, il, im, ib, is
  integer(kind=i4b) :: j, jlm, jl, jm, jb, js
  real(kind=r8b)    :: phi1, phi2, den, maglen, magdir(3), vxc, bxc, exc
  real(kind=r8b)    :: rho(4,nll), vxclm(nrmax,lmmax), bxclm(nrmax,3,lmmax)
  real(kind=r8b)    :: exc2(1), vxc2(2), fpirho(2), qcore, mcore
  complex(kind=c8b) :: work(nrmax)

! loop over radial points
  nr = nrpts(ia); vxclm = 0.d0; bxclm = 0.d0; kxc = 0.d0
! Test
  work = nrc(:,ia)
  qcore = real(radint(nr,work(1:nr),drmesh(1:nr,ia)))
  work = mrc(:,ia)
  mcore = real(radint(nr,work(1:nr),drmesh(1:nr,ia)))
  write(iodb,'("ia=",i4," qcore, mcore=",2f12.8)') ia, qcore, mcore
! Test
  open(file='xctest.dat',unit=iofile)
  do ir=1,nr
!   ++++++++++++++++++++++++++++++++++++++++++
!   loop over components of the density matrix
!   ++++++++++++++++++++++++++++++++++++++++++
    rho = 0.d0
    do j=1,nlmsba(ia)
      i3 = i2lmsb(:,j,ia)
      jb = i3(1); jlm = i3(2); js = i3(3)
      i2 = i2lm(:,jlm)
      jm = i2(1); jl = i2(2)
      phi1 = phiref(ir,jb,jl,js,ia)
      do i=1,nlmsba(ia)
        i3 = i2lmsb(:,i,ia)
        ib = i3(1); ilm = i3(2); is = i3(3)
        i2 = i2lm(:,ilm)
        im = i2(1); il = i2(2)
        phi2 = phiref(ir,ib,il,is,ia)
!       add up to the angular mesh for the density
!       go from spin labels to cartesian labels
        do k=1,4
          rho(k,:) = rho(k,:) + ds2c(k,is,js)*gfsum(i,j)*phi1*phi2*ylm(:,ilm)*ylm(:,jlm)
        end do
      end do
    end do
!    write(iodb,'("get_xc: rho=",4es16.8)') rho(:,1)*rmesh(ir,ia)**2
!   ++++++++++++++++++++++++++++++++++++++++++
!   Add the core charge and magnetization densities
    rho(3,:) = rho(3,:) + mrc(ir,ia)
    rho(4,:) = rho(4,:) + nrc(ir,ia)
!   Divide by r^2; the 4pi was taken care of in vxc_vwn
    rho = rho/rmesh(ir,ia)**2
!   Expand the xc potentials in spherical harmonics
    do ill=1,nll
      den = rho(4,ill)
      maglen = sqrt(dot_product(rho(1:3,ill),rho(1:3,ill)))
      magdir = rho(1:3,ill)/maglen
      fpirho(1) = den
      fpirho(2) = maglen
      call vxc_vwn(den,maglen,vxc,bxc,exc)
!      call vosko(exc2,fpirho,vxc2,1,1)
!      write(iodb,'("get_xc: diff vosko=",6es16.8)') vxc-bxc, vxc+bxc, exc, vxc2, exc2
      vxclm(ir,1:lmmax) = vxclm(ir,1:lmmax) + vxc*wll(ill)*ylm(ill,1:lmmax)
      do k=1,3
        bxclm(ir,k,1:lmmax) = bxclm(ir,k,1:lmmax) + magdir(k)*bxc*wll(ill)*ylm(ill,1:lmmax)
      end do
      kxc(ir) = kxc(ir) + wll(ill)*bxc/maglen
    end do
    kxc(ir) = kxc(ir)/rmesh(ir,ia)**2
    write(iofile,'(100es16.8)') rmesh(ir,ia), vr(ir,ia), br(ir,ia), vxclm(ir,1), bxclm(ir,:,1), kxc(ir)
  end do
  close(iofile)
! All done!
  end subroutine get_xc
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  subroutine vxc_vwn(den,mag,vxc,bxc,exc)
! xc potential from VWN parametrization
! densities with 4pi and without r^2
! check definition of rs

  implicit none

! Charge and magnetization densities
  real(kind=r8b), intent(in)  :: den, mag
! xc-potentials, charge and magnetic
  real(kind=r8b), intent(out) :: vxc, bxc, exc
! Parameters of the VWN parametrization (taken from the vosko.f routine)
  real(kind=r8b), parameter :: tol = 1.d-10
  real(kind=r8b), parameter :: fourpi = 16.d0*atan(1.d0)
  real(kind=r8b), parameter :: ap = 0.0621814d0, xp0 = -0.10498d0, bp = 3.72744d0
  real(kind=r8b), parameter :: cp = 12.9352d0, qp = 6.1519908d0
  real(kind=r8b), parameter :: cp1= 1.2117833d0, cp2= 1.1435257d0, cp3 = -0.031167608d0
  real(kind=r8b), parameter :: af = 0.0310907d0, xf0=-0.32500d0, bf = 7.06042d0
  real(kind=r8b), parameter :: cf = 18.0578d0, qf = 4.7309269d0
  real(kind=r8b), parameter :: cf1 = 2.9847935d0, cf2= 2.7100059d0, cf3= -0.1446006d0
! --------------------------------------------------------------------------------------
  real(kind=r8b) :: rs, zeta
  real(kind=r8b) :: fzeta, dfzeta, cbrt1, cbrt2, zeta3, zeta4, x, x2, x3, x4, xpx, xfx
  real(kind=r8b) :: beta, dbeta, atanp, atanf, ecp, ecf, ec
  real(kind=r8b) :: tp1, tf1, ucp, ucf, uc0, uc10, uc20, duc, duc1, duc2, uc1, uc2
  

! cut the densities for low values (vacuum)  
  if (den < tol) then
    vxc = 0.d0
    bxc = 0.d0
    return
  end if

! density and relative spin polarization
!  rs = (3.d0/(fourpi*den))**(1.d0/3.d0)
  rs = (3.d0/den)**(1.d0/3.d0)
!  rs = (3.d0/den)**(1.d0/3.d0)
  x = sqrt(rs); x2 = x*x; x3 = x*x2; x4 = x2*x2
  zeta = mag/den
  zeta3 = zeta**3
  zeta4 = zeta**4 - 1.d0

! exchange enhancement
  fzeta = ((1.d0 + zeta)**(4.d0/3.d0) + (1.d0 - zeta)**(4.d0/3.d0) - 2.d0)/(2.d0**(4.d0/3.d0) - 2.d0)
  cbrt1 = (1.d0 + zeta)**(1.d0/3.d0)
  cbrt2 = (1.d0 - zeta)**(1.d0/3.d0)
  dfzeta = (4.d0/3.d0)*(cbrt1 - cbrt2)/(2.d0**(4.d0/3.d0) - 2.d0)

! build up the correlation energy
  xpx = x2 + bp*x + cp
  xfx = x2 + bf*x + cf
  beta = 1.d0/(2.74208d0 + 3.182d0*x + 0.09873d0*x2 + 0.18268d0*x3)
  dbeta = -(0.27402d0*x + 0.09873d0 + 1.591d0/x)*beta**2
  atanp = atan(qp/(2.d0*x + bp))
  atanf = atan(qf/(2.d0*x + bf))
  ecp = ap*(log(x2/xpx) + cp1*atanp - cp3*(log((x - xp0)**2/xpx) + cp2*atanp))
  ecf = af*(log(x2/xfx) + cf1*atanf - cf3*(log((x - xf0)**2/xfx) + cf2*atanf))
  ec = ecp + fzeta*(ecf - ecp)*(1.d0 + zeta4*beta)
  exc = ec - 0.9163306d0/rs - 0.2381735d0/rs*fzeta

! build up the correlation potential
  tp1 = (x2 + bp*x)/xpx
  tf1 = (x2 + bf*x)/xfx
  ucp = ecp - ap/3.d0*(1.d0 - tp1 - cp3*(x/(x - xp0) - tp1 - xp0*x/xpx))
  ucf = ecf - af/3.d0*(1.d0 - tf1 - cf3*(x/(x - xf0) - tf1 - xf0*x/xfx))
  uc0 = ucp + (ucf - ucp)*fzeta
  uc10 = uc0 - (ecf - ecp)*(zeta - 1.d0)*dfzeta
  uc20 = uc0 - (ecf - ecp)*(zeta + 1.d0)*dfzeta
  duc = (ucf - ucp)*beta*zeta4*fzeta + (ecf - ecp)*(-rs/3.d0)*dbeta*zeta4*fzeta
  duc1 = duc - (ecf - ecp)*beta*(zeta - 1.d0)*(4.d0*zeta3*fzeta + zeta4*dfzeta)
  duc2 = duc - (ecf - ecp)*beta*(zeta + 1.d0)*(4.d0*zeta3*fzeta + zeta4*dfzeta)
  uc1 = uc10 + duc1
  uc2 = uc20 + duc2

! xc potentials
  vxc = 0.5d0*(uc1 + uc2) - 0.5d0*(cbrt1 + cbrt2)*1.221774d0/rs
  bxc = 0.5d0*(uc1 - uc2) - 0.5d0*(cbrt1 - cbrt2)*1.221774d0/rs

! All done!
  end subroutine vxc_vwn
!----------------------------------------------------------------------


  SUBROUTINE VOSKO(EXC,FPIRHO,VXC,IJEND,IJD)
!-----------------------------------------------------------------------
! calculate the spin-polarized exchange-correlation potential
! and the spin-polarized exchange-correlation energy from
! ceperley-alder ( parametrization of vosko, wilk and nusair )
!                                        ( m. manninen )
! use as input the density generated on an angular mesh (see
! subroutine vxclm) . fpirho(.,1) contains the charge density
! times 4 pi and fpirho(.,2) the spin density times 4 pi .
! then the ex.-cor. potential and the ex.-cor. energy on those
! mesh points is calculated .
! the spin-down potential is stored in vxc(.,1) .
!
!                              b.drittler    june 1987
!-----------------------------------------------------------------------
! .. Scalar Arguments ..
  INTEGER IJD,IJEND
! ..
! .. Array Arguments ..
  REAL*8 EXC(*),FPIRHO(IJD,2),VXC(IJD,2)
! ..
! .. Local Scalars ..
  REAL*8 AF,AP,ATNF,ATNP,BETA,BF,BP,CBRT1,CBRT2,CF,CF1, &
                   CF2,CF3,CP,CP1,CP2,CP3,DBETA,DFS,DUC,DUC1,DUC2, &
                   EC,ECF,ECP,FS,ONTHRD,QF,QP,RS,S,S4,SMAG,TF1,TP1, &
                   UC0,UC1,UC10,UC2,UC20,UCF,UCP,X,XF0,XFX,XP0,XPX
  INTEGER IJ
! ..
! .. Intrinsic Functions ..
  INTRINSIC ABS,DATAN,LOG,SQRT
! ..
! .. Save statement ..
  SAVE AP,XP0,BP,CP,QP,CP1,CP2,CP3,AF,XF0,BF,CF,QF,CF1,CF2,CF3
! ..
! .. Data statements ..
  DATA AP,XP0,BP,CP,QP,CP1,CP2,CP3/0.0621814D0,-0.10498D0,3.72744D0, &
       12.9352D0,6.1519908D0,1.2117833D0,1.1435257D0,-0.031167608D0/
  DATA AF,XF0,BF,CF,QF,CF1,CF2,CF3/0.0310907D0,-0.32500D0,7.06042D0, &
       18.0578D0,4.7309269D0,2.9847935D0,2.7100059D0,-0.1446006D0/
! ..
!
  ONTHRD = 1.0D0/3.0D0
!
!---> loop over the angular mesh points
!
  DO IJ = 1,IJEND
    FPIRHO(IJ,1) = MAX(1.0D-10,FPIRHO(IJ,1))
    SMAG = SIGN(1.0D0,FPIRHO(IJ,2))
    FPIRHO(IJ,2) = SMAG*MIN(FPIRHO(IJ,1)-1.0D-10,ABS(FPIRHO(IJ,2)))
    RS = (3.D0/FPIRHO(IJ,1))**ONTHRD
    S = FPIRHO(IJ,2)/FPIRHO(IJ,1)
    X = SQRT(RS)
    XPX = X*X + BP*X + CP
    XFX = X*X + BF*X + CF
    S4 = S**4 - 1.D0
    CBRT1 = (1.D0+S)** (1.D0/3.D0)
    CBRT2 = (1.D0-S)** (1.D0/3.D0)
    FS = ((1.D0+S)** (4.D0/3.D0)+ (1.D0-S)** (4.D0/3.D0)-2.D0)/  &
         (2.D0** (4.D0/3.D0)-2.D0)
    BETA = 1.D0/ (2.74208D0+3.182D0*X+0.09873D0*X*X+0.18268D0*X**3)
    DFS = 4.D0/3.D0* (CBRT1-CBRT2)/ (2.D0** (4.D0/3.D0)-2.D0)
    DBETA = - (0.27402D0*X+0.09873D0+1.591D0/X)*BETA**2
    ATNP = DATAN(QP/ (2.D0*X+BP))
    ATNF = DATAN(QF/ (2.D0*X+BF))
    ECP = AP* (LOG(X*X/XPX)+CP1*ATNP-  &
          CP3* (LOG((X-XP0)**2/XPX)+CP2*ATNP))
    ECF = AF* (LOG(X*X/XFX)+CF1*ATNF-  &
          CF3* (LOG((X-XF0)**2/XFX)+CF2*ATNF))
    EC = ECP + FS* (ECF-ECP)* (1.D0+S4*BETA)
!
!---> calculate ex.-cor. energy
!
    EXC(IJ) = EC - 0.9163306D0/RS - 0.2381735D0/RS*FS
    TP1 = (X*X+BP*X)/XPX
    TF1 = (X*X+BF*X)/XFX
    UCP = ECP - AP/3.D0* (1.D0-TP1-CP3* (X/ (X-XP0)-TP1-XP0*X/XPX))
    UCF = ECF - AF/3.D0* (1.D0-TF1-CF3* (X/ (X-XF0)-TF1-XF0*X/XFX))
    UC0 = UCP + (UCF-UCP)*FS
    UC10 = UC0 - (ECF-ECP)* (S-1.D0)*DFS
    UC20 = UC0 - (ECF-ECP)* (S+1.D0)*DFS
    DUC = (UCF-UCP)*BETA*S4*FS + (ECF-ECP)* (-RS/3.D0)*DBETA*S4*FS
    DUC1 = DUC - (ECF-ECP)*BETA* (S-1.D0)* (4.D0*S**3*FS+S4*DFS)
    DUC2 = DUC - (ECF-ECP)*BETA* (S+1.D0)* (4.D0*S**3*FS+S4*DFS)
    UC1 = UC10 + DUC1
    UC2 = UC20 + DUC2
!
!---> calculate exc.-cor. potential
!
    VXC(IJ,2) = UC1 - 1.221774D0/RS*CBRT1
    VXC(IJ,1) = UC2 - 1.221774D0/RS*CBRT2
    IF (ABS(FPIRHO(IJ,1)).LE.1.0D-10) THEN
       VXC(IJ,1) = 0.0D0
       VXC(IJ,2) = 0.0D0
    END IF
    
  END DO
END SUBROUTINE VOSKO


end module projection
