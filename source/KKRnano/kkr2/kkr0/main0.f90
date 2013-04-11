  program MAIN0

! Explanation of most variables follows below

!     ALAT                     : lattice constant (in a.u.)
!     ABASIS,BBASIS,CBASIS,    : scaling factors for rbasis
!     E1,E2,                   : energies needed in EMESHT
!     HFIELD                   : external magnetic field, for
!                              : initial potential shift in
!                              : spin polarised case
!     TK,                      : temperature
!     VCONST,                  : potential shift
!     BRAVAIS(3,3),            : bravais lattice vectors
!     RECBV(3,3),              : reciprocal basis vectors
!     RMTREF(NREFD),           : muffin-tin radius of reference system
!     RBASIS(3,NAEZD),         : position of atoms in the unit cell
!                              : in units of bravais vectors
!     RCLS(3,NACLSD,NCLSD),    : real space position of atom in cluster
!     RR(3,0:NRD)              : set of real space vectors (in a.u.)
!     VBC(2),                  : potential constants
!     WG(LASSLD),              : integr. weights for Legendre polynomials
!     YRG(LASSLD,0:LASSLD)     : spherical harmonics (GAUNT2)
!     ZAT(NAEZD)               : nuclear charge
!     INTERVX,INTERVY,INTERVZ, : number of intervals in x,y,z-direction
!                              : for k-net in IB of the BZ
!     ICST,                    : number of Born approximation
!     IEND,                    : number of nonzero gaunt coeffizients
!     IFILE,                   : unit specifier for potential card
!     IPE,IPF,IPFE,            : not real used, IPFE should be 0
!     KHFELD,                  : 0,1: no / yes external magnetic field
!     KVREL,                   : 0,1 : non / scalar relat. calculation
!     LMAX,                    : maximum l component in
!                              : wave function expansion
!     LPOT,                    : maximum l component in
!                              : potential expansion
!     NAEZ,                    : number of atoms in unit cell
!     NCLS,                    : number of reference clusters
!     NPNT1,NPNT2,NPNT3,       : number of E points (EMESHT)
!     NPOL,                    : number of Matsubara Pols (EMESHT)
!     NR,                      : number of real space vectors rr
!     NREF,                    : number of diff. ref. potentials
!     NSPIN,                   : counter for spin directions
!     IGUESS                   : 0,1 : no / yes (sc) initial guess, set
!                              : IGUESSD to 1 in inc.p if needed
!     BCP                      : 0,1 : no / yes bc-preconditioning, set
!                              : BCPD to 1 in inc.p if needed
!     QBOUND                   : exit condition for self-consistent
!                              : iteration
!     QMRBOUND                 : exit condition for QMR iterations
!     SCFSTEPS                 : number of scf iterations
!     INIPOL(NAEZD),           : initial spin polarisation

!     ATOM(NACLSD,NAEZD),      : atom at site in cluster
!     CLS(NAEZD),              : cluster around atom
!     NACLS(NCLSD),            : number of atoms in cluster
!     EZOA(NACLSD,NAEZD),      : EZ of atom at site in cluster
!     RMT(NAEZD)               : muffin-tin radius of true system
!     RMTNEW(NAEZD)            : adapted muffin-tin radius
!     RWS(NAEZD)               : Wigner Seitz radius
!     ICST                     : the regular non spherical wavefunctions, the
!                              : alpha matrix and the t-matrix in the ICST-th. born approx
!     IMT(NAEZD),              : r point at MT radius
!     IPAN(NAEZD),             : number of panels in non-MT-region
!     IRC(NAEZD),              : r point for potential cutting
!     IRCUT(0:IPAND,NAEZD),    : r points of panel borders
!     IRMIN(NAEZD),            : max r for spherical treatment
!     IRNS(NAEZD)              : number r points for non spher. treatm.
!     IRWS(NAEZD),             : r point at WS radius
!     LMSP(NAEZD,LMXSPD)       : 0,1 : non/-vanishing lm=(l,m) component
!                              : of non-spherical potential
!     LLMSP(NAEZD,NFUND)       : lm=(l,m) of 'nfund'th nonvanishing
!                              : component of non-spherical pot.
!     NTCELL(NAEZD),           : index for WS cell
!     REFPOT(NAEZD)            : ref. pot. card  at position

!     A(NAEZD),B(NAEZD)        : contants for exponential r mesh
!     R(IRMD,NAEZD)            : radial mesh ( in units a Bohr)
!     DRDI(IRMD,NAEZD)         : derivative dr/di
!     THETAS(IRID,NFUND,NCELLD): shape function
!                              :         ( 0 outer space
!                              : THETA = (
!                              :         ( 1 inside WS cell
!                              : in spherical harmonics expansion
!     LCORE(20,NPOTD)          : angular momentum of core states
!     NCORE(NPOTD)             : number of core states
! ----------------------------------------------------------------------

    use Config_Reader
    use InputParamsNew_mod
    use DimParams_mod
    use Main2Arrays_mod

    implicit none

!   double precision KIND constant
    integer, parameter :: DP = kind(1.0D0)

!     .. Parameters ..

    integer :: NSYMAXD
    integer :: MAXMSHD
!
!     Maximal number of Brillouin zone symmetries, 48 is largest
!     possible number
    parameter (NSYMAXD=48)
!     Maximal number of k-meshes used
    parameter (MAXMSHD=8)
!!     ..
!!     .. Local Scalars ..
!
!    double precision :: RMAX
!    double precision :: GMAX
!    double precision :: RCUTJIJ
!    double precision :: RCUTTRC
!
!!     .. KKR calculation options ..
    double precision :: VCONST
!
!    integer :: ICST
!    integer :: IRM    !     eq IRMD ?  totally USELESS
!    integer :: KHFELD
!    integer :: KPRE
!    integer :: KTE
!    integer :: KVMAD
!    integer :: KVREL
!    integer :: KXC
!    integer :: LMAX
!    integer :: NSPIN
!    integer :: KFORCE
!
!    logical :: JIJ
!    logical :: LDAU
!
!!     .. KKR calculation options options derived from others ..
!    integer :: LMPOT
!    integer :: LPOT
!    integer :: NSRA
!
!!     ..
!!     .. Local Arrays ..
    double precision :: VBC(2)
!
!!     .. Energy Mesh ..
!    double precision :: E1
    double precision :: E2
    double precision :: EFERMI
!    double precision :: TK
!
!    integer :: NPNT1
!    integer :: NPNT2
!    integer :: NPNT3
!    integer :: NPOL
    integer :: IELAST
!    integer :: MAXMESH
!
    complex(kind=DP), dimension(:), allocatable :: EZ
    complex(kind=DP), dimension(:), allocatable :: WEZ
    integer, dimension(:), allocatable :: KMESH
!
!!     .. Lattice ..
!    integer :: NAEZ
!    integer :: NR
!
!    double precision :: ALAT
    double precision :: VOLUME0
!
!    double precision :: BRAVAIS(3,3)
    double precision :: RECBV(3,3)
!
!    double precision, dimension(:,:), allocatable :: RBASIS
!    double precision, dimension(:),   allocatable :: RMTREF
!    double precision, dimension(:,:), allocatable :: RR
!    double precision, dimension(:),   allocatable :: ZAT
!
!!     .. Lattice aux. ..
!    double precision :: ABASIS
!    double precision :: BBASIS
!    double precision :: CBASIS
!
    integer, dimension(:),   allocatable :: NTCELL
!
!!     .. reference clusters
!    double precision :: RCUTXY
!    double precision :: RCUTZ
!
!    integer :: NCLS
!    integer :: NREF
!
!    double precision, dimension(:,:,:), allocatable :: RCLS
!
!    integer, dimension(:,:), allocatable :: ATOM
!    integer, dimension(:),   allocatable :: CLS
!    integer, dimension(:,:), allocatable :: EZOA
!    integer, dimension(:),   allocatable :: NACLS
!!     NUMN0(i) gives number of ref. cluster atoms around atom i
!    integer, dimension(:),   allocatable :: NUMN0
!!     INDN0(i,j) gives the atom index (from basis) of cluster atom j
!!     around central atom i
!    integer, dimension(:,:), allocatable :: INDN0
!    integer, dimension(:),   allocatable :: REFPOT
!
!!     .. reference system ..
!    double precision, dimension(:), allocatable :: VREF
!
!!     .. Brillouin zone ..
!!     kpoints in each direction
!    integer :: INTERVX
!    integer :: INTERVY
!    integer :: INTERVZ
!!     number of symmetries = number of symmetry matrices
!    integer :: NSYMAT
!
!!         symmetry matrices
!    complex(kind=DP), dimension(:,:,:), allocatable :: DSYMLL
!    integer, dimension(:), allocatable :: ISYMINDEX
!
!!     .. Self-consistency parameters ..
!    double precision :: FCM
!    double precision :: MIXING
!    double precision :: QBOUND
!
!    integer :: SCFSTEPS
!    integer :: IMIX
!!     has no effect, use ITDBRYD in inc.p instead
!    integer :: ITDBRY
!
!!     .. Iterative solver options ..
!    double precision :: QMRBOUND
!
    integer :: IGUESS
    integer :: BCP
!
!!     .. auxillary variables, not passed to kkr2
    double precision :: PI
!    double precision :: EREF
    double precision :: HFIELD
!
!    integer :: I1
    integer :: IE
!
    ! error code for allocations
    integer :: ierror
!
!    logical :: EVREF
!    logical :: LCARTESIAN
!
    character(len=40) POTENTIAL_FILENAME
    character(len=40) SHAPEFUN_FILENAME
!
    complex(kind=DP), dimension(:), allocatable :: DEZ ! needed for EMESHT
!
!    integer ::   LMMAXD,LMPOTD,LMXSPD

! ------------- parameters derived from others or calculated
!    integer ::   LPOTD
    integer ::   IEMXD
     integer, parameter :: KREL = 0

! ------------- parameters from global.conf (former inc.p, inc.cls)

!    integer :: LMAXD
!    integer :: NSPIND
!    integer :: NAEZD
!    integer :: IRNSD
!    integer :: IRMD
!    integer :: NREFD
!    integer :: NRD
!    integer :: IRID
!    integer :: NFUND
!    integer :: NCELLD
!    integer :: NACLSD
!    integer :: NCLSD
!    integer :: IPAND
!    integer :: NXIJD
!    integer :: KPOIBZ
    integer :: EKMD
!    integer :: IGUESSD
!    integer :: BCPD
!    integer :: NMAXD
!    integer :: ISHLD
!    integer :: LLY
!
!    integer :: SMPID
!    integer :: EMPID
!    integer :: NTHRDS
!    integer :: XDIM
!    integer :: YDIM
!    integer :: ZDIM
!    integer :: NATBLD
!    integer :: ITDBRYD
!
!    integer :: num_atom_procs

    type (InputParamsNew) :: input
    type (DimParams)      :: dims
    type (Main2Arrays)    :: arrays

    !-------- unused dummys
    !integer :: TRC

! ------------ end of declarations ---------------------------------

    call createDimParamsFromConf(dims)

    ierror = getInputParamsNewValues("input.conf", input)
    if (ierror /= 0) stop

    if (input%NPOL /= 0) then
      IEMXD = input%NPOL + input%NPNT1 + input%NPNT2 + input%NPNT3
    else
      ! DOS-calculation
      IEMXD = input%NPNT2
    end if

    dims%IEMXD = IEMXD

    ! important: determine IEMXD before creating arrays
    call createMain2Arrays(arrays, dims)

!    LPOTD = 2*LMAXD
!    LMMAXD= (LMAXD+1)**2
!    LMPOTD= (LPOTD+1)**2
!    LMXSPD= (2*LPOTD+1)**2
!    LPOT = LPOTD
!    LMPOT = LMPOTD
!

!-----------------------------------------------------------------------------
! Array allocations BEGIN 1
!-----------------------------------------------------------------------------
!
!    ! Lattice
!    allocate(RBASIS(3,NAEZD), stat=ierror)
!    allocate(RMTREF(NREFD), stat=ierror)
!    allocate(RR(3,0:NRD), stat=ierror)
!    allocate(ZAT(NAEZD), stat=ierror)
!
    allocate(NTCELL(dims%NAEZ), stat=ierror)
!
!    ! Reference Cluster
!    allocate(RCLS(3,NACLSD,NCLSD), stat=ierror)
!    allocate(ATOM(NACLSD,NAEZD), stat=ierror)
!    allocate(CLS(NAEZD), stat=ierror)
!    allocate(EZOA(NACLSD,NAEZD), stat=ierror)
!    allocate(NACLS(NCLSD), stat=ierror)
!    allocate(NUMN0(NAEZD), stat=ierror)
!    allocate(INDN0(NAEZD,NACLSD), stat=ierror)
!    allocate(REFPOT(NAEZD), stat=ierror)
!
!    ! Reference system
!    allocate(VREF(NAEZD), stat=ierror)
!
!    !         symmetry matrices
!    allocate(DSYMLL(LMMAXD,LMMAXD,NSYMAXD), stat=ierror)
!    allocate(ISYMINDEX(NSYMAXD), stat=ierror)

!-----------------------------------------------------------------------------
! Array allocations END 1
!-----------------------------------------------------------------------------

!    call RINPUT99(BRAVAIS,ALAT,RBASIS,ABASIS,BBASIS,CBASIS,CLS,NCLS, &
!                  E1,E2,TK,NPOL,NPNT1,NPNT2,NPNT3, &
!                  SCFSTEPS,IMIX,MIXING,QBOUND,FCM, &
!                  ITDBRY,NTCELL,NAEZ,IRM,ZAT, &
!                  NREF,ICST,IFILE,IPE,IPF,IPFE, &
!                  KHFELD,KPRE,KTE, &
!                  KVMAD,KVREL,KXC,LMAX,LMPOT,LPOT, &
!                  NSPIN,REFPOT, &
!                  INTERVX,INTERVY,INTERVZ, &
!                  HFIELD,VBC,VCONST, &
!                  POTENTIAL_FILENAME,SHAPEFUN_FILENAME, &
!                  RCUTZ,RCUTXY,RCUTJIJ,JIJ,RCUTTRC, &
!                  LDAU, &
!                  RMTREF,KFORCE, &
!                  IGUESS,BCP,QMRBOUND,LCARTESIAN,RMAX,GMAX, &
!                  LMAXD, IRNSD, TRC, LPOTD, NSPIND, &
!                  IRMD, NAEZD)

     call RINPUTNEW99(arrays%RBASIS,arrays%CLS,arrays%NCLS, NTCELL,&
                      arrays%NAEZ,arrays%ZAT,arrays%NREF, &
                      arrays%REFPOT,arrays%RMTREF)

! unnecessary parameters, read in for compatibility: IRM, KHFELD, ...

!     in case of a LDA+U calculation - read file 'ldauinfo'
!     and write 'wldau.unf', if it does not exist already
!    if (LDAU) then
!      call ldauinfo_read(LMAXD, NSPIND, ZAT, NAEZD)
!    end if

!===================================================================
! next short section used to search for file 'VREF'
! if present - DP-value given in that file will be used as EREF
!===================================================================
!    inquire(file='VREF',exist=EVREF)
!    if (EVREF) then
!      open(87,file='VREF',form='formatted')
!      read(87,*) EREF
!      close(87)
!    endif

    !     The default value for the repulsive reference potential is 8
    !     This value is used if file VREF does not exist
    !     (VREF should contain only one double precision value used as EREF)
    !     in future: move as parameter to inputfile
!    do I1 = 1,NAEZD
!      VREF(I1) = 8.D0
!      if (EVREF) VREF(I1) = EREF
!    end do

     arrays%VREF = 8.D0
!===================================================================
!===================================================================

!    NSRA = 1
!    if (KVREL >= 1) NSRA = 2

!    call TESTDIM(NSPIN,NAEZ,LMAX,NREF,LMAXD,NREFD,NSPIND)

    SHAPEFUN_FILENAME = 'shapefun'
    POTENTIAL_FILENAME = 'potential'

    open (19,file=SHAPEFUN_FILENAME,status='old',form='formatted')
    open (13,file=POTENTIAL_FILENAME,status='old',form='formatted')

    HFIELD = 0.0d0
    VCONST = 0.0d0
    VBC = 0.0d0

    ! read starting potential and shapefunctions
    call STARTB1_wrapper(input%alat, 13,6,9,0,0, &
                 HFIELD,VCONST, &
                 dims%LPOT,dims%NSPIND,NTCELL, &
                 EFERMI, VBC, arrays%ZAT, &
                 dims%IPAND, dims%IRID, dims%NFUND, dims%IRMD, dims%NCELLD, &
                 dims%NAEZ, dims%IRNSD)

    close(13)
    close(19)

! ----------------------------------------------------------------------
! update Fermi energy, adjust energy window according to running options

    IELAST = IEMXD

    ! Energy mesh
    allocate(EZ(IEMXD), stat=ierror)
    allocate(WEZ(IEMXD), stat=ierror)

    !   auxillary
    allocate(DEZ(IEMXD), stat=ierror)

! for non-DOS calculation upper energy bound corresponds to Fermi energy
! BAD: input%Emax is changed
    if ( input%NPOL /= 0 ) input%Emax = EFERMI

! --> set up energy contour
    PI = 4.0D0*ATAN(1.0D0)

    call EMESHT(EZ,DEZ,IELAST,input%Emin,input%Emax,EFERMI,input%tempr, &
    input%NPOL,input%NPNT1,input%NPNT2,input%NPNT3,IEMXD)
    do IE = 1,IELAST
      WEZ(IE) = -2.D0/PI*DEZ(IE)
    end do

! ================================================ deal with the lattice
    ! only for informative purposes
    arrays%BRAVAIS(:,1) = input%bravais_a
    arrays%BRAVAIS(:,2) = input%bravais_b
    arrays%BRAVAIS(:,3) = input%bravais_c

    call LATTIX99(input%ALAT,arrays%BRAVAIS,RECBV,VOLUME0, .true.)

! --> now generate the real-space lattice vectors for the
!     cluster generation
!
    call RRGEN(arrays%BRAVAIS,arrays%RR,arrays%NR,arrays%NRD)
    write(*,*)

    call SCALEVEC(arrays%RBASIS,input%basisscale(1), &
                  input%basisscale(2),input%basisscale(3), &
                  arrays%NAEZ,arrays%BRAVAIS,input%CARTESIAN)
! ======================================================================

! ======================================================================

!   initialise arrays for clsgen99 to garbage values
     arrays%EZOA = -1
     arrays%ATOM = -1
     arrays%NACLS = -1
     arrays%CLS = -1
     arrays%NUMN0 = -1
     arrays%INDN0 = -1

    call CLSGEN99(arrays%NAEZ,arrays%RR,arrays%NR,arrays%RBASIS,arrays%CLS,arrays%NACLS,arrays%REFPOT,arrays%ATOM, &
                  arrays%EZOA, &
                  arrays%RCLS,input%rclust,input%rclust, &
                  arrays%NUMN0,arrays%INDN0, &
                  arrays%NRD, arrays%NCLSD, arrays%NACLSD)

! xcpl test dimensions for Jij-calculation ..
!    call CLSJIJ0(NAEZ,RR,NR,RBASIS,RCUTJIJ,JIJ,NRD,NXIJD)
! xcpl .


! ======================================================================
!     setting up kpoints
! ======================================================================
    call BZKINT0(arrays%NAEZ, &
                 arrays%RBASIS,arrays%BRAVAIS,RECBV, &
                 arrays%NSYMAT,arrays%ISYMINDEX, &
                 arrays%DSYMLL, &
                 input%bzdivide(1),input%bzdivide(2),input%bzdivide(3), &
                 IELAST,EZ,arrays%KMESH,arrays%MAXMESH,MAXMSHD, &
                 arrays%LMAXD, IEMXD, KREL, arrays%KPOIBZ, EKMD)
    ! after return from bzkint0, EKMD contains the right value

    dims%EKMD = EKMD

! ======================================================================
! ======================================================================


! ======================================================================
! =        write out information for the other program parts           =
! ======================================================================

    IGUESS = dims%IGUESSD
    BCP = dims%BCPD

!    if (BCPD == 1 .and. NATBLD*XDIM*YDIM*ZDIM /= NAEZD) then
!      write(*,*) "ERROR: When BCPD==1 then NATBLD*XDIM*YDIM*ZDIM has to be equal to NAEZD."
!      stop
!    endif

!     Conversion of RMAX and GMAX to units of ALAT
    input%RMAX = input%RMAX*input%ALAT
    input%GMAX = input%GMAX/input%ALAT

    call TESTDIMLAT(input%ALAT,arrays%BRAVAIS,RECBV,input%RMAX,input%GMAX, &
                    dims%NMAXD, dims%ISHLD)
!

    call writeDimParams(dims)
    ierror = writeInputParamsNewToFile('input.unf', input)

!    call write_dimension_parameters( &
!        LMAXD, &
!        NSPIND, &
!        NAEZD, &
!        IRNSD, &
!        IRMD, &
!        NREFD, &
!        NRD, &
!        IRID, &
!        NFUND, &
!        NCELLD, &
!        NACLSD, &
!        NCLSD, &
!        IPAND, &
!        NXIJD, &
!        KPOIBZ, &
!        IGUESSD, &
!        BCPD, &
!        NMAXD, &
!        ISHLD, &
!        LLY, &
!        SMPID, &
!        EMPID, &
!        NTHRDS, &
!        XDIM, &
!        YDIM, &
!        ZDIM, &
!        NATBLD, &
!        ITDBRYD, &
!        IEMXD, &
!        EKMD, &
!        num_atom_procs)
!
    call WUNFILES_NEW(input%NPOL,input%NPNT1,input%NPNT2, &
    input%NPNT3,IELAST,input%tempr,input%Emin,input%Emax,EZ,WEZ, &
    arrays%BRAVAIS,input%RMAX,input%GMAX, &
    EFERMI, &
    input%SCFSTEPS, &
    input%NSRA,arrays%NREF, &
    arrays%NCLS,input%ICST,input%ALAT,arrays%ZAT, &
    arrays%REFPOT,arrays%RMTREF,arrays%VREF, &
    arrays%ATOM,arrays%CLS,arrays%RCLS,arrays%NACLS, &
    arrays%RBASIS,arrays%RR,arrays%EZOA, &
    arrays%KMESH,arrays%MAXMESH,arrays%NSYMAT, &
    arrays%DSYMLL, &
    input%IMIX,input%MIXING,input%FCM,input%KPRE, &
    input%KTE,input%KXC, &
    input%KFORCE, &
    dims%IEMXD,dims%NAEZ, &
    dims%LMMAXD,dims%NREFD, &
    dims%NACLSD,dims%NCLSD,dims%NRD, &
    dims%NSYMAXD, &
    arrays%NUMN0,arrays%INDN0, &
    IGUESS,BCP,input%QMRBOUND, &
    arrays%NR,input%RCUTJIJ,input%JIJ,input%LDAU,arrays%ISYMINDEX)

! ======================================================================

! deallocations

    ! Energy mesh
    deallocate(EZ, stat=ierror)
    deallocate(WEZ, stat=ierror)
!    deallocate(KMESH, stat=ierror)
!
!    ! Lattice
!    deallocate(RBASIS, stat=ierror)
!    deallocate(RMTREF, stat=ierror)
!    deallocate(RR, stat=ierror)
!    deallocate(ZAT, stat=ierror)
!
    deallocate(NTCELL, stat=ierror)
!
!    ! Reference Cluster
!    deallocate(RCLS, stat=ierror)
!    deallocate(ATOM, stat=ierror)
!    deallocate(CLS, stat=ierror)
!    deallocate(EZOA, stat=ierror)
!    deallocate(NACLS, stat=ierror)
!    deallocate(NUMN0, stat=ierror)
!    deallocate(INDN0, stat=ierror)
!
!    ! Reference system
!    deallocate(REFPOT, stat=ierror)
!    deallocate(VREF, stat=ierror)
!
!    ! Symmetry matrices
!    deallocate(DSYMLL, stat=ierror)
!    deallocate(ISYMINDEX, stat=ierror)
!
    !   auxillary
    deallocate(DEZ, stat=ierror)

  end program

