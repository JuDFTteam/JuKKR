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
!     CLEB(NCLEB,2),           : GAUNT coefficients (GAUNT)
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
!     ICLEB(NCLEB,3),          : pointer array
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
!     JEND(LMPOTD,             : pointer array for icleb()
!     LOFLM(LM2D),             : l of lm=(l,m) (GAUNT)
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

!     the inc.p wrapper module is used to include the inc.p and inc.cls
!     parameters - this works also for free form source files
    use inc_p_wrapper_module
    implicit none

!   double precision KIND constant
    integer, parameter :: DP = kind(1.0D0)

!      INCLUDE 'inc.p'
!      INCLUDE 'inc.cls'
!     .. Parameters ..
    integer ::   LASSLD,LMMAXD,LMPOTD,LMXSPD,LM2D, &
    MAXMSHD,NPOTD,NSYMAXD
    parameter (LMMAXD= (LMAXD+1)**2)
    parameter (NPOTD= NSPIND*NAEZD)
    parameter (LMPOTD= (LPOTD+1)**2)
    parameter (LMXSPD= (2*LPOTD+1)**2)
    parameter (LASSLD=4*LMAXD)
    parameter (LM2D= (2*LMAXD+1)**2)
!     Maximal number of Brillouin zone symmetries, 48 is largest
!     possible number
    parameter (NSYMAXD=48)
!     Maximal number of k-meshes used
    parameter (MAXMSHD=8)

!     inc.p, inc.cls parameters
!     LMAXD
!     NAEZD
!     LPOTD
!     LMAXD
!     NREFD
!     IRID
!     BCPD
!     NACLSD
!     NCLEB
!     IRMD
!     IEXMD
!     NGSHD
!     IGUESSD
!     IPAND
!     ISHLD
!     IRNSD
!     KPOIBZ
!     KREL
!     NFUND
!     NATRCD
!     NCLSD
!     NMAXD
!     NRD
!     NSPIND
!     NUTRCD
!     NXIJD


!     ..
!     .. Local Scalars ..

    double precision :: RMAX
    double precision :: GMAX
    double precision :: RCUTJIJ
    double precision :: RCUTTRC

!     .. KKR calculation options ..
    double precision :: VCONST

    integer :: ICST
!     eq IRMD ?
    integer :: IRM
    integer :: ISHIFT
    integer :: KHFELD
    integer :: KPRE
    integer :: KTE
    integer :: KVMAD
    integer :: KVREL
    integer :: KXC
    integer :: LMAX
    integer :: NSPIN
    integer :: KFORCE

    logical :: JIJ
    logical :: LDAU

!     .. KKR calculation options options derived from others ..
    integer :: LMPOT
    integer :: LPOT
    integer :: NSRA

!     ..
!     .. Local Arrays ..
    double precision :: VBC(2)

!     .. Energy Mesh ..
    double precision :: E1
    double precision :: E2
    double precision :: EFERMI
    double precision :: TK

    integer :: NPNT1
    integer :: NPNT2
    integer :: NPNT3
    integer :: NPOL
    integer :: IELAST
    integer :: MAXMESH

    complex(kind=DP) :: EZ(IEMXD)
    complex(kind=DP) :: WEZ(IEMXD)
    integer :: KMESH(IEMXD)

!     .. Lattice ..
    integer :: NAEZ
    integer :: NR

    double precision :: ALAT
    double precision :: VOLUME0

    double precision :: BRAVAIS(3,3)
    double precision :: RECBV(3,3)

    double precision :: RBASIS(3,NAEZD)
    double precision :: RMT(NAEZD)
    double precision :: RMTREF(NREFD)
    double precision :: RR(3,0:NRD)
    double precision :: RWS(NAEZD)
    double precision :: ZAT(NAEZD)

!     .. Lattice aux. ..
    double precision :: ABASIS
    double precision :: BBASIS
    double precision :: CBASIS


!     .. shape functions ..
    double precision :: THETAS(IRID,NFUND,NCELLD)
    double precision :: GSH(NGSHD)
    integer :: IFUNM(LMXSPD,NAEZD)
    integer :: IPAN(NAEZD)
    integer :: LLMSP(NFUND,NAEZD)
    integer :: LMSP(LMXSPD,NAEZD)
    integer :: ILM(NGSHD,3)
    integer :: NFU(NAEZD)
    integer :: NTCELL(NAEZD)
    integer :: IMAXSH(0:LMPOTD)

!     .. reference clusters
    double precision :: RCUTXY
    double precision :: RCUTZ

    integer :: NCLS
    integer :: NREF

    double precision :: RCLS(3,NACLSD,NCLSD)
    integer :: ATOM(NACLSD,NAEZD)
    integer :: CLS(NAEZD)
    integer :: EZOA(NACLSD,NAEZD)
    integer :: NACLS(NCLSD)
!     NUMN0(i) gives number of ref. cluster atoms around atom i
    integer :: NUMN0(NAEZD)
!     INDN0(i,j) gives the atom index (from basis) of cluster atom j
!     around central atom i
    integer :: INDN0(NAEZD,NACLSD)
    integer :: REFPOT(NAEZD)

!     .. reference system ..
    double precision :: VREF(NAEZD)


!     .. radial mesh(es) ..
    double precision :: A(NAEZD)
    double precision :: B(NAEZD)
    double precision :: DRDI(IRMD,NAEZD)
    double precision :: R(IRMD,NAEZD)
    integer :: IMT(NAEZD)
    integer :: IRC(NAEZD)
    integer :: IRCUT(0:IPAND,NAEZD)
    integer :: IRMIN(NAEZD)
    integer :: IRNS(NAEZD)
    integer :: IRWS(NAEZD)

!     .. Brillouin zone ..
!     kpoints in each direction
    integer :: INTERVX
    integer :: INTERVY
    integer :: INTERVZ
!     number of symmetries = number of symmetry matrices
    integer :: NSYMAT

!         symmetry matrices
    complex(kind=DP) DSYMLL(LMMAXD,LMMAXD,NSYMAXD)
    integer :: ISYMINDEX(NSYMAXD)

!     .. Spherical harmonics ..
    double precision :: WG(LASSLD) ! not passed to kkr2
    double precision :: YRG(LASSLD,0:LASSLD,0:LASSLD) ! not passed to kkr2

!     .. Clebsch-Gordon coefficients and spherical harmonics
    integer :: IEND

    double precision :: CLEB(NCLEB,2)

    integer :: ICLEB(NCLEB,3)
    integer :: JEND(LMPOTD,0:LMAXD,0:LMAXD)
    integer :: LOFLM(LM2D) ! gives l from LM index

!     .. core states ..
    integer :: ITITLE(20,NPOTD)
    integer :: LCORE(20,NPOTD)
    integer :: NCORE(NPOTD)

!     .. Self-consistency parameters ..
    double precision :: FCM
    double precision :: MIXING
    double precision :: QBOUND

    integer :: SCFSTEPS
    integer :: IMIX
!     has no effect, use ITDBRYD in inc.p instead
    integer :: ITDBRY

!     .. Iterative solver options ..
    double precision :: QMRBOUND

    integer :: IGUESS
    integer :: BCP

!     .. auxillary variables, not passed to kkr2
    double precision :: PI
    double precision :: EREF
    double precision :: HFIELD
    double precision :: E2IN

    integer :: I1
    integer :: IPE
    integer :: IPF
    integer :: IPFE
    integer :: IFILE
    integer :: IE

    logical :: EVREF
    logical :: LCARTESIAN

    character(len=40) I13
    character(len=40) I19

!     needed for EMESHT
    complex(kind=DP) DEZ(IEMXD)

    double precision :: RMTNEW(NAEZD)

    integer :: INIPOL(NAEZD)

! ------------ end of declarations ---------------------------------


!===================================================================
! next short section used to search for file 'VREF'
! if present - DP-value given in that file will be used as EREF
!===================================================================
    inquire(file='VREF',exist=EVREF)
    if (EVREF) then
      open(87,file='VREF',form='formatted')
      read(87,*) EREF
      close(87)
    endif

    !     The default value for the repulsive reference potential is 8
    !     This value is used if file VREF does not exist
    !     (VREF should contain only one double precision value used as EREF)
    !     in future: move as parameter to inputfile
    do I1 = 1,NAEZD
      VREF(I1) = 8.D0
      if (EVREF) VREF(I1) = EREF
    end do
!===================================================================
!===================================================================
    PI = 4.0D0*ATAN(1.0D0)
    EFERMI = 0.0d0


    call RINPUT99(BRAVAIS,ALAT,RBASIS,ABASIS,BBASIS,CBASIS,CLS,NCLS, &
                  E1,E2,TK,NPOL,NPNT1,NPNT2,NPNT3, &
                  SCFSTEPS,IMIX,MIXING,QBOUND,FCM, &
                  ITDBRY,IRNS,NTCELL,NAEZ,IRM,ZAT, &
                  NREF,ICST,IFILE,IPE,IPF,IPFE, &
                  KHFELD,KPRE,KTE, &
                  KVMAD,KVREL,KXC,LMAX,LMPOT,LPOT, &
                  NSPIN,REFPOT, &
                  ISHIFT,INTERVX,INTERVY,INTERVZ, &
                  HFIELD,VBC,VCONST,INIPOL, &
                  I13,I19, &
                  RCUTZ,RCUTXY,RCUTJIJ,JIJ,RCUTTRC, &
                  LDAU, &
                  RMTREF,KFORCE, &
                  IGUESS,BCP,QMRBOUND,LCARTESIAN,RMAX,GMAX, &
                  LMAXD, IRNSD, TRC, LPOTD, NSPIND, &
                  IRMD, NAEZD)

!     in case of a LDA+U calculation - read file 'ldauinfo'
!     and write 'wldau.unf', if it does not exist already
    if (LDAU) then
      call ldauinfo_read(LMAXD, NSPIND, ZAT, NAEZD)
    end if


!     IF(NAEZD.GT.100) IPE=0


    NSRA = 1
    if (KVREL >= 1) NSRA = 2

    call TESTDIM(NSPIN,NAEZ,LMAX,IRM,NREF,IRNS, &
    LMAXD,IRMD,IRNSD,NREFD,NSPIND)

    open (19,file=I19,status='old',form='formatted')
    open (IFILE,file=I13,status='old',form='formatted')

    E2IN = E2

    call STARTB1(IFILE,IPF,IPFE,IPE,KHFELD, &
                 1,NAEZ, &
                 RMTNEW,RMT,ITITLE,HFIELD,IMT,IRC,VCONST, &
                 IRNS,LPOT,NSPIN,IRMIN,NTCELL,IRCUT,IPAN, &
                 THETAS,IFUNM,NFU,LLMSP,LMSP,E2IN, &
                 VBC,RWS,LCORE,NCORE,DRDI, &
                 R,ZAT,A,B,IRWS,INIPOL,1, &
!                New after inc.p replace
                 IPAND, IRID, NFUND, IRMD, NCELLD, &
                 NAEZD, IRNSD)

    close(IFILE)
    close(19)

! ----------------------------------------------------------------------
! update Fermi energy, adjust energy window according to running options

    if ( NPOL == 0 ) EFERMI = E2IN
!     a test if E2IN is different from E2 after call to STARTB1,
!     before it was =E2
    if ( DABS(E2IN-E2) > 1D-10 .and. NPOL /= 0 ) E2 = E2IN

! --> set up energy contour

    call EMESHT(EZ,DEZ,IELAST,E1,E2,E2IN,TK, &
    NPOL,NPNT1,NPNT2,NPNT3,IEMXD)
    do IE = 1,IELAST
      WEZ(IE) = -2.D0/PI*DEZ(IE)
    end do


    call GAUNT2(WG,YRG,LMAX)
    call GAUNT(LMAX,LPOT,WG,YRG,CLEB,LOFLM,ICLEB,IEND,JEND,NCLEB)
!      OPEN(56,FILE='gaunt',FORM='formatted')
!      WRITE(56,FMT='(I10)') IEND
!      DO I = 1,IEND
!      WRITE(56,FMT='(3I5,1P,D25.17)')
!     +            ICLEB(I,1),ICLEB(I,2),ICLEB(I,3),CLEB(I,1)
!      END DO
!      CLOSE(56)

! --> setup of GAUNT coefficients C(l,m;l',m';l'',m'') for all
!     nonvanishing (l'',m'')-components of the shape functions THETAS

    call SHAPE(LPOT,NAEZ,GSH,ILM,IMAXSH,LMSP,NTCELL,WG,YRG,LMAX,NGSHD)
!      OPEN(56,FILE='gaunt_shape',FORM='formatted')
!      WRITE(56,FMT='(I10)') IMAXSH(LMPOTD)
!      DO I = 1,IMAXSH(LMPOTD)
!      WRITE(56,FMT='(3I5,1P,D25.17)')
!     +            ILM(I,1),ILM(I,2),ILM(I,3),GSH(I)
!      END DO
!      CLOSE(56)

! ================================================ deal with the lattice
    call LATTIX99(ALAT,BRAVAIS, &
                  RECBV,VOLUME0,RR,NR,NRD)

    call SCALEVEC(RBASIS,ABASIS,BBASIS,CBASIS, &
                  NAEZ,BRAVAIS,LCARTESIAN)
! ======================================================================

! ======================================================================



    call CLSGEN99(NAEZ,RR,NR,RBASIS,CLS,NACLS,REFPOT,ATOM, &
                  EZOA, &
                  RCLS,RCUTZ,RCUTXY, &
                  NUMN0,INDN0, &
                  NRD, NCLSD, NACLSD)

! xcpl test dimensions for Jij-calculation ..
    call CLSJIJ0(NAEZ,RR,NR,RBASIS,RCUTJIJ,JIJ,NRD,NXIJD)
    ! xcpl .

    ! C2 ===================================================================
    ! C2  generate cluster 2 in order to truncate GLLKE, GLLH, etc.
    ! C2 ===================================================================

    if (TRC == 1) then
      call CLSGEN_TRC(NAEZ,RR,NR,RBASIS, &
      RCUTTRC,ALAT, &
      NRD, NATRCD, NUTRCD)
    endif
! C2 ===================================================================
! C2 ===================================================================



! ======================================================================
!     setting up kpoints
! ======================================================================
    call BZKINT0(NAEZ, &
                 RBASIS,BRAVAIS,RECBV, &
                 NSYMAT,ISYMINDEX, &
                 DSYMLL, &
                 INTERVX,INTERVY,INTERVZ, &
                 IELAST,EZ,KMESH,MAXMESH,MAXMSHD, &
                 LMAX, IEMXD, KREL, KPOIBZ, EKMD, IGUESSD)
! ======================================================================
! ======================================================================


! ======================================================================
! =        write out information for the other program parts           =
! ======================================================================

    if (IGUESS > IGUESSD) then
      stop ' ERROR: to activate initial guess, set IGUESSD=1'
    endif
    if (BCP > BCPD) then
      stop ' ERROR: to activate bc-preconditioning, set BCPD=1'
    endif

!     Conversion of RMAX and GMAX to units of ALAT
    RMAX = RMAX*ALAT
    GMAX = GMAX/ALAT

    call TESTDIMLAT(ALAT,BRAVAIS,RECBV,RMAX,GMAX, &
                    NMAXD, ISHLD)


    call WUNFILES(NPOL,NPNT1,NPNT2,NPNT3,IELAST,TK,E1,E2,EZ,WEZ, &
    BRAVAIS,RECBV,VOLUME0,RMAX,GMAX, &
    EFERMI,VBC, &
    SCFSTEPS,LCORE,NCORE, &
    NSRA,NAEZ,NREF,NSPIN,LMAX, &
    NCLS,ICST,IPAN,IRCUT,ALAT,ZAT,R,DRDI, &
    REFPOT,RMTREF,VREF,IEND,JEND,CLEB,ICLEB, &
    ATOM,CLS,RCLS,NACLS,LOFLM, &
    RBASIS,RR,EZOA, &
    KMESH,MAXMESH,NSYMAT, &
    DSYMLL, &
    A,B,IFUNM,ITITLE, &
    LMSP,NTCELL,THETAS, &
    LPOT,LMPOT, &
    IMIX,MIXING,QBOUND,FCM,KPRE,KTE, &
    KVMAD,KXC,ISHIFT,KFORCE, &
    LLMSP,IMT,IRC,IRMIN,IRNS,IRWS,RWS,RMT,NFU, &
    GSH,ILM,IMAXSH, &
    IEMXD,IRMD,LMPOTD,NPOTD,NAEZD, &
    LMMAXD,IPAND,NREFD,LMAXD, &
    NCLEB,NACLSD,NCLSD,LM2D,NRD, &
    NSYMAXD,IRID,NFUND, &
    NCELLD,LMXSPD,NGSHD,NUMN0,INDN0, &
    IGUESS,BCP,QMRBOUND, &
    NR,RCUTJIJ,JIJ,LDAU,ISYMINDEX)
! ======================================================================

  end program

