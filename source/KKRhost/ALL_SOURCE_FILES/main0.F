      module MOD_MAIN0 
      
      use mod_wunfiles
      use mod_types, only: t_imp
#ifdef CPP_TIMING
      use mod_timing
#endif
      use mod_create_newmesh
      
      implicit none
      
      contains
          
      subroutine main0()      
C
C *********************************************************************
C * For KREL = 1 (relativistic mode)                                  *
C *                                                                   *
C *  NPOTD = 2 * NATYPD                                               *
C *  LMMAXD = 2 * (LMAXD+1)^2                                         *
C *  NSPIND = 1                                                       *
C *                                                                   *
C *********************************************************************
C
C Explanation of most variables follows below
C
C     ALAT                     : lattice constant (in a.u.)
C     ABASIS,BBASIS,CBASIS,    : scaling factors for rbasis
C     E1,E2,                   : energies needed in EMESHT
C     HFIELD                   : external magnetic field, for
C                              : initial potential shift in
C                              : spin polarised case
C     TK,                      : temperature
C     VCONST,                  : potential shift
C     BRAVAIS(3,3),            : bravais lattice vectors
C     RECBV(3,3),              : reciprocal basis vectors
C     MTFAC(NATYPD),           : scaling factor for radius MT
C     CLEB(NCLEB,2),           : GAUNT coefficients (GAUNT)
C     RMTREF(NREFD),           : muffin-tin radius of reference system
C     RBASIS(3,NAEZD+NEMBD),   : position of atoms in the unit cell
C                              : in units of bravais vectors
C     RCLS(3,NACLSD,NCLSD),    : real space position of atom in cluster
C     RR(3,0:NRD)              : set of real space vectors (in a.u.)
C     VBC(2),                  : potential constants
C     WG(LASSLD),              : integr. weights for Legendre polynomials
C     YRG(LASSLD,0:LASSLD)     : spherical harmonics (GAUNT2)
C     ZAT(NATYPD)              : nuclear charge
C     INTERVX,INTERVY,INTERVZ, : number of intervals in x,y,z-direction
C                              : for k-net in IB of the BZ
C     ICC,                     : enables the calculation of off-diagonal
C                              : elements of the GF.
C              = 0             : for SCF/DOS etc calc.
C              = 1             : writes out Gij for a cluster; cluster 
C                              : data is expected to be in the file I25
C              = -1            : writes out Gij for customised pairs
C                              : i,j that need to be set up elsewhere
C     NSHELL(0:NSHELD)         : index of atoms/pairs per shell 
C                              : (ij-pairs); nshell(0) = number of shells
C     NSH1(1..NSHELL(0))       : corresponding index of the sites I/J  
C     NSH2(1..NSHELL(0))       : (NSH1/2) in the unit cell in a shell
C                              :
C     IJTABCALC                : linear pointer, specifying whether
C                              : the block (i,j) has to be calculated
C                              : needs set up for ICC=-1, 
C                              : not used for ICC=1
C     IJTABSH                  : linear pointer, assigns pair (i,j) to
C                              : a shell in the array GS(*,*,*,NSHELD)
C     IJTABSYM                 : linear pointer, assigns pair (i,j) to
C                              : the rotation bringing GS into Gij
C     ISH,JSH                  : 
C     NOFGIJ                   : number of GF pairs IJ to be calculated
C                              : as determined from IJTABCALC<>0
C     IOFGIJ,JOFGIJ            : linear pointers, similar to NSH1/NSH2
C                              : but giving the actual index of sites
C                              : I,J = 1,NATOMIMP in the cluster
C     ICST,                    : number of Born approximation
C     IEND,                    : number of nonzero gaunt coeffizients
C     IFILE,                   : unit specifier for potential card
C     INS,                     : 0 (MT), 1(ASA), 2(Full Potential)
C     INSREF,                  : INS for reference pot. (usual 0)
C     IPE,IPF,IPFE,            : not real used, IPFE should be 0
C     KSCOEF,                  : 0,1: read shell structure from
C                              :      file 25
C     KHFELD,                  : 0,1: no / yes external magnetic field
C     KSHAPE,                  : exact treatment of WS cell
C     KVREL,                   : 0,1 : non / scalar relat. calculation
C     KWS,                     : 0 (MT), 1(ASA)
C     LMAX,                    : maximum l component in
C                              : wave function expansion
C     LPOT,                    : maximum l component in
C                              : potential expansion
C     NAEZ,                    : number of atoms in unit cell
C     NATYP,                   : number of kinds of atoms in unit cell
C     NCLS,                    : number of reference clusters
C     NEMB,                    : number of 'embedding' positions
C     NINEQ,                   : number of ineq. positions in  unit cell
C     NLAYER,                  : number of principal layer
C     NPNT1,NPNT2,NPNT3,       : number of E points (EMESHT)
C     NPOL,                    : number of Matsubara Pols (EMESHT)
C     NR,                      : number of real space vectors rr
C     NREF,                    : number of diff. ref. potentials
C     NSPIN,                   : counter for spin directions
C     SCFSTEPS                 : number of scf iterations
C     INIPOL(NATYPD),          : initial spin polarisation
C     IXIPOL(NATYPD),          : constraint of spin pol.
C     KAOEZ(NATYPD,NAEZD+NEMBD): kind of atom at site in elem. cell
C
C     ATOM(NACLSD,NAEZD),      : atom at site in cluster
C     CLS(NAEZD),              : cluster around atomic sites
C     NACLS(NCLSD),            : number of atoms in cluster
C     EZOA(NACLSD,NAEZD),      : EZ of atom at site in cluster
C     ICLEB(NCLEB,4),          : pointer array
C     RMT(NATYPD)              : muffin-tin radius of true system
C     RMTNEW(NATYPD)           : adapted muffin-tin radius 
C     RWS(NATYPD)              : Wigner Seitz radius
C     IMT(NATYPD),             : r point at MT radius
C     IPAN(NATYPD),            : number of panels in non-MT-region
C     IRC(NATYPD),             : r point for potential cutting
C     IRCUT(0:IPAND,NATYPD),   : r points of panel borders
C     IRMIN(NATYPD),           : max r for spherical treatment
C     IRNS(NATYPD)             : number r points for non spher. treatm.
C     IRWS(NATYPD),            : r point at WS radius
!     FPRADIUS(NATYPD)         : r point at which full-potential treatment starts
C     LMSP(NATYPD,LMXSPD)      : 0,1 : non/-vanishing lm=(l,m) component
C                              : of non-spherical potential
C     LLMSP(NATYPD,NFUND)      : lm=(l,m) of 'nfund'th nonvanishing
C                              : component of non-spherical pot.
C     JEND(LMPOTD,             : pointer array for icleb()
C     LOFLM(LM2D),             : l of lm=(l,m) (GAUNT)
C     NTCELL(NATYPD),          : index for WS cell
C     REFPOT(NATYPD+NEMBD)     : ref. pot. card  at position
C     A(NATYPD),B(NATYPD)      : contants for exponential r mesh
C     R(IRMD,NATYPD)           : radial mesh ( in units a Bohr)
C     DRDI(IRMD,NATYPD)        : derivative dr/di
C     THETAS(IRID,NFUND,NCELLD): shape function
C                              :         ( 0 outer space
C                              : THETA = (
C                              :         ( 1 inside WS cell
C                              : in spherical harmonics expansion
C     ECORE(20,NPOTD)          : core energies
C     LCORE(20,NPOTD)          : angular momentum of core states
C     NCORE(NPOTD)             : number of core states
C ----------------------------------------------------------------------
C  nlbasis              : number of basis layers of left 
C                         host (repeated units)
C  nleft                : number of repeated basis for left host 
C                         to get converged  electrostatic potentials
C  tleft(3,nlbasis)     : vectors of the basis for the left host
C  zperleft(3)          : vector to define how to repeat the basis 
C                         of the left host
C
C  nrbasis              : number of basis layers of right 
C                         host (repeated units)
C  nright               : number of repeated basis for right host 
C                         to get converged electrostatic potentials
C  tright(3,nlbasis)    : vectors of the basis for the right host
C  zperight(3)          : vector to define how to repeat the basis 
C                         of the right host
C  cmomhost
C  (lmpotd,1..nlbasis)   : charge moments of each atom of the left host
C  (lmpotd,nlbasis+1,..) : charge moments of each atom of the right host
C ======================================================================
C RELativistic mode 
C ----------------------------------------------------------------------
C Matrices to change from NON-RELATIVISTIC to RELATIVISTIC 
C representations or back
C
C  RREL(LMMAXD,LMMAXD)   : non-relat. REAL spher. harm.  >   (kappa,mue)
C                        : (kappa,mue)  > non-relat. REAL spher. harm.
C  CREL(LMMAXD,LMMAXD)   : non-relat. CMPLX. spher. harm. > (kappa,mue)
C                        : (kappa,mue)  > non-relat. CMPLX. spher. harm.
C  RC(LMMAXD,LMMAXD)     : NREL REAL spher. harm. >  CMPLX. spher. harm.
C                        : NREL CMPLX. spher. harm. > REAL spher. harm.
C to use the above matrices, the non-relat. representations have to 
C include the  spin index (i.e. should have the same dimension 
C LMMAXD = 2*(LMAXD+1)**2 as the relativistic matrices)
C
C alternatively, SRREL can be used combined with the pointer arrays 
C  IRREL/NRREL 
C  SRREL(2,2,LMMAXD)
C  IRREL(2,2,LMMAXD),NRREL(2,LMMAXD)
C     
C ----------------------------------------------------------------------
C  DROTQ(LMMAXD,LMMAXD,NAEZD) : rotation matrices to change between 
C                             : LOCAL/GLOBAL frame of reference for
C                             : magnetisation <> Oz or noncollinearity
C  SYMUNITARY(NSYMAXD)        : unitary/antiunitary symmetry flag
C     
C ----------------------------------------------------------------------
C internally used potential and mesh variables for the relativistic
C routines
C  VTREL(IRMD*KREL+(1-KREL),NATYPD)     : potential (spherical part)
C  BTREL(IRMD*KREL+(1-KREL),NATYPD)     : magnetic field
C  RMREL(IRMD*KREL+(1-KREL),NATYPD)     : radial mesh
C  DRDIREL(IRMD*KREL+(1-KREL),NATYPD)   : derivative of radial mesh
C  R2DRDIREL(IRMD*KREL+(1-KREL),NATYPD) : r**2 * drdi
C  JWSREL(NATYPD)                       : index of the WS radius
C  IRSHIFT(NATYPD)                      : shift of the REL radial mesh 
C                                       : with respect no NREL
C  ZREL(NATYPD)                         : atomic number (cast integer)
C 
C Note: IRMD-dimension depends on KREL to save space
C ======================================================================VINS()

      use mod_types
#ifdef CPP_OMP
      use omp_lib        ! necessary for omp functions
#endif
#ifdef CPP_MPI
      use mpi
      use mod_mympi, only: nranks
#endif
      use mod_version
      use mod_version_info
      use mod_md5sums

      IMPLICIT NONE
      INCLUDE 'inc.p'
C     .. Parameters ..
C parameter nembd1 avoids zero sized arrays.(2.1.01 R.Zeller)
      INTEGER NEMBD1
      PARAMETER (NEMBD1=NEMBD+1)
      INTEGER NPOTD
      PARAMETER (NPOTD= (2*(KREL+KORBIT) + 
     +           (1-(KREL+KORBIT))*NSPIND)*NATYPD)
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (KREL+KORBIT+1)*(LMAXD+1)**2)
      INTEGER NSPINDD
      PARAMETER (NSPINDD=NSPIND-KORBIT)
      INTEGER LMGF0D
      PARAMETER (LMGF0D= (LMAXD+1)**2)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER LMXSPD
      PARAMETER (LMXSPD= (2*LPOTD+1)**2)
      INTEGER LASSLD
      PARAMETER (LASSLD=4*LMAXD)
      INTEGER MMAXD
      PARAMETER (MMAXD = 2*LMAXD+1)
      INTEGER LM2D
      PARAMETER (LM2D= (2*LMAXD+1)**2)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
      INTEGER NOFGIJD
      PARAMETER (NOFGIJD = NATOMIMPD*NATOMIMPD+1)
      INTEGER NSYMAXD
      PARAMETER (NSYMAXD=48)
      INTEGER MAXMSHD
      PARAMETER (MAXMSHD=30)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ABASIS,ALAT,ALATNEW,BBASIS,C,CBASIS,
     +                 E1,E2,E2IN,EFERMI,ESHIFT,FCM,
     +                 HFIELD,MIXING,PI,QBOUND,RCUTXY,RCUTZ,
     +                 TK,VCONST,VOLUME0,LAMBDA_XC,RMAX,GMAX
      INTEGER I,I1,ICC,ICST,IE,IELAST,IEND,IFILE,IGF,
     +        IMIX,INS,INSREF,INTERVX,INTERVY,INTERVZ,INVMOD,
     +        IPE,IPF,IPFE,IRM,ISHIFT,
     +        ITDBRY,KCOR,KEFG,KFROZN,KHFELD,KHYP,KPRE,
     +        KSHAPE,KTE,KVMAD,KVREL,KWS,KXC,LM,LMAX,LMMAX,
     +        LMPOT,LPOT,MAXMESH,NAEZ,NATOMIMP,NVIRT,NATYP,
     +        NCLS,NEMB,NINEQ,NLAYER,NLBASIS,NLEFT,
     +        NPNT1,NPNT2,NPNT3,NPOL,NR,NRBASIS,NREF,NRIGHT,
     +        NSPIN,NSRA,NSYMAT,NOFGIJ,NQCALC,
     +        NPAN_LOG,NPAN_EQ,NCHEB,NPAN_TOT(NATYPD),
     +        NPAN_LOGNEW(NATYPD),NPAN_EQNEW(NATYPD)
      INTEGER NS, LRECABMAD, IREC
      DOUBLE PRECISION RPAN_INTERVALL(0:NTOTD,NATYPD),
     +                 RNEW(NTOTD*(NCHEBD+1),NATYPD)
      INTEGER          IPAN_INTERVALL(0:NTOTD,NATYPD)
      INTEGER IESEMICORE,NPOLSEMI,N1SEMI,N2SEMI,N3SEMI,IDOSEMICORE
      DOUBLE PRECISION FSEMICORE,EBOTSEMI,EMUSEMI,TKSEMI,R_LOG
      LOGICAL LINIPOL,LINTERFACE,LRHOSYM
      LOGICAL LCARTESIAN
      CHARACTER*40 I12,I13,I19,I25,I40
      CHARACTER*10 SOLVER
      DOUBLE PRECISION SOCSCL(KREL*LMAXD+1,KREL*NATYPD+(1-KREL))
      DOUBLE PRECISION SOCSCALE(NATYPD)
      DOUBLE PRECISION CSCL(KREL*LMAXD+1,KREL*NATYPD+(1-KREL))
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX DEZ(IEMXD),EZ(IEMXD),WEZ(IEMXD)
      DOUBLE COMPLEX DSYMLL(LMMAXD,LMMAXD,NSYMAXD),
     +               DSYMLL1(LMMAXD,LMMAXD,NSYMAXD)
      double complex, allocatable :: LEFTTINVLL(:,:,:,:,:), 
     +                               RIGHTTINVLL(:,:,:,:,:)
      DOUBLE PRECISION A(NATYPD),B(NATYPD),BRAVAIS(3,3),CLEB(NCLEB,2),
     +                 CMOMHOST(LMPOTD,NEMBD1),DRDI(IRMD,NATYPD),
     +                 DROR(IRMD,NATYPD),ECORE(20,NPOTD),
     +                 MTFAC(NATYPD),R(IRMD,NATYPD),RATOM(3,NSHELD),
     +                 RBASIS(3,NAEZD+NEMBD),RCLS(3,NACLSD,NCLSD),
     +                 RCLSIMP(3,NATOMIMPD),RECBV(3,3),
     +                 RMT(NATYPD),RMTNEW(NATYPD),FPRADIUS(NATYPD),
     +                 RMTREF(NREFD),RMTREFAT(NAEZD+NEMBD),RR(3,0:NRD),
     +                 RROT(48,3,NSHELD),RS(IRMD,0:LMAXD,NATYPD),
     +                 RSYMAT(64,3,3),RWS(NATYPD),S(0:LMAXD,NATYPD),
     +                 TLEFT(3,NEMBD1),
     +                 TRIGHT(3,NEMBD1),VBC(2),
     +                 VISP(IRMD,NPOTD),
     +                 VREF(NREFD),WG(LASSLD),
     +                 YRG(LASSLD,0:LASSLD,0:LASSLD),
     +                 ZAT(NATYPD),ZPERIGHT(3),ZPERLEFT(3)
      double precision, allocatable :: vins(:,:,:)
      INTEGER ATOM(NACLSD,NAEZD+NEMBD),ATOMIMP(NATOMIMPD),
     &        CLS(NAEZD+NEMBD),EZOA(NACLSD,NAEZD+NEMBD),
     +        ICHECK(NAEZD/NPRINCD,NAEZD/NPRINCD),
     +        ICLEB(NCLEB,4),IFUNM(NATYPD,LMXSPD),IFUNM1(LMXSPD,NATYPD),
     +        IMT(NATYPD),INIPOL(NATYPD),IPAN(NATYPD),IRC(NATYPD),
     +        IRCUT(0:IPAND,NATYPD),IRMIN(NATYPD),IRNS(NATYPD),
     +        IRWS(NATYPD),ISYMINDEX(NSYMAXD),ITITLE(20,NPOTD),
     +        IXIPOL(NATYPD),JEND(LMPOTD,0:LMAXD,0:LMAXD),
     +        KAOEZ(NATYPD,NAEZD+NEMBD),KFG(4,NATYPD),LCORE(20,NPOTD),
     +        LLMSP(NATYPD,NFUND),LMSP(NATYPD,LMXSPD),
     +        LMSP1(LMXSPD,NATYPD),LMXC(NATYPD),LOFLM(LM2D),
     +        NACLS(NCLSD),NCORE(NPOTD),NFU(NATYPD),
     +        NSH1(NSHELD),NSH2(NSHELD),NSHELL(0:NSHELD),
     +        NTCELL(NATYPD),REFPOT(NAEZD+NEMBD)
      INTEGER IJTABCALC(NOFGIJD),IJTABSYM(NOFGIJD),IJTABSH(NOFGIJD),
     +        IOFGIJ(NOFGIJD),JOFGIJ(NOFGIJD),IJTABCALC_I(NOFGIJD),
     +        IQCALC(NAEZD)
!      +        ISH(NSHELD,NOFGIJD),JSH(NSHELD,NOFGIJD),IQCALC(NAEZD)
      INTEGER, ALLOCATABLE :: ISH(:,:),JSH(:,:)
      INTEGER KMESH(IEMXD)
      INTEGER ITSCF,SCFSTEPS,KFORCE
      LOGICAL VACFLAG(2)
      CHARACTER*124 TXC(5)
      INTEGER ILM(NGSHD,3),IMAXSH(0:LMPOTD)
      DOUBLE PRECISION GSH(NGSHD)
      DOUBLE PRECISION TOLRDIF
C
C=======================================================================
C     CPA variables/ magnetisation angles -- description see RINPUT13
C
      INTEGER KMROT
      DOUBLE PRECISION QMTET(NAEZD),QMPHI(NAEZD)
C
      INTEGER NCPA,ICPA(NAEZD) 
      INTEGER ITCPAMAX,NOQ(NAEZD)
      INTEGER IQAT(NATYPD)
C
C=======================================================================
C     ITERMDIR running option introduced Apr 2003 -- Munich
C              (H. Ebert + V. Popescu) allows a self-consistent
C              determination of the magnetic configuration in REL mode
C
      DOUBLE PRECISION QMGAM(NAEZD)
      DOUBLE PRECISION QMGAMTAB(NAEZD,3),QMPHITAB(NAEZD,3),
     &                 QMTETTAB(NAEZD,3)
C=======================================================================
C
C     changes for impurity 20/02/2004 -- v.popescu according to 
C                                        n.papanikolaou VINS()
C
      INTEGER HOSTIMP(0:NATYPD)
      DOUBLE PRECISION CPATOL,CONC(NATYPD)
C-----------------------------------------------------------------------
      DOUBLE COMPLEX CREL(LMMAXD,LMMAXD),RC(LMMAXD,LMMAXD),
     &           RREL(LMMAXD,LMMAXD),SRREL(2,2,LMMAXD)
      INTEGER IRREL(2,2,LMMAXD),NRREL(2,LMMAXD)
      INTEGER J
      DOUBLE COMPLEX DROTQ(LMMAXD,LMMAXD,NAEZD)
      LOGICAL SYMUNITARY(NSYMAXD),PARA
      DOUBLE PRECISION FACT(0:100)
      DOUBLE PRECISION VTREL(IRMD*KREL+(1-KREL),NATYPD)
      DOUBLE PRECISION BTREL(IRMD*KREL+(1-KREL),NATYPD)
      DOUBLE PRECISION DRDIREL(IRMD*KREL+(1-KREL),NATYPD),
     &                 R2DRDIREL(IRMD*KREL+(1-KREL),NATYPD),
     &                 RMREL(IRMD*KREL+(1-KREL),NATYPD)
      INTEGER IRSHIFT(NATYPD),JWSREL(NATYPD),ZREL(NATYPD)
      INTEGER ISVATOM,NVATOM
      DOUBLE PRECISION ZATTEMP
C     .. Dummy variables needed only in IMPURITY program 
      DOUBLE PRECISION, allocatable :: THESME(:,:,:)

C
C-----------------------------------------------------------------------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C LDA+U LDA+U LDA+U     ph. mavropoulos according to Munich SPR-KKR
C                       h. ebert
C      input: 
C            UEFF, JEFF : input U,J parameters for each atom
C            EREFLDAU(1..NATYPD) : the energies of the projector's wave
C                                  functions (REAL)
C            LOPT(1..NATYPD): angular momentum QNUM for the atoms on
C                             which LDA+U should be applied (-1 to
C                             switch it OFF)
C      iteration index ITRUNLDAU
C      integer flag perform LDA+U IDOLDAU
C      integer flag LDA+U arrays available KREADLDAU
C      NTLDAU - number of atoms on which LDA+U is applied (<=NATYP)
C      arrays: ULDAU - calculated Coulomb matrix elements (EREFLDAU) 
C              WLDAU - potential matrix 
C              ITLDAU - integer pointer connecting the NTLDAU atoms to
C                       their corresponding index in the unit cell
C
      INTEGER IDOLDAU,ITRUNLDAU,KREADLDAU,NTLDAU
      INTEGER LOPT(NATYPD),ITLDAU(NATYPD)
      DOUBLE PRECISION EREFLDAU(NATYPD),UEFF(NATYPD),JEFF(NATYPD)
C     ..
C     .. distinguish between spin-dependent and spin-independent
C     .. quantities
       DOUBLE PRECISION ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD) 
      DOUBLE PRECISION WLDAU(MMAXD,MMAXD,NSPIND,NATYPD) ! Spin-Dependent
      DOUBLE COMPLEX PHILDAU(IRMD,NATYPD) 
C LDA+U LDA+U LDA+U
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Lloyds formula
      INTEGER LLY ! LLY <> 0 : apply Lloyds formula
      DOUBLE COMPLEX DELTAE  ! Energy difference for numerical derivative

C
C ruess: IVSHIFT test option
      INTEGER IVSHIFT
      
      !allocations:
      double precision, allocatable :: THETAS(:,:,:), THETASNEW(:,:,:)
      integer :: ierr
      
#ifdef CPP_OMP
!     .. OMP ..
      integer nth, ith          ! total number of threads, thread number
#endif

C     ..
C     .. External Functions ..
      LOGICAL OPT,TEST
      EXTERNAL OPT,TEST
C     ..
C     .. External Subroutines ..
      EXTERNAL BZKINT0,CINIT,CLSGEN_TB,DECIOPT,EPATHTB,GAUNT,GAUNT2,
     +         GFMASK,LATTIX99,RINIT,RINPUT13,SCALEVEC,
     +         STARTB1,STARTLDAU,TESTDIM,SHAPE
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,DABS,DBLE,DIMAG,LOG,MAX,SQRT
C     ..
C     ..


      allocate(THETAS(IRID,NFUND,NCELLD),THESME(IRID,NFUND,NCELLD),
     +         THETASNEW(NTOTD*(NCHEBD+1),NFUND,NCELLD),
     +         VINS(IRMIND:IRMD,LMPOTD,NSPOTD), stat=ierr)
      if(ierr/=0) stop '[main0] Error allocating thetas etc.'
      allocate(LEFTTINVLL(LMMAXD,LMMAXD,NEMBD1,NSPINDD,IEMXD),
     +    RIGHTTINVLL(LMMAXD,LMMAXD,NEMBD1,NSPINDD,IEMXD), stat=ierr)
      if(ierr/=0) stop '[main0] Error allocating lefttinvll etc.'

      
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!Write version info:
      WRITE(*,2004) 'Screened Korringa-Kohn-Rostoker ',
     &              'Electronic Structure Code',
     &              'for Bulk and Interfaces',
     &              'Juelich-Munich 2001 - 2016',
     &             ' Code version: ',trim(version(1)),
     &             ' Compile options:',(trim(version(i)),i=2,4),
     &             ' serial number for files:', serialnr
      WRITE(1337,2004) 'Screened Korringa-Kohn-Rostoker ',
     &              'Electronic Structure Code',
     &              'for Bulk and Interfaces',
     &              'Juelich-Munich 2001 - 2016',
     &             ' Code version: ',trim(version(1)),
     &             ' Compile options:',(trim(version(i)),i=2,4),
     &             ' serial number for files: ', serialnr
 2004  FORMAT( /79(1H*)/10X2A/27XA/10XA//2A/2(A,X)2(/A)//2A/79(1H*)/ )

#if defined(CPP_OMP) .OR. defined(CPP_HYBRID)
!$omp parallel shared(nth) private(ith)
      nth = omp_get_num_threads()
      ith = omp_get_thread_num()
      if(ith==0) then
      WRITE(*,'(1X,A,I5//79(1H*)/)')
     &              'Number of OpenMP threads used:',nth
      WRITE(1337,'(1X,A,I5//79(1H*)/)')
     &              'Number of OpenMP threads used:',nth
      endif
!$omp end parallel
#endif

#ifdef CPP_MPI
      WRITE(*,'(1X,A,I5//79(1H*)/)')
     &              'Number of MPI ranks used:',nranks
      WRITE(1337,'(1X,A,I5//79(1H*)/)')
     &              'Number of MPI ranks used:',nranks
#endif
!End write version info
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

!Consistency check
      IF ( (KREL.LT.0) .OR. (KREL.GT.1) )
     &     STOP ' set KREL=0/1 (non/fully) relativistic mode in inc.p'
      IF ( (KREL.EQ.1) .AND. (NSPIND.EQ.2) ) 
     &     STOP ' set NSPIND = 1 for KREL = 1 in inc.p'
C

      PI = 4.0D0*ATAN(1.0D0)
      EFERMI = 0.0d0
      CALL RINIT(NAEZD,QMGAM)
      ITSCF = 0  ! initialise SCF
      NATOMIMP = 0
      IDOLDAU = 0
      NOFGIJ = 0
      ! allocate and initialize testc and optc in t_params for run and test options
      ALLOCATE(t_params%TESTC(32),t_params%OPTC(32), stat=ierr) !CHARACTER*8
      if(ierr/=0) stop '[main0] Error allocating testc etc.'
      t_params%OPTC(1:32) =  '        '
      t_params%TESTC(1:32) = '        '
C
      CALL RINPUT13(ALAT,RBASIS,ABASIS,BBASIS,CBASIS,CLS,NCLS,
     +              E1,E2,TK,NPOL,NPNT1,NPNT2,NPNT3,
     +              EBOTSEMI,EMUSEMI,TKSEMI,NPOLSEMI,N1SEMI,N2SEMI,
     +              N3SEMI,FSEMICORE,ESHIFT,
     +              SCFSTEPS,IMIX,MIXING,QBOUND,FCM,
     +              ITDBRY,IRNS,FPRADIUS,NTCELL,NAEZ,NEMB,KAOEZ,IRM,ZAT,
     +              NINEQ,NREF,NREFD,ICST,IFILE,IGF,INS,INSREF,IPE,IPF,
     +              IPFE,KCOR,KEFG,KFROZN,KHFELD,KHYP,KPRE,KSHAPE,KTE,
     +              KFG,KVMAD,KVREL,KWS,KXC,LAMBDA_XC,LMAX,LMMAX,LMPOT,
     +              LPOT,NATYP,NSPIN,LMXC,TXC,ICC,REFPOT,
     +              ISHIFT,INTERVX,INTERVY,INTERVZ,
     +              HFIELD,MTFAC,VBC,VCONST,LINIPOL,INIPOL,
     +              IXIPOL,LRHOSYM,
     +              I12,I13,I19,I25,I40,
     +              NLBASIS,NRBASIS,NLEFT,NRIGHT,ZPERLEFT,
     +              ZPERIGHT,TLEFT,TRIGHT,LINTERFACE,RCUTZ,RCUTXY,
     +              RMTREF,RMTREFAT,KFORCE,
     &              KMROT,QMTET,QMPHI,NCPA,ICPA,ITCPAMAX,CPATOL,NOQ,
     &              IQAT,CONC,SOLVER,SOCSCL,CSCL,KREL,SOCSCALE,
     &              LOPT,UEFF,JEFF,EREFLDAU,KREADLDAU,
     &              LMAXD,LPOTD,NSPIND,NAEZD,NATYPD,NEMBD,NPRINCD,
     &              IRMD,IRNSD,NPAN_LOG,NPAN_EQ,NCHEB,R_LOG,IVSHIFT,
     &              TOLRDIF,LLY,DELTAE,
     &              LCARTESIAN,BRAVAIS,RMAX,GMAX)
c
C ================================================ deal with the lattice

      CALL LATTIX99(LINTERFACE,ALAT,NATYP,NAEZ,CONC,RWS,BRAVAIS,
     &              RECBV,VOLUME0,RR,NR,NRD,NATYPD)


      CALL SCALEVEC(LCARTESIAN,RBASIS,ABASIS,BBASIS,CBASIS,
     +              NLBASIS,NRBASIS,NLEFT,NRIGHT,ZPERLEFT,ZPERIGHT,
     +              TLEFT,TRIGHT,LINTERFACE,NAEZ,NEMB,BRAVAIS,KAOEZ,NOQ,
     &              NAEZD,NATYPD,NEMBD)
      ! After SCALEVEC all basis positions are in cartesian coords.

      NVIRT = 0
      IF ( OPT('VIRATOMS') ) THEN
         WRITE(1337,*) 'Calling ADDVIRATOMS'
         CALL ADDVIRATOMS14(
     &  LINTERFACE,NVIRT,NAEZ,NAEZD,NATYPD,NEMB,NEMBD,
     &  RBASIS,.TRUE.,BRAVAIS,NCLS,NINEQ,REFPOT,KAOEZ,NOQ,NREF,RMTREFAT,
     &  I25) 
      ENDIF



      CALL CLSGEN_TB(NAEZ,NEMB,NVIRT,RR,NR,RBASIS,KAOEZ,ZAT,CLS,NCLS,
     &               NACLS,ATOM,EZOA, 
     &               NLBASIS,NRBASIS,NLEFT,NRIGHT,ZPERLEFT,ZPERIGHT,
     &               TLEFT,TRIGHT,RMTREF,RMTREFAT,VREF,
     &               REFPOT,NREF,RCLS,RCUTZ,RCUTXY,LINTERFACE,ALAT,
     &               NAEZD,NATYPD,NEMBD,NPRINCD,NRD,NACLSD,NCLSD,NREFD)

! Now the clusters, reference potentials and muffin-tin radii have been set.
C ......................................................................
Consistency check 
C
      IF ( (KREL.EQ.1) .AND. (INS.NE.0) ) THEN
         WRITE(6,*)
     &        ' FULL-POTENTIAL RELATIVISTIC mode not implemented '
         STOP ' set INS = 0 in the input'
      END IF
C     
      IF ( KVREL.LE.1 ) THEN
         IF ( KREL.EQ.1 ) STOP
     &        ' KVREL <= 1 in input, but relativistic program used'
      ELSE
         IF ( KREL.EQ.0 ) STOP
     &        ' KVREL > 1 in input, but non-relativistic program used'
      END IF
C ......................................................................
C
      E2IN = E2
      NSRA = 1
      IF (KVREL.GE.1) NSRA = 2
      IF (KVREL.GE.2) NSRA = 3
C
      CALL TESTDIM(NSPIN,NAEZ,NEMB,NATYP,LMAX,IRM,INS,INSREF,NREF,
     &             IRNS,NCLS,NLAYER,
     &             KREL,LMAXD,NSPIND,NAEZD,NATYPD,NREFD,NCLSD,
     &             NEMBD,NPRINCD,KNOSPH,IRMD,IRNSD,KORBIT)
C
      IF ( INS.GT.0 )    OPEN (19,FILE=I19,STATUS='old',
     &                         FORM='formatted')
      IF ( IFILE.EQ.13 ) OPEN (IFILE,FILE=I13,STATUS='old',
     &                         FORM='formatted')
      IF ( ICC.GT.0 )    OPEN (25,FILE=I25,STATUS='unknown',
     &                         FORM='formatted')
C
      CALL STARTB1(IFILE,1337,1337,IPE,KVREL,KWS,LMAX,1,NATYP,
     +             ALATNEW,RMTNEW,RMT,ITITLE,IMT,IRC,VCONST,INS,
     +             IRNS,FPRADIUS,LPOT,NSPIN,VINS,IRMIN,KSHAPE,NTCELL,
     +             IRCUT,IPAN,THETAS,IFUNM,NFU,LLMSP,LMSP,E2IN,VBC,C,
     +             DROR,RS,S,VISP,RWS,ECORE,LCORE,NCORE,DRDI,R,ZAT,A,B,
     +             IRWS,1,LMPOTD,IRMIND,IRMD,LMXSPD,IPAND,IRID,
     +             IRNSD,LMAXD,NATYPD,NCELLD,NFUND,NSPOTD,IVSHIFT)
     
     
      ! find md5sums for potential and shapefunction
      call get_md5sums(INS, I13, I19)
      write(1337,'(A,A)') 
     &'Doing calculation with potential: ',md5sum_potential
      if(INS>0) then
         write(1337,'(A,A)') 
     &   'Doing calculation with shapefun: ',md5sum_shapefun
      end if


      IF ( TEST('Vspher  ') ) THEN
         WRITE(1337,*) 'TEST OPTION Vspher,', 
     &        'keeping only spherical component of potential.' 
         VINS(IRMIND:IRMD,2:LMPOTD,1:NSPOTD) = 0.D0
      ENDIF


      IF (OPT('zeropot ').OR.TEST('zeropot ')) THEN
        WRITE(1337,*) 'Using OPT zeropot, setting potential to zero.'
        WRITE(1337,*) 
     &           'Using OPT zeropot, setting nuclear charge to zero.'
        VINS(IRMIND:IRMD,1:LMPOTD,1:NSPOTD) = 0.D0
        VISP(1:IRMD,1:NPOTD) = 0.D0
        ZAT(1:NATYPD) = 0.D0
      ENDIF
C
      DO I1 = 1,NATYPD
         DO LM = 1,LMXSPD
            IFUNM1(LM,I1) = IFUNM(I1,LM)
            LMSP1(LM,I1) = LMSP(I1,LM)
         END DO
      END DO
C ----------------------------------------------------------------------
C update Fermi energy, adjust energy window according to running options
C
      IF ( NPOL.EQ.0 ) EFERMI = E2IN
      IF ( OPT('GF-EF   ') .OR. OPT('DOS-EF  ') ) THEN
         E1 = E2IN
         IF ( OPT('GF-EF   ') ) THEN
            WRITE (1337,FMT=9070)
         ELSE
            WRITE (1337,FMT=9080)
         END IF
      END IF

      IF ( DABS(E2IN-E2).GT.1D-10 .AND. NPOL.NE.0 ) E2 = E2IN
C ----------------------------------------------------------------------
      IF (OPT('GENPOT  ')) THEN
          REWIND(3)
          CALL GENERALPOT(3,1,NATYP,NSPIN,ZAT,ALAT,RMT,RMTNEW,RWS,
     +         R,DRDI,VISP,IRWS,A,B,INS,IRNS,
     +         LPOT,VINS,QBOUND,IRC,KSHAPE,E2IN,VBC,ECORE,
     +         LCORE,NCORE,LMPOTD,IRMD,IRMIND)
          CLOSE(3)
      END IF
C ----------------------------------------------------------------------
C --> Apply external magnetic field
!     from startb1 moved here
      IF (KHFELD.EQ.1) THEN
c
c---> maybe apply a magnetic field
c
         call BSHIFT_NS(VISP,VINS,NATYP,NSPIN,
     +   IRCUT,IRC,IRMIN,NTCELL,IMAXSH,ILM,IFUNM,LMSP,LMPOT,GSH,
     +   THETAS,THESME,R,KSHAPE,HFIELD,INIPOL)
      END IF
      if ( TEST('vpotout ') ) then !ruess
        open(unit=54633163,file='test_vpotout_bshift')
        do i1=1,natyp*nspin
            write(54633163,*) '# visp of atom ',i1
            write(54633163,'(50000E)') visp(:,i1)
        end do !iatom
        do i1=1,natyp*nspin
            write(54633163,*) '# vins of atom ',i1
            write(54633163,'(50000E)') vins(:,:,i1)
        end do !iatom
        close(54633163)
      end if
C ----------------------------------------------------------------------
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C -->  deal with the potential in the RELATIVISTIC CASE
C
      PARA = .TRUE.
      IF (KREL+KORBIT.EQ.1) THEN
C---------------------------------------------------------------
         IF (NSPIN.EQ.1) THEN
C
C --> for paramagnetic (NSPIN=1) input potential fill up also the 
C     V(DOWN), ECORE(DOWN), LCORE(DOWN), NCORE(DOWN) and ITITLE(DOWN)
C     arrays (needed)
C
            DO I = NATYP,1,-1
               J = 2*I - 1
               CALL DCOPY(IRMD,VISP(1,I),1,VISP(1,J),1)
               CALL DCOPY(IRMD,VISP(1,J),1,VISP(1,J+1),1)
C
               CALL DCOPY(20,ECORE(1,I),1,ECORE(1,J),1)
               CALL DCOPY(20,ECORE(1,J),1,ECORE(1,J+1),1)
C
               NCORE(J) = NCORE(I)
               NCORE(J+1) = NCORE(J)
C
               DO I1=1,20
                  LCORE(I1,J) = LCORE(I1,I)
                  LCORE(I1,J+1) = LCORE(I1,J)
                  ITITLE(I1,J ) = ITITLE(I1,I)
                  ITITLE(I1,J+1) = ITITLE(I1,J)
               END DO
            END DO
C---------------------------------------------------------------
         ELSE !NSPIN.EQ.1
C---------------------------------------------------------------
C --> check whether, although NSPIN=2 at input, the system is
C     paramagnetic (useful for symmetry cosiderations)
C
            DO I = 1,2*NATYP-1,2
               DO J=1,IRMD
                  IF (ABS(VISP(J,I)-VISP(J,I+1)).GT.1D-5) 
     &                 PARA = .FALSE.
               END DO
            END DO
            IF (PARA) THEN
               DO I=1,2*NATYP-1,2
                  CALL DCOPY(IRMD,VISP(1,I),1,VISP(1,I+1),1)
               END DO
            END IF
            
         END IF !NSPIN.EQ.1
C---------------------------------------------------------------
C
C --> finally, convert input potential to the internal relativistic
C     form VTREL,BTREL. Set up auxiliary arrays (needed in the REL 
C     routines) ZREL, JWSREL, RMREL, DRDIREL, R2DRDIREL, IRSHIFT
C
         IF(KREL.eq.1) ! call this only if relativisitic solver is used
     &   CALL RELPOTCVT(1,VISP,ZAT,R,DRDI,IRCUT,
     &        VTREL,BTREL,ZREL,RMREL,JWSREL,DRDIREL,R2DRDIREL,IRSHIFT,
     &        IPAND,IRMD,NPOTD,NATYPD)
      END IF !KREL+KORBIT.EQ.1
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C --> set up energy contour
C
      IDOSEMICORE = 0
      IF ( OPT('SEMICORE') ) IDOSEMICORE = 1
C
      CALL EPATHTB(EZ,DEZ,E2IN,IELAST,IESEMICORE,IDOSEMICORE,
     &             E1,E2,TK,NPOL,NPNT1,NPNT2,NPNT3,EBOTSEMI,EMUSEMI,
     &             TKSEMI,NPOLSEMI,N1SEMI,N2SEMI,N3SEMI,IEMXD)
      DO IE = 1,IELAST
        WEZ(IE) = -2.D0/PI*DEZ(IE)
        IF ( IE.LE.IESEMICORE ) WEZ(IE) = WEZ(IE)*FSEMICORE
      END DO
C
C --> update energy contour for Fermi-surface generation    ! fswrt=fermi-surface write
      IF (OPT('FERMIOUT'))THEN                              ! fswrt
        IF(AIMAG(EZ(1))>0d0) STOP 'E has imaginary part'    ! fswrt
        IELAST=3                                            ! fswrt
        EZ(2) = EZ(1) + CMPLX(1.0D-03,0.0D0)                ! fswrt
        EZ(3) = EZ(1) - CMPLX(1.0D-03,0.0D0)                ! fswrt
      END IF                                                ! fswrt
C
C --> update the value of NSPIN to be consistent with REL mode
C
      IF (KREL.EQ.1) NSPIN = 1
C
      CALL GAUNT2(WG,YRG,4*LMAXD)
      CALL GAUNT(LMAX,LPOT,WG,YRG,CLEB,LOFLM,ICLEB,IEND,JEND,
     &           NCLEB,LMAXD,LMGF0D,LMPOTD)
C
C --> set up of GAUNT coefficients C(l,m;l',m';l'',m'') for all 
C     nonvanishing (l'',m'')-components of the shape functions THETAS
C
      IF (KSHAPE.NE.0) 
     +     CALL SHAPE(LPOT,NATYP,GSH,ILM,IMAXSH,LMSP,NTCELL,WG,YRG,
     +                LASSLD,LMPOTD,NATYPD,NGSHD)

C

C ----------------------------------------------------------------------
C --> calculate Madelung constants (needed only for SCF calculations)
C ----------------------------------------------------------------------
cfivos      IF ( SCFSTEPS.GT.1 .OR. ICC.GT.0 ) THEN
      IF (NPOL.NE.0) THEN  ! No madelung calculation in case of DOS.
         !OPEN(99,FILE='madelinfo.txt')

! Use option 'ewald2d' if the madelung summation is to be carried out in
! single-slab mode, otherwise it is carried out in repeated (periodic)
! slab mode.
! Reason: the 2d-mode gives wrong results sometimes [e.g. in diamond
! structure (110)].
         IF (LINTERFACE.AND.( OPT('ewald2d ').OR.OPT('DECIMATE'))) THEN ! ewald2d
            WRITE(*,*) 'Calling MADELUNG2D'
C -------------------------------------------------------------- 2D case
            CALL MADELUNG2D(LPOT,YRG,WG,NAEZ,ALAT,VOLUME0,
     &                         BRAVAIS,RECBV,RBASIS,RMAX,GMAX,
     &                         NLBASIS,NLEFT,ZPERLEFT,TLEFT,
     &                         NRBASIS,NRIGHT,ZPERIGHT,TRIGHT,
     &                         LMXSPD,LASSLD,LPOTD,LMPOTD,
     &                         NMAXD,ISHLD,NEMBD1,WLENGTH)
            WRITE(*,*) 'Exited MADELUNG2D'
         ELSE
C -------------------------------------------------------------- 3D case
            IF ( LINTERFACE )
     &         CALL GETBR3(                                         ! ewald2d
     >                     NEMBD1,NLBASIS,ALAT,TLEFT,NRBASIS,TRIGHT,
     X                     BRAVAIS,RECBV,VOLUME0)

            WRITE(*,*) 'Calling MADELUNG3D'
            CALL MADELUNG3D(LPOT,YRG,WG,NAEZ,ALAT,VOLUME0,
     &                         BRAVAIS,RECBV,RBASIS,RMAX,GMAX,
     &                         NAEZD,LMXSPD,LASSLD,LPOTD,LMPOTD,
     &                         NMAXD,ISHLD,NEMBD,WLENGTH)
            WRITE(*,*) 'Exited MADELUNG3D'
         END IF

         !CLOSE(99)
      ELSE !NPOL==0
      ! write dummy files

         !DOUBLE PRECISION AVMAD(LMPOTD,LMPOTD),BVMAD(LMPOTD)
         LRECABMAD = WLENGTH*2*LMPOTD*LMPOTD + WLENGTH*2*LMPOTD
         OPEN (69,ACCESS='direct',RECL=LRECABMAD,
     &         FILE='abvmad.unformatted', FORM='unformatted')
          DO I = 1,NAEZ
             DO J = 1,NAEZ
                IREC = J + NAEZ*(I-1)
                WRITE (69,REC=IREC) 0.0d0, 0.0d0!AVMAD,BVMAD
             END DO
          END DO
          CLOSE(69)

      ENDIF !NPOL==0
C ----------------------------------------------------------------------
cfivos      END IF
C ======================================================================
C
C
C ======================================== set up I,J pairs for ICC = -1
C
      IF ( ICC.LT.0 ) 
     &     CALL SETGIJTAB(LINTERFACE,ICC,NAEZ,IQAT,RBASIS,BRAVAIS,
     &                    NATOMIMP,ATOMIMP,RCLSIMP,
     &                    NOFGIJ,IJTABCALC,IOFGIJ,JOFGIJ,NQCALC,IQCALC,
     &                    NATOMIMPD,IJTABCALC_I)
C
C ======================================================================
C
      DSYMLL=(0d0,0d0)
      DSYMLL1=(0d0,0d0)
      
      ! allocate ISH and JSH index arrays that might be quite big
      ALLOCATE(ISH(NSHELD,NOFGIJD),JSH(NSHELD,NOFGIJD), stat=ierr)
      if(ierr/=0) stop '[main0] Error allocating ISH, JSH'

      CALL BZKINT0(NSHELL,NAEZ,NATYP,NOQ,
     +             RBASIS,KAOEZ,ICC,BRAVAIS,RECBV,ATOMIMP,
     +             RSYMAT,ISYMINDEX,NSYMAT,I25,NATOMIMP,
     +             NSH1,NSH2,RCLSIMP,RATOM,
     +             IJTABSYM,IJTABSH,IJTABCALC,
     +             IOFGIJ,JOFGIJ,NOFGIJ,ISH,JSH,
     +             RROT,DSYMLL1,PARA,QMTET,QMPHI,SYMUNITARY,
     +             HOSTIMP,INTERVX,INTERVY,INTERVZ,
     +             IELAST,EZ,KMESH,MAXMESH,MAXMSHD,
     +             NSYMAXD,KREL+KORBIT,LMAXD,LMMAXD,KPOIBZ,NAEZD,NATYPD,
     +             NATOMIMPD,NSHELD,NEMBD)
C
C ======================================================================
C
      IF ( OPT('KKRFLEX ') ) THEN

      CALL WRITEHOSTSTRUCTURE(BRAVAIS,NATYP,RBASIS,NAEZD,NEMBD)

        OPEN (58,FILE='kkrflex_atominfo',FORM='FORMATTED')
        call version_print_header(58,
     &           '; '//md5sum_potential//'; '//md5sum_shapefun)
        NVATOM=0
        DO I = 1,NATOMIMP
          IF (KAOEZ(1,ATOMIMP(I))==-1) NVATOM=NVATOM+1
        END DO
        WRITE (58,'(500A)') '#NATOM   NTOTATOM' 
        WRITE (58,*) NATOMIMP,NATOMIMP-NVATOM 
        WRITE (58,'(500A)') '#Impurity positions x,y,z|Core Charge|Virtual 
     +                Atom?|Remove Atom?|LMAX'
        DO I = 1,NATOMIMP
          IF (KAOEZ(1,ATOMIMP(I))==-1) THEN
            ZATTEMP=0.D0
            ISVATOM=1
            NVATOM=NVATOM+1
          ELSE
            ISVATOM=0
            ZATTEMP=ZAT(KAOEZ(1,ATOMIMP(I)))
          END IF
          WRITE (58,'(3F,F6.2,3I)') (RCLSIMP(J,I),J=1,3), ZATTEMP,
     +                  ISVATOM,0,LMAXD
        END DO
        CLOSE (58)
      END IF
C
C ======================================================================
C
C
! fivos write out nshell and nsh1,nsh2 into standard output and in file shells.dat
      IF (ICC.NE.0 .and. .not.OPT('KKRFLEX ')) THEN
         OPEN(58,FILE='shells.dat')
         write(1337,*) 'Writing out shells (also in shells.dat):' ! fivos
         write(1337,*) 'itype,jtype,iat,jat,r(iat),r(jat)' ! fivos
         write(1337,*) NSHELL(0), 'NSHELL(0)' ! fivos
         write(58,*) NSHELL(0), 'NSHELL(0)' ! fivos
         do i1 = 1,NSHELL(0)    ! fivos
            write(1337,*) i1,NSHELL(i1), 
     &                 'No. of shell, No. of atoms in shell' ! fivos
            write(58,*) i1,NSHELL(i1), 
     &                 'No. of shell, No. of atoms in shell' ! fivos
            do lm = 1,NSHELL(i1) ! fivos
               write(1337,*) 'ish(i1,lm)',ish(i1,lm)
               if(ISH(i1,lm)>0 .and. JSH(i1,lm)>0) then   !fix bernd
                 write(1337,8614) NSH1(i1),NSH2(i1) ! fivos
     &              ,ISH(i1,lm),JSH(i1,lm) ! fivos
     &              ,(RCLSIMP(i,ISH(i1,lm)),i=1,3) ! fivos
     &              ,(RCLSIMP(i,JSH(i1,lm)),i=1,3) ! fivos
                 write(58,8614) NSH1(i1),NSH2(i1) ! fivos
     &              ,ISH(i1,lm),JSH(i1,lm) ! fivos
     &              ,(RCLSIMP(i,ISH(i1,lm)),i=1,3) ! fivos
     &              ,(RCLSIMP(i,JSH(i1,lm)),i=1,3) ! fivos
               else!fix bernd
                 write(1337,8615) NSH1(i1),NSH2(i1) ! fix bernd
     &              ,ISH(i1,lm),JSH(i1,lm)          ! fix bernd
                 write(58,8615) NSH1(i1),NSH2(i1)   ! fix bernd
     &              ,ISH(i1,lm),JSH(i1,lm)          ! fix bernd
               end if!fix bernd
 8614          format(4i5,6f16.6) ! fivos
 8615          format(4i5) ! fix bernd
            enddo               ! fivos
         enddo                  ! fivos
         write(1337,*) '###################'
         CLOSE(58)
      ENDIF
! end fivos

      CALL GFMASK(LINTERFACE,ICHECK,ICC,INVMOD,NSH1,NSH2,NAEZ,
     +            NSHELL(0),NAEZD,NPRINCD)
C
C======================================================================
C set up transformation matrices between REL/NREL representations
C 
      IF((KREL+KORBIT).EQ.1) CALL DRVBASTRANS(RC,CREL,RREL,SRREL,NRREL,
     &                                 IRREL,LMAXD+1,LMMAXD,2*(LMAXD+1),
     &                                 LMMAXD+2*(LMAXD+1),MMAXD,
     &                                 2*(LMAXD+1)*MMAXD)
      IF (OPT('NEWSOSOL')) THEN
       DO NS=1,NSYMAT
        CALL CHANGEREP(DSYMLL1(1,1,NS),'REL>RLM',DSYMLL(1,1,NS),LMMAXD,
     &                 LMMAXD,RC,CREL,RREL,'DSYMLL',0)
       ENDDO
c       DSYMLL(:,:,:)=DSYMLL1(:,:,:)
      ELSE
       DSYMLL(:,:,:)=DSYMLL1(:,:,:)
      ENDIF
C                
C======================================================================
C  for the case that the magnetisation is rotated with respect to
C  the (001)-direction (KMROT<>0) calculate the rotation matrices
C  to switch between the CRYSTAL and LOCAL frames of reference
C
      CALL CINIT(LMMAXD*LMMAXD*NAEZD,DROTQ)

      IF (KMROT.NE.0) THEN
         FACT(0) = 1.0D0
         DO I = 1,100
            FACT(I) = FACT(I-1)*DBLE(I)
         END DO
C
         DO I1=1,NAEZ
            CALL CALCROTMAT(MMAXD,(KREL+KORBIT)*3,
     &                      QMPHI(I1),QMTET(I1),0.0D0,DROTQ(1,1,I1),
     &                      FACT,LMMAXD)
         END DO
      END IF
C =========================================== treat decimation I/O cases
C
      IF ( OPT('deci-pot') ) 
     & CALL OUTPOTHOST(ALAT,INS,KREL+KORBIT,KMROT,NSPIN,NAEZ,NATYP,E2IN,
     &                 BRAVAIS,RBASIS,QMTET,QMPHI,NOQ,KAOEZ,IQAT,
     &                 ZAT,CONC,IPAN,IRCUT,SOLVER,SOCSCL,CSCL,
     &                 IRWS,RMTNEW,RWS,R,DRDI,VISP,
     &                 IRSHIFT,RMREL,DRDIREL,VTREL,BTREL,
     &                 LMAXD,NATYPD,NAEZD,IPAND,IRMD)
C
      IF ( OPT('deci-out') ) 
     &  CALL OUTTMATHOST(ALAT,INS,KREL+KORBIT,KMROT,NSPIN,NAEZ,LMMAX,
     &                   BRAVAIS,RBASIS,QMTET,QMPHI,E2IN,TK,
     &                   NPOL,NPNT1,NPNT2,NPNT3)
C
      IF ( OPT('DECIMATE') ) 
     &  CALL DECIOPT(ALAT,INS,KREL+KORBIT,KVREL,KMROT,NSPIN,NAEZ,LMMAX,
     &               BRAVAIS,TK,NPOL,NPNT1,NPNT2,NPNT3,
     &               EZ,IELAST,KAOEZ,
     &               LEFTTINVLL,RIGHTTINVLL,VACFLAG,NLBASIS,NRBASIS,
     &               CMOMHOST,VREF,RMTREF,NREF,REFPOT(NAEZ),
     &               LMAXD,LMGF0D,LMMAXD,LM2D,NEMBD1,IEMXD,
     &               NSPINDD,LMPOTD,NATYPD,IRMD,IPAND)
C
C ======================================================================
C
C ======================================================================
C ITERMDIR  -- initialise 
C
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
C
C ======================================================================
C LDA+U -- initialise 
C
      IF (OPT('LDA+U   ')) then
           CALL STARTLDAU(ITRUNLDAU,IDOLDAU,KREADLDAU,LOPT,UEFF,JEFF,
     &                    EREFLDAU,NATYP,NSPIN,WLDAU,ULDAU,PHILDAU,IRWS,
     &                    NTLDAU,ITLDAU,IRMD,NATYPD,NSPIND,MMAXD)
      END IF
C
C LDA+U
C ======================================================================
C
C ======================================================================
C =        write out information for the other program parts           =
C ======================================================================
C
C new solver for full-potential, spin-orbit, initialise
       IF (OPT('NEWSOSOL')) THEN

        CALL CREATE_NEWMESH(NATYPD,LMAXD,LPOTD,IRMD,IRNSD,IPAND,
     +                      IRID,NTOTD,NFUND,NCHEBD,NTOTD*(NCHEBD+1),
     +                      NSPIN,R,IRMIN,IPAN,IRCUT,
     +                      R_LOG,NPAN_LOG,NPAN_EQ,NCHEB,
     +                      NPAN_LOGNEW,NPAN_EQNEW,
     +                      NPAN_TOT,RNEW,RPAN_INTERVALL,IPAN_INTERVALL,
     +                      NCELLD,NTCELL,THETAS,THETASNEW)

      ENDIF

      CALL WUNFILES(NPOL,NPNT1,NPNT2,NPNT3,IELAST,TK,E1,E2,EZ,WEZ,
     &              EFERMI,NPOLSEMI,N1SEMI,N2SEMI,N3SEMI,IESEMICORE,
     &              TKSEMI,EBOTSEMI,EMUSEMI,FSEMICORE,
     &              VINS,VISP,VBC,VTREL,BTREL,RMREL,
     &              DRDIREL,R2DRDIREL,ZREL,JWSREL,IRSHIFT,
     &              ITSCF,SCFSTEPS,CMOMHOST,ECORE,LCORE,NCORE,
     &              QMTET,QMPHI,QMPHITAB,QMTETTAB,QMGAMTAB,DROTQ,
     &              NSRA,INS,NATYP,NAEZ,NINEQ,NREF,NSPIN,LMAX,
     &              NCLS,ICST,IPAN,IRCUT,ALAT,ZAT,R,DRDI,
     &              REFPOT,RMTREF,VREF,IEND,JEND,CLEB,ICLEB,
     &              ATOM,CLS,RCLS,NACLS,LOFLM,SOLVER,SOCSCL,CSCL,
     &              ICC,IGF,NLBASIS,NRBASIS,NCPA,ICPA,ITCPAMAX,
     &              CPATOL,RBASIS,RR,EZOA,NSHELL,NSH1,NSH2,
     &              IJTABCALC,IJTABCALC_I,ISH,JSH,IJTABSYM,IJTABSH,
     &              NOFGIJ,NQCALC,IQCALC,KMROT,KAOEZ,IQAT,NOQ,CONC,
     &              KMESH,MAXMESH,NSYMAT,SYMUNITARY,RROT,
     &              DSYMLL,INVMOD,ICHECK,
     &              NATOMIMP,RATOM,ATOMIMP,
     &              RC,CREL,RREL,SRREL,NRREL,IRREL,
     &              LEFTTINVLL,RIGHTTINVLL,VACFLAG,
     &              A,B,IFUNM,IFUNM1,INTERVX,INTERVY,INTERVZ,ITITLE,
     &              LMSP1,NTCELL,THETAS,
     &              LPOT,LMPOT,NRIGHT,NLEFT,LINTERFACE,
     &              IMIX,MIXING,QBOUND,FCM,ITDBRY,IRNS,KPRE,KSHAPE,KTE,
     &              KVMAD,KXC,LAMBDA_XC,TXC,ISHIFT,IXIPOL,LRHOSYM,
     &              KFORCE,LMSP,LLMSP,RMT,RMTNEW,RWS,IMT,IRC,IRMIN,IRWS,
     &              NFU,HOSTIMP,GSH,ILM,IMAXSH,IDOLDAU,ITRUNLDAU,NTLDAU,
     &              LOPT,ITLDAU,UEFF,JEFF,EREFLDAU,ULDAU,WLDAU,PHILDAU,
     &              IEMXD,IRMIND,IRMD,LMPOTD,NSPOTD,NPOTD,NATYPD,
     &              NEMBD1,LMMAXD,NAEZD,IPAND,NAEZD+NEMBD,NREFD,LMAXD,
     &              NCLEB,NACLSD,NCLSD,LM2D,LMAXD+1,MMAXD,NRD,NSHELD,
     &              NSYMAXD,NAEZD/NPRINCD,NATOMIMPD,NOFGIJD,
     &              NSPIND,IRID,NFUND,NCELLD,LMXSPD,NGSHD,KREL,NTOTD,
     &              NCHEBD,NPAN_LOG,NPAN_EQ,NPAN_LOGNEW,NPAN_EQNEW,
     &              NCHEB,R_LOG,NPAN_TOT,RNEW,RPAN_INTERVALL,
     &              IPAN_INTERVALL,NSPINDD,THETASNEW,SOCSCALE,TOLRDIF,
     &              LLY,DELTAE,RCLSIMP)


      IF (OPT('FERMIOUT'))THEN                                        ! fswrt
        CALL WRITE_TBKKR_FILES(LMAX,NEMB,NCLS,NATYP,NAEZ,IELAST,INS,  ! fswrt
     +                         ALAT,BRAVAIS,RECBV,RBASIS,CLS,NACLS,   ! fswrt
     +                         RCLS,EZOA,ATOM,RR,NSPIN)               ! fswrt
      END IF                                                          ! fswrt


      IF (OPT('GREENIMP')) THEN                                        ! GREENIMP
         ! fill array dimensions and allocate arrays in t_imp          ! GREENIMP
         call init_params_t_imp(t_imp,IPAND,NATYPD,IRMD,IRID,NFUND,    ! GREENIMP
     &                          NSPIN,IRMIND,LMPOTD)                   ! GREENIMP
         call init_t_imp(t_inc,t_imp)                                  ! GREENIMP
                                                                       ! GREENIMP
         ! next read impurity potential and shapefunction              ! GREENIMP
         CALL READIMPPOT(NATOMIMP,INS,1337,0,0,2,NSPIN,LPOT,           ! GREENIMP
     &                   t_imp%IPANIMP,t_imp%THETASIMP,t_imp%IRCUTIMP, ! GREENIMP
     &                   t_imp%IRWSIMP,KHFELD,HFIELD,t_imp%VINSIMP,    ! GREENIMP
     &                   t_imp%VISPIMP,t_imp%IRMINIMP,                 ! GREENIMP
     &                   t_imp%RIMP,t_imp%ZIMP)                        ! GREENIMP
      END IF                                                           ! GREENIMP


      IF (ISHIFT.EQ.2) THEN                                           ! fxf
         OPEN (67,FILE='vmtzero',FORM='formatted')                    ! fxf
         WRITE (67,9090) VBC(1)                                       ! fxf
         CLOSE(67)                                                    ! fxf
 9090    FORMAT(D20.12)                                               ! fxf
      END IF                                                          ! fxf

C Check for inputcard consistency in case of qdos option
      IF (OPT('qdos    ')) THEN
         write(1337,*)
         write(1337,*) '     < QDOS > : consistency check '
         IF ((NPOL.NE.0).AND.(NPNT1.EQ.0).AND.(NPNT3.EQ.0)) THEN
            STOP 'For qdos calculation change enery contour to dos path'
         ENDIF
         IF (TK.GT.50.d0) write(*,*) 'WARNING:  high energy
     & smearing due to high value of TEMPR for energy contour
     & integration could not be of advantage. Consider changeing
     & ''TEMPR'' to lower value'
         IF (TK.GT.50.d0) write(1337,*) 'WARNING:  high energy
     & smearing due to high value of TEMPR for energy contour
     & integration could not be of advantage. Consider changeing
     & ''TEMPR'' to lower value'
         write(1337,*) '       QDOS: consistecy check complete'
      ENDIF


C ======================================================================
C
      WRITE (1337,'(79(1H=),/,31X,"< KKR0 finished >",/,79(1H=),/)')
 9070 FORMAT (5X,'INFO:  Output of cluster Green function at E Fermi')
 9080 FORMAT (5X,'INFO:  Determination of DOS at E Fermi')
 
      ! deallocate unused arrays
      deallocate(thetas, thetasnew, thesme, vins, stat=ierr)
      if(ierr/=0) stop '[main0] Error deallocating thetas etc.'
      deallocate(LEFTTINVLL, RIGHTTINVLL, stat=ierr)
      if(ierr/=0) stop '[main0] Error deallocating lefttinvll etc.'
      deallocate(ISH, JSH, stat=ierr)
      if(ierr/=0) stop '[main0] Error deallocating ISH and JSH'

      END SUBROUTINE !MAIN0
      
      
      
      
      
      SUBROUTINE BSHIFT_NS(VISP,VINS,NATYP,NSPIN,
     +     IRCUT,IRC,IRMIN,NTCELL,IMAXSH,ILM,IFUNM,LMSP,LMPOT,GSH,
     +     THETAS,THESME,RMESH,KSHAPE,HFIELD,INIPOL)
      implicit none
c Adds a constant (=VSHIFT) to the potentials of atoms
c
c Parameters:
      include 'inc.p'
      INTEGER NPOTD,LMPOTD,LMXSPD,IRMIND
      PARAMETER (NPOTD=NSPIND*NATYPD,LMPOTD= (LPOTD+1)**2
     &     ,LMXSPD= (2*LPOTD+1)**2,IRMIND=IRMD-IRNSD)
c Input
      INTEGER KSHAPE,LMPOT,NATYP,NSPIN
      INTEGER IRCUT(0:IPAND,NATYPD),IRC(NATYPD),NTCELL(NATYPD)
     &     ,IMAXSH(0:LMPOTD),ILM(NGSHD,3),IFUNM(NATYPD,LMXSPD)
     &     ,LMSP(NATYPD,LMXSPD),IRMIN(NATYPD),INIPOL(NATYPD)
      DOUBLE PRECISION GSH(NGSHD),THETAS(IRID,NFUND,NCELLD)
     &     ,RMESH(IRMD,NATYPD)
      DOUBLE PRECISION THESME(IRID,NFUND,NCELLD)
      DOUBLE PRECISION VSHIFT,HFIELD


c Input/Output:
      DOUBLE PRECISION VISP(IRMD,NPOTD),VINS(IRMIND:IRMD,LMPOTD,NSPOTD)

c Inside
      DOUBLE PRECISION PSHIFTLMR(IRMD,LMPOTD),PSHIFTR(IRMD)
      INTEGER ISPIN,IH,IPOT,IR,LM,IMT1,IRC1,IRMIN1
      DOUBLE PRECISION RFPI
      
      RFPI = SQRT(16.0D0*ATAN(1.0D0))


      DO IH = 1,NATYP
      
      
         IMT1 = IRCUT(1,IH)
         IRC1 = IRC(IH)
         IRMIN1 = IRMIN(IH)

         DO ISPIN = 1,NSPIN
         
            ! shift potential spin dependent
            VSHIFT = -DBLE(2*ISPIN-3)*HFIELD*INIPOL(IH)

            WRITE (1337,*) 'SHIFTING OF THE POTENTIALS OF ATOM',IH,
     &           'spin',ispin,' BY', VSHIFT, 'RY.'
            IPOT = NSPIN * (IH-1) + ISPIN

            CALL RINIT(IRMD*LMPOTD,PSHIFTLMR)
            CALL RINIT(IRMD,PSHIFTR)
            DO IR = 1,IRC1
               PSHIFTLMR(IR,1) = VSHIFT
            ENDDO

            IF (KSHAPE.EQ.0) THEN ! ASA

               DO IR = 1,IRC1
                  VISP(IR,IPOT) = VISP(IR,IPOT) + PSHIFTLMR(IR,1)
               END DO

            ELSE                ! Full-potential

! 
               CALL CONVOL(IMT1,IRC1,NTCELL(IH),
     &              IMAXSH(LMPOT),ILM,IFUNM,LMPOT,GSH,
     &              THETAS,THESME,0.d0,RFPI,
     &              RMESH(1,IH),PSHIFTLMR,PSHIFTR,LMSP)


               DO IR = 1,IRC1
                  VISP(IR,IPOT) = VISP(IR,IPOT) + PSHIFTLMR(IR,1)
               ENDDO


               DO LM = 2,LMPOT
                  DO IR = IRMIN1,IRC1
                VINS(IR,LM,IPOT)=VINS(IR,LM,IPOT)+PSHIFTLMR(IR,LM)*RFPI
                  ENDDO
               ENDDO

            END IF              ! (KSHAPE.EQ.0)


         END DO

      END DO



      END SUBROUTINE !bshift_ns
      
      
      END MODULE !MOD_MAIN0
      
