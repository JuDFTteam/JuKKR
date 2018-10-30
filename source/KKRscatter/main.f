c     MAIN PROGRAM 
c
c     This program enables us to investigate the scattering process as well as 
c     the spin-dependent transport using Boltzmann approach
c
c     reconstructed by Swantje
c     reconstructed by Long 17.09.13
c     updated and commented by Long 13.05.2016
c ************************************************************************
c     
c     Options OPT implemented by Long, updated by Bernd, Guillaume
c
c
c     FERMIOUT     :   writedown data for FS calculation 
c                  :   also need option OPERATOR for spin expectation value
c
c     NO-BREAK     :   write green functions for the host directly
c                      implemented by Bearnd Zimmermann
c
c     BREAK-1 or 2 :   write data to the files and then calculate
c                      green functions for the host in parallel
c
c     WRTGREEN     :   calculate green functions for the host
c                      required fullBZ (or SETGROUPnoinvers) in TEST option
c
c     GREENIMP     :   calculate green functions/tmatrix for the impurity
c
c
c     NEWSOSOL     :   new solver for spin-orbit coupling 
c                      implemented following David'code by Long
c 
c     OPERATOR     :   calculate the torque and spin flux
c                      implemented by Guillaume Geranton
c
******************************************************************************* 
c     Additional options 
c
c     FERMI-VL     :   the fermi-velocity for the k-points given in kpoints_FS 
c                      writedown wavefuntcion coefficients without rotation
c                      writedown spin expectation value 
c
c     SOFIELD      :   calculate spin-orbit field by perturbation
c                      for Dyakonov-Perel mechanism
c
c     BAND-STR     :   calculate the band structure 
c
c     ROTCOEFF     :   rotate wavefunction coefficients with Gauge condition
c
c     LIFETIME     :   calculate lifetime 
c                     
c     HALLCOND     :   conductivity using Boltzmann approach
c                      
******************************************************************************* 
c

      PROGRAM TBKKR
      IMPLICIT NONE
      include 'inc.p' 
      include 'inc.cls'
C 
      INTEGER LMMAXD,LMPOTD
      PARAMETER (LMMAXD= (LMAXD+1)**2,LMPOTD= (LPOTD+1)**2)
      INTEGER LMXSPD
      PARAMETER (LMXSPD= (2*LPOTD+1)**2)
      INTEGER NPOTD
      PARAMETER (NPOTD=NSPIND*NATYPD)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
      INTEGER LASSLD
      PARAMETER (LASSLD=4*LMAXD)
      INTEGER LM2D
      PARAMETER (LM2D = (2*LMAXD+1)**2)
      INTEGER NCPAIRD
      PARAMETER(NCPAIRD=100)
      INTEGER NKPOID
      PARAMETER(NKPOID=9000)
      DOUBLE PRECISION ONE,ONEM,TWO
      PARAMETER (ONE = 1.D0,ONEM = -1.D0,TWO=2.D0)
      DOUBLE COMPLEX CONE,CONEM,CZERO,CI
      PARAMETER (CONE  = ( 1.0D0,0.0D0),
     +           CONEM = (-1.0D0,0.0D0),
     +           CZERO = ( 0.0D0,0.0D0),
     +           CI    = ( 0.0D0,1.0D0))
C     ..
C     .. Local Scalars ..
c
      LOGICAL       
     +     COMPLX,
     +     LINIPOL,LSTART,LRHOSYM,
     +     OPT,
     +     TEST,
     +     TRANS_TMAT
C
      CHARACTER*40         I13,I19,I25,I40,I12
      CHARACTER*9          I16,I17,I18,I21,I22,I23
C
C     .. DOUBLE PRECISION ....
C
      DOUBLE PRECISION 
     +     ALATC,BLATC,CLATC,       ! lattice constants (in a.u.)
     +     C,                       ! speed of light
     +     ABASIS,BBASIS,CBASIS,    ! scaling factors for rbasis
     +     E3,
     +     E1,E2IN,E2,              ! energies needed in EMESHT
     +     EFERMI,                  ! Fermi energy (for scaling DOS)
     +     HFIELD,                  ! external magnetic field, for 
                                    ! initial potential shift in 
                                    ! spin polarised case
     +     STIME,STIME0,STIME1,     ! real time
     +     STIMEK,STIMEREF,TMTIME,TMTIME1
      DOUBLE PRECISION 
     +     ETIME,ETIME0,ETIME1,
     +     TK,                      ! temperature
     +     KB,                      ! Boltzmann constant
     +     VCONST,                  ! potential shift
c
     +     ALAT,ALATNEW,EIN2,ESHIFT,FCM,MIX,MIXING,MIXING0,QBOUND,
     +     CHRGNT,DENEF,E2SHIFT,EFCTOR,
     +     CHRGOLD,DENOLD,DENI,
     +     F,FF,PI,RTEMP,RMSAVM,RMSAVQ,RMSAV1,SIGN,X,
     +     RFOURIER                  ! cutting radius for Fourier-transf.
C
C     .. DOUBLE PRECISION ARRAYS ....
C
      DOUBLE PRECISION
     +     BRAVAIS(3,3),            ! bravais lattice vectors
     +     RECBV(3,3),              ! reciprocal basis vectors 
     +     MTFAC(NATYPD),           ! scaling factor for radius MT
c
     +     A(NATYPD),B(NATYPD),     ! cosntants for exponential r mesh
     +     CLEB(NCLEB,2),           ! GAUNT coefficients (GAUNT)
     +     DRDI(IRMD,NATYPD),       ! derivative dr/di
     +     DROR(IRMD,NATYPD),       ! logarithmic derivative (1/r)*(dr/di)
     +     ECORE(20,NPOTD),         ! core states
     +     R(IRMD,NATYPD),          ! radial r mesh (in units a Bohr)
     +     RBASIS(3,NAEZD+NEMBD),   ! position of atoms in the unit cell
                                    ! in units of bravais vectors
     +     RCLS(3,NACLSD,NCLSD),    ! real space position of atom in cluster
     +     RMT(NATYPD),             ! Muffin-Tin-radius
     +     RMTNEW(NATYPD),          ! adapted MT radius
     +     RMTREF(NREFD),           ! adapted MT radius for reference pot
     +     RR(3,0:NRD),             ! set of real space vectors (in a.u.)
     +     RS(IRMD,0:LMAXD,NATYPD), ! r mesh for relat. calc.
     +     RWS(NATYPD),             ! Wigner Seitz radius
     +     S(0:LMAXD,NATYPD),
     +     VREF(NREFD)
      DOUBLE PRECISION
     +     THETAS(IRID,NFUND,NCELLD), ! shape function
                                      !         ( 0 outer space
                                      ! THETA = (
                                      !         ( 1 inside WS cell
                                      ! in spherical harmonics expansion
     +     VBC(2),                  ! potential constants
     +     VINS(IRMIND:IRMD,        ! nonsperical potential
     +                 LMPOTD,NSPOTD),       
     +     VISP(IRMD,NPOTD),        ! spherical input potential
     +     dVdr(IRMD,NPOTD),        ! radial derivative of the spherical input potential
c     +     VONS(IRMD,LMPOTD,NPOTD), ! output potential
     +     WG(LASSLD),              ! integr. weights for Legendre pol.
                                    ! (GAUNT2)
     +     YRG(LASSLD,0:LASSLD,
     +                0:LASSLD),    ! spherical harmonics (GAUNT2)
     +     Z(NATYPD)               ! nuclues charge
c
      DOUBLE PRECISION
     +     GSH(NGSHD),
     +     RHOC(IRMD,NPOTD),        ! core charge density
     +     RIJ(IJD,3),
     +     RHYPF(150,NPOTD),HYPSUM(10,NPOTD),  ! for hyperfine fields
c    p. zahn, 10.3.99,
c     +     WORK(IRMD,LMPOTD,NATYPD),
     +     WORK,
     +     WTYR(IJD,LMPOTD),
     +     YR(IJD,LMPOTD),
     +     W((LMMAXD+1)**2*NATYPD)
c
c   new madelung variables
       INTEGER NGMAX,NRMAX,NSHLG,NSHLR,IS1
       DOUBLE PRECISION VOLUME0,RCUTZ,RCUTXY,TEMP(3)

 
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DOUBLE PRECISION SMALL
      DOUBLE COMPLEX DELTA,T
      DOUBLE COMPLEX   LSM(LMMAXD*2,LMMAXD*2),
     +               PHASE_SHIFT(0:LMAXD,NATYPD)
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     .. INTEGER ....
      integer
     +     INTERVX,INTERVY,INTERVZ, ! number of intervals in x-,y-,z-direction
                                    ! for k-net in IB of the BZ
     +     ICC,                     ! center of cluster for output of GF
c                                   ! for DOS calc. icc = 0 
     +     ICLS,                    ! counter for cluster
     +     ICST,                    ! number of Born approximation
     +     IELAST,                  ! number of complex energy values
     +     IEND,                    ! number of nonzero gaunt coeffizients
     +     IFILE,                   ! unit specifier for potential card
     +     IINFO,IINFO1,            ! info's about r-mesh in CALRMT
     +     IN,                      ! I/O channel for system potential
     +     SINN,SOUT,RIN,ROUT,      ! in/out of GF, WF, T-matrices
     +     INS,                     ! 0 (ASA), 1(Full Potential)
     +     INSREF,                  ! INS for reference pot. (usual 0)
     +     IOWFCT,IOWREF,           ! unit specifier for wave functions
     +     IOPREF,IOPSYS,           !                    potentials
     +     IOTREF,IOTSYS,           !                    T matrices
     +     IOGREF,                  !                    GF(ref sys.)
     +     IPE,IPF,IPFE,            ! not real used, IPFE should be 0
     +     IGF,IMIX,IPOTOU,IPRCOR,
     +     IRM,IRNUMX,
     +     KSCOEF,                  ! 0,1: read shell structure from
     +     ISHIFT,ITCCOR,ITCLST,ITDBRY,
     +     KHFELD,                  ! 0,1: no / yes external magnetic field
     +     KSHAPE,                  ! exact treatment of WS cell
     +     KVREL,                   ! 0,1 : non / scalar relat. calculation
     +     KWS,                     ! 0 (MT), 1(ASA)
     +     KMT,                     ! scaling of RMT with MTFAC
c                                   ! 0: RMT from crystal structure 
c                                   ! 1: RMT from crystal structure
c                                   !    scaled by RMTFAC
c                                   ! 2: RMT = RMTFAC
c                                   ! 3: RMT from ref. pot. card
     +     KCOR,KEFG,KFROZN,
     +     KHYP,KHYPO,KPRE,KTE,KVMAD,KXC,KFORCE
      INTEGER
     +     LATT,                    ! 0: scu  1: fcc  2: bcc
c                                   ! 3: tet  4: bct  5: hex
c                                   !                 8: ort
     +     LMAX,                    ! maximum l component in 
c                                   ! wave function expansion
     +     LPOT,                    ! maximum l component in 
c                                   ! potential expansion
     +     LMMAX,LMPOT,MD,
     +     MAXMESH,                 ! number of k-mesh sets
     +     MMIN,MMAX,               ! min. max. eigenvalue
     +     M2,                      ! maximum number of bands
c
     +     NMESH,                   ! number of k set 
     +     NAEZ,                    ! number of atoms in unit cell
     +     NATYP,                   ! number of kinds of atoms in unit cell
     +     NCLS,                    ! number of reference clusters
     +     NCL_IMP,                ! number of atoms in the cluster
     +                              !      taken into account in the lifetime calculation
     +     IMPLAYER,                ! number of the layer in which the
                                    !      impurity is located   
     +     ILAYERS,                   ! number of the layers without Vc
     +     NEMB,                    ! number of 'embedding' positions
     +     NEMBZ,                   ! inequiv. 'embedding' positions 
     +     NINEQ,                   ! number of ineq. positions in  unit cell
     +     NLAYER,                  ! number of principal layer
     +     NPNT1,NPNT2,NPNT3,       ! number of E points (EMESHT)
     +     NPOL,                    ! number of Matsubara Pols (EMESHT)
     +     NR,                      ! number of real space vectors rr
     +     NREF,                    ! number of diff. ref. potentials
     +     NSPIN,                   ! counter for spin directions
     +     NSPO,                    ! counter for spin directions for
                                    !                  spin-orbit coupling 
     +     NSPOH,                   ! counter for spin-orbit coupling of
     +                              !                             the host 
     +     IRMINSO,                 ! number of points for which SP-O
     +     NSTEPS,                  ! number of iterations
     +     NZ,                      ! number of atoms at centers of inversion
     +     I,J,IK,                  ! counters
     +     IPOT,ISPIN,IS,IE1,NKTOT,
     +     I1,IC,IE,II,INV,I2,IET,
     +     IA,IATYP,IDIM,IH,IJEND,IP,IPARA,
     +     IR,IRC1,IRMIN1,IT,ITC,IWRIT,
     +     L,LM,LM1,LM2,LM1SO,LM2SO,
     +     N,N1,N2,NITMIX,NMIX,IOS,ISP,
     +     RF,IHANDLE,IHANDLE1,IHANDLES,NKDOS,NS,NS1,IATOM,IEPRINT,
     +     INATYP,IIRM,ILMPOT,ILMC,IWF,ILC,IMC,IJ,IENT,ENTG,
     +     NLMAX,NKMMAX,NMUEMAX,NKMPMAX,NKMAX,LINMAX
      INTEGER NPAN_LOG,NPAN_EQ,NCHEB
      DOUBLE PRECISION R_LOG,THETA,PHI
C     
C     .. INTEGER ARRAYS ....
C
      INTEGER
     +     EQINV(NAEZD),            ! site equiv. by invers. symmetry
     +     INIPOL(NATYPD),          ! initial spin polarisation
     +     IXIPOL(NATYPD),          ! constraint of spin pol.
     +     KAOEZ(NAEZD+NEMBD),      ! kind of atom at site in elem. cell
c
c     ..  cluster arrays
c
     +     ATOM(NACLSD,NAEZD),     ! atom at site in cluster
     +     CLS(NATYPD),             ! cluster around atom
     +     NACLS(NCLSD),            ! number of atoms in cluster
     +     EZOA(NACLSD,NAEZD),      ! EZ of atom at site in cluster
c
     +     ICLEB(NCLEB,4),          ! pointer array
     +     IFUNM(NATYPD,LMXSPD),
     +     IMT(NATYPD),             ! r point at MT radius
     +     IPAN(NATYPD),            ! number of panels in non-MT-region
     +     IRC(NATYPD),             ! r point for potential cutting
     +     IRCUT(0:IPAND,NATYPD),   ! r points of panel borders
     +     IRMIN(NATYPD),           ! max r for spherical treatment
     +     IRNS(NATYPD),            ! number r points for non spher. treatm.
     +     IRWS(NATYPD),            ! r point at WS radius
     +     ITITLE(20,NPOTD)         ! title line in potential card
      INTEGER
     +     JEND(LMPOTD,             ! pointer array for icleb()
     +             0:LMAXD,0:LMAXD),
     +     KFG(4,NATYPD),
     +     LCORE(20,NPOTD),         ! angular momentum of core states
     +     LLMSP(NATYPD,NFUND),     ! lm=(l,m) of 'nfund'th nonvanishing
                                    ! component of non-spherical pot.
     +     LMSP(NATYPD,LMXSPD),     ! 0,1 : non/-vanishing lm=(l,m) component
                                    ! of non-spherical potential
     +     LOFLM(LM2D),             ! l of lm=(l,m) (GAUNT)
     +     LMXC(NATYPD),
     +     NCORE(NPOTD),            ! number of core states
     +     NFU(NATYPD),
     +     NSHELL(0:NSHELD),        ! number of atoms in shell
c                                   ! nshell(0) = number of shells
     +     NTCELL(NATYPD),          ! index for WS cell 
     +     NTCELLR(NREFD),          ! index for WS cell of ref. pot.
     +     REFPOT(NATYPD+NEMBD),    ! ref. pot. card  at position
c     
     +     ILM(NGSHD,3),
     +     IMAXSH(0:LMPOTD),ATOMIMP(NATOMIMPD),
     &     ICOUPLMAT(NAEZD,NAEZD)

      INTEGER IATCONDL(NCPAIRD),IATCONDR(NCPAIRD),NCONDPAIR,IEGFOUT ! for conductivity
C
C     .. DOUBLE COMPLEX ....
c     
      DOUBLE COMPLEX DF,E,EK,CFAC,CFACL
c
C
C     .. DOUBLE COMPLEX ARRAYS ....
c     
      DOUBLE COMPLEX 
c
     +     ALPHA(0:LMAXD,NSPD,NATYPD),
     +     AR(NSPD*LMMAXD,NSPD*LMMAXD,NATYPD),
     +     DEZ(IEMXD),              ! weights for complex energy integration
     +     EZ(IEMXD),               ! energy points for compl.integration
c     +     PNSSO(NSPD*LMMAXD,NSPD*LMMAXD,IRMD,2,NATYPD),
c     +     GINP(NACLSD*LMMAXD,LMMAXD,NCLSD,3), ! cluster GF(ref syst.)
     +     PZI(IRMD,0:LMAXD,NATYPD,NSPIND),
     +     QZI(IRMD,0:LMAXD,NSPIND),
     +     SZI(IRMD,0:LMAXD,NSPIND),
     +     FZI(IRMD,0:LMAXD,NATYPD,NSPIND),
     +     WEZ(IEMXD)               ! modified weights for energy integration
                                    ! = -2.0/pi*DEZ
      INTEGER NATOMIMP
      DOUBLE PRECISION RCLSIMP(3,NATOMIMPD)
      INTEGER, ALLOCATABLE :: IPANIMP(:),IRCUTIMP(:,:),
     +        IRMINIMP(:),IRWSIMP(:)
      DOUBLE PRECISION, ALLOCATABLE :: RIMP(:,:),ZIMP(:),
     +        THETASIMP(:,:,:),
     +        VISPIMP(:,:),
     +        VINSIMP(:,:,:)
      DOUBLE COMPLEX, ALLOCATABLE ::  DTMTRX(:,:)
      INTEGER IHOST,HOSTIMP(NATYPD)
      DOUBLE COMPLEX 
     +     TMATLL(NSPD*LMMAXD,NSPD*LMMAXD,NAEZD),
     +     TMATLL1(NSPD*LMMAXD,NSPD*LMMAXD,NATYPD),
     +     TMATLLE(NSPD*LMMAXD,NSPD*LMMAXD,NATYPD,3),
     +     TSTAR(NSPD*LMMAXD,NSPD*LMMAXD),
     +     TREFLL(LMMAXD,LMMAXD,NREFD)
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c    Lines with a !new1 in the end are added after 
c    20.9.99
      INTEGER NLAY,                               !new1 
     +        NLBASIS,NRBASIS,                    !new1 
     +        NLEFT,NRIGHT,IER                    !new1
      DOUBLE PRECISION                            !new1 
     +     ZPERLEFT(3),ZPERIGHT(3),               !new1
     +     TLEFT(3,NEMBD),TRIGHT(3,NEMBD)         !new1
      double precision tx,ty,tz
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      DOUBLE COMPLEX BESSJW1(0:LMAXD+1),BESSYW1(0:LMAXD+1),
     &               HANKWS1(0:LMAXD+1),BESSJW2(0:LMAXD+1),
     &               BESSYW2(0:LMAXD+1),HANKWS2(0:LMAXD+1)

ccccccccccccccccccc
      LOGICAL STRCON,LINTERFACE,LCARTESIAN,VACFLAG(2)        !new1
      CHARACTER*50 HOME
      CHARACTER*80 UIO
      DOUBLE PRECISION  E_FS_IN,DELTA_K_FS,VOLBZ,VSHIFT

      INTEGER       :: NKSUM,NK,IROTMAX

      DOUBLE COMPLEX  DELTALAMDAK(NKPOID,3),DELTALAMDADELE(NKPOID,3),
     +                PMUP(NKPOID),PMDOWN(NKPOID),SPERTURB(NKPOID,2,4)
      DOUBLE COMPLEX, ALLOCATABLE :: GINP(:,:,:)   ! reference Green func.
c      DOUBLE COMPLEX, ALLOCATABLE :: GINP(:,:,:,:)  
      DOUBLE COMPLEX, ALLOCATABLE :: PNSSO(:,:,:,:,:) ! right redial wavefunction
c     Left regular solutions are used to represent the Green function that we need for Kubo
      DOUBLE COMPLEX, ALLOCATABLE :: PNSSO_LEFT(:,:,:,:,:) 
c
c
C     .. Common blocks ..
c
      COMMON /CMACH/STRCON
      COMMON /RFOUR/ RFOURIER
c
c     .. EXTERNAL STATEMENT
c
      EXTERNAL   
     +     CLSGEN99,DELTAMAT,
     +     RTEST,ROPT,TEST,OPT,
     +     RINPUT1,TMATRX01,
     +     TMWRIT,WFWRIT,
     +     BRYDBM,CHEACC,CINIT,CONVOL,CORELB,CSEND,DAXPY,
     +     DCOPY,ECOUB,EMESHT,EPOTINB,ESPCB,ESPVB,ETOTB1,
     +     GAUNT,GAUNT2,KLOOPZ1,MIXSTR,MTZERO,RCSTOP,RHOLMB01,
     +     RHONSB01,RHOTOTB,RINIT,RITES,SHAPE,SPHERE,
     +     STARTB1,TMREAD,VINTRAS,VMADEL,VXCLM,WFREAD
c
C     .. Intrinsic Functions ..
      INTRINSIC ACOS,DABS,DATAN,DIMAG,DMAX1,
     +          DMIN1,DSIGN,DSQRT,LOG,MAX,REAL,SIN,ZABS,ZEXP
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DCLOCK,VPOT,BFIELD,SCALEB
      DOUBLE COMPLEX EXPIDL
      INTEGER MAPBLOCK,MYNODE,NUMNODES
      EXTERNAL DCLOCK,EXPIDL
C     ..
      CHARACTER*24 TXC(5)
C     ..
C     .. Data statements ..
c
      DATA 
     +     SMALL    / 1.0D-8       /,
     +     KB       / 0.6333659D-5 /
c
      DATA I16 /'tmat.ref '/,
     +     I17 /'gll.ref  '/,
     +     I18 /'wfct.ref '/,
     +     I21 /'tmat.sys '/,
     +     I22 /'gll.sys  '/
c
C     .. Save statement ..
c
      SAVE
c      ALLOCATE(GINP(NACLSD*LMMAXD,LMMAXD,NCLSD,3))
      ALLOCATE(GINP(NACLSD*LMMAXD,LMMAXD,NCLSD))
      ALLOCATE(PNSSO(NSPD*LMMAXD,NSPD*LMMAXD,IRMD,2,NATYPD))
      ALLOCATE(PNSSO_LEFT(NSPD*LMMAXD,NSPD*LMMAXD,IRMD,2,NATYPD))
      
c
c ---------------------------------------------------CONSTANTS -----------
      PI = 4.0D0*DATAN(1.0D0)
      EFERMI = 0.0d0
      NITMIX = 0
      GINP=0d0
      STIME0 = DCLOCK()

c ------------------------------------------------------------------------
c read in data from inputcard
c ------------------------------------------------------------------------

       CALL RINPUT99(ALAT,RBASIS,ABASIS,BBASIS,CBASIS,
     &           ALATC,BLATC,CLATC,LATT,CLS,NCLS,
     &           E1,E2,TK,NPOL,NPNT1,NPNT2,NPNT3,ESHIFT,
     +           DELTA_K_FS,E_FS_IN,                
     &           ITCLST,NSTEPS,IMIX,MIXING,QBOUND,FCM,ITDBRY,
     &           IRNS,RMTREF,NTCELL,NAEZ,NEMB,KAOEZ,EQINV,IRM,Z,
     &           NINEQ,NREF,NTCELLR,
     &           ICST,IFILE,IGF,INS,INSREF,IPE,IPF,IPFE,
     &           KCOR,KEFG,KFROZN,KHFELD,KHYP,KPRE,KSHAPE,KTE,
     &           KFG,KVMAD,KVREL,KWS,KXC,LMAX,LMMAX,LMPOT,LPOT, 
     &           NATYP,NSPIN,NSPO,NSPOH,NCL_IMP,IMPLAYER,ILAYERS,
     &           LMXC,TXC,KSCOEF,ICC,REFPOT,
     &           IPOTOU,IPRCOR,IRNUMX,ISHIFT,ITCCOR,
     &           MD,INTERVX,INTERVY,INTERVZ,
     &           HFIELD,COMPLX,
     &           KMT,MTFAC,VBC,VCONST,LINIPOL,INIPOL,IXIPOL,LRHOSYM,
     &           MMIN,MMAX,SINN,SOUT,RIN,ROUT,M2,I12,I13,I19,I25,I40,
     &           NLBASIS,NRBASIS,NLEFT,NRIGHT,ZPERLEFT,ZPERIGHT,    
     &           TLEFT,TRIGHT,LINTERFACE,RCUTZ,RCUTXY,
     &           NPAN_LOG,NPAN_EQ,NCHEB,R_LOG)  



c ------------------------------------------------------------------------
c cross check read in data (can be commented out)
c ------------------------------------------------------------------------

       CALL TESTDIM(NSPIN,NSPO,NSPOH,NAEZ,NEMB,NATYP,LMAX,IRM,INS,
     &              INSREF,NREF,M2,IRNS,NCLS,NLAYER,LATT)     

c ------------------------------------------------------------------------
c read in impurity coordinate in scoef file
c ------------------------------------------------------------------------

      IF (INS.GT.0) OPEN (19,FILE=I19,status='old',FORM='formatted')

      OPEN (25,FILE=I25,status='unknown',FORM='formatted')
      IHOST=0
      IF ((OPT('WRTGREEN').OR.OPT('GREENIMP'))
     &     .AND..NOT.OPT('KUBO    ')) THEN
c      REWIND(25)
      READ(25,*) NATOMIMP
      DO I=1,NATOMIMP
       READ(25,*) (RCLSIMP(J,I),J=1,3),ATOMIMP(I)
      ENDDO
      DO 125 I=1,NATYPD
       DO J=1,NATOMIMP
        IF (ATOMIMP(J).EQ.I) THEN
         IHOST=IHOST+1 
         HOSTIMP(IHOST)=ATOMIMP(J)
         GO TO 125
        ENDIF
       ENDDO
  125  CONTINUE
      ENDIF

      IF (OPT('WRTGREEN')) THEN
       OPEN (UNIT=58,FILE='green_host',FORM='FORMATTED')
      ENDIF

      IF (OPT('BAND-STR')) THEN
       OPEN (61,FILE='bandstruct.data',FORM='formatted',
     +       STATUS='unknown')
       OPEN (62,FILE='reallarge.data',FORM='formatted',
     +       STATUS='unknown')
      ENDIF
 

c
c ---> dependencies of input parameters
c
      IF (TEST('c/a     ')) CLATC = CLATC*ALATC
      ALAT = ALATC
      MD = intervx

      IF (ICC.NE.0 .AND. IGF.EQ.0) THEN
        IGF = 1
      END IF

      IF (ICC.EQ.0 .AND. IGF.NE.0) THEN
        ICC = -1
      END IF


      IF(ICC.GT.NINEQ) THEN
        write(6,*) 'The choice of ICC (=',ICC,') is not possible.'
        write(6,*) 'ICC is changed to the equivalent atom ',
     +       EQINV(ICC),'.'
        ICC = EQINV(ICC)
      ENDIF

c ---> end of dimension and input data check
c
      write(6,2100) 
c
c ------------------------------------------------------------------------
c
c determines the direct and reciprocal space basis vectors
c and a number of direct space vectors RR(3,0:NRD)
c
c reads the bravais vectors directly from the input card and returns
c the same things as LATTIX: BRAVAIS, RECBV, RR, NR.
c BRAVAIS in units of au/alat
c ------------------------------------------------------------------------

      CALL LATTIX99(ALATC,BRAVAIS,RECBV,RR,NR)

c ------------------------------------------------------------------------
c normalization of basis vectors
c
c in the inputcard RBASIS are the basis vectors in the
c units of the bravais vectors and au/alat 
c
c now RBASIS will be the basis vectors in units of au/alat 
c in (xyz) reference
c ------------------------------------------------------------------------

      CALL SCALEVEC(RBASIS,ABASIS,BBASIS,CBASIS,NLBASIS,
     &              NRBASIS,NLEFT,NRIGHT,ZPERLEFT,ZPERIGHT,            
     &              TLEFT,TRIGHT,LINTERFACE,NAEZ,NEMB,BRAVAIS,KAOEZ)

c ------------------------------------------------------------------------
c determine the configuration of the clusters in the lattice
c ------------------------------------------------------------------------

      CALL CLSGEN99(NAEZ,RR,NR,RBASIS,KAOEZ,Z,CLS,NACLS,REFPOT,ATOM,
     &              EZOA,NLBASIS,NRBASIS,NLEFT,NRIGHT,ZPERLEFT,
     &              ZPERIGHT,TLEFT,TRIGHT,RCLS,RCUTZ,RCUTXY,LINTERFACE,
     &              ALAT,NAEZD,NATYPD,NEMBD,NPRINCD,NRD,NACLSD,NCLSD) 

      IF (OPT('full inv')) THEN
      WRITE(6,*) '>>>>>>>>>>> ================ <<<<<<<<<<<'
      WRITE(6,*) '>>    F U L L   I N V E R S I O N     <<'
      WRITE(6,*) '>>>>>>>>>>> ================ <<<<<<<<<<<'
      END IF
      call flush(6)

      IOPSYS = 13   ! potential file

      OPEN(IFILE,FILE=I13,STATUS='OLD',FORM='FORMATTED')

      NSHELL(0) = 1
      LSTART = .true.

c ------------------------------------------------------------------------
c calculate the Gaunt's coefficients
c ------------------------------------------------------------------------
        
      CALL GAUNT2(WG,YRG)            
      CALL GAUNT(LMAX,LPOT,WG,YRG,CLEB,LOFLM,ICLEB,IEND,JEND) 

c ------------------------------------------------------------------------
c calculate the spherical harmonics 
c ------------------------------------------------------------------------

      IF (KXC.LT.3) THEN 
          CALL SPHERE_NOGGA(LPOT,YR,WTYR,RIJ,IJEND,IJD,W)
      ELSE
          CALL SPHERE(LPOT,YR,WTYR,RIJ,IJEND,IJD)
      END IF

      IN=IFILE

      E2IN = E2

      STIME0 = DCLOCK()

      IF (IPE.EQ.1) IPF = IPFE

      IINFO = 1

c ------------------------------------------------------------------------
c read in potential and shape function in case of full-potential 
c ------------------------------------------------------------------------

      call STARTB1(IN,IPF,IPFE,IPE,KVREL,KWS,KHFELD,LMAX,
     +       1,NATYP,ALATNEW,RMTNEW,RMT,    
     +       ITITLE,HFIELD,IMT,IRC,VCONST,
     +       INS,IRNS,LPOT,NSPIN,VINS,IRMIN,KSHAPE,NTCELL,
     +       IRCUT,IPAN,THETAS,IFUNM,NFU,LLMSP,LMSP,E2IN,
     +       VBC,C,DROR,RS,S,VISP,RWS,ECORE,LCORE,NCORE,DRDI,
     +       R,Z,A,B,IRWS,INIPOL,IINFO)

c ------------------------------------------------------------------------
c setup of GAUNT coefficients C(l,m;l',m';l'',m'') for all 
c nonvanishing (l'',m'')-components of the shape functions THETAS
c ------------------------------------------------------------------------

      IF (KSHAPE.NE.0) 
     +     CALL SHAPE(LPOT,NATYP,GSH,ILM,IMAXSH,LMSP,NTCELL,WG,YRG)

c ------------------------------------------------------------------------
c add potential if needed
c ------------------------------------------------------------------------

c constant potential shift by hand 
c only up to muffin-tin sphere
c      DO I1=1,NPOTD
c       IF (I1.EQ.IMPLAYER) THEN
c        DO IR=1,IRCUT(1,I1)
c         DO ISPIN=1,NSPIN
c          IPOT = NSPIN*(I1-1)+ISPIN
c          VISP(IR,IPOT)=VISP(IR,IPOT)+0.01
c         ENDDO
c        ENDDO
c       ENDIF
c      ENDDO
c      VCONST = 0.D0

c for full-potential, convoluted with shape function
c      VSHIFT=-0.18
c      CALL POTENSHIFT1(VISP,VINS,NATYP,NSPIN,IRCUT,IRC,IRMIN,NTCELL,
c     +                 IMAXSH,ILM,IFUNM,LMSP,LMPOT,GSH,THETAS,R,KSHAPE,
c     +                 VSHIFT)

      IF (NPOL.EQ.0) EFERMI = E2IN
c
c --->  test E2IN and EMU (EMESHT)
c
      IF (DABS(E2IN-E2).GT.1d-10.AND.NPOL.NE.0) THEN
        write(6,*) 'EFERMI(file = ',in,' ) .NE. E2(emesht)',
     +         ' .AND. NPOL.NE.0 <<<<<<'
        write(6,FMT='('' E2IN = '',f14.6,''  E2 = '',f14.6)') 
     +         E2IN,E2
c
        E2 = E2IN 
c
        write(6,FMT='('' FERMI ENERGY = '',f12.5)') E2
c     
c --->  doesn't read TREFLL, GREFLL, TMATLL1 and WF's from files
c     
        SINN = 0
        RIN = 0
        write(6,FMT='('' SINN, RIN     = '',2I6)') SINN,RIN
          
      END IF                      ! (DABS(E2IN-E2).GT.1d-10.AND.NPOL.NE.0)

c ------------------------------------------------------------------------
c core charge
c ------------------------------------------------------------------------

          DO 310 IP = 1,NATYP
            CALL CORELB(IPF,ITCCOR,KHYPO,KCOR-1,IPRCOR,IRNUMX,IP,IRM,
     +                  NSPIN,RHOC,IMT,KSHAPE,VISP,RWS,ECORE,DRDI,R,Z,A,
     +                  B,IRWS,LMXC,KFG,RHYPF,HYPSUM)
 310      CONTINUE

c ------------------------------------------------------------------------
c write potential to the file fort.11
c important if the potential is modified
c ------------------------------------------------------------------------

          REWIND (11)
          CALL RITES(11,1,NATYP,NSPIN,Z,ALATC,RMT,RMTNEW,RWS,
     +               ITITLE,R,DRDI,VISP,IRWS,A,B,TXC,KXC,INS,IRNS,
     +               LPOT,VINS,QBOUND,IRC,KSHAPE,E2IN,VBC,ECORE,
     +               LCORE,NCORE)
          REWIND (11)

      IFILE = 11

c ------------------------------------------------------------------------
c creates a mesh in the complex energy plane
c ------------------------------------------------------------------------

      CALL EMESHT(EZ,DEZ,IELAST,E1,E2,TK,NPOL,NPNT1,NPNT2,NPNT3)

      WRITE (6,FMT=9200) IELAST
      WRITE (6,FMT=9210) TK
      WRITE (6,FMT=9220) NPOL,NPNT1,NPNT2,NPNT3

      IF (IELAST.GT.IEMXD) then
        write(6,*) 'Please, change the parameter iemxd in',
     +         ' inc.p to a value greater equal ',ielast
        CALL RCSTOP('IEMXD   ')
      ENDIF

c ------------------------------------------------------------------------
c scales the weight WEZ for the energy integration and determine
c number of necessary k-mesh sets
c ------------------------------------------------------------------------

      DO IE = 1,IELAST
        IF (TEST('e-net   ')) 
     +       WRITE (6,FMT='(2(f14.5,f10.5))') EZ(IE),DEZ(IE)
        WEZ(IE) = -2.d0/PI*DEZ(IE)
      END DO

      MAXMESH = 1
      IF (.NOT.TEST('fix mesh')) THEN
        DO IE = 1,IELAST
          IF (DIMAG(EZ(IE)).GT.1.D-8) THEN
            N=1+log(DIMAG(EZ(IE))/DIMAG(EZ(IELAST)))/log(2.0d0)
            MAXMESH=MAX(MAXMESH,N)
          END IF
        END DO
        IF (NPOL.EQ.0) MAXMESH = 1 
        write (6,*) 'MAXMESH : ',maxmesh
      END IF                      

c ------------------------------------------------------------------------
c     |                                |
c     |   BEGIN DO LOOP OVER ENERGIES  |
c     |                                |
c ------------------------------------------------------------------------

      STIME = DCLOCK()
      TMTIME = 0.0d0

      IF (OPT('FERMI-VL')) THEN
       IELAST=3
       EZ(2)=EZ(1)+0.001
       EZ(3)=EZ(1)-0.001
      ENDIF 

      IF (OPT('FERMIOUT').AND..NOT.OPT('ONLY_EF ')) THEN
       IELAST=3
       EZ(2)=EZ(1)+CMPLX(1.0D-03,0.0D0)
       EZ(3)=EZ(1)-CMPLX(1.0D-03,0.0D0)
      ELSEIF (OPT('FERMIOUT')) THEN
       IELAST=1
      ENDIF 

c ------------------------------------------------------------------------
c read in impurities potential and shape function in case of FP
c require shapefun_imp, potential_imp
c ------------------------------------------------------------------------

      IF (OPT('GREENIMP')) THEN
       IF (IELAST.EQ.3) OPEN(UNIT=59,FILE='GMATLL_GES',FORM='FORMATTED')
       IF (IELAST.EQ.3) OPEN(UNIT=60,FILE='green_host',FORM='FORMATTED')
        ALLOCATE (IPANIMP(NATOMIMP))
        ALLOCATE (IRCUTIMP(0:IPAND,NATOMIMP))
        ALLOCATE (IRMINIMP(NATOMIMP))
        ALLOCATE (IRWSIMP(NATOMIMP))
        ALLOCATE (RIMP(IRMD,NATOMIMP))
        ALLOCATE (ZIMP(NATOMIMP))
        ALLOCATE (THETASIMP(IRID,NFUND,NATOMIMP))
        ALLOCATE (VISPIMP(IRMD,NATOMIMP*NSPIN))
        ALLOCATE (VINSIMP(IRMIND:IRMD,LMPOTD,NATOMIMP*NSPIN))
        CALL READIMPPOT(NATOMIMP,INS,IPF,IPFE,IPE,KWS,NSPIN,LPOT,
     &                 IPANIMP,THETASIMP,IRCUTIMP,IRWSIMP,KHFELD,
     &                 HFIELD,VINSIMP,VISPIMP,IRMINIMP,RIMP,ZIMP)

c ------------------------------------------------------------------------
c scale magnetic potential of impurity
c ------------------------------------------------------------------------

c      SCALEB=0.8
c      DO I1=1,NATOMIMP
c       ISPIN=1
c       IPOT = NSPIN* (I1-1) + ISPIN
c       IF (mod(IPOT,2).EQ.1) THEN
c        VPOT=0d0
c        BFIELD=0d0
c        DO IR=1,IRMD
c         VPOT=0.5*(VISPIMP(IR,IPOT)+VISPIMP(IR,IPOT+1))
c         BFIELD=SCALEB*0.5*(VISPIMP(IR,IPOT)-VISPIMP(IR,IPOT+1))
c         VISPIMP(IR,IPOT)=VPOT+BFIELD
c         VISPIMP(IR,IPOT+1)=VPOT-BFIELD
c        ENDDO
c        VPOT=0d0
c        BFIELD=0d0
c        DO IR=IRMIND,IRMD
c         DO LM1=1,LMPOTD
c          VPOT=0.5*(VINSIMP(IR,LM1,IPOT)+VINSIMP(IR,LM1,IPOT+1)) 
c         BFIELD=SCALEB*0.5*(VINSIMP(IR,LM1,IPOT)-VINSIMP(IR,LM1,IPOT+1))
c          VINSIMP(IR,LM1,IPOT)=VPOT+BFIELD
c          VINSIMP(IR,LM1,IPOT+1)=VPOT-BFIELD
c         ENDDO
c        ENDDO
c       ENDIF
c       ENDDO
c
      ENDIF ! IF ((OPT'GREENIMP'))

      DO 360 IE = 1,IELAST

        IF(TEST('ie      ')) then
          write(6,9170) IE,EZ(IE)
        END IF

        E = EZ(IE)
        DF = WEZ(IE)/REAL(NSPIN)

        NMESH = 1
        IF (DIMAG(EZ(IE)).GT.1.D-10) THEN
          IF (.NOT.TEST('fix mesh'))
     +      NMESH = 1+
     +      log(DIMAG(EZ(IE))/DIMAG(EZ(IELAST)))/log(2.0d0)
        END IF 

        STIMEREF = DCLOCK()
              
        IINFO = 0
        IF (IE.EQ.1) 
     +      IINFO = 1
        IINFO1 = IINFO

        E = EZ(IE)
                
c ------------------------------------------------------------------------
c calculate t-matrix of reference system
c ------------------------------------------------------------------------

        DO 20 I1 = 1,NREF
         VREF(I1)=8.d0
         CALL CALCTREF(E,VREF(I1),RMTref(I1),LMAX,LM1,
     &                 TREFLL(1,1,I1),LMAXD+1,LMMAXD)
         WRITE(6,*) "after CALCTREF" 
 20     CONTINUE          
              
c ------------------------------------------------------------------------
c calculate GF of reference system on clusters
c first search for atom at center
c ------------------------------------------------------------------------

        DO 90 ICLS=1,NCLS

          I1 = 1
          IC = 0

          DO WHILE (IC.EQ.0 .AND. I1.LE.NINEQ)
            IF (CLS(I1).EQ.ICLS) IC = I1
            I1 = I1 + 1
          END DO
          IF (I1.GT.NINEQ) THEN
            WRITE(*,*) 'I1.GT.NINEQ'
            WRITE(*,*) 'I1=',I1
            WRITE(*,*) 'NINEQ=',NINEQ
          END IF

          IF (IC.EQ.0) STOP 'Error in CLS(*) array in main'
          IF(TEST('flow    ')) 
     +         write(6,*) 'CLUSTER ',ICLS,' at ATOM ',IC

          CALL GLL95(LSTART,E,CLEB(1,2),ICLEB,LOFLM,IEND,
     +             TREFLL,ATOM(1,IC),KAOEZ,REFPOT,
     +             RCLS(1,1,ICLS),NACLS(ICLS),ALATC,0,
     +             GINP(1,1,ICLS))
c     +             GINP(1,1,ICLS,IE))
            
 90     CONTINUE              ! ICLS=1,NCLS

c ------------------------------------------------------------------------
c read potential and shape function from fort.11
c ------------------------------------------------------------------------

        IF (IE.EQ.1)
     +        WRITE (6,FMT=9080)  DCLOCK()-STIMEREF
              
        IN = IFILE 
        REWIND (IN)

        IINFO = 0
        IF (IE.EQ.1 .AND. TEST('pot head') ) IINFO = 1

        call STARTB1(IN,IPF,IPFE,IPE,KVREL,KWS,KHFELD,LMAX,
     +             1,NATYP,ALATNEW,RMTNEW,RMT, 
     +             ITITLE,HFIELD,IMT,IRC,VCONST,
     +             INS,IRNS,LPOT,NSPIN,VINS,IRMIN,KSHAPE,NTCELL,
     +             IRCUT,IPAN,THETAS,IFUNM,NFU,LLMSP,LMSP,E2IN,
     +             VBC,C,DROR,RS,S,VISP,RWS,ECORE,LCORE,NCORE,DRDI,
     +             R,Z,A,B,IRWS,INIPOL,IINFO)

        E = EZ(IE)

        IF(TEST('flow    '))  write(6,*) ' >>> after STARTB1'

        TMTIME1 = DCLOCK()    


c ------------------------------------------------------------------------
c calculate wavefunction and t-matrix of real system
c ------------------------------------------------------------------------

        IF (OPT('NEWSOSOL').OR.OPT('WRTGREEN').OR.
     &      OPT('GREENIMP').OR.OPT('BAND-STR')) THEN
         OPEN(UNIT=10,FILE='nonco_angle.dat',form='formatted')
        ENDIF

        DO 120 I1 = 1,NATYP

          IF(TEST('flow    '))  write(6,*) 'IATOM =',I1
          ISPIN=1
          IPOT = NSPIN* (I1-1) + ISPIN

          IRMINSO=1 

         IF (MOD(IRMINSO-1,2) .EQ. 1) STOP "IRMINSO must be not even"
         
c          IF (OPT('NEWSOSOL').OR.OPT('WRTGREEN').OR.
c     &        OPT('GREENIMP').OR.OPT('BAND-STR')) THEN

         IF (OPT('WRTGREEN').OR.OPT('BAND-STR').OR.
     &        OPT('FERMI-VL').OR.OPT('FERMIOUT')) THEN


c for new solver
          IF (OPT('NEWSOSOL')) THEN      
           READ(10,*) THETA,PHI
           THETA=THETA/360.0D0*2.0D0*PI
           PHI=PHI/360.0D0*2.0D0*PI
           IF (I1.EQ.1) WRITE(*,*) "Before TMAT_NEWSOLVER"
           WRITE(6,*) 'Atom',I1,IPOT

           CALL TMAT_NEWSOLVER(INS,NSPIN,LMAX,R(1,I1),Z(I1),E,KVREL,
     +           IRWS(I1),IPAN(I1),IRCUT(0,I1),IRMIN(I1),C,CLEB(1,1),
     +           ICLEB,IEND,NPAN_LOG,NPAN_EQ,NCHEB,R_LOG,
     +           VINS(IRMIND:IRMD,1:LMPOTD,IPOT:IPOT+NSPIN-1),
     +           VISP(1:IRMD,IPOT:IPOT+NSPIN-1),
     +           PNSSO(1,1,1,1,I1),PNSSO_LEFT(1,1,1,1,I1),
     +           TMATLL1(1,1,I1),THETA,PHI,I1,
     +           IMPLAYER,TSTAR,DRDI(1,I1))

          IF (I1.EQ.NATYP) WRITE(*,*) "After TMAT_NEWSOLVER"

          ELSE

c for old solver

          IF (I1.EQ.1) WRITE(*,*) "Before TMAT"

           CALL TMATRX01_SO(.FALSE.,ALPHA(:,:,I1),AR(:,:,I1),
     +           DRDI(1,I1),E,ICST,INS,NSPO,NSPOH,NSPIN,LMAX,
     +           PZI(1,1,I1,1),QZI,
     +           FZI(1,1,I1,1),SZI,
     +           PNSSO(1,1,1,1,I1),
     +           TMATLL1(1,1,I1),R(1,I1),
     +           VINS(IRMIND:IRMD,1:LMPOTD,IPOT:IPOT+NSPIN-1),
     +           VISP(1:IRMD,IPOT:IPOT+NSPIN-1),Z(I1),IRWS(I1),IPAN(I1),
     +           IRCUT(0,I1),IRMIN(I1),IRMINSO,KVREL,C,DROR(1,I1),
     +           RS(1,0,I1),S(0,I1),CLEB(1,1),LOFLM,ICLEB,IEND,
     +           IE,IELAST,EZ(IELAST),.TRUE.,PHASE_SHIFT(:,I1),LSM,I1,
     +           IMPLAYER,TSTAR)

          IF (I1.EQ.NATYP) WRITE(*,*) "After TMAT"

          ENDIF
         ENDIF
  100    CONTINUE          
c ------------------------------------------------------------------------
c transforms the t-matrix from the basis of real spherical harmonics
c to the relativistic kappa-mu-basis, in which the tmatrix is diagonal
c this is for full-relativistic ASA
c ------------------------------------------------------------------------
c
c          TRANS_TMAT=.FALSE.
cc          TRANS_TMAT=.TRUE.
c          IF (TRANS_TMAT == .TRUE. .AND. NSPO==2) THEN
c            CALL TRANSFORM_TMAT(I1,LMAX,LMMAXD,TMATLL1(1,1,I1),LSM,
c     +                      PNSSO(:,:,484,1,I1),IE,IELAST,EZ(IE),.TRUE.)
c          END IF
c
c ------------------------------------------------------------------------
c  check whether the AR transformation transforms the complex radial
c  wavefunctions to real wavefunctions
c  U_L (r)=sum _L' AR^-1 _LL' R_L' (r)
c ------------------------------------------------------------------------
c
c          IF (I1 == 3 ) THEN
c             WRITE(123,"(I5)") I1
c             CALL CHECK_REAL_RL(YR,PNSSO(:,:,IRMIND:IRMD,:,I1),
c     +                         I1,ALPHA(:,:,I1),AR(:,:,I1),NSPD,NSPIND,
c     +                         LOFLM,LM2D,
c     +                         LMAX,LMMAXD,LMPOTD,IRMIND,IRMD,IJEND,R)
c            STOP "after check AR " 
c          END IF

 120    CONTINUE              ! I1 = 1,NATYP

        TMTIME = TMTIME + DCLOCK() - TMTIME1

c ------------------------------------------------------------------------
c determine deltat = tmat(sys) - tmat(ref) = TMATLLSO1 - TMATLL
c ------------------------------------------------------------------------

        DO 430 I1 = 1,NAEZ

          INV = EQINV(I1)
          CFAC = CONE
          IF (I1.NE.INV) CFAC = CONEM

          INV = KAOEZ(INV)
          RF  = REFPOT(INV)

          DO LM1 = 1,LMMAXD
            DO LM2 = 1,LMMAXD
              CFACL = CFAC**(LOFLM(LM1)+LOFLM(LM2))
              TMATLL(LM2,LM1,I1) = CFACL*
     +           (TMATLL1(LM2,LM1,INV)-TREFLL(LM2,LM1,RF))
              IF (NSPD .NE. 1 ) THEN
                TMATLL(LM2+LMMAXD,LM1+LMMAXD,I1) = CFACL*
     +           (TMATLL1(LM2+LMMAXD,LM1+LMMAXD,INV)-TREFLL(LM2,LM1,RF))
              END IF
            END DO
            IF (NSPD .NE. 1 ) THEN
              DO LM2 = LMMAXD+1,2*LMMAXD
                CFACL = CFAC**(LOFLM(LM1)+LOFLM(LM2-LMMAXD))
                TMATLL(LM2,LM1,I1) = CFACL*TMATLL1(LM2,LM1,INV)
              END DO
            END IF
          END DO     

          IF (NSPD .NE. 1 ) THEN
            DO LM1 = LMMAXD+1,2*LMMAXD
              DO LM2 = 1,LMMAXD
                CFACL = CFAC**(LOFLM(LM1-LMMAXD)+LOFLM(LM2))
                  TMATLL(LM2,LM1,I1) = CFACL*TMATLL1(LM2,LM1,INV)
              END DO
            END DO
          END IF

 430    CONTINUE                ! I1 = 1,NAEZ

c ------------------------------------------------------------------------
c determine the spin operator, torque operator and spin flux operator
c ------------------------------------------------------------------------

        IF (IE.EQ.1) THEN
             WRITE(*,*) 'Computing spin operator'
             CALL NORMCOEFF_SO(IRMINSO,IRCUT,LMAX,
     +                       LMMAX,PNSSO,THETAS,NTCELL,
     +                       IFUNM,IPAN,LMSP,KVREL,CLEB,ICLEB,IEND,DRDI,
     +                       IRWS,LINTERFACE,ISP,IMPLAYER,NSPOH)
         IF (OPT('OPERATOR')) THEN
             WRITE(*,*) 'Computing torq operator'
             CALL NORMCOEFF_SO_TORQ(IRMINSO,IRCUT,LMAX,
     +                    LMMAX,PNSSO,THETAS,NTCELL,
     +                    IFUNM,IPAN,LMSP,KVREL,CLEB,ICLEB,IEND,DRDI,
     +                    IRWS,VISP,NSPIN,VINS,IRMIN)
             WRITE(*,*) 'Computing spinflux operator'
             CALL NORMCOEFF_SO_SPINFLUX(IRMINSO,IRCUT,LMAX,
     +                    LMMAX,PNSSO,THETAS,NTCELL,
     +                    IFUNM,IPAN,LMSP,KVREL,CLEB,ICLEB,IEND,DRDI,
     +                    IRWS)
         ENDIF!OPT('OPERATOR')
        ENDIF

        CLOSE(10)

c ------------------------------------------------------------------------
c calculate GF for the host and write down to file name green_host
c ------------------------------------------------------------------------

        IF (OPT('WRTGREEN')) THEN
         CALL WRITEGREEN(ALAT,NSPIN,MAXMESH,IE,IELAST,E,NSHELL,
     +        MD,INTERVY,INTERVZ,BLATC/ALATC,CLATC/ALATC,NAEZ,NATYP,CLS,
     +        EQINV,NACLS,RR,RBASIS,EZOA,ATOM,RCLS,KAOEZ,ICC,BRAVAIS,
     +        RECBV,LPOT,YR,WTYR,RIJ,IJEND,IATCONDL,IATCONDR,NCONDPAIR,
     +        LINTERFACE,INS,GINP,TMATLL,NATOMIMP,ATOMIMP,
c     +        LINTERFACE,INS,GINP(1,1,1,IE),TMATLL,NATOMIMP,ATOMIMP,
     +        RCLSIMP,IHOST)
        ENDIF

c ------------------------------------------------------------------------
c calculate GF,t-matrix and delta-matrix for impurity
c if IE = 1, calculate t- and delta-matrix, write to the file: DTMTRX
c if IE = 3, calculate GF, write to the file: GMATLL_GES
c ------------------------------------------------------------------------

        IF (OPT('GREENIMP')) THEN

         ALLOCATE(DTMTRX(NSPD*LMMAXD*NATOMIMP,NSPD*LMMAXD*NATOMIMP))
         DTMTRX=CZERO

         CALL TMATIMP_NEWSOLVER(INS,NSPIN,LMAX,R,Z,IELAST,E,
     +        KVREL,IRWS,IPAN,IRCUT,IRMIN,C,CLEB(1,1),ICLEB,IEND,
     +        NPAN_LOG,NPAN_EQ,NCHEB,R_LOG,VINS,VISP,NATOMIMP,RCLSIMP,
     +        ATOMIMP,IHOST,HOSTIMP,RIMP,ZIMP,IRWSIMP,
     +        IPANIMP,IRCUTIMP,IRMINIMP,VINSIMP,VISPIMP,DTMTRX)

         IF (IELAST.EQ.3) THEN
          CALL GREENIMP(NATOMIMP,DTMTRX,E)
         ENDIF

         DEALLOCATE(DTMTRX)    

        ENDIF
c
        STIMEK = DCLOCK()
           IF (KVREL.GE.1) THEN
              EK = SQRT(E+E*E/(C*C))
           ELSE
              EK = SQRT(E)
           END IF

c ------------------------------------------------------------------------
c write down data need to run FS-code (Bernd's code)
c ------------------------------------------------------------------------


      IF (OPT('FERMIOUT')) THEN

       WRITE(6,*) 'write down for FS calculation'

        CALL WRITETBKKR_DATA(ie,lmaxd, lmmaxd, lmax, nspd, nspo,
     +                         nspoh, nrd, nembd, nemb, nclsd, ncls,
     +                         natypd, natyp, naezd, naez, naclsd,
     +                         ielast, ins, ALAT, BRAVAIS, RECBV,
     +                         RBASIS, CLS, NACLS, EQINV, EZOA, ATOM,
     +                         KAOEZ, RCLS, RR, EZ(IE), TMATLL, 
     +                         GINP)
c     +                         GINP(1,1,1,IE))
      ENDIF


c ------------------------------------------------------------------------
c calculate spin expectation value, scattering probability, 
c conductivity by Boltzmann
c ------------------------------------------------------------------------

c      WRITE(6,*) 'before fermi-vl'
c      IF (OPT('FERMI-VL')) THEN
c        CALL FERMIVL(LSTART,MD,TMATLL,
c     +               EZ,DF,E2,ALAT,NSPIN,NSPOH,MAXMESH,NMESH,NSPD*LMMAX,
c     +               IE,IELAST,IGF,KSCOEF,NSHELL,INTERVY,INTERVZ,
c     +               BLATC/ALATC,CLATC/ALATC,NAEZ,NATYP,CLS,EQINV,NACLS,
c     +               RR,RBASIS,EZOA,ATOM,RCLS,KAOEZ,LATT,ICC,GINP,
c     +               BRAVAIS,RECBV,LPOT,YR,WTYR,RIJ,IJEND,
c     &               ATOMIMP,IATCONDL,IATCONDR,NCONDPAIR,IEGFOUT,
c     &               VOLBZ,LINTERFACE,INS,IMPLAYER,
c     &               DELTALAMDAK,TMATLLE,DELTALAMDADELE,TSTAR,
c     &               PMUP,PMDOWN,SPERTURB)
c      ENDIF

c ------------------------------------------------------------------------
c calculate bandstructure 
c ------------------------------------------------------------------------
c
c      IF (OPT('BAND-STR')) THEN
c       CALL BANDSTRSO(EZ,NSPD*LMMAX,NSPO,TMATLL,IE,ALAT,NAEZ,
c     +               CLS,EQINV,NACLS,RR,EZOA,ATOM,KAOEZ,
c     +               GINP,RCLS)
cc     +               GINP(1,1,1,IE),RCLS)
c      ENDIF

cc rotate coefficients along the given quantization axis
c        IF (OPT('ROTCOEFF')) THEN
c          CALL ROTATECOEFF_SO(LMMAX,LINTERFACE,IMPLAYER,ILAYERS)
c        ENDIF

cc calculate lifetime
c      IF (OPT('LIFETIME')) THEN
c
c        WRITE (6,*) "before LIFETIME"
c
c        CALL CALCULATE_LIFETIME(LMAXD,NSPO,NSPOH,NSPIN,
c     +       NCL_IMP,ALATC,LINTERFACE,NATYP,IMPLAYER,RR,
c     +       RBASIS,NR,LATT,MD,INTERVY,INTERVZ,BLATC/ALATC,CLATC/ALATC)
c
c      END IF    

 360  CONTINUE                  ! IE = 1,IELAST

      DEALLOCATE(GINP)
      DEALLOCATE(PNSSO_LEFT)
      DEALLOCATE(PNSSO)
      IF (OPT('GREENIMP')) THEN
       DEALLOCATE (IPANIMP,IRCUTIMP,IRMINIMP,IRWSIMP)
       DEALLOCATE (RIMP,ZIMP)
       DEALLOCATE (THETASIMP,VISPIMP,VINSIMP)
       CLOSE(59)
       CLOSE(60)
      ENDIF
      IF (OPT('WRTGREEN')) THEN
       CLOSE(58)
      ENDIF
      IF (OPT('BAND-STR')) THEN
       CLOSE(61)
       CLOSE(62)
      ENDIF

 2100 FORMAT(79(1H-))
 9080 format (' TIME IN GLL95             : ',f9.2)
 9170 FORMAT(' IE    = ',I3,', E = (',1p,d11.3,' ,',d11.3,')')
 9200 FORMAT(' The calculations are done with',I4,' points.')
 9210 FORMAT(' The temperature used is',F14.6,' K.')
 9220 FORMAT(I4,' poles and',I4,' and',I4,' and',I4,
     +     ' points on the contours are used.')

      WRITE(6,*) 'End of KKR'

      END










