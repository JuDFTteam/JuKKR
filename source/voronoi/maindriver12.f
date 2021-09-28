      PROGRAM KKRGEOMETRY
      use mod_version, only: version1, version2, version3, version4
      use mod_version_info, only: serialnr
      use mod_version_info, only: construct_serialnr
      use mod_version_info, only: version_print_header
      implicit none
      include 'inc.geometry' 
      INTEGER IBMAXD
      PARAMETER (IBMAXD=(LMAXD1+1)*(LMAXD1+1))
      REAL*8    ONE,ONEM,TWO
      PARAMETER (ONE = 1.D0,ONEM = -1.D0,TWO=2.D0)
      DOUBLE COMPLEX CONE,CONEM,CZERO,CI
      PARAMETER (CONE  = ( 1.0D0,0.0D0),
     +           CONEM = (-1.0D0,0.0D0),
     +           CZERO = ( 0.0D0,0.0D0),
     +           CI    = ( 0.0D0,1.0D0))
c#@# KKRcodes: VORONOI KKRhost KKRimp
c#@# KKRtags: geometry initialization input-output potential
c *****************************************************************
c * Program description and small help.
c * This is a utility of the tb-kkr and impurity programs. The 
c * purpose to construct the potential files possibly shape
c * functions and visualize the lattice using an external ray-tracer
c * this version is using povray.
c * As input the tb-kkr "inputcard" is used with some extra parameters
c * In case of impurity calculations an extra file with the atomic   
c * positions is neaded.
c * I describe some of the problems the program handles: 
c *
c * If you have no potentials just set the parameters you want
c * and leave blanc space as potential file 
c *
c *
c * 1. ASA potentials for some lattice, start from potentials 
c * in the GENERAL MESH or start from scratch
C *
C * i.  Define the lattice in the inputcard
C * ii. Put old potential as input and use ASA parameters in the
C *     inputcard
C * The program calculates asa spheres, etc and produces the potential
C * in the correct format to start the tb-kkr
C *
C * 2. I nead FP for some lattice 
C *
C * i Define lattice and give parameters for full potential
C * ii. Use some old potential or start from scratch
C * iii. Define the "muffin-tin-ization" in case of shifting
c *      the atoms later. This must be given in a.u. for each atom
C * The program calculates the first neighbours for each atom, makes
C * a Voronoi construction (posible weights) and then constructs the
C * shape functions. It sets the shape functions to the given mt-radius
c * and continues by constructing the potentials in the obtrained 
c * radial mesh.
c *
c * 3. Potentials for impurity calculation
c *  i. prepare "inputcard" for all parameters, and use the option
c *     "IMPURITY" make a file with the impurity atomic positions
c *     filename:  "impurity.atoms"
c *     in the format given in subroutine 'readimpatoms_kkrflex'
c !
c *
c * The program creates two files : lattice.pov and voronoi.pov
c * you nead to put the file "povray.ini" in your path and run
c * povray with :     povray +I lattice.pov
c * 
c * Rasmol with connectivity info also avaliable 3.3.2002, protein 
c * database format 
c * visualize file: lattice.pdb 
c *
c * the camera positions etc are in the end of the files.
c * help on povray : http://www.povray.org  
c *                                                    ver.  10.2000
c *
c *  parameters due to change to do different things:
c *    bbox    :  data statement changes the bounding box
c *               for drawing atoms with povray all atoms
c *               in the box are written out in file lattive.pov
c *    npoi    : data number of points for the shape function
c *              usualy set to 125 and can be changed 
c *    dlt = 0.05 controls the accuracy of angular integration
c *            for producing the shape functions can be also changed
c *    nrad : number of points for the muffin-tinizing 
c *           used in sub mtmesh
c *                                                  planed : smeared shapes
c *
c * KNOWN PROBLEM : Sometimes the voronoi construction fails.
c * This usualy means that the coordinates are not accurate enought
c * Try to change the weights by 0.001 or so and run again 
c *           
c * 
c *******************************************************************     
      REAL*8    
     +     ALATC,BLATC,CLATC,       ! lattice constants (in a.u.)
     +     ABASIS,BBASIS,CBASIS,    ! scaling factors for rbasis
     +     TOLVDIST,                ! Max. tolerance for distance of two vertices
     +     TOLAREA,                 ! Max. tolerance for area of polygon face
     +     TOLEULER,                ! Used in calculation of Euler angles, subr. EULER
     &     TOLHS,                   ! Tolerance for halfspace routine
     +     VOLUC,                   ! volume of unit cell in units of alat**3
     +     EFSET                    ! wished Fermi level of generated potential
C
C     .. REAL*8    ARRAYS ....
C
      REAL*8   
     +     ATWGHT(NATYPD),
     +     BRAVAIS(3,3),            ! bravais lattice vectors
     +     RECBV(3,3),              ! reciprocal basis vectors 
c
     +     A(NATYPD),B(NATYPD),     ! contants for exponential r mesh
     +     R(IRMD,NATYPD),          ! radial r mesh (in units a Bohr)
     +     RBASIS(3,NAEZD+NEMBD),   ! position of atoms in the unit cell
                                    ! in units of bravais vectors
!     +     RBASIS1(3,NAEZD+NEMBD),  ! pos. of atoms in atomic units
     +     RCLS(3,NACLSD,NCLSD),    ! real space position of atom in cluster
     +     RMT(NATYPD),             ! Muffin-Tin-radius
     +     RMTNEW(NATYPD),          ! adapted MT radius
     +     RMTREF(NAEZD+NEMBD),     ! MT-radius of ref. system
     +     RR(3,0:NRD),             ! set of real space vectors (in a.u.)
     +     RWS(NATYPD),             ! Wigner Seitz radius
     +     ZATOM(NTOTD)             ! Nuclear charge
      INTEGER ICC,                  ! center of cluster for output of GF
     +     ICLS, NAEZ,              ! number of atoms in unit cell
     +     NATYP,                   ! number of kinds of atoms in unit cell
     +     NCLS,                    ! number of reference clusters
     +     NEMB,                    ! number of 'embedding' positions
     &     NUMCELL,                 ! Number of inequiv. cells (disregarding shift)
     &     NUMSHAPE,                ! Number of different shapefunctions (cell + shift)
     &     NUMCONSTRUCTED           ! Number of different constructed shapefuncts.
      INTEGER
     +     EQINV(NAEZD),            ! site equiv. by invers. symmetry
     +     KAOEZ(NAEZD+NEMBD),      ! kind of atom at site in elem. cell
     &     SHAPEREF(NTOTD),         ! representative atom having equiv. shape (cell + shift)
     &     CELLREF(NTOTD),          ! repr. atom having equiv. cell (without shift)
     &     IDSHAPE(NTOTD),          ! which shape is assigned to a given atom
     &     SHAPEATOM(NSHAPED),      ! representative atom for each shape
     &     SITEAT(NATYPD)           ! Site of each atom in unit cell for CPA calculations (1.LE.SITEAT(I).LE.NAEZ, I=1,..,NATYP)
      LOGICAL LSHIFT(NTOTD),        ! true if the expansion is arround shifted position
     &        LCONSTRUCTED(NSHAPED),! true if particular shape has been constructed
     &        LMTREF,               ! true if MT-radius of ref. system was read in from input.
     &        LCARTESIAN,           ! cartesian or internal coordinates of basis
     &        LCARTESIMP            ! cartesian or internal coordinates of impurities

      REAL*8 
     &     RMTHLF(NAEZD+NEMBD),     ! Touching muffin-tin from bisection to nearest neighbor
     &     RMT0_ALL(NTOTD),         ! Muffin-tin radius per atom
     &     ROUT_ALL(NTOTD),         ! Outer cell-radius per atom
     &     DISTNN(NAEZD+NIMPD),     ! Distance from cell center to nearest-neighbor cell center (2*RMTHLF)
     &     VOLUME_ALL(NTOTD),       ! Volume per atom
     &     A3_ALL(NFACED,NTOTD),    ! A3,B3,C3,D3: Defining the faces per atom
     &     B3_ALL(NFACED,NTOTD),
     &     C3_ALL(NFACED,NTOTD),
     &     D3_ALL(NFACED,NTOTD),
     &     XVERT_ALL(NVERTD,NFACED,NTOTD), ! X,Y,Z coords per face per atom
     &     YVERT_ALL(NVERTD,NFACED,NTOTD),
     &     ZVERT_ALL(NVERTD,NFACED,NTOTD), 
     &     AOUT_ALL(NATYPD)         ! Parameter A for exponential mesh of output-pot.

      INTEGER NFACE_ALL(NTOTD),        ! Number of cell-faces per atom
     &        NVERT_ALL(NFACED,NTOTD)  ! Number of vertices per face per atom
c
c     ..  cluster arrays
c
      INTEGER
     +     ATOM(NACLSD,NTOTD),      ! atom at site in cluster
     +     CLS(NTOTD),              ! cluster around atom
     +     COCLS(NCLSD),            ! center of cluster (kaez )
     +     NACLS(NCLSD),            ! number of atoms in cluster
     +     EZOA(NACLSD,NAEZD),      ! EZ of atom at site in cluster
     +     IMT(NATYPD),             ! r point at MT radius
     +     IPAN(NATYPD),            ! number of panels in non-MT-region
     +     IRC(NATYPD),             ! r point for potential cutting
     +     IRCUT(0:IPAND,NATYPD),   ! r points of panel borders
     +     IRMIN(NATYPD),           ! max r for spherical treatment
     +     IRNS(NATYPD),            ! number r points for non spher. treatm.
     +     IRWS(NATYPD),            ! r point at WS radius
     &     LMXC(NATYPD),KFG(4,NATYPD),NTCELL(NATYPD),
     &     IREFPOT(NAEZD+NEMBD)      
c     
      INTEGER NLAY,NLBASIS,NRBASIS,                     
     +        NLEFT,NRIGHT,IER,NTOTAL,N,II1,II2,IVEC   
      REAL*8                                
     +     ZPERLEFT(3),ZPERIGHT(3),               
     +     TLEFT(3,NEMBD),TRIGHT(3,NEMBD),RMTCORE(NTOTD)
c
c Impurity:
      INTEGER NUMIMP,NKILLATOM  ! Number of impurity sites and sites to be killed
      INTEGER KILLATOM(NIMPD)   ! (0/1) Host atoms to be removed in impurity calculation
      INTEGER CLSIMP(NIMPD)     ! Cluster-type of each impurity
      INTEGER NACLSIMP(NACLSD)
      INTEGER ATOMIMP(NACLSD,NIMPD)  ! Atom-type of each position in an imp-cluster
      INTEGER NCLSIMP
      REAL*8 RCLSIMP(3,NACLSD,NACLSD)                  ! Impurity-cluster positions
      REAL*8 RIMPURITY(3,NIMPD),RKILL(3,NIMPD) ! Coords of imp. sites and "killed" sites
      REAL*8 DXIMP(NIMPD),DYIMP(NIMPD),DZIMP(NIMPD) ! Impurity shift
      REAL*8 ZIMP(NIMPD)        ! Impurity atomic number
      REAL*8 RMTIMP(NIMPD)      ! Impurity core-radius
      REAL*8 RMTHLFIMP(NIMPD)   ! Touching muffin-tin from bisection to nearest neighbor
c End Impurity
c
      INTEGER NVEC,NFACELIM ! Number of vectors defining faces, maximum allowed number of faces
      REAL*8  RVEC(3,NFACED)
c
c  Shape functions
c
      INTEGER NCELL
      INTEGER NPAN_ALL(NSHAPED),MESHN_ALL(NSHAPED),NFUN_ALL(NSHAPED)
      INTEGER NM_ALL(NPAND,NSHAPED),LMIFUN_ALL(IBMAXD,NSHAPED)
      REAL*8    XRN_ALL(IRID,NSHAPED),DRN_ALL(IRID,NSHAPED),! Radial mesh and weight in shape-region
     &                 THETAS(IRID,IBMAXD),
     &                 THETAS_ALL(IRID,IBMAXD,NSHAPED)  ! Radial shape functions
      REAL*8    SCALE_ALL


      INTEGER NFACE,NPANEL
      INTEGER NVERTMAX,NFACEMAX ! Max. no. of vertices and faces among all cells.
      INTEGER NVERT(NFACED)
      INTEGER NFACECL(NCLSD),NVERTCL(NFACED,NCLSD)
      REAL*8    VOLUME
      REAL*8    A3(NFACED),B3(NFACED),C3(NFACED),D3(NFACED)
      REAL*8    A3CL(NFACED,NCLSD),B3CL(NFACED,NCLSD),
     &                 C3CL(NFACED,NCLSD),D3CL(NFACED,NCLSD)
      REAL*8    XVERT(NVERTD,NFACED),YVERT(NVERTD,NFACED),
     &       ZVERT(NVERTD,NFACED)
      REAL*8    XVERTCL(NVERTD,NFACED,NCLSD),
     &          YVERTCL(NVERTD,NFACED,NCLSD),
     &          ZVERTCL(NVERTD,NFACED,NCLSD),
     &          VOLUMECL(NSHAPED),RWSCL(NSHAPED),RMTCL(NSHAPED)
      REAL*8    DX(NTOTD),DY(NTOTD),DZ(NTOTD)
      REAL*8    ROUT,RTEST,DLT,CRAD,RX,RY,RZ,MTRADIUS,VTOT,
     &                 SHAPESHIFT(3,NTOTD)
      CHARACTER*256 UIO
      INTEGER NATOMS,NSITES,NSHAPE ! Number of atoms, sites, shapes
      INTEGER LMAX,KEYPAN,NPOI,NA,IAT,JAT,ICL,N1A,I2,II,ISITE
      INTEGER NBEGIN,ISITEBEGIN
      INTEGER IFACE,ISHAPE,JV,CELLREFI,IVERT
      INTEGER NM(NPAND),SHAPECL(NATYPD)
      INTEGER IL,I,J,IFILE,NREF,INS,KWS,LPOT,LMPOT,NSPIN,KSHAPE
      INTEGER NR,KMT,NINEQ,LMMAX,I1,IC0,IX,IIMP,ICLSIMP
      LOGICAL LINTERFACE,LJELL,MAKESHAPE,VOROPLOT,LCOMPARE,LSKIP
      REAL*8    RCUTZ,RCUTXY,QBOUND,SHIFT(3),DIFF,RMT0,DXI,DYI,DZI
      REAL*8    WEIGHT0,RMTCOREI 
      REAL*8    BBOX(3),IMPSIZE(NIMPD),WEIGHT(NFACED)
      REAL*8    SIZEFAC(-NLEMBD*NEMBD:NTOTD)
      INTEGER IRM,KXC,ICLUSTER,LMAXSHAPE,NRAD,NMIN,NSMALL
      CHARACTER*124 TXC(3)
      LOGICAL OPT,TEST
      CHARACTER*3 ELEM_NAME
      CHARACTER*40 I13
      REAL*8 DISTPLANE,PI
      EXTERNAL DISTPLANE

      REAL*8 CRT(NPAND)  ! Critical points (end-points of panels)
      REAL*8 DENPT !"Optimal" density of points between RMT and outer radius (eg 100/a_B)
      REAL*8 STARTFP ! Defining the radius to start full-potential treatment
      REAL*8 FPRADIUS(NTOTD) ! RMT * STARTFP
      REAL*8 BOUT
      INTEGER NPAN,NMESH ! Number of panels, suggested number of mesh-points in outer region
      INTEGER NMT ! Number of points in MT-radius
      INTEGER IP

      LOGICAL POT_EXISTS           ! set to true if a potential file of the  2014/06 Benedikt
                                   ! name as given in inputcard exists       2014/06 Benedikt

      SAVE
      DATA NM/NPAND*25/
c
c  Data paremeters due to change for npoi check dimensions in inc.geometry
c  for description look at main header description
c
c     -----------------------------------------------------------------------
      DATA BBOX/2.0d0,2.0d0,3.0d0/
      DATA DLT/0.05d0/  ! Parameter for theta-integration (Gauss-Legendre rule). Usually 0.05
      DATA NPOI/125/    ! Total number of shapefunction points
      DATA NRAD/10/     ! Muffintinization points
      DATA NMIN/7/      ! Minimum number of points in panel
      DATA NSMALL/10000/ ! A large number to start (See subr. divpanels)
      DATA NMT/349/     ! Number of points in MT-radius
      DATA STARTFP/0.2D0/ ! Radius to start full-potential treatment = STARTFP * RMT
      DATA DENPT/100.D0/  ! "Optimal" density of points between RMT and outer radius (eg 100/a_B)
      DATA TOLVDIST/1.D-10/ ! Max. tolerance for distance of two vertices
      DATA TOLHS/1.D-16/   ! Used in soutine HALFSPACE
      DATA TOLAREA/1.D-10/  ! Max. tolerance for area of polygon face
      DATA TOLEULER/1.D-10/ ! Used in calculation of Euler angles, subr. EULER
c     -----------------------------------------------------------------------     


      ! Find serial number and print version
      call construct_serialnr()
      WRITE(*,'(A)')       '##########################################'
      WRITE(*,'(A)')       'This is the Voronoi program'
      WRITE(*,'(A,A)')     'Code version: ',trim(version1)
      WRITE(*,'(A,A,A,A)') 'Compile options:', trim(version2), 
     +                                         trim(version3), 
     +                                         trim(version4)
      WRITE(*,'(A,A)')     'serial number for files:', serialnr
      WRITE(*,'(A)')       '##########################################'
      WRITE(*,'(A)')

      IF (IRMD.LT.IRNSD) STOP 'Error: IRMD.LT.IRNSD'
      IF (IRNSD.LT.IRID) STOP 'Error: IRNSD.LT.IRID'
      PI = 4.D0*DATAN(1.D0)
      QBOUND = 1.D-07
      TXC(1) = ' Morruzi,Janak,Williams  #serial: ' // serialnr
      TXC(2) = ' von Barth,Hedin         #serial: ' // serialnr
      TXC(3) = ' Vosko,Wilk,Nusair       #serial: ' // serialnr
c     
c     Read Lattice form inputcard
c
      WRITE(6,2000)
      WRITE(6,2001)
      WRITE(6,2002)
      WRITE(6,2000)
     
      OPEN(15,file='shapefun',status='unknown')

      DX(:) = 0.D0      
      DY(:) = 0.D0     
      DZ(:) = 0.D0     

      CALL READINPUT(BRAVAIS,LCARTESIAN,RBASIS,ABASIS,BBASIS,CBASIS,
     &     DX,DY,DZ,
     &     ALATC,BLATC,CLATC,
     &     IRNS,NAEZ,NEMB,KAOEZ,IRM,ZATOM,SITEAT,
     &     INS,KSHAPE,
     &     LMAX,LMMAX,LPOT, 
     &     NATYP,NSPIN,
     &     NMIN,NRAD,NSMALL,
     &     TOLHS,TOLVDIST,TOLAREA,
     &     KXC,TXC,
     &     I13,
     &     NLBASIS,NRBASIS,NLEFT,NRIGHT,ZPERLEFT,ZPERIGHT,    
     &     TLEFT,TRIGHT,LINTERFACE,RCUTZ,RCUTXY,RMTCORE,
     &     LMTREF,RMTREF,SIZEFAC,NFACELIM, EFSET, AOUT_ALL)



      CALL LATTIX12(LINTERFACE,ALATC,BRAVAIS,RECBV,RR,NR,VOLUC)   

c     
c     Scale 
c     
      CALL SCALEVEC2000(LCARTESIAN,RBASIS,ABASIS,BBASIS,CBASIS,
     &     NLBASIS,NRBASIS,NLEFT,NRIGHT,ZPERLEFT,ZPERIGHT,            
     &     TLEFT,TRIGHT,LINTERFACE,NAEZ,NEMB,BRAVAIS,KAOEZ)
c
c     Rationalise basis vectors
      CALL RATIONALBASIS(
     >     BRAVAIS,RECBV,NAEZ+NEMB,
     X     RBASIS)

      IF (LINTERFACE) THEN
        CALL RATIONALBASIS(
     >       BRAVAIS,RECBV,NLBASIS,
     X       TLEFT)
        CALL RATIONALBASIS(
     >       BRAVAIS,RECBV,NRBASIS,
     X       TRIGHT)
      ENDIF
c
      CALL CLSGEN_VORONOI(NATYP,NAEZ,NEMB,RR,NR,RBASIS,
     &        KAOEZ,ZATOM,CLS,NCLS,
     &        NACLS,ATOM,EZOA,
     &        NLBASIS,NRBASIS,NLEFT,NRIGHT,ZPERLEFT,
     &        ZPERIGHT,TLEFT,TRIGHT,
     &        RCLS,RMTHLF,RCUTZ,RCUTXY,LINTERFACE,
     &        ALATC)

      DISTNN(1:NAEZ) = 2.D0*RMTHLF(1:NAEZ)

      NSITES = NAEZ
      NATOMS = NATYP

!------------------------------------
! Begin impurity part
      NUMIMP = 0 ! Initialize
      IF (OPT('IMPURITY')) THEN 

      stop 'not working properly'

      ! For impurity calculation, the host clusters are used as input, 
      ! while impurity clusters are generated as output.
!         CALL READIMPATOMS12(
!     >        ALATC,LCARTESIAN,
!     <        NUMIMP,RIMPURITY,NKILLATOM,RKILL,DXIMP,DYIMP,DZIMP,
!     <        RMTIMP,IMPSIZE,ZIMP,LCARTESIMP)
!         CALL SCALEVECIMP(
!     >        NUMIMP,NKILLATOM,BRAVAIS,LINTERFACE,LCARTESIMP,
!     X        RIMPURITY,RKILL,DXIMP,DYIMP,DZIMP)
!         CALL CLSGENIMP12(
!     >        NUMIMP,RIMPURITY,NKILLATOM,RKILL,
!     >        CLS,NCLS,NACLS,ATOM,RCLS,
!     >        BRAVAIS,RECBV,NAEZ,RBASIS,RCUTZ,RCUTXY,
!     <        CLSIMP,NCLSIMP,NACLSIMP,ATOMIMP,RCLSIMP,RMTHLFIMP)
!         
!         DO IIMP=1,NUMIMP
!            DO IX=1,3
!               RIMPURITY(IX,IIMP) = RIMPURITY(IX,IIMP)*ALATC
!            END DO
!         END DO  

         ! From here on continue as if the impurities were host-sites
         ! and just re-name the arrays.
         
         ! Indexing: Crystal-sites are stored from 1 to NAEZ, 
         ! impurity sites from NAEZ+1 to NAEZ+NUMIMP.
         ! This indexing helps to avoid multiple calculation of
         ! shapes later.
         ! This was already anticipated in subroutine CLSGENIMP_KKRFLEX
         ! when creating the array ATOMIMP.
         NSITES = NAEZ + NUMIMP
         NATOMS = NATYP + NUMIMP  ! In cpa, NATYP>NAEZ
         IF (NSITES.GT.NTOTD) THEN
            WRITE(*,*) 'Found ',NSITES,' atoms, but NTOTD=',NTOTD
            STOP 'NTOTD'
         ENDIF
         IF (NCLS + NCLSIMP.GT.NCLSD) THEN
            WRITE(*,*) 'Found ', NCLS + NCLSIMP,
     &                 ' clusters, but NCLSD=',NCLSD
            STOP 'NCLSD'
         ENDIF
      
         DO IIMP = 1,NUMIMP
            ISITE = NAEZ + IIMP
            SITEAT(NATYP+IIMP) = ISITE          ! ZATOM is atom-dependent, not site dependent (CPA)
            CLS(ISITE)        = CLSIMP(IIMP)
            ATOM(:,ISITE)     = ATOMIMP(:,IIMP)
            DX(ISITE)         = DXIMP(IIMP)
            DY(ISITE)         = DYIMP(IIMP)
            DZ(ISITE)         = DZIMP(IIMP)
            SIZEFAC(ISITE)    = IMPSIZE(IIMP)
            ZATOM(NATYP+IIMP) = ZIMP(IIMP)       ! ZATOM is atom-dependent, not site dependent (CPA)
            RMTCORE(ISITE)    = RMTIMP(IIMP)
            DISTNN(ISITE)     = 2.D0*RMTHLFIMP(IIMP)
            ! IRNS for impurities is initialized arbitrarily 
            ! but can be read-in if necessary. It is reset
            ! later: see "Check consistency of param. IRNS"
            IRNS(ISITE) = 208
         END DO
         DO ICLSIMP = 1,NCLSIMP
            ICLS = NCLS + ICLSIMP
            RCLS(:,:,ICLS) = RCLSIMP(:,:,ICLSIMP)
            NACLS(ICLS)    = NACLSIMP(ICLSIMP)
         ENDDO
         NCLS = NCLS + NCLSIMP

      END IF       ! IF (OPT('IMPURITY'))
! End impurity part
!------------------------------------


!------------------------------------
! Set weights in array sizefac according to touching muffin tin
      IF (OPT('FINDSIZE').OR.OPT('findsize')) THEN

         SIZEFAC(1:NAEZ) = RMTHLF(1:NAEZ)**2

         IF (LINTERFACE) THEN
            II = 0
            DO IL = 1,NLEFT
               DO IAT = 1,NEMB
                  II = II - 1
                  SIZEFAC(II) = RMTHLF(NAEZ+IAT)**2
               ENDDO
            ENDDO
         ENDIF

         DO IIMP = 1,NUMIMP
            SIZEFAC(NAEZ+IIMP) =  RMTHLFIMP(IIMP)**2
         ENDDO

      ENDIF
!------------------------------------


c The following parameters can be made atom-dependent and e.g. read in.
      DO IAT = 1,NSITES
         if (AOUT_ALL(IAT)<0.d0) then
           ! set to default value if negative value is found
           AOUT_ALL(IAT) = 0.025D0  ! Parameter A for exponential mesh of output-pot.
         end if
         IRWS(IAT) = IRM          ! Index of outmost point.
         IMT(IAT) = NMT
      ENDDO

c     
c     Create geometry info for Voronoi construction
c
c -----------------------------------------------------------------

c Find all inequivalent Voronoi cells

      DO 20 IAT = 1,NSITES  ! Loop over all atoms to find their cells

c        Prepare cluster for this voronoi cell. 
c        Exclude first vector, which is always zero by cluster construction.         
c        Store in array RVEC.

         ICLUSTER = CLS(IAT)
         ! Upper limit of allowed number of faces can be smaller than dimension for speedup
         NVEC = MIN(NACLS(ICLUSTER) - 1,NFACELIM)  ! No. of cluster atoms excl. central
         WRITE(6,*) 'Preparing neighbours of Site:',IAT
         WRITE(6,*) 'Number of considered neighbours is:',NFACELIM
         IF (NVEC.GT.NFACED) THEN
            WRITE(6,*) 'Increase NFACED ',NVEC,NFACED
            STOP
         END IF
         DO IVEC = 2,NVEC + 1           
            RVEC(1:3,IVEC-1) = RCLS(1:3,IVEC,ICLUSTER) 
         END DO
         WRITE(6,*) 'Max. neighbour distance is:',
     &            DSQRT(RVEC(1,NVEC)**2+RVEC(2,NVEC)**2+RVEC(3,NVEC)**2)


c        Cluster mapping done, now map weights.
         DO IVEC = 1,NVEC 
            II = ATOM(IVEC+1,IAT)   ! +1 because of exclusion of central atom
                                    ! II<0 for the embedded positions around the slab
            IF (II.GT.NSITES) II = 0  ! should not happen

c           In case of slab or decimation you have to take care of the boundary
c           atoms. Give all atoms outside the system weight 1.0.
c           Therefore, sizefac(0) = 1.0 is defined earlier.
            WEIGHT(IVEC) = SIZEFAC(II)     

         END DO
         WEIGHT0 = SIZEFAC(IAT)
     

         WRITE(6,*) 'Entering VORONOI12 for atom=',IAT            
         CALL VORONOI12(
     >    NVEC,RVEC,NVERTD,NFACED,WEIGHT0,WEIGHT,TOLVDIST,TOLAREA,TOLHS,
     <    RMT0,ROUT,VOLUME,NFACE,A3,B3,C3,D3,NVERT,XVERT,YVERT,ZVERT)

c        Now store results in atom-dependent array.
         RMT0_ALL(IAT) = RMT0
         ROUT_ALL(IAT) = ROUT
         VOLUME_ALL(IAT) = VOLUME
         NFACE_ALL(IAT) = NFACE
         A3_ALL(:,IAT) = A3(:)
         B3_ALL(:,IAT) = B3(:)
         C3_ALL(:,IAT) = C3(:)
         D3_ALL(:,IAT) = D3(:)
         NVERT_ALL(:,IAT) = NVERT(:)
         XVERT_ALL(:,:,IAT) = XVERT(:,:)
         YVERT_ALL(:,:,IAT) = YVERT(:,:)
         ZVERT_ALL(:,:,IAT) = ZVERT(:,:)


 20   ENDDO                 ! DO 20 IAT = 1,NSITES 


c-------------------------------------------------------------------------------
c Compare all Voronoi cells to find representative ones.
c If all faces are the same, then the cells are equivalent (except center-shift).

      NUMCELL = 1     ! Cell of first atom is always accepted.
      CELLREF(1) = 1
      
      DO 30 IAT = 2,NSITES
         NFACE = NFACE_ALL(IAT)
         A3(:) = A3_ALL(:,IAT)
         B3(:) = B3_ALL(:,IAT)
         C3(:) = C3_ALL(:,IAT)
         D3(:) = D3_ALL(:,IAT)


         LCOMPARE = .FALSE.
         DO JAT = 1,IAT - 1 ! Loop over all previous atoms
            IF (NFACE.EQ.NFACE_ALL(JAT)) THEN
               DIFF = 0
               DO IFACE = 1,NFACE
                  DIFF = DIFF + DABS(A3(IFACE)-A3_ALL(IFACE,JAT)) + 
     &                          DABS(B3(IFACE)-B3_ALL(IFACE,JAT)) +
     &                          DABS(C3(IFACE)-C3_ALL(IFACE,JAT)) + 
     &                          DABS(D3(IFACE)-D3_ALL(IFACE,JAT)) 
               ENDDO
               IF (DIFF.LT.1D-5) LCOMPARE=.TRUE. ! Cells are equivalent
            ENDIF

            IF (LCOMPARE) THEN
               CELLREF(IAT) = JAT
               GOTO 30 ! Jump loop
            ENDIF

         ENDDO ! JAT = 1,IAT - 1

         IF (.NOT.LCOMPARE) THEN
            NUMCELL = NUMCELL + 1
            CELLREF(IAT) = IAT
         ENDIF

 30   ENDDO  ! IAT = 2,NSITES
      WRITE(*,*) 'Found ',NUMCELL,' inequivalent cells.'

c Now the number of inequivalent cells is NUMCELL.
c Each atom I has the same cell as atom CELLREF(I).

c------------------------------------------------------------------------------
c Now find all different shape functions. 
c Consider two shape functions equivalent if the Voronoi cells are
c equivalent (same CELLREF), if the shift DX,DY,DZ is the same,
c and if the chosen core radius is the same.

      DO IAT = 1,NSITES 
         IF (DABS(DX(IAT))+DABS(DY(IAT))+DABS(DZ(IAT)).GE.1.D-5) THEN
            LSHIFT(IAT) = .TRUE.
         ELSE
            LSHIFT(IAT) = .FALSE.
         ENDIF
      END DO 

      NUMSHAPE = 1
      SHAPEREF(1) = 1
      SHAPEATOM(1) = 1
      IDSHAPE(1) = 1
      
      DO 40 IAT = 2,NSITES

         DXI = DX(IAT)
         DYI = DY(IAT)
         DZI = DZ(IAT)
         RMTCOREI = RMTCORE(IAT)
         CELLREFI = CELLREF(IAT)

         DO JAT = 1,IAT - 1
            ! If atoms have equivalent cells, compare center shift and core-radius
            IF (CELLREFI.EQ.CELLREF(JAT)) THEN    
               DIFF = DABS(DXI-DX(JAT)) + DABS(DYI-DY(JAT)) + 
     &                DABS(DZI-DZ(JAT)) + DABS(RMTCOREI - RMTCORE(JAT))
               IF (DIFF.LE.1.D-5) THEN
                  SHAPEREF(IAT) = SHAPEREF(JAT)
                  IDSHAPE(IAT) = IDSHAPE(JAT)
                  GOTO 40 ! Jump loop
               ENDIF
            ENDIF
         ENDDO
c        If this point is reached, then no equivalent shape was found.         
         NUMSHAPE = NUMSHAPE + 1
         SHAPEREF(IAT) = IAT         ! Representative atom of particular atom concerning shapes
         IDSHAPE(IAT) = NUMSHAPE     ! Shape index for particular atom
         SHAPEATOM(NUMSHAPE) = IAT   ! Representative atom of shape
 40   ENDDO

      WRITE(*,*) 'Found ',NUMSHAPE,' inequivalent shapes.'
      DO IAT = 1,NSITES
         WRITE(*,FMT='(10(A10I5))') '%    Atom ',IAT,
     &        '  cellref ',CELLREF(IAT),' shaperef ',SHAPEREF(IAT),
     &        '  idshape ',IDSHAPE(IAT)
      ENDDO

c------------------------------------------------------------------------------
c Define RWS and RMT in case of ASA. Otherwise these are recalculated.
      DO ISHAPE = 1,NUMSHAPE
         IAT = SHAPEATOM(ISHAPE)
         VOLUMECL(ISHAPE) = VOLUME_ALL(IAT)
         RWSCL(ISHAPE) = 
     &                ALATC*(3.D0*VOLUMECL(ISHAPE)/4.D0/PI)**(1.D0/3.D0)

         MTRADIUS = MIN( RMT0_ALL(IAT)*ALATC , RMTCORE(IAT) ) ! in atomic units
         RMTCL(ISHAPE) = MTRADIUS
      ENDDO
c------------------------------------------------------------------------------


      ISITEBEGIN = 1
      NBEGIN = 1
      IF (OPT('IMPURITY')) THEN 
         ISITEBEGIN = NAEZ + 1
         NBEGIN = NATYP + 1  ! In cpa, NATYP>NAEZ
      ENDIF

c**********************************************************************************
cINS-INS-INS-INS-INS-INS-INS-INS-INS-INS-INS-INS-INS-INS-INS-INS-INS-INS-INS-INS-INS-
      IF (INS.GT.0) THEN

c Loop to construct different shapes is from 1 to NSITES=NAEZ for host calculation
c or from NAEZ+1 to NSITES=NAEZ+NUMIMP for impurity calculation.


      LCONSTRUCTED(1:NUMSHAPE) = .FALSE. ! Initialize. Shapes not yet constructed.
      NUMCONSTRUCTED = 0

      DO 50 ISHAPE = 1,NUMSHAPE

         ! In case of impurity calculation: check if the shape is only relevant
         ! for host atoms but not for impurity atoms. If so, then skip construction.
         LSKIP = .TRUE.
         DO IAT = ISITEBEGIN,NSITES
            IF (IDSHAPE(IAT).EQ.ISHAPE) LSKIP = .FALSE.
         ENDDO
         IF (LSKIP) GOTO 50     ! Avoid calculation

         IAT = SHAPEATOM(ISHAPE)

         WRITE(*,*) 'Constructing shape',ISHAPE,' with rep. atom',IAT

         NFACE = NFACE_ALL(IAT)
         A3(:) = A3_ALL(:,IAT)
         B3(:) = B3_ALL(:,IAT)
         C3(:) = C3_ALL(:,IAT)
         D3(:) = D3_ALL(:,IAT)
         NVERT(:) = NVERT_ALL(:,IAT)
         XVERT(:,:) = XVERT_ALL(:,:,IAT)
         YVERT(:,:) = YVERT_ALL(:,:,IAT)
         ZVERT(:,:) = ZVERT_ALL(:,:,IAT)

c------------------------------------------------------------------------------
C     Perform shift of polyhedron center by DX,DY,DZ.
         IF (LSHIFT(IAT)) THEN
            DO IFACE = 1,NFACE
               D3(IFACE) = D3(IFACE) - 
     &         A3(IFACE)*DX(IAT) - B3(IFACE)*DY(IAT) - C3(IFACE)*DZ(IAT)
               DO JV = 1,NVERT(IFACE)                  
                  XVERT(JV,IFACE) = XVERT(JV,IFACE) - DX(IAT)
                  YVERT(JV,IFACE) = YVERT(JV,IFACE) - DY(IAT)
                  ZVERT(JV,IFACE) = ZVERT(JV,IFACE) - DZ(IAT) 
               END DO
            END DO
            ! Update touching radius rmt0 and external radius rmax round new center.
            RMT0 = 1.D10
            ROUT = 0.D0
            DO IFACE = 1,NFACE
               DIFF = DISTPLANE(A3(IFACE),B3(IFACE),C3(IFACE),D3(IFACE))
               IF (DIFF.LT.RMT0) RMT0 = DIFF
               DO JV = 1,NVERT(IFACE) 
                  DIFF = DSQRT( XVERT(JV,IFACE)**2 + YVERT(JV,IFACE)**2 
     &                                             + ZVERT(JV,IFACE)**2)
                  IF (DIFF.GT.ROUT) ROUT = DIFF
               ENDDO
            ENDDO
            ! redefine inscribed and outer radius after shift
            RMT0_ALL(IAT) = RMT0  
            ROUT_ALL(IAT) = ROUT
            ! This should be the same as calculated by sub. shape as XRN_ALL(1,ISHAPE)
         END IF

c------------------------------------------------------------------------------
c Find number of panels without entering shape routine
         CALL FINDPANELS(
     >     NFACE,A3,B3,C3,D3,NVERT,XVERT,YVERT,ZVERT,TOLEULER,TOLVDIST,
     <     NPAN,CRT)
c Suggest number of points between MT radius and outer radius based on panels
c to prepare for mesh-calculation in subr. shape.
         CALL SUGGESTPTS(
     >                   NPAN,CRT,DENPT,ALATC,NMIN,
     <                   NMESH)
         WRITE(*,FMT=
     &   '("No. of panels by sub FINDPANELS for shape:",I5,":",I5)')
     &      ISHAPE,  NPAN
         WRITE(*,FMT='(10F24.16)') (CRT(IP),IP=1,NPAN)
         WRITE(*,*) 
     &   'Suggested number of points before muffintinization:',NMESH
         NMESH = MAX(NMESH,NPOI)
         IRWS(IAT) = IMT(IAT) + NMESH + NRAD
c Some checks on the dimensions:
         IF (NMESH+NRAD.GT.IRID) THEN
            WRITE(*,*) 'INCREASE IRID IN inc.geometry',IRID
            STOP 'IRID'
         ENDIF
         IF (IRWS(IAT).GT.IRMD) THEN
            WRITE(*,*) 'INCREASE IRMD IN inc.geometry',IRMD,IRWS(IAT)
            STOP 'IRMD'
         ENDIF

c------------------------------------------------------------------------------
c Calculate lm-expansion of shape functions
         LMAXSHAPE = 4*LMAX
         KEYPAN = 0
         DO IP=1,NPAND
            NM(IP) = 0   
         END DO

         MTRADIUS = MIN( RMT0_ALL(IAT)*ALATC , RMTCORE(IAT) ) ! in atomic units


         IF (.NOT.OPT('SIMULASA').AND.KSHAPE.GT.0) THEN


            WRITE(*,*) 'Calculating Shape:',ISHAPE
            CALL SHAPE(NMESH,A3,B3,C3,D3,
     &           TOLVDIST,TOLEULER,NMIN,
     &           NVERT,XVERT,YVERT,ZVERT,NFACE,LMAXSHAPE,
     &           DLT,KEYPAN,NM,ICLUSTER,NCELL,SCALE_ALL,NPAN_ALL(ISHAPE)
     &          ,MESHN_ALL(ISHAPE),NM_ALL(1,ISHAPE),XRN_ALL(1,ISHAPE),
     &           DRN_ALL(1,ISHAPE),NFUN_ALL(ISHAPE),LMIFUN_ALL(1,ISHAPE)
     &          ,THETAS )
         ELSE
            CALL FAKESHAPE(
     >           VOLUME_ALL(IAT),MTRADIUS,ALATC,NMESH,
     <           NPAN_ALL(ISHAPE),MESHN_ALL(ISHAPE),NM_ALL(1,ISHAPE),
     <           NFUN_ALL(ISHAPE),XRN_ALL(1,ISHAPE),DRN_ALL(1,ISHAPE),
     <           LMIFUN_ALL(1,ISHAPE),THETAS )
                 LCONSTRUCTED(ISHAPE) = .TRUE.
!                 NUMCONSTRUCTED = NUMCONSTRUCTED + 1

         ENDIF
c------------------------------------------------------------------------------

c The first mesh point of the shape is the mt-radius, save it
         RMTCL(ISHAPE) = XRN_ALL(1,ISHAPE)*ALATC  

         WRITE(6,1020) NFACE
         WRITE(6,1000) VOLUMECL(ISHAPE)
         WRITE(6,1010) RWSCL(ISHAPE)       
         WRITE(6,1030) RMTCL(ISHAPE)

         IF (MTRADIUS.GT.RMTCL(ISHAPE)*0.99d0) THEN
            WRITE(6,*) 'MTRADIUS,RMTCL',MTRADIUS,RMTCL(ISHAPE)
            MTRADIUS = 0.99D0*RMTCL(ISHAPE)
            WRITE(6,*) 'RMT OF SHAPE',ISHAPE,' REDUCED TO ',MTRADIUS
         ENDIF
            
c------------------------------------------------------------------------------
c Insert extra panel between core radius and touching-MT radius
         IF ( (MTRADIUS.LT.RMTCL(ISHAPE)) .AND. (MTRADIUS.GT.0.D0)) THEN
            WRITE(6,1050)
            WRITE(6,1070) MTRADIUS
            WRITE(6,1080)
     &           (RMTCL(ISHAPE)-MTRADIUS)/RMTCL(ISHAPE)*100.0D0
c     
c           ! Now change MT and write out the shape function
c     
            CALL MTMESH(NRAD,NPAN_ALL(ISHAPE),MESHN_ALL(ISHAPE),
     &           NM_ALL(1,ISHAPE),XRN_ALL(1,ISHAPE),DRN_ALL(1,ISHAPE),
     &           NFUN_ALL(ISHAPE),THETAS,LMIFUN_ALL(1,ISHAPE),
     &           MTRADIUS/ALATC)     
            THETAS_ALL(:,:,ISHAPE) = THETAS(:,:)
  
c           ! Redefine
            RMTCL(ISHAPE) = XRN_ALL(1   ,ISHAPE)*ALATC
            IF (.NOT.OPT('SIMULASA')) THEN
              RWSCL(ISHAPE) = XRN_ALL(NMESH+NRAD,ISHAPE)*ALATC
            ELSE
              RWSCL(ISHAPE) = XRN_ALL(MESHN_ALL(ISHAPE),ISHAPE)*ALATC
            END IF
            WRITE(6,*) 'rmt  = ',RMTCL(ISHAPE)
            WRITE(6,*) 'rmax = ',RWSCL(ISHAPE)
            LCONSTRUCTED(ISHAPE) = .TRUE.
            NUMCONSTRUCTED = NUMCONSTRUCTED + 1

         ELSE  

            WRITE(6,1090)
            WRITE(6,1070) RMTCORE(IAT)
            LCONSTRUCTED(ISHAPE) = .FALSE.

         END IF
c        END IF

c     
c     Now the shape has an extra panel, and starts at mtradius. 
c     
         WRITE(6,*) 'Shapes construction finished.'

         WRITE(6,*) '<<<<<<<<<<<<<<<<<<<<<<--->>>>>>>>>>>>>>>>>>>>>'

 50   ENDDO

c------------------------------------------------------------------------------

C Divide large panels into smaller ones, if NSMALL is smaller
C than the number of points in the largest panel.
      DO ISHAPE = 1,NUMSHAPE
         CALL DIVPANELS(
     >        NFUN_ALL(ISHAPE),NMIN,NSMALL,     
     X        THETAS_ALL(1,1,ISHAPE),XRN_ALL(1,ISHAPE),DRN_ALL(1,ISHAPE)
     X        ,NPAN_ALL(ISHAPE),MESHN_ALL(ISHAPE),NM_ALL(1,ISHAPE))
         write(*,*) 'meshn',meshn_all(ishape)
      ENDDO
c------------------------------------------------------------------------------

c Shapes constructed, now write out

      IF (OPT('WRITEALL')) THEN ! Write shape for every atom
         WRITE(15,FMT='(I5)') NSITES - ISITEBEGIN + 1

         IF (OPT('IMPURITY')) THEN  ! write old-style scaling factor, not used any more
            WRITE(15,FMT='(E20.12)') 1.D0      
         ELSE
            DO I1 = ISITEBEGIN,NSITES,4 ! write old-style scaling factor, not used any more
               WRITE(15,FMT='(4E20.12)') 1.D0,1.D0,1.D0,1.D0 
            ENDDO
         ENDIF
         DO IAT = ISITEBEGIN,NSITES
            ISHAPE = IDSHAPE(IAT)
            IF (LCONSTRUCTED(ISHAPE)) THEN
               WRITE(*,*) 'Writing out shape',ISHAPE,' for atom ',IAT
               CALL WRITESHAPE(NPAN_ALL(ISHAPE),MESHN_ALL(ISHAPE),
     &         NM_ALL(1,ISHAPE),XRN_ALL(1,ISHAPE),DRN_ALL(1,ISHAPE),
     &         NFUN_ALL(ISHAPE),LMIFUN_ALL(1,ISHAPE),
     &         THETAS_ALL(1,1,ISHAPE),ISHAPE)
            ELSE
               WRITE(*,*) 'Skipping unconstructed shape',ISHAPE,
     &                    ' for atom ',IAT
            ENDIF
         ENDDO
      ELSE ! write shape for representative atoms only.
         WRITE(15,FMT='(I5)') NUMCONSTRUCTED
         DO ISHAPE = 1,NUMCONSTRUCTED,4
            WRITE(15,FMT='(4E20.12)') 1.D0,1.D0,1.D0,1.D0
         ENDDO
         DO ISHAPE = 1,NUMSHAPE
            IF (LCONSTRUCTED(ISHAPE)) THEN
               WRITE(*,*) 'Writing out shape',ISHAPE
               CALL WRITESHAPE(NPAN_ALL(ISHAPE),MESHN_ALL(ISHAPE),
     &         NM_ALL(1,ISHAPE),XRN_ALL(1,ISHAPE),DRN_ALL(1,ISHAPE),
     &         NFUN_ALL(ISHAPE),LMIFUN_ALL(1,ISHAPE),
     &         THETAS_ALL(1,1,ISHAPE),ISHAPE)
            ELSE
               WRITE(*,*) 'Skipping unconstructed shape',ISHAPE
            ENDIF
         ENDDO
      ENDIF !  (OPT('IMPURITY'))



      END IF                    ! (INS.GT.0)
cINS-INS-INS-INS-INS-INS-INS-INS-INS-INS-INS-INS-INS-INS-INS-INS-INS-INS-INS-INS-INS-
c**********************************************************************************


c     
c     Print general information about lattice
c     
      WRITE(6,*) '*******************************************'
      WRITE(6,*) ' Analyzing lattice finised ...     '
      WRITE(6,*) ' Number of SITES..............',NSITES
      WRITE(6,*) ' Number of CLUSTERS ..........',NCLS
      WRITE(6,*) ' Number of CELLS..............',NUMCELL
      WRITE(6,*) ' Number of SHAPES.............',NUMSHAPE
      WRITE(6,*) ' Number of constructed shapes.',NUMCONSTRUCTED

      DO IAT = ISITEBEGIN,NSITES
         WRITE(6,901) IDSHAPE(IAT),CLS(IAT),DX(IAT),DY(IAT),DZ(IAT),IAT
      END DO
         
      WRITE(6,*) ' ****** Volume for each atom ******' 
      IF (INS.GT.0) 
     &     WRITE(6,*) ' RMAX is the maximum radius of the shape, 
     & RMT is updated'
      IF (INS.EQ.0) WRITE(6,*) 'RMAX is the Atomic Sphere radius'
      VTOT = 0.D0
         
      DO IAT = ISITEBEGIN,NSITES
         WRITE(6,1150) IAT,VOLUMECL(IDSHAPE(IAT)),
     &        RMTCL(IDSHAPE(IAT)),RWSCL(IDSHAPE(IAT))
         VTOT = VTOT + VOLUMECL(IDSHAPE(IAT))
      END DO
         
      WRITE(6,1160) VTOT,VTOT*ALATC*ALATC*ALATC
      WRITE(6,1170) VOLUC
      WRITE(6,*)
      WRITE(6,*) '+++++++++++ Note about volumes +++++++++++++++++'    
      WRITE(6,*) 
     &     ' Volume is probably different if you are in 2d mode '
      WRITE(6,*) ' In case of 3d mode and volume inconsistency '
      WRITE(6,*) ' try increasing the cluster atoms by changing '
      WRITE(6,*) ' RCLUSTZ and RCLUSTXY  in the inputcard'
      WRITE(6,*) '++++++++++++++++++++++++++++++++++++++++'
      WRITE(6,*)
      WRITE(6,*)
c     
c     ----------------------------------------------------------------------
c     
c         goto  88    ! skip potentials
      WRITE(6,*) 'Preparing potentials ................. '
      IF (I13(1:7).NE.'       ') THEN
                                                                        !    2014/06 Benedikt
         POT_EXISTS = .FALSE.                                           !    2014/06 Benedikt
         INQUIRE(FILE=TRIM(ADJUSTL(I13)),EXIST=POT_EXISTS)              !    2014/06 Benedikt
         IF (.NOT. POT_EXISTS) THEN                                     !    2014/06 Benedikt
           WRITE(6,*)                                                   !    2014/06 Benedikt
     &        ' No potential file with the name ', TRIM(ADJUSTL(I13)),  !    2014/06 Benedikt
     &        ' (as given in inputcard) found.'                         !    2014/06 Benedikt
           WRITE(6,*)                                                   !    2014/06 Benedikt
     &        ' I will use jellium potentials instead.'                 !    2014/06 Benedikt
                                                                        !    2014/06 Benedikt
           LJELL = .TRUE.                                               !    2014/06 Benedikt
                                                                        !    2014/06 Benedikt
         ELSE                                                           !    2014/06 Benedikt
           WRITE(6,1190) I13
           WRITE(6,*) ' If potentials do not match program will stop '
           WRITE(6,*) ' Expect to find the following potentials '
            
           DO IAT = ISITEBEGIN,NSITES
              CALL ELEMENTDATABASE(ZATOM(IAT),ELEM_NAME)
              IF (NSPIN.EQ.2) THEN
                 WRITE(6,1200) ELEM_NAME
                 WRITE(6,1210) ELEM_NAME
              ELSE 
                 WRITE(6,1220) ELEM_NAME
              END IF
           END DO
            
           LJELL = .FALSE.

         END IF                                                         !    2014/06 Benedikt
            
      ELSE
         WRITE(6,*) 
     &      ' No potential file specified I will use jellium potentials'
         LJELL = .TRUE.   
      END IF
         
      IF (INS.EQ.0) THEN 
         WRITE(6,*) ' Spherical Potential Output '
      ELSE
         WRITE(6,*) ' Full  Potential Output '
      END IF
c     
c     
c     
c Define param. IRNS and IMT. Check consistency. 
c If unreasonable, set it equal to core-radius point.
      DO IAT = 1,NSITES
        WRITE(*,*) 'MESHN_ALL(IDSHAPE(IAT))',
     &   IAT,IDSHAPE(IAT),MESHN_ALL(IDSHAPE(IAT))
         IRWS(IAT) = IMT(IAT) + MESHN_ALL(IDSHAPE(IAT))
         IF (IRWS(IAT).GT.IRMD) THEN
            WRITE(*,*) 'INCREASE IRMD IN inc.geometry',IRMD,IRWS(IAT)
            STOP 'IRMD'
         ENDIF

c         IMT(IAT) = IRWS(IAT) - MESHN_ALL(IDSHAPE(IAT))

         FPRADIUS(IAT) = STARTFP * RMT0_ALL(IAT)  ! Not yet scaled with latt. parameter
         BOUT = RMT0_ALL(IAT) / 
     &        ( EXP( AOUT_ALL(IAT)*DFLOAT(IMT(IAT)-1) ) - 1.0D0 )
         IRMIN(IAT) = NINT( DLOG(FPRADIUS(IAT)/BOUT)/AOUT_ALL(IAT) ) + 1
         IRNS(IAT) = IRWS(IAT) - IRMIN(IAT)
         IF (IRNS(IAT).LT.MESHN_ALL(IDSHAPE(IAT))) 
     &                  IRNS(IAT) = MESHN_ALL(IDSHAPE(IAT))
         FPRADIUS(IAT) = FPRADIUS(IAT) * ALATC ! Scale with latt. parameter

      ENDDO


      IF (.NOT.LJELL) THEN
         CALL GENPOTSTART(NSPIN,11,I13,INS,ISITEBEGIN,NSITES,ZATOM,
     &        SITEAT,KSHAPE,IDSHAPE,VOLUMECL,LPOT,AOUT_ALL,RWSCL,RMTCL,
     &        RMTCORE,MESHN_ALL,XRN_ALL,DRN_ALL,
     &        IRWS,IRNS,ALATC,
     &        QBOUND,KXC,TXC,LJELL)
      ENDIF
c     Using jellium potentials from dir:  jellpot
c
      IF (LJELL) THEN 
         CALL JELLSTART(NSPIN,11,I13,INS,NBEGIN,NATOMS,ZATOM,SITEAT,
     &        KSHAPE,IDSHAPE,VOLUMECL,LPOT,AOUT_ALL,RWSCL,RMTCL,
     &        RMTCORE,MESHN_ALL,XRN_ALL,DRN_ALL,THETAS_ALL,LMIFUN_ALL,
     &        NFUN_ALL,IRWS,IRNS,ALATC,
     &        QBOUND,KXC,TXC, EFSET)   
      END IF
      CLOSE(11)
c     
c     Everything finished now make some nice pictures !
c     
c     -----------------------------------------------------------
c     This is for drawing with povray
c     -----------------------------------------------------------
 88   CONTINUE                  ! skiping potentials


c Write out cell info.

      OPEN(69,FILE='vertices.dat',FORM='FORMATTED')
      call version_print_header(69)
      DO IAT = 1,NSITES
         IF (CELLREF(IAT).EQ.IAT) THEN
            WRITE(69,FMT='("# Cell of representative atom:",I5)') IAT
            DO IFACE = 1,NFACE_ALL(IAT)
               WRITE(69,FMT='("# Face:",I5)') IFACE
               DO IVERT = 1,NVERT_ALL(IFACE,IAT)
                  WRITE(69,FMT='(3F10.5)') XVERT_ALL(IVERT,IFACE,IAT),
     &                                     YVERT_ALL(IVERT,IFACE,IAT),
     &                                     ZVERT_ALL(IVERT,IFACE,IAT)
               ENDDO
               WRITE(69,FMT='(3F10.5)') XVERT_ALL(1,IFACE,IAT),
     &                                  YVERT_ALL(1,IFACE,IAT),
     &                                  ZVERT_ALL(1,IFACE,IAT)
               WRITE(69,FMT='(A2)') '  '
               WRITE(69,FMT='(A2)') '  '
            ENDDO
         ENDIF
      ENDDO
      CLOSE(69)

c Write out inner and outer radii per atom.
      OPEN(69,FILE='radii.dat',FORM='FORMATTED')
      call version_print_header(69)
      WRITE(69,FMT='(2A52)') 
     &   '# Radii after shift of cell centers (units of alat),',
     &   ' ratio Rmt0/Rout, NN-dist. before shift, ratio      '
      WRITE(69,FMT='(2A48)') 
     &   '# IAT    Rmt0           Rout            Ratio(%)',
     &   '   dist(NN)      Rout/dist(NN) (%)              '
      DO IAT = 1,NSITES
         WRITE(69,FMT='(I5,2F15.10,F12.2,F15.10,F12.2)') 
     &   IAT,RMT0_ALL(IAT),ROUT_ALL(IAT),100*RMT0_ALL(IAT)/ROUT_ALL(IAT)
     &   ,DISTNN(IAT),100*ROUT_ALL(IAT)/DISTNN(IAT)
      ENDDO
      CLOSE(69)

c Write out cell details to be used as input by KKR program that
c includes shape calculation.
      ! First find highest number of faces and vertices 
      ! to assist array allocation in other program
      NVERTMAX = 0
      NFACEMAX = 0
      DO IAT = ISITEBEGIN,NSITES
         IF (NFACE_ALL(IAT).GT.NFACEMAX) NFACEMAX = NFACE_ALL(IAT)
         DO IFACE = 1,NFACE_ALL(IAT)
            IF (NVERT_ALL(IFACE,IAT).GT.NVERTMAX) 
     &                             NVERTMAX = NVERT_ALL(IFACE,IAT)
         ENDDO
      ENDDO

      OPEN(69,FILE='cellinfo.dat',FORM='FORMATTED')
      call version_print_header(69)
      WRITE(69,FMT='("# Faces and vertices for all atoms")')
      WRITE(69,FMT='("# Number of atoms:")')
      WRITE(69,FMT='(2I6)') NSITES - ISITEBEGIN + 1
      WRITE(69,FMT='("# Largest number of faces and vertices:")')
      WRITE(69,FMT='(2I6)') NFACEMAX,NVERTMAX

      DO IAT = ISITEBEGIN,NSITES
         WRITE(69,FMT='(I6,"  Atom: Z=",F6.2)') 
     &                      IAT - ISITEBEGIN + 1,ZATOM(IAT)
         WRITE(69,FMT='(2I6,"  Faces")') NFACE_ALL(IAT)
         DO IFACE = 1,NFACE_ALL(IAT)
            WRITE(69,FMT='(I6,4E16.8)') IFACE,A3_ALL(IFACE,IAT),
     &      B3_ALL(IFACE,IAT),C3_ALL(IFACE,IAT),D3_ALL(IFACE,IAT)
            WRITE(69,FMT='(I6,"  Vertices")') NVERT_ALL(IFACE,IAT)
            DO IVERT = 1,NVERT_ALL(IFACE,IAT)
               WRITE(69,FMT='(I6,4E16.8)') IVERT,
     &         XVERT_ALL(IVERT,IFACE,IAT),YVERT_ALL(IVERT,IFACE,IAT),
     &         ZVERT_ALL(IVERT,IFACE,IAT)
            ENDDO
         ENDDO
      ENDDO
      CLOSE(69)

         
c Find clusters for use in TB-KKR input. These can be different from the
c Voronoi clusters because the reference-potential type should be compared.
c The reference-potential type is assumed to be uniquely given by the radius
c RMTREF that can be >= 0. Subroutine clsgen_tb is used. It is similar to 
c clsgen_voronoi but it accounts for RMTREF.
c If RMTREF is not given in the input it is set to depend on the values
c of the touching rmt (see earlier in this routine).
c First return rbasis to the shifted positions (they were pulled back to the
c unshifted in sub. readinput).

      RBASIS(1,1:NAEZ) = RBASIS(1,1:NAEZ) + DX(1:NAEZ)
      RBASIS(2,1:NAEZ) = RBASIS(2,1:NAEZ) + DY(1:NAEZ)
      RBASIS(3,1:NAEZ) = RBASIS(3,1:NAEZ) + DZ(1:NAEZ)
      CALL CLSGEN_TB(NAEZ,NEMB,RR,NR,RBASIS,
     &        KAOEZ,ZATOM,CLS,NCLS,
     &        NACLS,ATOM,EZOA,
     &        NLBASIS,NRBASIS,NLEFT,NRIGHT,ZPERLEFT,
     &        ZPERIGHT,TLEFT,TRIGHT,LMTREF,RMTREF,IREFPOT,
     &        RCLS,RCUTZ,RCUTXY,LINTERFACE,
     &        ALATC)


!---------------------------------------------------------------
!-- Write out some results for the kkr input
      OPEN(69,FILE='atominfo.txt',FORM='FORMATTED')
      call version_print_header(69)

      WRITE(69,FMT='(A8)') 'ATOMINFO'
      WRITE(69,2003)
      DO IAT = 1,NAEZ
         WRITE(69,FMT='(F5.2,5I2,3I7,F6.0,I8,2F12.7)') 
     &        ZATOM(IAT),1,3,3,0,0,CLS(IAT),IREFPOT(IAT),
     &        IDSHAPE(IAT),1.D0,IRNS(IAT),RMTREF(IAT),FPRADIUS(IAT)
      ENDDO
      WRITE(69,*) 'NCLS=',MAXVAL(CLS(1:NAEZ))

      IF (LINTERFACE) THEN
         WRITE(69,FMT='(A10,43X,A3,A12)') 
     &        'LEFTBASIS ','CLS','   LEFTMTREF'
         DO JAT = 1,NLBASIS
            IAT = NAEZ + JAT
            WRITE(69,FMT='(3E17.8,I5,F12.5)') 
     &            TLEFT(1:3,JAT),CLS(IAT),RMTREF(IAT)
         ENDDO
         WRITE(69,FMT='(A10,43X,A3,A12)') 
     &        'RIGHBASIS ','CLS','   RIGHMTREF'
         DO JAT = 1,NRBASIS
            IAT = NAEZ + NLBASIS + JAT
            WRITE(69,FMT='(3E17.8,I5,F12.5)') 
     &             TRIGHT(1:3,JAT),CLS(IAT),RMTREF(IAT)
         ENDDO
      ENDIF

      WRITE(69,*) ' '
      WRITE(69,FMT='(A8)') '<SHAPE> '
      IF (NATYP.EQ.NAEZ) THEN
         WRITE(69,FMT='(I6)')  (IDSHAPE(IAT),IAT=1,NAEZ)
      ELSE
         WRITE(69,FMT='(2I6)')  
     &        (SITEAT(IAT),IDSHAPE(SITEAT(IAT)),IAT=1,NATYP) ! CPA, NATYP>NAEZ
      ENDIF
         

      CLOSE(69)
!---------------------------------------------------------------



      VOROPLOT = .FALSE. 
      IF (VOROPLOT) THEN
         STOP 'Arrays for povray output not correct in version Jan 2012'
         WRITE(6,*) ' NOW PREPARING VORONOI PICTURES '
        
         OPEN(UNIT=12,STATUS='unknown',FILE='voronoi.pov')
         WRITE(12,201)  
         
         DO IAT = ISITEBEGIN, NSITES

            ICL = IDSHAPE(IAT)  !CLS(IAT)
            RX = RBASIS(1,IAT)
            RY = RBASIS(2,IAT)
            RZ = RBASIS(3,IAT)
            NFACE = NFACECL(ICL) 
c     
c     
C     The equation        A1*x + A2*y + A3*z = A4
c     is transformed to   A1*x'+ A2*y'+ A3*z'= A4 -(A1*dx+A2*dy+A3dz)
c     
c     
c     A4 = A4 - A1*NEW0(1) - A2*NEW0(2) - A3* NEW0(3)
c
            DO I=1,NFACE
               NVERT(I) = NVERTCL(I,ICL)
               DO J=1,NVERT(I)   
                  XVERT(J,I) = (XVERTCL(J,I,ICL) + RX)
                  YVERT(J,I) = (YVERTCL(J,I,ICL) + RY)
                  ZVERT(J,I) = (ZVERTCL(J,I,ICL) + RZ)
               ENDDO
            ENDDO     
            DO I=1,NFACE
               
               WRITE(12,99) NVERT(I)+1
               DO J=1,NVERT(I)
                  WRITE(12,100) XVERT(J,I),YVERT(J,I),ZVERT(J,I)
               ENDDO
               WRITE(12,101) XVERT(1,I),YVERT(1,I),ZVERT(1,I)
               WRITE(12,102)
            ENDDO
            CRAD = 0.01D0
            DO I=1,NFACE
               NA = NVERT(I)
               DO J=1,NA-1
                  WRITE(12,500)
                  WRITE(12,100) XVERT(J,I),YVERT(J,I),ZVERT(J,I)
                  WRITE(12,100) XVERT(J+1,I),YVERT(J+1,I),ZVERT(J+1,I)
                  WRITE(12,502) CRAD
                  WRITE(12,503)
               ENDDO
               WRITE(12,500)
               WRITE(12,100) XVERT(NA,I),YVERT(NA,I),ZVERT(NA,I)
               WRITE(12,100) XVERT(1,I),YVERT(1,I),ZVERT(1,I)
               WRITE(12,502) CRAD
               WRITE(12,503)
            ENDDO

         ENDDO     
   
         WRITE(12,200)
         WRITE(12,202)
         WRITE(12,203)
      ENDIF

c     -------------------------------------------------------
 500     format('cylinder { ')
 502     format(F10.5)
 503     format(' texture { pigment { color Red }}}')
         
 99      format('polygon { ' /
     &        I5, ',')
 100     format('< ',F10.5,',',F10.5,',',F10.5,' > ,')
 101     format('< ',F10.5,',',F10.5,',',F10.5,' > ')
 102     format('pigment { color rgb <1, 0, 0> transmit .95 } }')
 200     format('camera {' /
     &        '   location < 4 , 3, 2 > ' /
     &        '   look_at < 0, 0, 1 > }')
 201     format('#include "colors.inc" ' /
     &        ' background {color White }')
 202     format(' light_source { < 5, 3, -10 > color White }')
 203     format(' light_source { < 4, 8, 10 > color Red }')
 900     format( ' ****** NEW POLYHEDRON AND SHAPE FUNCTION ******'//)
 901     format('Shape...',I4,' has cluster..',I4,' extra shift..',3F8.4
     &        ,' for atom..',I4)
 1000    format('Volume ....(alat^3).......',F14.8)
 1010    format('ASA Sphere (a.u.).........',F14.8)
 1030    format('MT  Sphere (a.u.).........',F14.8)
 1020    format('Faces for polyhedron......',I8)
 1050    format(/'Your shape function will be updated '/
     &        'to match to the smaller Muffin Tin sphere.')
 1070    format('New RMT (a.u.)............',F14.8)
 1080    format('Percentage Change ........',F8.2)
 1090    format('The Shape Function was not updated ' /
     &        'Your Muffin-Tin Sphere was too big or zero' /
     &        'or you are doing ASA calculation' )
 1112    format(' ASA Radius  set to ........',F12.8/
     &        '***************************************'/)
 1150    format
     & (' Atom ..',I5,' Volume(alat^3)  : ',F12.8,
     &        '   RMT : ',F8.5,'   RMAX : ',F8.5) 
 1160    format(/'Total volume (alat^3)    ...',F16.8/
     &          'Total Volume (a.u.^3)   ....',F16.8)
 1170    format(/'Bravais cross prod.(alat^3):',F16.8)
 1190    format('Potential file will be used ........ ',A40/
     &        '**** MUST BE IN GENERAL FORMAT ****'/) 
 1200    format(A3,'  POTENTIAL SPIN DOWN ')
 1210    format(A3,'  POTENTIAL SPIN UP   ' ) 
 1220    format(A3,'  POTENTIAL ')
 2000    format
     &('**********************************************************')
 2001    format
     &('*         SCREENED KKR POTENTIAL PREPARATION  UTILITY    *')
 2002    format
     &('*                                                        *')
 2003    FORMAT(
     & 'ZAT   LMXC  KFG   <CLS> <REFPOT> <NTC>  FAC  <IRNS> <RMTREF>',
     & '   <FPRADIUS>')
    

         END 
      
      

