c ************************************************************************
      SUBROUTINE CLSGENIMP12(
     >        NUMIMP,RIMPURITY,NKILLATOM,RKILL,
     >        CLS,NCLS,NACLS,ATOM,RCLS,
     >        BRAVAIS,RECBV,NAEZ,RBASIS,RCUTZ,RCUTXY,
     <        CLSIMP,NCLSIMP,NACLSIMP,ATOMIMP,RCLSIMP,RMTHLFIMP)
      implicit none
c#@# KKRtags: VORONOI KKRimp geometry
c ************************************************************************
c This subroutine is used to create the clusters around each atom 
c where in order to prepare the shape functions and produce povray 
c input files. (Based on clsgen99.f)
c
c STRATEGY : 
c Calculate the cluster of each atom by the lattice
c parameters avaliable. Sort the atoms in a unique way :big r, big z, big y
c compare the positions with the previous clusters to see if there is 
c a difference. If not keep only previous clusters and make indexing if
c a new cluster is found then check dimensions and continue for the new
c atom.  
c
c
      include 'inc.geometry'    
c
c
c     .. arguments
c
c
      
      INTEGER
     +     CLS(*),                  ! sort of cluster around atom
     +     NACLS(*),                ! number of atoms in cluster
     +     ATOM(NACLSD,*)           ! Index to atom in elem/cell at site in cluster
c Input:
      INTEGER NAEZ                    ! number of host-basis atoms
      INTEGER NUMIMP                  ! Number of impurity atoms
      INTEGER NKILLATOM               ! Number of host-sites to be "killed"
      INTEGER NCLS                    ! number of different host clusters
c 
      REAL*8 BRAVAIS(3,3),RECBV(3,3)  ! Bravais vectors of real and reciprocal lattice
      REAL*8 RBASIS(3,*)              ! Positions of host-basis atoms
      REAL*8 RIMPURITY(3,*)           ! positions of impurity sites
      REAL*8 RCLS(3,NACLSD,*)         ! real space position of atom in host cluster
      REAL*8 RKILL(3,*)               ! Positions of sites to be "killed"
      REAL*8 RCUTXY,RCUTZ             ! Cutoff radius of cluster
      
c
c Output:
      INTEGER NCLSIMP                     ! Number of inequivalent clusters
      INTEGER NACLSIMP(*)                 ! Number of sites in each cluster
      INTEGER CLSIMP(*)                   ! Type of cluster of each atom
      INTEGER ATOMIMP(NACLSD,*)           ! Index to atom in elem/cell at site in cluster
      REAL*8 RCLSIMP(3,NACLSD,*)          ! real space position of atom in imp. cluster
      REAL*8 RMTHLFIMP(NIMPD)             ! Muffin tin radius by bisection to nearest neighbour
c
c Local:
      INTEGER IATOM(NACLSD),ATOM1(NACLSD) ! Atom-index in cluster 
      INTEGER ISORT(NACLSD)           ! Auxiliary array for sorting clusters
      INTEGER NIMPCLS ! Number of atoms in impurity cluster
      INTEGER IIMP,JIMP,IBASIS,JKILL,JVEC,IX,NB1,NB2,NB3,ICLS,JCLS,POS
     +     ,IA,IN,IB,IAT,I1
      INTEGER NEAREST                 ! Nearest found host-basis to certain imp-atom
      REAL*8 RG(3,NACLSD),RS1(NACLSD)             ! Auxiliary array for sorting clusters
      REAL*8 RSORT(NACLSD)            ! Parameter for sorting sites in cluster
      REAL*8 RLATT(3)                 ! Lattice vector
      REAL*8 RNEAREST(3)              ! Position of nearest host-atom to imp-atom
      REAL*8 RCURRENT(3)              ! Position of a site.
      REAL*8 DISTSQ,DISTMINSQ         ! Distance**2 between atoms
      REAL*8 TOLIMP                   ! If imp. atom is closer than this to host atom, 
                                      ! it is considered to replace it
      

      LOGICAL LACCEPT,LFOUND ! Accept-flag,Found-flag
c
      REAL*8 EPSSHL,DELTA(3),RCLS1(3,NACLSD)
      REAL*8 RCUT2,RCUTXY2,R2
c
      LOGICAL  L2DIM,CLUSTCOMP_VORONOI
      CHARACTER*256 UIO
      INTEGER IER
c
c
      LOGICAL TEST,OPT
      EXTERNAL DSORT,CLUSTCOMP_VORONOI,IOINPUT
      INTRINSIC MIN,SQRT
c
      DATA     EPSSHL   / 1.0D-4 /
      DATA     TOLIMP   / 1.D-6  / 


      WRITE(6,*) 
     &     '>>> CLSGENIMP_KKRFLEX: generation of cluster coordinates'

c  Cutoffs for max. cluster radius
      RCUTXY2 = (RCUTXY+EPSSHL)*(RCUTXY+EPSSHL)
      RCUT2   = (RCUTZ+EPSSHL)*(RCUTZ+EPSSHL)
      RCUT2 = MAX(RCUTXY2,RCUT2)

c Initialize number of non-equivalent impurity clusters.
      NCLSIMP = 0

      CALL IoInput('TOLIMP          ',UIO,1,7,IER)
      IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) TOLIMP
      ! Default is set by data-statement above.

c
c Strategy:
c Loop over impurities
c   1  Find a near crystal point to impurity (preferably the nearest)
c   2  Span cluster of all impurities
c   3  Span host-cluster around this crystal point
c   -  Exclude common points of host and impurity
c   -  Exclude "killed atoms".
c   4  Center cluster around impurity.
c End loop over impurities
c

      DO 999 IIMP = 1,NUMIMP   ! Loop over impurities

c 1  Find closest lattice point to impurity. 
c Use reciprocal lattice vectors (remember they are given in units 2 pi/a).
c   
      NB1 = NINT( RIMPURITY(1,IIMP)*RECBV(1,1) + 
     &            RIMPURITY(2,IIMP)*RECBV(2,1) + 
     &            RIMPURITY(3,IIMP)*RECBV(3,1)   )
      NB2 = NINT( RIMPURITY(1,IIMP)*RECBV(1,2) + 
     &            RIMPURITY(2,IIMP)*RECBV(2,2) + 
     &            RIMPURITY(3,IIMP)*RECBV(3,2)   )
      NB3 = NINT( RIMPURITY(1,IIMP)*RECBV(1,3) + 
     &            RIMPURITY(2,IIMP)*RECBV(2,3) + 
     &            RIMPURITY(3,IIMP)*RECBV(3,3)   )
      DO IX = 1,3
         RLATT(IX) = NB1 * BRAVAIS(IX,1) + 
     &               NB2 * BRAVAIS(IX,2) + NB3 * BRAVAIS(IX,3)
         DELTA(IX) = RIMPURITY(IX,IIMP) - RLATT(IX)
      ENDDO
c
c The closest lattice point is RLATTT = NB1 * Bravais1 + NB2 * Bravais2 + NB3 * Bravais3.
c The impurity position relative to this point is DELTA.
c Now find nearest basis vector relative to this lattice point.
c
      NEAREST = 1
      DISTMINSQ = 1.D10
      DO IBASIS = 1,NAEZ
         DISTSQ = (DELTA(1) - RBASIS(1,IBASIS))**2 + 
     &            (DELTA(2) - RBASIS(2,IBASIS))**2 +
     &            (DELTA(3) - RBASIS(3,IBASIS))**2
         IF (DISTSQ.LT.DISTMINSQ) THEN
            NEAREST = IBASIS
            DISTMINSQ = DISTSQ
         ENDIF
      ENDDO

      RNEAREST(:) = RLATT(:) + RBASIS(:,NEAREST)
c
c
c 2  Span impurity-cluster around impurity.
c
      NIMPCLS = 1 ! Initialize number of sites in cluster
      RCLS1(:,1) = RIMPURITY(:,IIMP) ! First vector is always the site itself.

      DO JIMP = 1,NUMIMP  ! Loop over (other) impurities
         LACCEPT = .TRUE.
         DISTSQ = (RIMPURITY(1,IIMP) - RIMPURITY(1,JIMP))**2 +
     &            (RIMPURITY(2,IIMP) - RIMPURITY(2,JIMP))**2 +
     &            (RIMPURITY(3,IIMP) - RIMPURITY(3,JIMP))**2
         IF (DISTSQ.GT.RCUT2) LACCEPT = .FALSE. ! Impurity site is considered to be too far away
         IF (IIMP.EQ.JIMP) LACCEPT = .FALSE.    ! Avoid self
         IF (LACCEPT) THEN
            NIMPCLS = NIMPCLS + 1                   ! Add to cluster one more site...
            RCLS1(:,NIMPCLS) = RIMPURITY(:,JIMP)    ! ... and its coordinates
            ATOM1(NIMPCLS) = NAEZ + JIMP            ! Index to store what kind of atoms are in the cluster.
                                                    ! (value between 1 and NAEZ refers to host)
         ENDIF
      ENDDO
c
c 3  Span host-cluster around impurity excluding common atoms with impurity-cluster and killed atoms.
c
      DO JVEC = 1,NACLS(CLS(NEAREST))    ! Loop over cluster-sites of nearest host-atom
                                         ! Remember, CLS(I) is the cluster-index of basis-atom I.
         RCURRENT(:) = RNEAREST(:) + RCLS(:,JVEC,CLS(NEAREST))  ! RCLS is the "reduced" position w.resp.
                                                                ! to cluster-center. For JVEC=1, RCLS=(0,0,0).
         LACCEPT = .TRUE.

C        Check if host atom is identical with any of the impurities
         DO JIMP = 1,NUMIMP       ! Loop over impurity atoms 
            DISTSQ = (RCURRENT(1) - RIMPURITY(1,JIMP))**2 +
     &               (RCURRENT(2) - RIMPURITY(2,JIMP))**2 +
     &               (RCURRENT(3) - RIMPURITY(3,JIMP))**2
            IF (DISTSQ.LT.TOLIMP) LACCEPT = .FALSE. ! Host site is considered identical with impurity site
         ENDDO

C        Check if host atom is identical with any of the killed atoms
         DO JKILL = 1,NKILLATOM    ! Loop over killed atoms 
            DISTSQ = (RCURRENT(1) - RKILL(1,JKILL))**2 +
     &               (RCURRENT(2) - RKILL(2,JKILL))**2 +
     &               (RCURRENT(3) - RKILL(3,JKILL))**2
            IF (DISTSQ.LT.1.D-8) LACCEPT = .FALSE. ! Host site is considered identical with killed site
         ENDDO

c        Check if host atom is too far from current impurity.
         DISTSQ = (RCURRENT(1) - RIMPURITY(1,IIMP))**2 +
     &            (RCURRENT(2) - RIMPURITY(2,IIMP))**2 +
     &            (RCURRENT(3) - RIMPURITY(3,IIMP))**2
         IF (DISTSQ.GT.RCUT2) LACCEPT = .FALSE. ! Host site is considered too far away


         IF (LACCEPT) THEN
            NIMPCLS = NIMPCLS + 1
            RCLS1(:,NIMPCLS) = RCURRENT(:)
            ATOM1(NIMPCLS) = ATOM(JVEC,NEAREST)    ! Index to store what kind of atoms are in the cluster.
         ENDIF

      ENDDO

c Now the impurity IIMP has a cluster of NIMPCLS sites around it, with positions at RCLS1.
c Store number of sites in array NACLSIMP.
      NACLSIMP(IIMP) = NIMPCLS
c 4   Center cluster around impurity (coordinate shift)
      DO JVEC = 1,NIMPCLS
         RCLS1(:,JVEC) = RCLS1(:,JVEC) - RIMPURITY(:,IIMP)
      ENDDO

c-------------------------------------------------------------------------------------
c sort sort sort sort sort sort sort sort sort sort sort sort sort sort sort sort sort 
c Sort cluster sites according to: radius, then z, then y, then x.
      DO IA = 1,NIMPCLS
         RSORT(IA) = SQRT(RCLS1(1,IA)**2+RCLS1(2,IA)**2+RCLS1(3,IA)**2)
         RSORT(IA) = 100000000.D0*RSORT(IA)+
     &                  100000.D0*RCLS1(3,IA)+
     &                     100.D0*RCLS1(2,IA)+
     &                      0.1D0*RCLS1(1,IA)
      ENDDO      
C     

      CALL DSORT(RSORT,ISORT,NIMPCLS,POS)
c     Rearange exchange IA with IB
c     Map temporarily to another array         
      DO IA = 1,NIMPCLS       
         RG(:,IA) = RCLS1(:,IA)
         IATOM(IA) = ATOM1(IA)
         RS1(IA) = RSORT(IA)
      END DO    
c     Now use correct order
      DO IA = 1,NIMPCLS
         IB = ISORT(IA)
         RCLS1(:,IA) = RG(:,IB)
         ATOM1(IA) = IATOM(IB)
         RSORT(IA) = RS1(IB)
      END DO


c sort sort sort sort sort sort sort sort sort sort sort sort sort sort sort sort sort  
c-------------------------------------------------------------------------------------
c
c Store array ATOM1. The I'th site of cluster around IIMP contains atom-type ATOMIMP(I,IIMP).
c Remember, the indexing contains first all host-basis atoms, then all impurity atoms.
      ATOMIMP(:,IIMP) = ATOM1(:)


c-------------------------------------------------------------------------------------
c    Now compare new cluster to previously found and keep only if different.

      LFOUND = .FALSE.
      DO ICLS = 1,NCLS                     ! First check host-clusters
         IF ( CLUSTCOMP_VORONOI(RCLS,ICLS,NACLS(ICLS),RCLS1,NIMPCLS) ) 
     &      THEN
            CLSIMP(IIMP) = ICLS            ! Impurity IIMP is mapped to host cluster index
            LFOUND = .TRUE.
         ENDIF
      ENDDO

      IF (.NOT.LFOUND) THEN                ! Now check impurity clusters
         DO ICLS = 1,NCLSIMP
            IF ( CLUSTCOMP_VORONOI(RCLSIMP,ICLS,
     &              NACLSIMP(ICLS),RCLS1,NIMPCLS) ) THEN  ! Cluster has been found before
               CLSIMP(IIMP) = ICLS + NCLS ! Impurity IIMP is mapped to imp. cluster index
               LFOUND = .TRUE.                
         ENDIF
      ENDDO
      ENDIF

      IF (.NOT.LFOUND) THEN                ! New cluster
         NCLSIMP = NCLSIMP + 1             ! Number of new clusters due to imp.
         CLSIMP(IIMP) = NCLSIMP + NCLS     ! Impurity IIMP is mapped to cluster index;
                                           ! plus NCLS to have a unified indexing of clusters
         NACLSIMP(NCLSIMP) = NIMPCLS       ! Number of atoms in this cluster
         RCLSIMP(:,:,NCLSIMP) = RCLS1(:,:) ! Positions in cluster
      END IF

c     All different clusters have been found
c-------------------------------------------------------------------------------------



 999  ENDDO                     ! IIMP = 1,NUMIMP     ! Loop over impurities


c Write out some info.
      WRITE(*,*) '>>>>>>> Info from clsgenimp_kkrflex'
      DO IIMP = 1,NUMIMP
         WRITE(*,*) 
     &        'Impurity ',IIMP,' has cluster ',CLSIMP(IIMP),
     &        ' with ',NACLSIMP(IIMP),' sites.'
      ENDDO


c Compare impurity clusters with host clusters for chek reasons.

      DO ICLS = 1,NCLSIMP       ! Loop over impurity clusters
         DO JCLS = 1,NCLS       ! Loop over host clusters
            IF ( CLUSTCOMP_VORONOI(RCLSIMP,ICLS,NACLSIMP(ICLS),
     &                     RCLS(1,1,JCLS),NACLS(JCLS))    ) 
     &           WRITE(*,*) 'Impurity cluster ',ICLS,
     &           ' is geometrically identical with host cluster ',JCLS
         ENDDO                  ! JCLS = 1,NCLS
      ENDDO                     ! ICLS = 1,NCLSIMP


c Now the cluster vectors are contained in the array RCLSIMP
c and the number of sites in the cluster in NACLSIMP.

! Find the muffin-tin of each atom by bisection of the distance to the
! nearest atom
      DO IIMP = 1,NUMIMP
         DISTMINSQ = 1.D100
         ICLS = CLSIMP(IIMP)
         DO JVEC = 2,NACLSIMP(ICLS)
            R2 = RCLSIMP(1,JVEC,ICLS)**2 + RCLSIMP(2,JVEC,ICLS)**2 + 
     &                                     RCLSIMP(3,JVEC,ICLS)**2
            IF (R2.LT.DISTMINSQ) DISTMINSQ = R2
         ENDDO
         RMTHLFIMP(IIMP) = DSQRT(DISTMINSQ)/2.D0
      ENDDO

      WRITE(6,*) ' Sub clsgen99  exiting <<<<<<<<<<<<<'

c ------------------------------------------------------------------------


c
c ------------------------------------------------------------------------

 1000   format(' cluster around atom     ',10I4/,
     +         ' number atoms in cluster ',10I4)
 1001   format(I4,2I5,3F8.2,I6,4f7.2)
 1002   format(' cocls : naez =',I3,' (x,y,z)= ',3F10.4)
 1010   format(12x,I6,3F10.4)
 1020   format('  Nr  naez kaoez     x       y       z',
     +         '  ezoa  RR(1)  RR(2)  RR(3)      R')
 1030   FORMAT(3I8)
 1040   FORMAT('Cu  ',3D24.16,2I8,F18.12)
 1050   FORMAT(3F12.7,'  scaling factor')
 1060   FORMAT(I4,3F12.7,'  center',/,(I4,3f12.7))
 1070   FORMAT(I4,3F12.7,'  center of gravity')
 1080   FORMAT('contains ',I4,'  atoms.')
 9005   FORMAT(I4)
 9010   FORMAT(('# Z     ',20F4.0))
 9020   FORMAT(('# KAOEZ ',20I4))
 9030   FORMAT(F12.7,6x,'ALAT')
 9040   FORMAT('> cluster ',I4,' at atom ',I4,
     +       ' of type ',I4,'.')
 9050   FORMAT('      **** COUPLING MATRIX ****',I3,' x ',I3)
 9060   FORMAT(I4,1X,200I1)
 9080   FORMAT(I6,3F15.6)
 9090   FORMAT('The Number of layers is each Principal Layer = ',
     &                I5)
 8000   format('sphere { ')
 8010   format ('<',F10.5,',',F10.5,',',F10.5,'> ,')
 8011   format ('<',F10.5,',',F10.5,',',F10.5,'> ')
 8020   format (F10.5)
 8030   format(' texture { ')
 8040   format('  pigment { color ',A24,' } } } ')
 8050   format('  pigment { color ',A24,'  } ')
c
 8500  format('cylinder { ')
 8502  format(F10.5)
 8503  format(' texture { pigment { color Red }}}')
c
 8099   format('polygon { ' /
     &    I5 ,',')
 8102  format('pigment { color rgb <1, 0, 0> transmit .95 } }')
 8200  format('camera {' /
     & '   location < 4 , 3, -4 > ' /
     & '   look_at < 0, 0, -1 > }')
 8201  format('#include "colors.inc" ' /
     &      ' background {color White }')
 8202  format(' light_source { < 5, 3, -10 > color White }')
 8203  format(' light_source { < 4, 8, 10 > color Red }')
 8600  format('finish{ ambient 0 '/
     &        '        diffuse 0 '/
     &        '        reflection 0.25 '/
     &        '        roughness 0.001 }}')
c             'interior{ior 1.33} }' )
      RETURN
      END


