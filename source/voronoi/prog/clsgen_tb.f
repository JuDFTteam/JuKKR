c ************************************************************************
      SUBROUTINE CLSGEN_TB(NAEZ,NEMB,RR,NR,RBASIS,
     &                   KAOEZ,ZAT,CLS,NCLS,
     &                   NACLS,ATOM,EZOA, 
     &                   NLBASIS,NRBASIS,NLEFT,NRIGHT,ZPERLEFT,ZPERIGHT,
     &                   TLEFT,TRIGHT,LMTREF,RMTREF,IREFPOT,
     &                   RCLS,RCUT,RCUTXY,L2DIM,
     &                   ALAT)
      use mod_version_info
      implicit none
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
      REAL*8       ALAT         ! lattice constant A
      REAL*8       RCUT,RCUTXY
      REAL*8      
     +     RBASIS(3,*),             ! pos. of basis atoms in EZ
     +     RCLS(3,NACLSD,*),        ! real space position of atom in cluster
     +     RR(3,0:NRD),             ! set of lattice vectors
     +     ZAT(*),                  ! nucleus charge
     +     RMTREF(*)
!     +     RWS(*),
!     &     BBOX(3)                  ! bounding box for povray plots
c
      INTEGER
     +     KMT,                     ! scaling of RMT with MTFAC (not used)
c                                   ! 0: RMT from crystal structure 
c                                   ! 1: RMT - " -  scaled by RMTFAC
c                                   ! 2: RMT = RMTFAC
c                                   ! 3: RMT from ref. pot. card
     +     NAEZ,                    ! number of atoms in EZ
     +     NEMB,                    ! number of embedding postions
     +     NCLS,                    ! number of diff. clusters
     +     NINEQ,                   ! number of nonequivalent atomic 
                                    ! positions in EZ
     +     NR                       ! number of lattice vectors RR
c
      INTEGER
     +     CLS(*),                  ! type of cluster around atom
     +     KAOEZ(*),                ! type of atom at position in EZ
     +     NACLS(*),                ! number of atoms in cluster
     +     ATOM(NACLSD,*),          ! index to atom in elem/cell at site in cluster
     +     EZOA (NACLSD,*)          ! index to bravais lattice  at site in cluster
c
c     .. locals
c
      INTEGER 
     +     AJ,C,ILAY,J,N1,INUM,ISUM,IR,INEI,ISITE,JSITE,IAT1,
     +     NA,NUMBER,MAXNUMBER,NC,NPRIN,ITEST1,ITEST,IX,
     +     POS,IA,IN,IB,II,JATOM,ICU,IC,IAT,I0,I1,ICLUSTER
      INTEGER IATOM(NACLSD),IEZOA(NACLSD),
     +     ISORT(NACLSD),ICOUPLMAT(NAEZD,NAEZD),
     &     IREP(NCLSD) ! representative atom of cluster (inverse of CLS)
      INTEGER IREFPOT(NAEZD+NEMBD),NREFPOT
c
      REAL*8        
     +     RAD,R1,R2,RABS,RD,T,EPSSHL,TOL,TOL2,DISTMIN,
     +     ASC(3),RCLS1(3,NACLSD),
     +     R0(3,20),RG(3,NACLSD),TMP(3),RSORT(NACLSD)
      INTEGER NLAY,                                
     +        NLBASIS,NRBASIS,                    
     +        NLEFT,NRIGHT           
      REAL*8                                   
     +        ZPERLEFT(3),ZPERIGHT(3),            
     +        TLEFT(3,*),TRIGHT(3,*)
      REAL*8       RCUT2,RCUTXY2,RXY2
c
      LOGICAL LFOUND
      LOGICAL  L2DIM,CLUSTCOMP_TB
c
c
      LOGICAL TEST,OPT,pov,LMTREF
      EXTERNAL DSORT,CLUSTCOMP_TB
      INTRINSIC MIN,SQRT
c
      DATA     EPSSHL   / 1.D-4 /
      DATA     TOL   / 1.D-7 /
      DATA     TOL2   / 1.D-7 /
c ------------------------------------------------------------------------
      WRITE(6,*) '>>> CLSGEN_TB: generation of cluster coordinates'
c This is generating the clusters which have a distance smaller
c than RCUT and RCUTXY in plane .
c The cluster atoms are ordered with radius and then z>y>x 
c The ordering allows an easy comparison of clusters
c The principal layer for each layer (atom in unit cell) is
c calculated also for each cluster and the maximum number
c is returned. Some dimension tests are also done      
      WRITE(6,*) 'RCUT = ',RCUT,' RCUTXY = ',RCUTXY
      IF (ABS(RCUTXY - RCUT).LT.1.D-4) THEN
          WRITE(6,*) 'Spherical Clusters are created'
c          LSPHER=.TRUE.
      END IF 
      OPEN(8,FILE='clusters',status='unknown')
      call version_print_header(8)
      WRITE(8,9005) NAEZ
      WRITE(8,9030) ALAT
      WRITE(8,9010) (ZAT(KAOEZ(IAT)),IAT=1,NAEZ)
      WRITE(8,9020) (KAOEZ(IAT),IAT=1,NAEZ)
      
      RCUTXY2 = RCUTXY**2
      RCUT2 = RCUT**2
c===============================================================
c Check the if the dimension NACLSD is enough before starting to
c assign clusters. Write out necessary dimension.
c Also find touching MT radius.
      DO 100 JATOM = 1,NAEZ + NEMB      ! loop in all atoms or layers

         MAXNUMBER = 0
         NUMBER = 0             ! counter for atoms in cluster
         DISTMIN = 1.D100       ! Large initial value for RMT**2
         DO NA = 1,NAEZ         ! loop in all atoms
            DO IR = 0,NR        ! loop in all bravais vectors    
               TMP(1:3) = RR(1:3,IR)+RBASIS(1:3,NA)-RBASIS(1:3,JATOM)
               RXY2 =  TMP(1)**2+TMP(2)**2
               R2   =  TMP(3)**2 + TMP(1)**2+TMP(2)**2
               IF ( (RXY2.LE.RCUTXY2).AND.(R2.LE.RCUT2) )  THEN
                  NUMBER = NUMBER + 1
               END IF
               IF (R2.GT.TOL2) DISTMIN = MIN(DISTMIN,R2)
            END DO              ! N loop in bravais
         END DO                 ! NA loop in NAEZ
  
c     
c     In the case of 2 dimensional case loop in the atoms outside.
c     
         IF (L2DIM) THEN
c     Somehow messy   
            DO IR = 0,NR
               DO ILAY = NLEFT,1,-1  ! loop in some layers on left side
                  DO I1 = NLBASIS,1,-1 ! loop in representative atoms on left side
                     TMP(1:3) = RR(1:3,IR) + TLEFT(1:3,I1) +
     &                    (ILAY-1)*ZPERLEFT(1:3) - RBASIS(1:3,JATOM)
                     RXY2 =  TMP(1)**2 + TMP(2)**2
                     R2  =   TMP(3)**2 + RXY2
                  
                     IF ((RXY2.LE.RCUTXY2).AND.(R2.LE.RCUT2)) THEN
                        NUMBER = NUMBER + 1
                     END IF
                     IF (R2.GT.TOL2) DISTMIN = MIN(DISTMIN,R2)
                  END DO 
               END DO
c
               DO ILAY = 1,NRIGHT
                  DO I1 = 1,NRBASIS
                     TMP(1:3) = RR(1:3,IR)+ TRIGHT(1:3,I1) + 
     &                    (ILAY-1)*ZPERIGHT(1:3) - RBASIS(1:3,JATOM)  
                     RXY2 = TMP(1)**2 + TMP(2)**2
                     R2  =  TMP(3)**2 + RXY2
                     IF ((RXY2.LE.RCUTXY2).AND.(R2.LE.RCUT2)) THEN
                        NUMBER = NUMBER + 1
                     END IF
                     IF (R2.GT.TOL2) DISTMIN = MIN(DISTMIN,R2)
                  END DO 
               END DO
c     
         
            END DO              ! loop in all bravais lattices
         END IF                 ! L2DIM Interface calculation

         DISTMIN = DSQRT(DISTMIN)/2.D0 ! Touching RMT
         ! Define MT-radius of TB-ref. potentials if not read-in from the input
         IF (.NOT.LMTREF) RMTREF(JATOM) = 
     &                    0.99D0*INT( DISTMIN*ALAT*100.D0)/100.D0 ! round up the decimals
c     
         MAXNUMBER = MAX(MAXNUMBER,NUMBER) ! Find largest cluster size

         WRITE(*,*) 'clsgen_tb: cluster size of site:',JATOM,':',NUMBER
         WRITE(*,*) 'clsgen_tb: Touching RMT of site:',JATOM,':',DISTMIN


 100  ENDDO

      IF (NUMBER.GT.NACLSD) THEN 
         WRITE(6,*) '(a) Increase the parameter NACLSD ',
     &        'to a value greater equal ',MAXNUMBER,'.'
         STOP 'clsgen_tb: Dimension error (a).'
      ENDIF

c===============================================================

c Find different types of ref. potential according to the rmtref.
      IREFPOT(1) = 1
      NREFPOT = 1
      DO ISITE = 2,NAEZ + NEMB
         LFOUND = .FALSE.
         DO JSITE = 1,ISITE - 1
            IF (ABS(RMTREF(ISITE)-RMTREF(JSITE)).LE.TOL) THEN
               IREFPOT(ISITE) = IREFPOT(JSITE)
               LFOUND = .TRUE.
            ENDIF
         ENDDO
         IF (.NOT.LFOUND) THEN
            NREFPOT = NREFPOT + 1
            IREFPOT(ISITE) = NREFPOT
         ENDIF
      ENDDO

      ICLUSTER = 1
      RCUTXY2 = (RCUTXY+EPSSHL)*(RCUTXY+EPSSHL)
      RCUT2   = (RCUT+EPSSHL)*(RCUT+EPSSHL)

      


c===============================================================
c
      DO 200 JATOM = 1,NAEZ + NEMB      ! loop in all atoms or layers
      
         CLS(JATOM) = 0   

         NUMBER = 0             ! counter for atoms in cluster
         DO NA = 1,NAEZ         ! loop in all atoms
            DO IR = 0,NR        ! loop in all bravais vectors    
               TMP(1:3) = RR(1:3,IR)+RBASIS(1:3,NA)-RBASIS(1:3,JATOM)
               RXY2 =  TMP(1)**2+TMP(2)**2
               R2   =  TMP(3)**2 + TMP(1)**2+TMP(2)**2
            
               IF ( (RXY2.LE.RCUTXY2).AND.(R2.LE.RCUT2) )  THEN
                  NUMBER = NUMBER + 1
                  ATOM(NUMBER,JATOM) = NA ! store the atom in elem cell
                  EZOA(NUMBER,JATOM) = IR ! store the bravais vector
                  RCLS1(1:3,NUMBER) = TMP(1:3)
               END IF
            END DO              ! N loop in bravais
         END DO                 ! NA loop in NAEZ
  
c     
c     In the case of 2 dimensional case loop in the atoms outside.
c     
         IF (L2DIM) THEN
c     Somehow messy
c     ATOM gives the kind of atom 
c   
            DO IR = 0,NR
               DO ILAY = NLEFT,1,-1  ! loop in some layers on left side
                  DO I1 = NLBASIS,1,-1 ! loop in representative atoms on left side
                     TMP(1:3) = RR(1:3,IR) + TLEFT(1:3,I1) +
     &                    (ILAY-1)*ZPERLEFT(1:3) - RBASIS(1:3,JATOM)
                     RXY2 =  TMP(1)**2 + TMP(2)**2
                     R2  =   TMP(3)**2 + RXY2
                  
                     IF ((RXY2.LE.RCUTXY2).AND.(R2.LE.RCUT2)) THEN
                        NUMBER = NUMBER + 1
                        ATOM(NUMBER,JATOM) = -NAEZ - I1 ! negative values are used in dlke1.f
                        EZOA(NUMBER,JATOM) = IR ! ILAY,I1 are negative
                        RCLS1(1:3,NUMBER) = TMP(1:3)
                     END IF
                  END DO 
               END DO
c     
c     
               DO ILAY = 1,NRIGHT
                  DO I1 = 1,NRBASIS
                     TMP(1:3) = RR(1:3,IR)+ TRIGHT(1:3,I1) + 
     &                    (ILAY-1)*ZPERIGHT(1:3) - RBASIS(1:3,JATOM)  
                     RXY2 = TMP(1)**2 + TMP(2)**2
                     R2  =  TMP(3)**2 + RXY2
                     IF ((RXY2.LE.RCUTXY2).AND.(R2.LE.RCUT2)) THEN
                        NUMBER = NUMBER + 1
                        ATOM(NUMBER,JATOM) = -NAEZ - NLBASIS - I1 
                        EZOA(NUMBER,JATOM) = IR
                        RCLS1(1:3,NUMBER) = TMP(1:3)
                     END IF
                  END DO 
               END DO
c     
         
            END DO              ! loop in all bravais lattices
         END IF                 ! L2DIM Interface calculation
c     
c     Now the atom JATOM has its cluster.  
c     Sort the atoms of the cluster in increasing order. 
c     First by distance, then by z, then by y, then by x.  
c     
         IF (NUMBER.GT.NACLSD) THEN ! should not hit here, this was checked earlier
            WRITE(6,*) '(b) Increase the parameter NACLSD ',
     &                 'to a value greater equal ',NUMBER,'.'
            STOP 'clsgen_tb: Dimension error (b).'
         END IF

         DO IA=1,NUMBER
           RSORT(IA) = DSQRT(RCLS1(1,IA)**2+
     &                       RCLS1(2,IA)**2+
     &                       RCLS1(3,IA)**2)

           RSORT(IA) = 1.D9*RSORT(IA)+
     &                 1.D6*RCLS1(3,IA)+
     &                 1.D3*RCLS1(2,IA)+
     &                 1.D0*RCLS1(1,IA)
         END DO      
c     
         CALL DSORT(RSORT,ISORT,NUMBER,POS)
c     Rearange exchange ia with ib
c     MAP temporarily to another array         
         DO IA=1,NUMBER       
            RG(1:3,IA) = RCLS1(1:3,IA)
            IATOM(IA) = ATOM(IA,JATOM)
            IEZOA(IA) = EZOA(IA,JATOM)
         END DO    
c Now use correct order
         DO IA =1,NUMBER
            IB = ISORT(IA)
            RCLS1(1:3,IA) = RG(1:3,IB)
            ATOM(IA,JATOM) = IATOM(IB)
            EZOA(IA,JATOM) = IEZOA(IB) 
         END DO
c     
c     Now the clusters have a unique sorting and can be compared with
c     each other Check if ICLUSTER was found previously       
c     
         DO ICU = 1,ICLUSTER-1
            N1 = NACLS(ICU)
            IAT1 = IREP(ICU)
            IF ( CLUSTCOMP_TB( RCLS,IREFPOT,ATOM,IAT1,
     &                         ICU,N1,RCLS1,NUMBER,JATOM) )
     &           CLS(JATOM) = ICU
            ! return true if found before
         END DO

         IF (CLS(JATOM).EQ.0) THEN ! no equivalent found, add new cluster
            NCLS = ICLUSTER
            IF (NCLS.GT.NCLSD) THEN
               WRITE(6,*) '(c) Increase the parameter NCLSD ',
     &              '  to a value greater equal ',ICLUSTER,' .'
               STOP 'clsgen_tb: Dimension error (c).' 
            END IF
            CLS(JATOM) = ICLUSTER
            NACLS(ICLUSTER) = NUMBER
            IREP(ICLUSTER) = JATOM ! cluster-class is represented by the cluster around jatom
            WRITE(6,FMT='(A27,I5,A7,I5)')
     &           'clsgen_voronoi: Cluster No.',ICLUSTER,' sites:',NUMBER
            DO IN = 1,NUMBER
               RCLS(1:3,IN,ICLUSTER) = RCLS1(1:3,IN)
               WRITE(6,800) JATOM,ATOM(IN,JATOM),EZOA(IN,JATOM),
     &                  (RCLS1(IX,IN),IX=1,3),
     &              SQRT(RCLS1(1,IN)**2+RCLS1(2,IN)**2+RCLS1(3,IN)**2)
 800           FORMAT(3I5,4F8.4)
            END DO   
            ICLUSTER = ICLUSTER + 1
         END IF 
c ******************************************************
 200  ENDDO                     ! JATOM = 1,NAEZ + NEMB

c Now all clusters of all atoms are found.
c

      WRITE(6,*) 'Clusters from clsgen_tb:'
      DO JATOM = 1,NAEZ
         WRITE(6,8000) JATOM,IREFPOT(JATOM),RMTREF(JATOM),
     &        CLS(JATOM),NACLS(CLS(JATOM))
      ENDDO
      IF (L2DIM) THEN
         WRITE(6,*) 'Clusters from clsgen_tb in outer region, left:'
         DO IA = 1,NLBASIS
            JATOM = NAEZ + IA
            WRITE(6,8000) JATOM,IREFPOT(JATOM),RMTREF(JATOM),
     &           CLS(JATOM),NACLS(CLS(JATOM))
         ENDDO
         
         WRITE(6,*) 'Clusters from clsgen_tb in outer region, right:'
         DO IA = 1,NRBASIS
            JATOM = NAEZ + NLBASIS + IA
            WRITE(6,8000) JATOM,IREFPOT(JATOM),RMTREF(JATOM),
     &           CLS(JATOM),NACLS(CLS(JATOM))
         ENDDO
      ENDIF
         




c
c ------------------------------------------------------------------------
       WRITE(6,*) ' Sub clsgen_tb  exiting <<<<<<<<<<<<<'
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
 8000   FORMAT('CLSGEN_TB: Atom',I5,' Refpot',I5,' Rmtref',F10.7,
     &         ' TB-cluster',I5,' Sites',I5)
 9005   FORMAT(I4)
 9010   FORMAT(('# Z     ',20F4.0))
 9020   FORMAT(('# KAOEZ ',20I4))
 9030   FORMAT(F12.7,6x,'ALAT')
 9040   FORMAT('> cluster ',I4,' at atom ',I4,
     +       ' of type ',I4,'.')
 9060   FORMAT(I4,1X,200I1)
 9080   FORMAT(I6,3F15.6)
 9090   FORMAT('The Number of layers is each Principal Layer = ',
     &                I5)


      RETURN
      END


