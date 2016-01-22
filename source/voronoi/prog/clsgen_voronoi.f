c ************************************************************************
      SUBROUTINE CLSGEN_VORONOI(NATYP,NAEZ,NEMB,RR,NR,RBASIS,
     &                   KAOEZ,ZAT,CLS,NCLS,
     &                   NACLS,ATOM,EZOA, 
     &                   NLBASIS,NRBASIS,NLEFT,NRIGHT,ZPERLEFT,ZPERIGHT,
     &                   TLEFT,TRIGHT,
     &                   RCLS,RMTHLF,RCUT,RCUTXY,L2DIM,
     &                   ALAT)
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
     &     RMTHLF(NAEZD+NEMBD)      ! touching rmt by bisection of dist. to nearest neighbour
!     +     RWS(*),
!     +     RMT(*),
!     &     BBOX(3)                  ! bounding box for povray plots
c
      INTEGER
     +     KMT,                     ! scaling of RMT with MTFAC (not used)
c                                   ! 0: RMT from crystal structure 
c                                   ! 1: RMT - " -  scaled by RMTFAC
c                                   ! 2: RMT = RMTFAC
c                                   ! 3: RMT from ref. pot. card
     +     NATYP,                   ! number of sorts of atoms
     +     NAEZ,                    ! number of atoms in EZ
     +     NEMB,                    ! number of embedding postions
     +     NCLS,                    ! number of diff. clusters
     +     NINEQ,                   ! number of nonequivalent atomic 
                                    ! positions in EZ
     +     NR                       ! number of lattice vectors RR
c
      INTEGER
     +     CLS(*),                  ! sort of cluster around atom
     +     KAOEZ(*),                ! sort of atom at position in EZ
     +     NACLS(*),                ! number of atoms in cluster
     +     ATOM (NACLSD,*),         ! index to atom in elem/cell at site in cluster
     +     EZOA (NACLSD,*)         ! index to bravais lattice  at site in cluster
c
c     .. locals
c
      INTEGER 
     +     AJ,C,ILAY,J,N1,INUM,ISUM,IR,INEI,
     +     NA,NUMBER,NC,NPRIN,ITEST1,ITEST,IX,
     +     POS,IA,IN,IB,II,JATOM,ICU,IC,IAT,I0,I1,ICLUSTER
      INTEGER IATOM(NACLSD),IEZOA(NACLSD),
     +     ISORT(NACLSD),ICOUPLMAT(NAEZD,NAEZD)
c
      REAL*8        
     +     RAD,R1,R2,RABS,RD,T,EPSSHL,DISTMIN,TOL2,
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
      LOGICAL  L2DIM,CLUSTCOMP_VORONOI
c
c
      LOGICAL TEST,OPT,pov
      EXTERNAL DSORT,CLUSTCOMP_VORONOI
      INTRINSIC MIN,SQRT
c
      DATA     EPSSHL   / 1.0D-4 /
      DATA     TOL2   / 1.D-7 /
c ------------------------------------------------------------------------
      WRITE(6,*) '>>> CLSGEN99: generation of cluster coordinates'
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
      WRITE(8,9005) NAEZ
      WRITE(8,9030) ALAT
      WRITE(8,9010) (ZAT(KAOEZ(IAT)),IAT=1,NAEZ)
      WRITE(8,9020) (KAOEZ(IAT),IAT=1,NAEZ)
      
      
      ICLUSTER = 1
      RCUTXY2 = (RCUTXY+EPSSHL)*(RCUTXY+EPSSHL)
      RCUT2   = (RCUT+EPSSHL)*(RCUT+EPSSHL)

      
c===============================================================
c Check the if the dimension NACLSD is enough before starting to
c assign clusters. Write out necessary dimension.
      DO 100 JATOM = 1,NAEZ+NEMB     ! loop in all atoms or layers

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
      
         DISTMIN = DSQRT(DISTMIN)/2.D0 ! Touching RMT by bisection of dist. to nearest neighbour
         RMTHLF(JATOM) = DISTMIN

         WRITE(*,*) 'clsgen_voronoi: Max. cluster size found: ',NUMBER
         IF (NUMBER.GT.NACLSD) THEN 
            WRITE(6,*) '(a) Increase the parameter NACLSD ',
     &                 'to a value greater equal ',NUMBER,'.'
            STOP 'clsgen: Dimension error (a).'
         ENDIF

 100  ENDDO ! JATOM = 1,NAEZ+NEMB
c===============================================================
c
      DO 200 JATOM = 1,NAEZ+NEMB       ! loop in all atoms or layers
      
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
                        ATOM(NUMBER,JATOM) = -NAEZ - I1 ! negative values are used outside slab
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
                        ATOM(NUMBER,JATOM) = -NAEZ - NLBASIS - I1 ! negative values are used outside slab
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
            STOP 'clsgen13: Dimension error (b).'
         END IF

         DO IA=1,NUMBER
           RSORT(IA) = SQRT(RCLS1(1,IA)**2+
     &                      RCLS1(2,IA)**2+
     &                      RCLS1(3,IA)**2)

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
            IF (CLUSTCOMP_VORONOI(RCLS,ICU,N1,RCLS1,NUMBER))  ! return true if found before 
     &           CLS(JATOM) = ICU 
         END DO

         IF (CLS(JATOM).EQ.0) THEN
            NCLS = ICLUSTER
            IF (ICLUSTER.GT.NCLSD) THEN
               WRITE(6,*) '(c) Increase the parameter NCLSD ',
     &              '  to a value greater equal ',ICLUSTER,' .'
               STOP 'clsgen: Dimension error (c).' 
            END IF
            CLS(JATOM) = ICLUSTER
            NACLS(ICLUSTER) = NUMBER
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
 200  ENDDO                     ! JATOM = 1,NAEZ+NEMB

c
c Now all clusters of all atoms are found print out
c and test the results...
c

      WRITE(6,*) 'Clusters from clsgen_voronoi:'
      DO JATOM = 1,NAEZ
         WRITE(6,8000) JATOM,RMTHLF(JATOM),CLS(JATOM),NACLS(CLS(JATOM))
      ENDDO
      IF (L2DIM) THEN
         WRITE(6,*) 'Clusters from clsgen_tb in outer region, left:'
         DO IA = 1,NLBASIS
            JATOM = NAEZ + IA
            WRITE(6,8000) 
     &           JATOM,RMTHLF(JATOM),CLS(JATOM),NACLS(CLS(JATOM))
         ENDDO
         
         WRITE(6,*) 'Clusters from clsgen_tb in outer region, right:'
         DO IA = 1,NRBASIS
            JATOM = NAEZ + NLBASIS + IA
            WRITE(6,8000) 
     &           JATOM,RMTHLF(JATOM),CLS(JATOM),NACLS(CLS(JATOM))
         ENDDO
      ENDIF

      ! The arrays CLS,NACLS,RCLS for atoms between [NAEZ+1,NAEZ+NEMB]
      ! are disregarded in the rest of the program. In case of impurity
      ! calculation, these array positions can be occupied by the impurity
      ! cluster data. However, RMTHLF is used also for the embedded atoms
      ! in order to construct "optimized" weights.

c
c ------------------------------------------------------------------------
      DO 22 JATOM = 1,NAEZ
          IC = CLS(JATOM)
          NUMBER = NACLS(IC)
            WRITE(8,FMT=1030) NUMBER
            WRITE(8,FMT=1030) JATOM,IC
            DO 105 INEI = 1,NUMBER
              RAD = SQRT(RCLS(1,INEI,IC)**2 + 
     &                   RCLS(2,INEI,IC)**2+RCLS(3,INEI,IC)**2)
              WRITE(8,1040) (RCLS(II,INEI,IC)*ALAT,II=1,3)      
 105        END DO
 22    CONTINUE

       WRITE(6,*) ' Sub clsgen_voronoi  exiting <<<<<<<<<<<<<'
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
 8000   FORMAT('CLSGEN_VORONOI: Atom',I5,' Rmthlf',F10.7,
     &   ' cluster',I5,' Sites',I5)

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


      RETURN
      END


