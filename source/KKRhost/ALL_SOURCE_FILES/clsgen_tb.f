! ************************************************************************
      SUBROUTINE CLSGEN_TB(NAEZ,NEMB,NVIRT,RR,NR,RBASIS,KAOEZ,ZAT,
     &                   CLS,NCLS,NACLS,ATOM,EZOA, 
     &                   NLBASIS,NRBASIS,NLEFT,NRIGHT,ZPERLEFT,ZPERIGHT,
     &                   TLEFT,TRIGHT,RMTREF,RMTREFAT,VREF,
     &                   IREFPOT,NREFPOT,RCLS,RCUT,RCUTXY,L2DIM,ALAT,
     &                NAEZD,NATYPD,NEMBD,NPRINCD,NRD,NACLSD,NCLSD,NREFD)
      use mod_version_info
      implicit none
! ************************************************************************
! This subroutine is used to create the clusters around each atom 
! (Based on clsgen99.f). Also the reference potential height and radius is set
! (vref and rmtref).
!
! STRATEGY : 
! Calculate the cluster of each atom by the lattice
! parameters avaliable. Sort the atoms in a unique way :big r, big z, big y
! compare the positions with the previous clusters to see if there is 
! a difference. If not keep only previous clusters and make indexing if
! a new cluster is found then check dimensions and continue for the new
! atom.  
!
!
!
!
!     .. arguments
!
      INTEGER NAEZD,NATYPD,NEMBD,NPRINCD,NRD,NACLSD,NCLSD,NREFD
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
!
      INTEGER
     +     NAEZ,                    ! number of atoms in EZ
     +     NEMB,                    ! number of embedding postions
     +     NCLS,                    ! number of diff. clusters
     +     NR,                      ! number of lattice vectors RR
     &     NLR,                     ! =NEMB in decimation, =0 in slab or bulk
     &     NVIRT,                   ! Number of virtual atoms
     &     NPRINC                   ! Calculated number of layers in a principal layer
!
      INTEGER
     +     CLS(*),                   ! type of cluster around atom
     +     KAOEZ(NATYPD,NAEZD+NEMBD),! type of atom at position in EZ
     +     NACLS(*),                 ! number of atoms in cluster
     +     ATOM(NACLSD,*),           ! index to atom in elem/cell at site in cluster
     +     EZOA (NACLSD,*)           ! index to bravais lattice  at site in cluster

!     .. locals

      INTEGER 
     +     ILAY,N1,IR,ISITE,JSITE,IAT1,
     +     NA,NUMBER,MAXNUMBER,!IX,
     +     POS,IA,IN,IB,II,JATOM,ICU,IC,IAT,I1,ICLUSTER,NCLSALL
      INTEGER IATOM(NACLSD),IEZOA(NACLSD),
     +     ISORT(NACLSD),ICOUPLMAT(NAEZD,NAEZD),
     &     IREP(NCLSD) ! representative atom of cluster (inverse of CLS)
      INTEGER IREFPOT(NAEZD+NEMBD),NREFPOT
      REAL*8 RMTREFAT(NAEZD+NEMBD),RMTREF1(NAEZD+NEMBD)
      REAL*8 VREFAT(NAEZD+NEMBD),VREF1(NAEZD+NEMBD),VREF(NREFD)


      REAL*8        
     +     R2,EPSSHL,TOL,TOL2,DISTMIN,
     +     RCLS1(3,NACLSD),
     +     RG(3,NACLSD),TMP(3),RSORT(NACLSD)
      INTEGER NLBASIS,NRBASIS,                    
     +        NLEFT,NRIGHT           
      REAL*8                                   
     +        ZPERLEFT(3),ZPERIGHT(3),            
     +        TLEFT(3,*),TRIGHT(3,*)
      REAL*8       RCUT2,RCUTXY2,RXY2,DIST

      LOGICAL LFOUND
      LOGICAL  L2DIM,CLUSTCOMP_TB


      EXTERNAL DSORT,CLUSTCOMP_TB
      INTRINSIC MIN,SQRT

      DATA     EPSSHL   / 1.D-4 /
      DATA     TOL   / 1.D-7 /
      DATA     TOL2   / 1.D-7 /

! ------------------------------------------------------------------------
      WRITE(1337,*) '>>> CLSGEN_TB: generation of cluster coordinates'
! This is generating the clusters which have a distance smaller
! than RCUT and RCUTXY in plane .
! The cluster atoms are ordered with radius and then z>y>x 
! The ordering allows an easy comparison of clusters
! The principal layer for each layer (atom in unit cell) is
! calculated also for each cluster and the maximum number
! is returned. Some dimension tests are also done      
      WRITE(1337,*) 'RCUT = ',RCUT,' RCUTXY = ',RCUTXY
      IF (ABS(RCUTXY - RCUT).LT.1.D-4) THEN
          WRITE(1337,*) 'Spherical Clusters are created'
!          LSPHER=.TRUE.
      END IF 
      OPEN(8,FILE='clusters',status='unknown')
      call version_print_header(8)
      WRITE(8,9005) NAEZ
      WRITE(8,9030) ALAT
      WRITE(8,9010) (ZAT(KAOEZ(1,IAT)),IAT=1,NAEZ-NVIRT)
      WRITE(8,9020) (KAOEZ(1,IAT),IAT=1,NAEZ)
      
      RCUTXY2 = RCUTXY**2
      RCUT2 = RCUT**2
      NLR = 0
      IF (L2DIM) NLR = NEMB
      VREFAT(:) = 8.D0      ! Set to 8 Rydbergs
      VREF1(:) = 8.D0      ! Set to 8 Rydbergs
      VREF(:) = 8.D0      ! Set to 8 Rydbergs

!===============================================================
! Check the if the dimension NACLSD is enough before starting to
! assign clusters. Write out necessary dimension.
! Also find touching MT radius.
      DO 100 JATOM = 1,NAEZ + NLR ! loop in all sites incl. left/right host in decimation case

         MAXNUMBER = 0
         NUMBER = 0             ! counter for atoms in cluster
         DISTMIN = 1.D100       ! Large initial value for RMT**2
         DO NA = 1,NAEZ         ! loop in all sites in unit cell
            IF (KAOEZ(1,NA).NE.-1) THEN ! Exclude virtual atoms from clusters (except clust. center)
               DO IR = 0,NR     ! loop in all bravais vectors    
                  TMP(1:3) = RR(1:3,IR)+RBASIS(1:3,NA)-RBASIS(1:3,JATOM)
                  RXY2 =  TMP(1)**2+TMP(2)**2
                  R2   =  TMP(3)**2 + TMP(1)**2+TMP(2)**2
                  IF ( (RXY2.LE.RCUTXY2).AND.(R2.LE.RCUT2) )  THEN
                     NUMBER = NUMBER + 1
                  END IF
                  IF (R2.GT.TOL2) DISTMIN = MIN(DISTMIN,R2)
               END DO           ! IR loop in bravais
            ENDIF               ! (KAOEZ(1,NA).NE.-1) 
         END DO                 ! NA loop in NAEZ
  
!     
!        In the case of 2 dimensional case loop in the atoms outside.
!     
         IF (L2DIM) THEN
!        Somehow messy (lionel messi?)
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
!     
         
            END DO              ! loop in all bravais lattices
         END IF                 ! L2DIM Interface calculation

         DISTMIN = DSQRT(DISTMIN)/2.D0 ! Touching RMT
         ! Define MT-radius of TB-ref. potentials if not read-in from the input
         IF (RMTREFAT(JATOM).LT.0.D0)  ! I.e., if not read in 
     &      RMTREFAT(JATOM) = INT( DISTMIN*ALAT*99.D0 )/100.D0 ! round up the decimals
         IF (KAOEZ(1,JATOM).EQ.-1) THEN ! Virtual atom
            RMTREFAT(JATOM) = 1.D-20 
            VREFAT(JATOM) = 0.D0
         ENDIF

!     
         MAXNUMBER = MAX(MAXNUMBER,NUMBER) ! Find largest cluster size

         WRITE(1337,*) 'clsgen_tb: cluster size of site:',JATOM,':',
     &                                                          NUMBER
         WRITE(1337,*) 'clsgen_tb: Touching RMT of site:',JATOM,':',
     &                                                          DISTMIN


 100  ENDDO

      IF (NUMBER.GT.NACLSD) THEN 
         WRITE(6,*) '(a) Increase the parameter NACLSD ',
     &        'to a value greater equal ',MAXNUMBER,'.'
         STOP 'clsgen_tb: Dimension error (a).'
      ENDIF

!===============================================================

! Find different types of ref. potential according to the rmtref and vref
      IREFPOT(1) = 1
      RMTREF1(1) = RMTREFAT(1)
      VREF1(1) = VREFAT(1)
      NREFPOT = 1
      DO ISITE = 2,NAEZ + NEMB
         LFOUND = .FALSE.
         DO JSITE = 1,ISITE - 1
            IF ( ABS(RMTREFAT(ISITE) - RMTREFAT(JSITE)) + 
     &           ABS(VREFAT(ISITE)   - VREFAT(JSITE)      ).LE.TOL) THEN
               IREFPOT(ISITE) = IREFPOT(JSITE)
               LFOUND = .TRUE.
            ENDIF
         ENDDO
         IF (.NOT.LFOUND) THEN
            NREFPOT = NREFPOT + 1
            IREFPOT(ISITE) = NREFPOT
            ! RMTREFAT goes over all sites, RMTREF1 only over all inequivalent ref. potentials.
            RMTREF1(NREFPOT) = RMTREFAT(ISITE) 
            VREF1(NREFPOT) = VREFAT(ISITE)
         ENDIF
      ENDDO

      IF (NREFPOT.GT.NREFD) THEN
         WRITE(*,*) 'clsgen_tb: NREFPOT.GT.NREFD:',NREFPOT,NREFD
         STOP 'clsgen_tb: NREFPOT.GT.NREFD'
      ENDIF
      ! Now that the dimension is known, copy to array RMTREF
      DO I1 = 1,NREFPOT
         IF (RMTREF(I1).LT.0.D0) RMTREF(I1) = RMTREF1(I1)
      ENDDO
      VREF(1:NREFPOT) = VREF1(1:NREFPOT)



      ICLUSTER = 1
      RCUTXY2 = (RCUTXY+EPSSHL)*(RCUTXY+EPSSHL)
      RCUT2   = (RCUT+EPSSHL)*(RCUT+EPSSHL)



!===============================================================

      DO 200 JATOM = 1,NAEZ + NLR   ! loop in all atoms or layers
      
         CLS(JATOM) = 0   

         NUMBER = 0                      ! counter for sites in cluster
         DO NA = 1,NAEZ                  ! loop in all sites
            IF (KAOEZ(1,NA).NE.-1) THEN  ! proceed only if the neighbour is not virtual atom
               DO IR = 0,NR              ! loop in all bravais vectors    
                  TMP(1:3) = RR(1:3,IR)+RBASIS(1:3,NA)-RBASIS(1:3,JATOM)
                  RXY2 =  TMP(1)**2+TMP(2)**2
                  R2   =  TMP(3)**2 + TMP(1)**2+TMP(2)**2
            
                  IF ( (RXY2.LE.RCUTXY2).AND.(R2.LE.RCUT2) )  THEN
                     NUMBER = NUMBER + 1
                     ATOM(NUMBER,JATOM) = NA ! store the atom in elem cell
                     EZOA(NUMBER,JATOM) = IR ! store the bravais vector
                     RCLS1(1:3,NUMBER) = TMP(1:3)
                  END IF
               END DO           ! IR loop in bravais
            ENDIF
         END DO                 ! NA loop in NAEZ
  
!     
!        In the case of 2 dimensional case loop in the atoms outside.
!     
         IF (L2DIM) THEN
!        Somehow messy (eh? lionel?)
!        Index ATOM gives the kind of atom 
!   
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
!     
!     
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
!     
         
            END DO              ! loop in all bravais lattices
         END IF                 ! L2DIM Interface calculation
!     
!     Now the atom JATOM has its cluster.  
!     Sort the atoms of the cluster in increasing order. 
!     First by distance, then by z, then by y, then by x.  
!     
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
!     
         CALL DSORT(RSORT,ISORT,NUMBER,POS)
!     Rearange exchange ia with ib
!     MAP temporarily to another array         
         DO IA=1,NUMBER       
            RG(1:3,IA) = RCLS1(1:3,IA)
            IATOM(IA) = ATOM(IA,JATOM)
            IEZOA(IA) = EZOA(IA,JATOM)
         END DO    
! Now use correct order
         DO IA =1,NUMBER
            IB = ISORT(IA)
            RCLS1(1:3,IA) = RG(1:3,IB)
            ATOM(IA,JATOM) = IATOM(IB)
            EZOA(IA,JATOM) = IEZOA(IB) 
         END DO
!     
!     Now the clusters have a unique sorting and can be compared with
!     each other Check if ICLUSTER was found previously       
!     
         DO ICU = 1,ICLUSTER-1
            N1 = NACLS(ICU)
            IAT1 = IREP(ICU)
            IF ( CLUSTCOMP_TB( RCLS,IREFPOT,ATOM,IAT1,
     &                         ICU,N1,RCLS1,NUMBER,JATOM,NACLSD) )
     &           CLS(JATOM) = ICU
            ! return true if found before
         END DO

         IF (CLS(JATOM).EQ.0) THEN ! no equivalent found, add new cluster
            NCLSALL = ICLUSTER ! incl. embedded atoms of left/right host
            IF (JATOM.LE.NAEZ) NCLS = ICLUSTER  ! Excl. embedded atoms of left/right host
            IF (NCLSALL.GT.NCLSD) THEN
               WRITE(6,*) '(c) Increase the parameter NCLSD ',
     &              '  to a value greater equal ',ICLUSTER,' .'
               STOP 'clsgen_tb: Dimension error (c).' 
            END IF
            CLS(JATOM) = ICLUSTER
            NACLS(ICLUSTER) = NUMBER
            IREP(ICLUSTER) = JATOM ! cluster-class is represented by the cluster around jatom
            DO IN = 1,NUMBER
               RCLS(1:3,IN,ICLUSTER) = RCLS1(1:3,IN)
!               WRITE(6,800) JATOM,ATOM(IN,JATOM),EZOA(IN,JATOM),
!     &                  (RCLS1(IX,IN),IX=1,3),
!     &              SQRT(RCLS1(1,IN)**2+RCLS1(2,IN)**2+RCLS1(3,IN)**2)
! 800           FORMAT(3I5,4F8.4)
            END DO   
            ICLUSTER = ICLUSTER + 1
         END IF 
! ******************************************************
 200  ENDDO                     ! JATOM = 1,NAEZ + NEMB

! Now all clusters of all atoms are found.


      WRITE(1337,*) 'Clusters from clsgen_tb:'
      DO JATOM = 1,NAEZ
         WRITE(1337,8000) JATOM,IREFPOT(JATOM),RMTREFAT(JATOM),
     &        VREFAT(JATOM),CLS(JATOM),NACLS(CLS(JATOM))
      ENDDO
      IF (L2DIM) THEN
         WRITE(1337,*) 'Clusters from clsgen_tb in outer region, left:'
         DO IA = 1,NLBASIS
            JATOM = NAEZ + IA
            WRITE(1337,8000) JATOM,IREFPOT(JATOM),RMTREFAT(JATOM),
     &           VREFAT(JATOM),CLS(JATOM),NACLS(CLS(JATOM))
         ENDDO
         
         WRITE(1337,*) 'Clusters from clsgen_tb in outer region, right:'
         DO IA = 1,NRBASIS
            JATOM = NAEZ + NLBASIS + IA
            WRITE(1337,8000) JATOM,IREFPOT(JATOM),RMTREFAT(JATOM),
     &           VREFAT(JATOM),CLS(JATOM),NACLS(CLS(JATOM))
         ENDDO
      ENDIF
         
      ! Write out clusters in file
      OPEN(8,FILE='clusters',STATUS='UNKNOWN')
      WRITE(8,9005) NAEZ
      WRITE(8,9030) ALAT
      WRITE(8,9010) (ZAT(KAOEZ(1,I1)),I1=1,NAEZ)
      WRITE(8,9020) (KAOEZ(1,I1),I1=1,NAEZ)
      DO JATOM = 1,NAEZ
         IC = CLS(JATOM)
         NUMBER = NACLS(IC)
         WRITE(8,FMT=1030) NUMBER
         WRITE(8,FMT=1030) JATOM,IC
         DO I1=1,NUMBER
            DIST = SQRT(
     &           RCLS(1,I1,IC)**2+RCLS(2,I1,IC)**2+RCLS(3,I1,IC)**2)
            WRITE(8,1041) (RCLS(II,I1,IC),II=1,3),ATOM(I1,JATOM),
     &           ZAT(ABS(ATOM(I1,JATOM))),DIST     
         END DO

      ENDDO

      ! Write out the coupling matrix
      WRITE(1337,*) 'Coupling matrix:'
      DO JATOM = 1,NAEZ
          DO IAT = 1,NAEZ
             ICOUPLMAT(JATOM,IAT) = 0
             DO I1=1,NUMBER 
                IF (ATOM(I1,JATOM).EQ.IAT) THEN
                   ICOUPLMAT(JATOM,IAT) = 1
                END IF
             END DO
          END DO
          WRITE(1337,9060) JATOM,(ICOUPLMAT(JATOM,IAT),IAT=1,NAEZ)          
       END DO

       IF (L2DIM) THEN
       ! Calculate number of layers in principal layer
       NPRINC = 1
       DO JATOM = 1,NAEZ  ! loop over rows
          DO IAT = 1,JATOM - 1 ! loop over columns before the diagonal
           IF (ICOUPLMAT(JATOM,IAT).EQ.1) NPRINC = MAX(NPRINC,JATOM-IAT)
          ENDDO
          DO IAT = JATOM + 1,NAEZ ! loop over columns after the diagonal
           IF (ICOUPLMAT(JATOM,IAT).EQ.1) NPRINC = MAX(NPRINC,IAT-JATOM)
          ENDDO
       ENDDO
       WRITE(1337,*)
     &      'CLSGEN_TB: Number of layers in a principal layer: NPRINC=',
     &       NPRINC
       ENDIF



! ------------------------------------------------------------------------
       WRITE(1337,*) ' Sub clsgen_tb  exiting <<<<<<<<<<<<<'
! ------------------------------------------------------------------------

 1000   format(' cluster around atom     ',10I4/,
     +         ' number atoms in cluster ',10I4)
 1001   format(I4,2I5,3F8.2,I6,4f7.2)
 1002   format(' cocls : naez =',I3,' (x,y,z)= ',3F10.4)
 1010   format(12x,I6,3F10.4)
 1020   format('  Nr  naez kaoez     x       y       z',
     +         '  ezoa  RR(1)  RR(2)  RR(3)      R')
 1030   FORMAT(3I8)
 1040   FORMAT('Cu  ',3D24.16,2I8,F18.12)
 1041   FORMAT(3E27.19,I5,F5.1,E17.9)
 1050   FORMAT(3F12.7,'  scaling factor')
 1060   FORMAT(I4,3F12.7,'  center',/,(I4,3f12.7))
 1070   FORMAT(I4,3F12.7,'  center of gravity')
 1080   FORMAT('contains ',I4,'  atoms.')
 8000   FORMAT('CLSGEN_TB: Atom',I5,' Refpot',I5,' Rmtref',F10.7,
     &   ' Vref',F10.7,' TB-cluster',I5,' Sites',I5)
 9005   FORMAT(I4)
 9010   FORMAT(('# ZAT   ',20F6.2))
 9020   FORMAT(('# KAOEZ ',20I4))
 9030   FORMAT(F12.7,6x,'ALAT')
 9040   FORMAT('> cluster ',I4,' at atom ',I4,
     +       ' of type ',I4,'.')
 9060   FORMAT(I4,1X,500I1)
 9080   FORMAT(I6,3F15.6)
 9090   FORMAT('The Number of layers is each Principal Layer = ',
     &                I5)


      RETURN
      END


