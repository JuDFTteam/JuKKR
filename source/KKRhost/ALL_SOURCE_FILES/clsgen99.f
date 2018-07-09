c **********************************************************************
      SUBROUTINE CLSGEN99(NAEZ,RR,NR,RBASIS,KAOEZ,Z,CLS,NACLS,REFPOT,
     &                    ATOM,EZOA,NLBASIS,NRBASIS,NLEFT,NRIGHT,
     &                    ZPERLEFT,ZPERIGHT,TLEFT,TRIGHT,RCLS,
     &                    RCUT,RCUTXY,L2DIM,ALAT,
     &                    NAEZD,NATYPD,NEMBD,NPRINCD,NRD,NACLSD,NCLSD)
      use mod_version_info
      IMPLICIT NONE
c **********************************************************************
c This subroutine is used to create the clusters around each atom 
c where repulsive potentials will be positioned.
c
c STRATEGY : 
c Calculate the cluster of each atom by the lattice
c parameters avaliable. Sort the atoms in a unique way :
C                        big r, big z, big y
c compare the positions with the previous clusters to see if there is 
c a difference. If not keep only previous clusters and make indexing if
c a new cluster is found then check dimensions and continue for the new
c atom.  
c
c Small bug in assigning clusters removed 29/04/2003 v.popescu 
C IATCLUS(NCLSD) is pointing to the first atomic site associated 
C with a tb-cluster
c
cccc      include 'inc.p'
C     ..
C     .. Arguments ..
      INTEGER NACLSD,NCLSD
      INTEGER NAEZD,NATYPD,NEMBD,NPRINCD,NRD
      DOUBLE PRECISION ALAT         ! lattice constant A
      DOUBLE PRECISION RCUT,RCUTXY
      DOUBLE PRECISION
     +     RBASIS(3,*),             ! pos. of basis atoms in EZ
     +     RCLS(3,NACLSD,*),        ! real space position of atom in cluster
     +     RR(3,0:NRD),             ! set of lattice vectors
     +     Z(*)                     ! nucleus charge
c
      INTEGER
     +     NAEZ,                    ! number of atoms in EZ
     +     NR                       ! number of lattice vectors RR
c
      INTEGER
     +     CLS(*),                  ! sort of cluster around atom
     +     REFPOT(*),              
     +     KAOEZ(NATYPD,
     +           NAEZD+NEMBD),      ! sort of atom at position in EZ
     +     NACLS(*),                ! number of atoms in cluster
     +     ATOM(NACLSD,*),          ! index to atom in elem/cell at site in cluster
     +     EZOA(NACLSD,*)           ! index to bravais lattice  at site in cluster
c
c     .. locals
c
      INTEGER 
     +     I,N1,INUM,ISUM,
     +     NA,NUMBER,N,NPRIN,ITEST1,ITEST,
     +     POS,IA,IN,IB,II,JATOM,ICU,IC,IAT,I0,I1,ICLUSTER
      INTEGER IATOM(NACLSD),IEZOA(NACLSD),
     +     ISORT(NACLSD),ICOUPLMAT(NAEZD,NAEZD)
      INTEGER IATCLS(NCLSD)
c
      DOUBLE PRECISION  
     +     R,R2,EPSSHL,
     +     RCLS1(3,NACLSD),
     +     RG(3,NACLSD),TMP(3),RSORT(NACLSD)
      INTEGER NLBASIS,NRBASIS,
     +        NLEFT,NRIGHT           
      DOUBLE PRECISION                             
     +        ZPERLEFT(3),ZPERIGHT(3),            
     +        TLEFT(3,NEMBD+1),TRIGHT(3,NEMBD+1)
      DOUBLE PRECISION RCUT2,RCUTXY2,RXY2 
c
      LOGICAL  L2DIM,CLUSTCOMP
c
c
      LOGICAL TEST,LSPHER
      EXTERNAL DSORT,CLUSTCOMP
      INTRINSIC MIN,SQRT
c
      DATA     EPSSHL   / 1.0D-4 /
c
c ------------------------------------------------------------------------
c This is generating the clusters which have a distance smaller
c than RCUT and RCUTXY in plane .
c The cluster atoms are ordered with radious and then z>y>x 
c The ordering allows an easy comparison of clusters
c The principal layer for each layer (atom in unit cell) is
c calculated also for each cluster and the maximum number
c is returned. Some dimension tests are also done      

C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (1337,'(79(1H=))')
      WRITE (1337,'(16X,A)') 
     &     'CLSGEN99: generation of TB-clusters coordinates'
      WRITE (1337,'(79(1H=))')
      WRITE (1337,*)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
      LSPHER = .FALSE.
      WRITE(1337,*) 'RCUT = ',rcut,' RCUTXY = ',rcutxy
      IF (ABS(rcutxy - rcut).LT.1.D-4) THEN
          WRITE(1337,*) 'Spherical Clusters are created'
          LSPHER = .TRUE.
      END IF 
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF (TEST('clusters')) THEN
         OPEN(8,FILE='clusters',STATUS='UNKNOWN')
         call version_print_header(8)
         WRITE(8,9005) NAEZ
         WRITE(8,9030) ALAT
         WRITE(8,9010) (Z(KAOEZ(1,I)),I=1,NAEZ)
         WRITE(8,9020) (KAOEZ(1,I),I=1,NAEZ)
      END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ICLUSTER = 1
      DO N = 1,NCLSD
         IATCLS(N) = 0
         NACLS(N) = 0
      END DO
      CALL RINIT(3*NACLSD*NCLSD,RCLS)
C
      RCUTXY2 = (RCUTXY+EPSSHL)*(RCUTXY+EPSSHL)
      RCUT2   = (RCUT+EPSSHL)*(RCUT+EPSSHL)
            
      DO 2 JATOM = 1,NAEZ       ! loop in all atoms or layers
         CLS(JATOM) = 0   
         NUMBER = 0             ! counter for atoms in cluster
         DO NA = 1,NAEZ  ! loop in all atoms
            DO N=0,NR    ! loop in all bravais vectors    
               DO I=1,3
                  TMP(I) = RR(I,N)+RBASIS(I,NA)-RBASIS(I,JATOM)
               END DO
               RXY2 =  TMP(1)**2+TMP(2)**2
               R2   =  TMP(3)**2 
               IF (LSPHER) R2 = R2 + RXY2

               IF ( (RXY2.LE.RCUTXY2).AND.(R2.LE.RCUT2) )  THEN
                  NUMBER = NUMBER + 1
                  IF (NUMBER.GT.NACLSD) THEN 
                     WRITE (6,*) 
     &                 ' ERROR: Dimension NACLSD in inc.cls too small',
     &                 NUMBER, NACLSD
                     STOP '   < CLSGEN99 >'
                  END IF
C
                  ATOM(NUMBER,JATOM) = NA ! store the atom in elem cell
                  EZOA(NUMBER,JATOM) = N ! store the bravais vector
                  DO I=1,3
                     RCLS1(I,NUMBER) = TMP(I)
                  END DO
               END IF
 800           format(3I5,4F8.4)
            END DO              ! N loop in bravais
            
         END DO                 ! NA loop in NAEZ
  
c     
c     In the case of 2 dimensional case loop in the atoms
c     outside.
c     
         IF (L2DIM) THEN
c     Somehow meshy
c     ATOM gives the kind of atom 
c   
            DO N=0,NR
               DO I=NLEFT,1,-1  ! loop in some layers on left side
                  DO I1=NLBASIS,1,-1 ! loop in representative atoms on left side
                  DO I0=1,3
                  TMP(I0) = RR(I0,N) + TLEFT(I0,i1) + (I-1)*ZPERLEFT(I0)
     &                               - RBASIS(I0,JATOM)
                  END DO
                  RXY2 =  TMP(1)**2+TMP(2)**2
                  R2   =  TMP(3)**2 
                  IF (LSPHER) R2 = R2 + RXY2
                  
                  IF ((RXY2.LE.RCUTXY2).AND.(R2.LE.RCUT2)) THEN
                     
                     NUMBER = NUMBER + 1
                     IF (NUMBER.GT.NACLSD) THEN 
                        WRITE (6,*) 
     &                  ' ERROR: Dimension NACLSD in inc.cls too small',
     &                       NUMBER, NACLSD
                        STOP '   < CLSGEN99 >'
                     END IF
C
                     ATOM(NUMBER,JATOM) = -NAEZ-I1 ! negative values are used in dlke1.f
                     EZOA(NUMBER,JATOM) = N ! I,I1 are negative
                     DO I0=1,3
                        RCLS1(I0,NUMBER) = TMP(I0)
                     END DO
                  END IF
                  
               END DO 
            END DO
c     
c
            DO I=1,NRIGHT
               DO I1=1,NRBASIS
               DO I0=1,3
                 TMP(I0) = RR(I0,N)+ TRIGHT(I0,i1) + (I-1)*ZPERIGHT(I0)
     &                             - RBASIS(I0,JATOM)  
               END DO
               RXY2 =  TMP(1)**2+TMP(2)**2
               R2  =  TMP(3)**2 + TMP(1)**2+TMP(2)**2
               IF ((RXY2.LE.RCUTXY2).AND.(R2.LE.RCUT2)) THEN
                  NUMBER = NUMBER + 1
                  IF (NUMBER.GT.NACLSD) THEN 
                     WRITE (6,*) 
     &                 ' ERROR: Dimension NACLSD in inc.cls too small',
     &                    NUMBER, NACLSD
                     STOP '   < CLSGEN99 >'
                  END IF
                  ATOM(NUMBER,JATOM) = -NAEZ-NLBASIS-I1 
                  EZOA(NUMBER,JATOM) = N
                  DO I0=1,3
                     RCLS1(I0,NUMBER) = TMP(I0)
                  END DO
               END IF
            END DO 
         END DO
c     
         
      END DO                    ! loop in all bravais lattices
      END IF                 ! L2DIM Interface calculation
c     
c     Now the atom JATOM Has it's cluster first 
c     sort the atoms of the cluster in increasing order. First by distance
c     Then by z then by y    
c     
         do ia=1,number
            rsort(ia) = SQRT(RCLS1(1,ia)**2+
     &                       RCLS1(2,ia)**2+
     &                       RCLS1(3,ia)**2)
            rsort(ia) = 100000000.d0*rsort(ia)+
     &                    10000.d0*RCLS1(3,ia)+
     &                   10.d0*RCLS1(2,ia)+
     &                  0.1d0*RCLS1(1,ia) 
         end do      
c     
         CALL DSORT(RSORT,ISORT,NUMBER,POS)
c     Rearange exchange ia with ib
c MAP temporarily to another array         
         do IA=1,NUMBER       
            do I=1,3
               RG(I,IA)    = RCLS1(I,IA)
            end do
            IATOM(IA) = ATOM(IA,JATOM)
            IEZOA(IA) = EZOA(IA,JATOM)
         end do    
c Now use correct order
         do IA =1,NUMBER
            IB = ISORT(IA)
             do I=1,3
               RCLS1(I,IA) = RG(I,IB)
            end do
            ATOM(IA,JATOM) = IATOM(IB)
            EZOA(IA,JATOM) = IEZOA(IB) 
         END DO
c     
c     Now the clusters have a unique sorting and can be compared with
c     each other Check if ICLUSTER was found previously       
c     
         DO ICU = 1,ICLUSTER-1
            
           N1 = NACLS(ICU)
c return true if found before
           IF( CLUSTCOMP(RCLS,REFPOT,ATOM,IATCLS(ICU),ICU,N1,RCLS1,
     &          NUMBER,JATOM,NACLSD)) CLS(JATOM) = ICU
         END DO
         IF (CLS(JATOM).EQ.0) THEN
            IF (ICLUSTER.GT.NCLSD) THEN
               write(6,*) 'Please, increase the parameter NCLSD in',
     &              ' inc.cls to a value greater equal ',ICLUSTER,' .'
               STOP 'Dimension error.' 
            END IF
            CLS(JATOM) = ICLUSTER
            NACLS(ICLUSTER) = NUMBER
            IATCLS(ICLUSTER) = JATOM
            DO IN = 1,NUMBER
               DO II=1,3
               RCLS(II,IN,ICLUSTER) = RCLS1(II,IN)
               END DO
               write(1337,800) jatom,atom(in,jatom),ezoa(in,jatom),
     &                  (rcls1(i,in),i=1,3),
     &      sqrt(rcls1(1,in)**2+rcls1(2,in)**2+rcls1(3,in)**2)
            END DO   
            ICLUSTER = ICLUSTER + 1
         END IF 
c ******************************************************
       write(1337,*) 'Atom ',JATOM,' has cluster ', CLS(JATOM),
     &            'with ',NUMBER,' sites'
 2    CONTINUE                      ! JATOM = 1,NAEZ
c
c Now all clusters of all atoms are found print out
c and test the results...
c
c
      DO 22 JATOM = 1,NAEZ
c     
c ------------------------------------------------------------------------
          IC = CLS(JATOM)
          NUMBER = NACLS(IC)
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          IF (TEST('clusters')) THEN
             WRITE(8,FMT=1030) NUMBER
             WRITE(8,FMT=1030) JATOM,IC
             DO I=1,NUMBER
                R = SQRT(
     &               RCLS(1,I,IC)**2+RCLS(2,I,IC)**2+RCLS(3,I,IC)**2)
                WRITE(8,1041) (RCLS(II,I,IC),II=1,3),ATOM(I,JATOM),
     &               Z(ABS(ATOM(I,JATOM))),R      
             END DO
          END IF
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Now print out the coupling matrix
c
c 
          DO IAT = 1,NAEZ
             ICOUPLMAT(JATOM,IAT) = 0
             DO I=1,NUMBER 
                IF (ATOM(I,JATOM).EQ.IAT) THEN
                   ICOUPLMAT(JATOM,IAT) = 1
                END IF
             END DO
          END DO

          WRITE(1337,9060) JATOM,(ICOUPLMAT(JATOM,IAT),IAT=1,NAEZ)          

 22       END DO ! Do loop in JATOM (second test loop)
!
! ----------------------------------------------------------------------
!     output coupling matrix to file                              ! GODFRIN
      open(file='couplings.dat',unit=123456,status='replace')     ! GODFRIN
      write(123456,'("# Couplings between atoms via gref")')      ! GODFRIN
      write(123456,'(i8)') naez                                   ! GODFRIN
      do ia=1,naez                                                ! GODFRIN
        write(123456,'(1000i1)') icouplmat(ia,1:naez)             ! GODFRIN
      end do                                                      ! GODFRIN
      close(123456)                                               ! GODFRIN
! ----------------------------------------------------------------------
c
c Now testing the shape of the dyson equation
c
          ITEST1 = ICOUPLMAT(NAEZ,1)+ICOUPLMAT(1,NAEZ)
          IF (ITEST1.NE.0) THEN 
          WRITE (1337,*) ' This is not a banded matrix ' 
          END IF
          nprin = 0
          do i=1,naez-1
             itest = icouplmat(1,i)
             itest1 = icouplmat(1,i+1)
             if (itest.eq.1.and.itest1.eq.0) then
                nprin = i-1
             end if
          end do
          if (nprin.eq.0) nprin = 1
          WRITE(1337,*) '***********************************************'
          WRITE(1337,*) '********** TESTING THE COUPLING MATRIX ********'
          WRITE(1337,*) '***********************************************'
          WRITE(1337,9090) NPRIN
          IF (NPRIN.NE.NPRINCD) THEN 
              WRITE(1337,*) 'Please change NPRINCD in your inc.p file'
              WRITE(1337,*) 'from ',NPRINCD, 'to ',NPRIN
              WRITE(1337,*) ' ******** RESULTS COULD BE WRONG ******** '
           END IF 
c Now check if you can divide the matrix correctly
          IF (MOD(NAEZ,NPRIN).NE.0) THEN
              WRITE(1337,*) ' Your matrix cannot be divided in '
              WRITE(1337,*) ' Principal layers. Use a number of layers '
              write(1337,*) ' which is multiple of ',NPRIN
          END IF 
c                       NPR 
c NL + 2*NPR*NL - 4* sum  {n}
c                      n=1 
          isum = 0
          do i=1,nprin
             isum = isum + i
          end do
          inum = NAEZ + 2*NPRIN*NAEZ - 2*ISUM
c
c
          isum = 0
          do i=1,naez
             do i0=1,naez
               isum = isum + icouplmat(i,i0)         
             end do
          end do

          IF (ISUM.EQ.INUM) THEN 
              WRITE(1337,*) ' Your matrix is BAND DIAGONAL' 
          ELSE 
              WRITE(1337,*) ' Your matrix is *NOT* BAND DIAGONAL',ISUM,INUM
          END IF  
          WRITE(1337,*) 
          WRITE(1337,*) ' Sub clsgen99  exiting <<<<<<<<<<<<<'
c ------------------------------------------------------------------------

 1030   FORMAT(3I8)
c 1041   FORMAT(3F15.8,I6,F7.1,F12.6)
 1041   FORMAT(3E27.19,I5,F5.1,E17.9)
 9005   FORMAT(I4)
 9010   FORMAT(('# Z     ',20F4.0))
 9020   FORMAT(('# KAOEZ ',20I4))
 9030   FORMAT(F12.7,6x,'ALAT')
 9060   FORMAT(I4,1X,200I1)
 9090   FORMAT('The Number of layers in each Principal Layer = ',
     &                I5)
      RETURN
      END
