c ************************************************************************
      SUBROUTINE CLSGEN99(NAEZ,RR,NR,RBASIS,
     &                   CLS,
     &                   NACLS,REFPOT,ATOM,EZOA, 
     &                   RCLS,RCUT,RCUTXY,
     &                   NUMN0,INDN0)
      IMPLICIT NONE
c ************************************************************************
c This subroutine is used to create the clusters around each atom 
c where repulsive potentials will be positioned.
c
c STRATEGY : 
c Calculate the cluster of each atom by the lattice
c parameters avaliable. Sort the atoms in a unique way :big r, big z, big y
c compare the positions with the previous clusters to see if there is 
c a difference. If not keep only previous clusters and make indexing if
c a new cluster is found then check dimensions and continue for the new
c atom.  
c
      INCLUDE 'inc.p'
      INCLUDE 'inc.cls'
c
c     .. arguments
c
      DOUBLE PRECISION RCUT,RCUTXY
      DOUBLE PRECISION
     +     RBASIS(3,*),             ! pos. of basis atoms in EZ
     +     RCLS(3,NACLSD,*),        ! real space position of atom in cluster
     +     RR(3,0:NRD)              ! set of lattice vectors
c
      INTEGER
     +     NAEZ,                    ! number of atoms in EZ
     +     NR                       ! number of lattice vectors RR
c
      INTEGER NUMN0(NAEZD),INDN0(NAEZD,NACLSD),
     +     CLS(*),                  ! sort of cluster around atom
     +     REFPOT(*),              
     +     NACLS(*),                ! number of atoms in cluster
     +     ATOM(NACLSD,*),          ! index to atom in elem/cell at site in cluster
     +     EZOA(NACLSD,*)           ! index to bravais lattice  at site in cluster
c
c     .. locals
c
      INTEGER 
     +     I,N1,
     +     NA,NUMBER,N,
     +     POS,IA,IN,IB,II,JATOM,ICU,IC,IAT,ICLUSTER
      INTEGER IATOM(NACLSD),IEZOA(NACLSD),
     +     ISORT(NACLSD)
      INTEGER IATCLS(NCLSD)
      LOGICAL*1 ICOUPLMAT(NAEZD)
c
      DOUBLE PRECISION  
     +     R2,EPSSHL,
     +     RCLS1(3,NACLSD),
     +     RG(3,NACLSD),TMP(3),RSORT(NACLSD)
      DOUBLE PRECISION RCUT2,RCUTXY2,RXY2 
c
      LOGICAL CLUSTCOMP
c
c
      EXTERNAL XSORT,CLUSTCOMP
      INTRINSIC MIN,SQRT
c
      DATA     EPSSHL   / 1.0D-4 /
c
c ------------------------------------------------------------------------
c This is generating the clusters which have a distance smaller
c than RCUT and RCUTXY in plane .
c The cluster atoms are ordered with radious and then z>y>x 
c The ordering allows an easy comparison of clusters.

C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (6,'(79(1H=))')
      WRITE (6,'(16X,A)') 
     &     'CLSGEN99: generation of TB-clusters coordinates'
      WRITE (6,'(79(1H=))')
      WRITE (6,*)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
      WRITE(6,*) 'RCUT = ',rcut,' RCUTXY = ',rcutxy
      IF (ABS(RCUTXY - RCUT).LT.1.D-4) THEN
          WRITE(6,*) 'Spherical Clusters are created'
      END IF 
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ICLUSTER = 1
      DO N = 1,NCLSD
         IATCLS(N) = 0
      END DO
      RCUTXY2 = (RCUTXY+EPSSHL)*(RCUTXY+EPSSHL)
      RCUT2   = (RCUT+EPSSHL)*(RCUT+EPSSHL)
            
      DO JATOM = 1,NAEZ       ! loop in all atoms
         CLS(JATOM) = 0   
         NUMBER = 0           ! counter for atoms in cluster
         DO NA = 1,NAEZ  ! loop in all atoms
            DO N=0,NR    ! loop in all bravais vectors    
               DO I=1,3
                  TMP(I) = RR(I,N)+RBASIS(I,NA)-RBASIS(I,JATOM)
               END DO
               RXY2 =  TMP(1)**2+TMP(2)**2
               R2   =  TMP(3)**2 + TMP(1)**2+TMP(2)**2

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
            END DO              ! N loop in bravais
            
         END DO                 ! NA loop in NAEZ
  
c     Now the atom JATOM Has it's cluster first 
c     sort the atoms of the cluster in increasing order. First by distance
c     Then by z then by y    
c     
         DO IA=1,NUMBER
            RSORT(IA) = SQRT(RCLS1(1,IA)**2+
     &                       RCLS1(2,IA)**2+
     &                       RCLS1(3,IA)**2)
            RSORT(IA) = 100000000.D0*RSORT(IA)+
     &                     100000.D0*RCLS1(3,IA)+
     &                        100.D0*RCLS1(2,IA)+
     &                         0.1D0*RCLS1(1,IA) 
         END DO
c     
         CALL XSORT(RSORT,ISORT,NUMBER,POS)
c     Rearange exchange ia with ib
c MAP temporarily to another array         
         DO IA=1,NUMBER       
            DO I=1,3
               RG(I,IA)    = RCLS1(I,IA)
            END DO
            IATOM(IA) = ATOM(IA,JATOM)
            IEZOA(IA) = EZOA(IA,JATOM)
         END DO    
c Now use correct order
         DO IA =1,NUMBER
            IB = ISORT(IA)
             DO I=1,3
               RCLS1(I,IA) = RG(I,IB)
            END DO
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
           IF( CLUSTCOMP(RCLS,REFPOT,ATOM,IATCLS(ICU),N1,RCLS1,
     &          NUMBER,JATOM)) CLS(JATOM) = ICU
         END DO
         IF (CLS(JATOM).EQ.0) THEN
            IF (ICLUSTER.GT.NCLSD) THEN
               WRITE(6,*) 'Please, increase the parameter NCLSD in',
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
            WRITE(6,'(3I5,4F8.4)') JATOM,ATOM(IN,JATOM),EZOA(IN,JATOM),
     &                  (RCLS1(I,IN),I=1,3),
     &      SQRT(RCLS1(1,IN)**2+RCLS1(2,IN)**2+RCLS1(3,IN)**2)
            END DO   
            ICLUSTER = ICLUSTER + 1
         END IF 
c ******************************************************
       WRITE(6,*) 'Atom ',JATOM,' has cluster ', CLS(JATOM),
     &            'with ',NUMBER,' sites'
      END DO
c
c Now all clusters of all atoms are found
c
       DO JATOM = 1,NAEZ
      
          IC = CLS(JATOM)
          NUMBER = NACLS(IC)
  
          DO IAT = 1,NAEZ
             ICOUPLMAT(IAT) = .FALSE.
             DO I=1,NUMBER 
                IF (ATOM(I,JATOM).EQ.IAT) THEN
                   ICOUPLMAT(IAT) = .TRUE.
                END IF
             END DO
          END DO

          NUMN0(JATOM) = 0
          DO IAT = 1,NAEZ
          IF(ICOUPLMAT(IAT)) THEN
          NUMN0(JATOM) = NUMN0(JATOM) + 1
          INDN0(JATOM,NUMN0(JATOM)) = IAT
c     WRITE(6,*) 'JATOM,NUMN0,INDN0',JATOM,NUMN0(JATOM),IAT
          ENDIF
          ENDDO

        END DO

          WRITE(6,*) 
          WRITE(6,*) ' Sub clsgen99  exiting <<<<<<<<<<<<<'

      RETURN
      END
