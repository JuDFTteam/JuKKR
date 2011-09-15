c ************************************************************************
      SUBROUTINE CLSGEN_TRC(
     >                      NAEZ,RR,NR,RBASIS,
     >                      RCUTTRC,ALAT,
     &                      NRD, NATRCD, NUTRCD)
C
      IMPLICIT NONE
c ************************************************************************
c This subroutine is used to create the clusters around each atom 
c where repulsive potentials will be positioned.
c
c STRATEGY : 
c Calculate the truncation-cluster of each atom by the lattice
c parameters avaliable. Sort the atoms in a unique way :big r, big z, big y
c
c
C      INCLUDE 'inc.p'
C      INCLUDE 'inc.cls'
c
c     .. input arguments
c

c     new after inc.p replace
      INTEGER NRD
      INTEGER NATRCD
      INTEGER NUTRCD

      DOUBLE PRECISION ALAT,               ! lattice constant A
     +                 RCUTTRC,            ! cutoff radius defining truncation zone
     +                 RBASIS(3,NAEZ),    ! pos. of basis atoms in EZ
     +                 RR(3,0:NRD)         ! set of lattice vectors
      INTEGER          NAEZ,               ! number of atoms in EZ
     +                 NR                  ! number of lattice vectors RR
c
c     .. output, written to 'trnc.unf'
c
      INTEGER          NUTRC,              ! number of inequivalent atoms in the cluster
     +                 INTRC(NATRCD),      ! pointer to atoms in the unit cell
     +                 NATRC,              ! number of atoms in cluster
     +                 ATTRC(NATRCD),      ! index to atom in elem/cell at site in cluster
     +                 EZTRC(NATRCD)       ! index to bravais lattice  at site in cluster
c
c     .. locals
c
      INTEGER          LRECTRC

      INTEGER          I,N1,
     +                 NA,NUMBER,N,
     +                 POS,IA,IN,IB,JATOM,IAT,
     +                 MXNATRC,MNNATRC,MXNUTRC,MNNUTRC
      INTEGER          IATOM(NATRCD),IEZOA(NATRCD),
     +                 ISORT(NATRCD)
      LOGICAL          ICOUPLMAT(NAEZ)
      DOUBLE PRECISION RCLS1(3,NATRCD),
     +                 RG(3,NATRCD),TMP(3),RSORT(NATRCD)
      DOUBLE PRECISION R,RCUTTOL,R2,EPSSHL
c
c
c
      EXTERNAL         XSORT
      INTRINSIC        MIN,SQRT
c
      DATA             EPSSHL   / 1.0D-4 /

c     initialise record length for trnc.unf file
      LRECTRC   = 4*(NATRCD*3+2)
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
     &     'CLSGEN_TRC: generation of TRC-clusters coordinates'
      WRITE (6,'(79(1H=))')
      WRITE (6,*)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
      RCUTTOL   = (RCUTTRC+EPSSHL)*(RCUTTRC+EPSSHL)
      MXNATRC   = 0
      MNNATRC   = 999999
      MXNUTRC   = 0
      MNNUTRC   = 999999
C
      OPEN (37,ACCESS='direct',RECL=LRECTRC,FILE='trnc.unf',
     +     FORM='unformatted')
C
C======================================================================
C loop in all atoms begin
C======================================================================

      DO JATOM = 1,NAEZ
         NUMBER = 0      ! counter for atoms in cluster
         DO NA = 1,NAEZ  ! loop in all atoms
            DO N=0,NR    ! loop in all bravais vectors    
               DO I=1,3
                  TMP(I) = RR(I,N)+RBASIS(I,NA)-RBASIS(I,JATOM)
               END DO
               R2   =  TMP(3)**2 + TMP(1)**2+TMP(2)**2
C
               IF (R2.LE.RCUTTOL)  THEN
                  NUMBER = NUMBER + 1
                  IF (NUMBER.GT.NATRCD) THEN 
                     WRITE (6,*) 
     &              ' ERROR: Dimension NATRCD in inc.cls too small',
     &                 NUMBER, NATRCD
                     STOP '   < CLSGEN_TRC >'
                  END IF
C
                  ATTRC(NUMBER) = NA ! store the atom in elem cell
                  EZTRC(NUMBER) = N ! store the bravais vector
                  DO I=1,3
                     RCLS1(I,NUMBER) = TMP(I)
                  END DO
               END IF
C
            END DO              ! N loop in bravais

         END DO                 ! NA loop in NAEZ

         NATRC = NUMBER
  
c     Now the atom JATOM Has it's cluster first 
c     sort the atoms of the cluster in increasing order. First by distance
c     Then by z then by y    
c     
         DO IA=1,NUMBER
            RSORT(IA) = SQRT(RCLS1(1,IA)**2+
     &                       RCLS1(2,IA)**2+
     &                       RCLS1(3,IA)**2)
            RSORT(IA) = 100000000.D0*RSORT(IA)+
     &                    10000.D0*RCLS1(3,IA)+
     &                   10.D0*RCLS1(2,IA)+
     &                  0.1D0*RCLS1(1,IA) 
         END DO      
c     
         CALL XSORT(RSORT,ISORT,NUMBER,POS)
c     Rearange exchange ia with ib
c MAP temporarily to another array         
         DO IA=1,NUMBER       
            DO I=1,3
               RG(I,IA)    = RCLS1(I,IA)
            END DO
            IATOM(IA) = ATTRC(IA)
            IEZOA(IA) = EZTRC(IA)
         END DO    
c Now use correct order
         DO IA =1,NUMBER
            IB = ISORT(IA)
             DO I=1,3
               RCLS1(I,IA) = RG(I,IB)
            END DO
            ATTRC(IA) = IATOM(IB)
            EZTRC(IA) = IEZOA(IB) 
         END DO

C
C define NUTRC and INTRC
C

         NUMBER = NATRC

         DO IAT = 1,NAEZ
           ICOUPLMAT(IAT) = .FALSE.
           DO I=1,NUMBER 
             IF (ATTRC(I).EQ.IAT) THEN
               ICOUPLMAT(IAT) = .TRUE.
             END IF
           END DO
         END DO

         NUTRC = 0
         DO IAT = 1,NAEZ
           IF(ICOUPLMAT(IAT)) THEN
             NUTRC = NUTRC + 1
             INTRC(NUTRC) = IAT
           ENDIF
         ENDDO
C
C write information on max and min truncation cluster size
C
         IF (MXNATRC.LT.NATRC) MXNATRC = NATRC
         IF (MNNATRC.GT.NATRC) MNNATRC = NATRC
         IF (MXNUTRC.LT.NUTRC) MXNUTRC = NUTRC
         IF (MNNUTRC.GT.NUTRC) MNNUTRC = NUTRC
C
C write information on individual truncation zones to trnc.unf
C

         WRITE(37,REC=JATOM) NATRC,ATTRC,EZTRC,NUTRC,INTRC

       END DO

C======================================================================
C loop in all atoms end
C======================================================================

      CLOSE(37)

      WRITE(6,*) 
      WRITE(6,*) ' No. of atoms in Truncation zone is in range [',
     +           MNNATRC,',',MXNATRC,']'
      WRITE(6,*) 
      WRITE(6,*) ' No. of ineq. atoms in Truncation zone is in range [',
     +           MNNUTRC,',',MXNUTRC,']'
      WRITE(6,*)
      IF (MXNUTRC.NE.NUTRCD) THEN
        WRITE(6,*) 'ERROR: set parameter NUTRCD in inc.p to ',MXNUTRC
        STOP '   < CLSGEN_TRC >'
      ENDIF
      WRITE(6,*) ' Sub clsgen_trc  exiting <<<<<<<<<<<<<'

      RETURN
      END
