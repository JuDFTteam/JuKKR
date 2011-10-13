c ************************************************************************
      SUBROUTINE CLSJIJ(
     >                  I1,NAEZ,RR,NR,RBASIS,RCUT,LMPIC,NSYMAT,
     >                  ISYMINDEX,
     <                  IXCP,NXCP,NXIJ,RXIJ,RXCCLS,ZKRXIJ)
c ************************************************************************
c This subroutine is used to create the clusters around each atom 
c where Jij's are calculated
c
c called by: main2
c
c STRATEGY :
c Calculate the cluster of each atom by the lattice
c parameters avaliable. Sort the atoms in a unique way :big r, big z, big y
c this routine is executed only for the assigned atom I1 of the 
c responsible processor
c                                                          A.Thiess 7/2009
c ************************************************************************
      IMPLICIT NONE
c
      INCLUDE 'inc.p'
c
c
c     .. array arguments
c
      DOUBLE PRECISION RBASIS(3,NAEZD),   ! pos. of basis atoms in EZ
     +                 RR(3,0:NRD),       ! set of lattice vectors
     +                 RXIJ(NXIJD),       ! interatomic distance Ri-Rj 
     +                 ZKRXIJ(48,3,NXIJD) ! enters the exp-factor of G in kkrmat01
      INTEGER          IXCP(NXIJD),       ! index to atom in elem/cell at site in cluster
     +                 NXCP(NXIJD),       ! index to bravais lattice  at site in cluster
     +                 ISYMINDEX(48)
c
c
c     .. scalar arguments
c
      DOUBLE PRECISION RCUT
      INTEGER          NXIJ,            ! number of atoms in cluster
     +                 NAEZ,            ! number of atoms in EZ
     +                 NR,              ! number of lattice vectors RR
     +                 I1,              ! processor calling this routines deals with atom I1
     +                 LMPIC,
     +                 NSYMAT
c
c     .. local arrays
c
      DOUBLE PRECISION RXCCLS(3,NXIJD),   ! real space pos of atom in cluster
     +                 TMP(3),
     +                 IRCLS(3,NXIJD),RSORT(NXIJD),
     +                 RMAT(64,3,3)
      INTEGER          IIXCP(NXIJD),INXCP(NXIJD),ISORT(NXIJD)
c
c     .. local scalars
c
      DOUBLE PRECISION EPSSHL,RCUT2,RTMP
      INTEGER          IAEZ,IB,ID,IR,IV,IX,POS
      CHARACTER*10     ROTNAME(64)
c
c
      EXTERNAL         XSORT,POINTGRP
      INTRINSIC        SQRT
c
      DATA             EPSSHL   / 1.0D-4 /
c
c ------------------------------------------------------------------------
c This is generating the clusters which have a distance smaller
c than RCUT and RCUTXY in plane .
c The cluster atoms are ordered with radious and then z>y>x 
c The ordering allows an easy comparison of clusters.
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       RCUT2   = (RCUT+EPSSHL)*(RCUT+EPSSHL)
C======================================================================
C loop in all atoms begin
C======================================================================

       NXIJ = 0           ! counter for atoms in cluster
       DO IAEZ = 1,NAEZ  ! loop in all atoms
         DO IR = 0, NR    ! loop in all bravais vectors    
           DO IV=1,3
             TMP(IV) = RR(IV,IR)+RBASIS(IV,IAEZ)-RBASIS(IV,I1)
           ENDDO
           RTMP   =  TMP(3)**2 + TMP(1)**2+TMP(2)**2
           IF (RTMP.LE.RCUT2)  THEN
             NXIJ = NXIJ + 1
             IF (NXIJ.GT.NXIJD) THEN 
               WRITE (6,*) 
     &         ' ERROR: Dimension NXIJD in inc.cls too small',
     &         NXIJ, NXIJD
               STOP '   < CLSJIJ >'
             ENDIF
c
             IXCP(NXIJ) = IAEZ       ! store the atom in elem cell
             NXCP(NXIJ) = IR         ! store the bravais vector
c
             DO IV=1,3
               RXCCLS(IV,NXIJ) = TMP(IV)
             ENDDO
           ENDIF
         ENDDO              ! IR loop in bravais

       ENDDO                 ! IAEZ loop in NAEZ
c
c     sort the atoms of the cluster in increasing order. First by distance
c     Then by z then by y
c
        DO IX=1,NXIJ
            RSORT(IX) = SQRT(RXCCLS(1,IX)**2+
     &                       RXCCLS(2,IX)**2+
     &                       RXCCLS(3,IX)**2)
            RSORT(IX) = 100000000.D0*RSORT(IX)+
     &                      10000.D0*RXCCLS(3,IX)+
     &                         10.D0*RXCCLS(2,IX)+
     &                         0.1D0*RXCCLS(1,IX)
        ENDDO
c
        CALL XSORT(RSORT,ISORT,NXIJ,POS)
c
c     Rearange exchange ia with ib
c MAP temporarily to another array
c
        DO IX=1,NXIJ
          DO IV=1,3
            IRCLS(IV,IX)    = RXCCLS(IV,IX)
          ENDDO
          IIXCP(IX) = IXCP(IX)
          INXCP(IX) = NXCP(IX)
        ENDDO
c
c Now use correct order
c
        CALL POINTGRP(RMAT,ROTNAME)
c
        DO IX =1,NXIJ
          IB = ISORT(IX)
          DO IV=1,3
            RXCCLS(IV,IX) = IRCLS(IV,IB)
          ENDDO
          IXCP(IX) = IIXCP(IB)
          NXCP(IX) = INXCP(IB) 
          RXIJ(IX) = 
     +    SQRT(RXCCLS(1,IX)**2+RXCCLS(2,IX)**2+RXCCLS(3,IX)**2) ! store interatomic distance
c
          DO ID = 1,NSYMAT
            DO IV = 1,3
              ZKRXIJ(ID,IV,IX) = RMAT(ISYMINDEX(ID),IV,1)*RXCCLS(1,IX) +
     +                           RMAT(ISYMINDEX(ID),IV,2)*RXCCLS(2,IX) +
     +                           RMAT(ISYMINDEX(ID),IV,3)*RXCCLS(3,IX) -
     +                           RBASIS(IV,IXCP(IX)) +
     +                           RBASIS(IV,I1)    !ART
            ENDDO
          ENDDO
c
        ENDDO
c
      END




