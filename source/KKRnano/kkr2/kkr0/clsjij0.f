c ************************************************************************
      SUBROUTINE CLSJIJ0(
     &                   NAEZ,RR,NR,RBASIS,RCUT,JIJ,NRD,NXIJD)
c ************************************************************************
c This subroutine is used check the cluster size around each atom 
c where Jij's are calculated
c all arguments are input arguments
c
c called by: main0
c                                                          A.Thiess 7/2009
c ************************************************************************
      IMPLICIT NONE
c
c      INCLUDE 'inc.p'
c
c
c     .. array arguments
c
      DOUBLE PRECISION RBASIS(3,NAEZ),   ! pos. of basis atoms in EZ
     +                 RR(3,0:NRD)        ! set of lattice vectors
c
c
c     .. scalar arguments
c
      INTEGER, INTENT(IN) :: NRD
      INTEGER, INTENT(IN) :: NXIJD
      DOUBLE PRECISION RCUT
      INTEGER          NAEZ,            ! number of atoms in EZ
     +                 NR               ! number of lattice vectors RR
      LOGICAL          JIJ
c
c     .. local arrays
c
      DOUBLE PRECISION TMP(3)
c
c     .. local scalars
c
      DOUBLE PRECISION EPSSHL,RCUT2,RTMP
      INTEGER          IAEZ,IR,IV,NXIJ,I1
c
c
      DATA             EPSSHL   / 1.0D-4 /
c
c ------------------------------------------------------------------------
c This is generating the clusters which have a distance smaller
c than RCUT2.
c ------------------------------------------------------------------------
       IF (JIJ) THEN
C
       RCUT2   = (RCUT+EPSSHL)*(RCUT+EPSSHL)
C======================================================================
C loop in all atoms begin
C======================================================================
       DO I1 = 1,NAEZ ! loop over all Jij-centers
C
         NXIJ = 0           ! counter for atoms in cluster
         DO IAEZ = 1,NAEZ   ! loop in all atoms
           DO IR = 0, NR    ! loop in all bravais vectors    
             DO IV=1,3
               TMP(IV) = RR(IV,IR)+RBASIS(IV,IAEZ)-RBASIS(IV,I1)
             ENDDO
             RTMP = TMP(3)**2 + TMP(1)**2+TMP(2)**2
             IF (RTMP.LE.RCUT2)  THEN
               NXIJ = NXIJ + 1
               IF (NXIJ.GT.NXIJD) THEN 
                 WRITE (6,*) 
     &           ' ERROR: Dimension NXIJD in inc.cls too small',
     &           NXIJ, NXIJD
                 STOP '   < CLSJIJ >'
               ENDIF
             ENDIF
           ENDDO              ! IR loop in bravais
         ENDDO                ! IAEZ loop in NAEZ
C
       ENDDO ! loop over all Jij-centers

      WRITE(6,'(79(1H=),/,15X,A)') 
     &     'CLSJIJ0: checking Jij-cluster size ........ OK'
      WRITE (6,'(79(1H=),/)')
C======================================================================
C======================================================================
       ENDIF
C
      END




