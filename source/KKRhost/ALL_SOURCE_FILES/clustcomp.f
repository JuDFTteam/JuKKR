       LOGICAL FUNCTION CLUSTCOMP_TB(
     >     RCLS,IREFPOT,ATOM,IAT1,IC1,N1,RCLS1,N2,IAT2,NACLSD)
      implicit none
c     This function returns true if cluster ic1 is equal to new cluster
c     RCLS        coordinates of all (already found) clusters
c     IC1         First cluster index
c     N1          Number of atoms in IC1 cluster
c     rcls1       coordinates of new cluster
c     n2          number of atoms in ic2 cluster
c     ATOM:       atom-type at a certain position in a cluster
c     IAT1, IAT2: central atoms of 1st,2nd cluster
c     IREFPOT:    Type of reference potential
      INTEGER NACLSD
      REAL*8  RCLS(3,NACLSD,*),RCLS1(3,NACLSD)
      INTEGER ATOM(NACLSD,*),IREFPOT(*)
      INTEGER IC1,N1,IC2,N2,IAT1,IAT2,N,I
      REAL*8  RD,TOL
      LOGICAL LREFLOG
      TOL = 1.D-5
      CLUSTCOMP_TB = .FALSE.
      LREFLOG = .TRUE.
      IF (N1.EQ.N2) THEN
         RD = 0.D0
         DO N=1,N1
            IF (IREFPOT(ABS(ATOM(N,IAT1))).NE.                 ! compare ref-potential types
     +          IREFPOT(ABS(ATOM(N,IAT2)))) LREFLOG = .FALSE.
            DO I=1,3
               RD = RD + ABS(RCLS(I,N,IC1) - RCLS1(I,N))  ! compare coordinares
            END DO
         END DO
         IF (ABS(RD).LT.TOL.AND.LREFLOG) CLUSTCOMP_TB = .TRUE.
      END IF
      END 
