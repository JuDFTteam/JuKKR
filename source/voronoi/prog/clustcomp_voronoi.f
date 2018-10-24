      LOGICAL FUNCTION CLUSTCOMP_VORONOI(RCLS,IC1,N1,RCLS1,N2)
      implicit none
c#@# KKRtags: VORONOI geometry
c  This func returns true if cluster number ic1
c  is equal to cluster ic2
c  RCLS        clusters coordinates
c  IC1         First cluster
c  N1          Number of atoms in IC1 cluster
c  IC2         Second cluster
c  rcls1       coordinates of second cluster
c  n2          number of atoms in ic2 cluster
c        include 'inc.p'
      include 'inc.geometry'
      REAL*8  RCLS(3,NACLSD,*),RCLS1(3,NACLSD)
      INTEGER IC1,N1,IC2,N2,N,I
      REAL*8  RD,TOL
      TOL = 1.D-5
      CLUSTCOMP_VORONOI = .FALSE.
      IF (N1.EQ.N2) THEN
         RD = 0.D0
         DO N=1,N1
            DO I=1,3
               RD = RD + ABS(RCLS(I,N,IC1) - RCLS1(I,N)) 
            END DO
         END DO
         IF (ABS(RD).LT.TOL) CLUSTCOMP_VORONOI = .TRUE.
      END IF
      END 
