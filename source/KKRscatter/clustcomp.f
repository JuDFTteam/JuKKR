      LOGICAL FUNCTION CLUSTCOMP(RCLS,REFPOT,ATOM,IATOM,
     &                           IC1,N1,RCLS1,N2,JATOM,NACLSD)
c  This function returns true if cluster number ic1
c  is equal to cluster ic2
c  RCLS        clusters coordinates
c  IC1         First cluster
c  N1          Number of atoms in IC1 cluster
c  IC2         Second cluster
c
C     .. Scalar Arguments ..
      INTEGER IATOM,IC1,JATOM,N1,N2,NACLSD
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION RCLS(3,NACLSD,*),RCLS1(3,NACLSD)
      INTEGER ATOM(NACLSD,*),REFPOT(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION R
      INTEGER I,N
      LOGICAL REFLOG
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      CLUSTCOMP = .FALSE.
      REFLOG = .TRUE.
      IF (N1.EQ.N2) THEN
        R = 0
        DO N = 1,N1
          IF (REFPOT(ABS(ATOM(N,IATOM))).NE.
     +        REFPOT(ABS(ATOM(N,JATOM)))) REFLOG = .FALSE.
          DO I = 1,3
            R = R + ABS(RCLS(I,N,IC1)-RCLS1(I,N))
          END DO
        END DO
        IF (ABS(R).LT.1.D-6 .AND. REFLOG) THEN
        CLUSTCOMP = .TRUE.
        ENDIF
      END IF
      END
