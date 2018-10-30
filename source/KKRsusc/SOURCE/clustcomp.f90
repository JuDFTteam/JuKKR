MODULE MOD_CLUSTCOMP
CONTAINS
      LOGICAL FUNCTION CLUSTCOMP(RCLS,RMT,ATOM,IATOM, &
                                 IC1,N1,RCLS1,N2,JATOM,NACLSD,NR,LMAXATOM)
!  This function returns true if cluster number ic1
!  is equal to cluster ic2
!  RCLS        clusters coordinates
!  IC1         First cluster
!  N1          Number of atoms in IC1 cluster
!  IC2         Second cluster
!
!     .. Scalar Arguments ..
      INTEGER IATOM,IC1,JATOM,N1,N2,NACLSD,NR
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION RCLS(3,NACLSD,NR),RCLS1(3,NACLSD)
      INTEGER ATOM(NACLSD,NR)
      DOUBLE PRECISION RMT(NR)
      INTEGER LMAXATOM(NR)

!            IF( CLUSTCOMP(RCLS,RMT,ATOM,IATCLS(ICU),ICU,N1,RCLS1, &
!                 NUMBER,IATOM,NACLSD)) CLS(IATOM) = ICU

!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION R
      INTEGER I,N
      LOGICAL REFLOG
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS
!     ..
      CLUSTCOMP = .FALSE.
      REFLOG = .TRUE.
      IF (N1.EQ.N2) THEN
        R = 0
        DO N = 1,N1
          IF (RMT(ABS(ATOM(N,IATOM))).NE. RMT(ABS(ATOM(N,JATOM))) .OR. &
              LMAXATOM(ABS(ATOM(N,IATOM))) .NE. LMAXATOM(ABS(ATOM(N,JATOM)))  ) REFLOG = .FALSE.
          DO I = 1,3
            R = R + ABS(RCLS(I,N,IC1)-RCLS1(I,N))
          END DO
        END DO
        IF (ABS(R).LT.1.D-6 .AND. REFLOG) CLUSTCOMP = .TRUE.
      END IF
      END FUNCTION
END MODULE MOD_CLUSTCOMP
