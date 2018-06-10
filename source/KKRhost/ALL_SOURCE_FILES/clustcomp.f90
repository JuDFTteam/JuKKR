LOGICAL FUNCTION clustcomp_tb(  &
        rcls,irefpot,atom,iat1,ic1,n1,rcls1,n2,iat2,naclsd)
!     This function returns true if cluster ic1 is equal to new cluster
!     RCLS        coordinates of all (already found) clusters
!     IC1         First cluster index
!     N1          Number of atoms in IC1 cluster
!     rcls1       coordinates of new cluster
!     n2          number of atoms in ic2 cluster
!     ATOM:       atom-type at a certain position in a cluster
!     IAT1, IAT2: central atoms of 1st,2nd cluster
!     IREFPOT:    Type of reference potential
IMPLICIT NONE
INTEGER NACLSD
REAL*8  RCLS(3,NACLSD,*),RCLS1(3,NACLSD)
INTEGER ATOM(NACLSD,*),IREFPOT(*)
INTEGER IC1,N1,N2,IAT1,IAT2,N,I
REAL*8  RD,TOL
LOGICAL LREFLOG

tol = 1.d-5
clustcomp_tb = .false.
lreflog = .true.
IF (n1 == n2) THEN
  rd = 0.d0
  DO n=1,n1
    IF (irefpot(ABS(atom(n,iat1))) /=  &               ! compare ref-potential types
        irefpot(ABS(atom(n,iat2)))) lreflog = .false.
    DO i=1,3
      rd = rd + ABS(rcls(i,n,ic1) - rcls1(i,n))  ! compare coordinares
    END DO
  END DO
  IF (ABS(rd) < tol.AND.lreflog) clustcomp_tb = .true.
END IF
END FUNCTION clustcomp_tb
