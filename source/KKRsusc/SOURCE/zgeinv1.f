      MODULE mod_zgeinv1
      CONTAINS
      SUBROUTINE ZGEINV1(A,U,AUX,IPIV,DIM)
c ************************************************************************
C   - inverts a general double complex matrix A,
c   - the result is return in U,
c   - input matrix A is returned unchanged,
c   - AUX is a auxiliary matrix,
c   - A,U and AUX are of dimension (DIM,DIM),
c ------------------------------------------------------------------------
!       USE mod_cinit
      IMPLICIT NONE
      INTEGER DIM,IPIV(:)
      DOUBLE COMPLEX A(DIM,DIM),AUX(DIM,DIM),U(DIM,DIM)
c
C     .. PARAMETER
C
      DOUBLE COMPLEX CONE
      PARAMETER(CONE=(1.D0,0.D0))
      DOUBLE COMPLEX CZERO
      PARAMETER(CZERO=(0.D0,0.D0))
C
      INTEGER LM1,INFO
      EXTERNAL ZCOPY,ZGETRS,ZGETRF
c ------------------------------------------------------------------------
!       CALL CINIT(DIM*DIM,U)
      U = CZERO
      DO 10 LM1=1,DIM
        U(LM1,LM1) = CONE
 10   END DO

      CALL ZCOPY(DIM*DIM,A,1,AUX,1)
      CALL ZGETRF(DIM,DIM,AUX,DIM,IPIV,INFO)
      CALL ZGETRS('N',DIM,DIM,AUX,DIM,IPIV,U,DIM,INFO)
      
      RETURN
      END SUBROUTINE
      END MODULE mod_zgeinv1