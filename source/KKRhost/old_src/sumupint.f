      SUBROUTINE SUMUPINT(SUM,VG,G,WG,VF,F,WF,N)
      IMPLICIT NONE
C
C
C Dummy arguments
C
      INTEGER N
      COMPLEX*16 SUM
      REAL*8 VF,VG
      COMPLEX*16 F(2,2),G(2,2)
      REAL*8 WF(2,2),WG(2,2)
C
C Local variables
C
      INTEGER I,J
C
C   ********************************************************************
C   *                                                                  *
C   ********************************************************************
C
      SUM = 0.0D0
      DO J = 1,N
         DO I = 1,N
            SUM = SUM + VG*G(I,J)*WG(I,J) + VF*F(I,J)*WF(I,J)
         END DO
      END DO
C
      END
