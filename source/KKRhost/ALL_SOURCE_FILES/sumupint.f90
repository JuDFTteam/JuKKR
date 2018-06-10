SUBROUTINE sumupint(sum,vg,g,wg,vf,f,wf,n)
!   ********************************************************************
!   *                                                                  *
!   ********************************************************************
IMPLICIT NONE


! Dummy arguments
INTEGER N
COMPLEX*16 SUM
REAL*8 VF,VG
COMPLEX*16 F(2,2),G(2,2)
REAL*8 WF(2,2),WG(2,2)

! Local variables
INTEGER I,J

sum = 0.0D0
DO j = 1,n
  DO i = 1,n
    sum = sum + vg*g(i,j)*wg(i,j) + vf*f(i,j)*wf(i,j)
  END DO
END DO

END SUBROUTINE sumupint
