SUBROUTINE rotate(t1,mode,t2,n,rot,nkmmax)
!   ********************************************************************
!   *                                                                  *
!   *   performs the rotation of the matrix  T1  using the rotation-   *
!   *   matrix  ROT, set up by <CALCROTMAT>                            *
!   *                                                                  *
!   *          T2 = ROT  * T1 * ROT+     IF  MODE = 'L->G'             *
!   *          T2 = ROT+ * T1 * ROT      IF  MODE = 'G->L'             *
!   *                                                                  *
!   *   see:     E.M. ROSE  ELEMENTARY THEORY OF ANGULAR MOMENTUM      *
!   *                                                                  *
!   * 01/11/00                                                         *
!   ********************************************************************
IMPLICIT NONE

! PARAMETER definitions
COMPLEX*16 C0,C1
PARAMETER (C0=(0.0D0,0.0D0),C1=(1.0D0,0.0D0))

! Dummy arguments
CHARACTER*4 MODE
INTEGER N,NKMMAX
COMPLEX*16 ROT(NKMMAX,NKMMAX),T1(NKMMAX,NKMMAX),T2(NKMMAX,NKMMAX)

! Local variables
CHARACTER*1 FL1,FL2
COMPLEX*16 W1(NKMMAX,NKMMAX)


IF ( mode == 'L->G' ) THEN
  fl1 = 'N'
  fl2 = 'C'
ELSE IF ( mode == 'G->L' ) THEN
  fl1 = 'C'
  fl2 = 'N'
ELSE
  WRITE (*,*) ' MODE = ',mode
  STOP 'in <ROTATE>  MODE not allowed'
END IF

CALL zgemm(fl1,'N',n,n,n,c1,rot,nkmmax,t1,nkmmax,c0,w1,nkmmax)

CALL zgemm('N',fl2,n,n,n,c1,w1,nkmmax,rot,nkmmax,c0,t2,nkmmax)

END SUBROUTINE rotate
