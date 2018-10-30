C*==rotate.f    processed by SPAG 6.05Rc at 14:31 on  1 Nov 2000
      SUBROUTINE ROTATE(T1,MODE,T2,N,ROT,NKMMAX)
C   ********************************************************************
C   *                                                                  *
C   *   performs the rotation of the matrix  T1  using the rotation-   *
C   *   matrix  ROT, set up by <CALCROTMAT>                            *
C   *                                                                  *
C   *          T2 = ROT  * T1 * ROT+     IF  MODE = 'L->G'             *
C   *          T2 = ROT+ * T1 * ROT      IF  MODE = 'G->L'             *
C   *                                                                  *
C   *   see:     E.M. ROSE  ELEMENTARY THEORY OF ANGULAR MOMENTUM      *
C   *                                                                  *
C   * 01/11/00                                                         *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      COMPLEX*16 C0,C1
      PARAMETER (C0=(0.0D0,0.0D0),C1=(1.0D0,0.0D0))
C
C Dummy arguments
C
      CHARACTER*4 MODE
      INTEGER N,NKMMAX
      COMPLEX*16 ROT(NKMMAX,NKMMAX),T1(NKMMAX,NKMMAX),T2(NKMMAX,NKMMAX)
C
C Local variables
C
      CHARACTER*1 FL1,FL2
      COMPLEX*16 W1(NKMMAX,NKMMAX)
C
C*** End of declarations rewritten by SPAG
C
      IF ( MODE.EQ.'L->G' ) THEN
         FL1 = 'N'
         FL2 = 'C'
      ELSE IF ( MODE.EQ.'G->L' ) THEN
         FL1 = 'C'
         FL2 = 'N'
      ELSE
         WRITE (*,*) ' MODE = ',MODE
         STOP 'in <ROTATE>  MODE not allowed'
      END IF
C
      CALL ZGEMM(FL1,'N',N,N,N,C1,ROT,NKMMAX,T1,NKMMAX,C0,W1,NKMMAX)
C
      CALL ZGEMM('N',FL2,N,N,N,C1,W1,NKMMAX,ROT,NKMMAX,C0,T2,NKMMAX)
C
      END
