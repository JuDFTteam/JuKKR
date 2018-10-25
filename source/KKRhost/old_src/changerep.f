C*==changerep.f    processed by SPAG 6.05Rc at 11:40 on 21 Jul 2001
      SUBROUTINE CHANGEREP(A,MODE,B,N,M,RC,CREL,RREL,TEXT,LTEXT)
C   ********************************************************************
C   *                                                                  *
C   *   change the representation of matrix A and store in B           *
C   *   according to MODE:                                             *
C   *                                                                  *
C   *   RLM>REL   non-relat. REAL spher. harm.  >   (kappa,mue)        *
C   *   REL>RLM   (kappa,mue)  > non-relat. REAL spher. harm.          *
C   *   CLM>REL   non-relat. CMPLX. spher. harm.  >   (kappa,mue)      *
C   *   REL>CLM   (kappa,mue)  > non-relat. CMPLX. spher. harm.        *
C   *   RLM>CLM   non-relat. REAL spher. harm.  >  CMPLX. spher. harm. *
C   *   CLM>RLM   non-relat. CMPLX. spher. harm.  >  REAL spher. harm. *
C   *                                                                  *
C   *   the non-relat. representations include the  spin index         *
C   *                                                                  *
C   *   for LTEXT > 0 the new matrix  B  is printed                    *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      COMPLEX*16 C1,C0
      PARAMETER (C1=(1.0D0,0.0D0),C0=(0.0D0,0.0D0))
C
C Dummy arguments
C
      INTEGER LTEXT,M,N
      CHARACTER*7 MODE
      CHARACTER*(*) TEXT
      COMPLEX*16 A(M,M),B(M,M),CREL(M,M),RC(M,M),RREL(M,M)
C
C Local variables
C
      INTEGER KEY
      COMPLEX*16 W1(M,M)
C
C*** End of declarations rewritten by SPAG
C
C---------------------- transform MAT from (kappa,mue) to REAL (l,ml,ms)
      IF ( MODE.EQ.'REL>RLM' ) THEN
         CALL ZGEMM('N','N',N,N,N,C1,RREL,M,A,M,C0,W1,M)
         CALL ZGEMM('N','C',N,N,N,C1,W1,M,RREL,M,C0,B,M)
         KEY = 2
      ELSE IF ( MODE.EQ.'RLM>REL' ) THEN
         CALL ZGEMM('C','N',N,N,N,C1,RREL,M,A,M,C0,W1,M)
         CALL ZGEMM('N','N',N,N,N,C1,W1,M,RREL,M,C0,B,M)
         KEY = 3
      ELSE IF ( MODE.EQ.'REL>CLM' ) THEN
         CALL ZGEMM('N','N',N,N,N,C1,CREL,M,A,M,C0,W1,M)
         CALL ZGEMM('N','C',N,N,N,C1,W1,M,CREL,M,C0,B,M)
         KEY = 2
      ELSE IF ( MODE.EQ.'CLM>REL' ) THEN
         CALL ZGEMM('C','N',N,N,N,C1,CREL,M,A,M,C0,W1,M)
         CALL ZGEMM('N','N',N,N,N,C1,W1,M,CREL,M,C0,B,M)
         KEY = 3
      ELSE IF ( MODE.EQ.'CLM>RLM' ) THEN
         CALL ZGEMM('N','N',N,N,N,C1,RC,M,A,M,C0,W1,M)
         CALL ZGEMM('N','C',N,N,N,C1,W1,M,RC,M,C0,B,M)
         KEY = 2
      ELSE IF ( MODE.EQ.'RLM>CLM' ) THEN
         CALL ZGEMM('C','N',N,N,N,C1,RC,M,A,M,C0,W1,M)
         CALL ZGEMM('N','N',N,N,N,C1,W1,M,RC,M,C0,B,M)
         KEY = 2
      ELSE
         WRITE (*,*) ' MODE = ',MODE
         STOP 'in <ROTATE>  MODE not allowed'
      END IF
C
      IF ( LTEXT.GT.0 ) CALL CMATSTR(TEXT,LTEXT,B,N,M,KEY,KEY,0,1D-8,6)
C     IF ( LTEXT.GT.0 ) CALL CMATSTR(TEXT,LTEXT,B,N,M,KEY,KEY,0,1D-12,6)
      END
