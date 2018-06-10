SUBROUTINE changerep(a,mode,b,n,m,rc,crel,rrel,text,ltext)
!   ********************************************************************
!   *                                                                  *
!   *   change the representation of matrix A and store in B           *
!   *   according to MODE:                                             *
!   *                                                                  *
!   *   RLM>REL   non-relat. REAL spher. harm.  >   (kappa,mue)        *
!   *   REL>RLM   (kappa,mue)  > non-relat. REAL spher. harm.          *
!   *   CLM>REL   non-relat. CMPLX. spher. harm.  >   (kappa,mue)      *
!   *   REL>CLM   (kappa,mue)  > non-relat. CMPLX. spher. harm.        *
!   *   RLM>CLM   non-relat. REAL spher. harm.  >  CMPLX. spher. harm. *
!   *   CLM>RLM   non-relat. CMPLX. spher. harm.  >  REAL spher. harm. *
!   *                                                                  *
!   *   the non-relat. representations include the  spin index         *
!   *                                                                  *
!   *   for LTEXT > 0 the new matrix  B  is printed                    *
!   *                                                                  *
!   ********************************************************************
IMPLICIT none

! PARAMETER definitions
COMPLEX*16 C1,C0
PARAMETER (C1=(1.0D0,0.0D0),C0=(0.0D0,0.0D0))

! Dummy arguments
INTEGER LTEXT,M,N
CHARACTER*7 MODE
CHARACTER*(*) TEXT
COMPLEX*16 A(M,M),B(M,M),CREL(M,M),RC(M,M),RREL(M,M)

! Local variables
INTEGER KEY
COMPLEX*16 W1(M,M)


!---------------------- transform MAT from (kappa,mue) to REAL (l,ml,ms)
IF ( mode == 'REL>RLM' ) THEN
  CALL zgemm('N','N',n,n,n,c1,rrel,m,a,m,c0,w1,m)
  CALL zgemm('N','C',n,n,n,c1,w1,m,rrel,m,c0,b,m)
  key = 2
ELSE IF ( mode == 'RLM>REL' ) THEN
  CALL zgemm('C','N',n,n,n,c1,rrel,m,a,m,c0,w1,m)
  CALL zgemm('N','N',n,n,n,c1,w1,m,rrel,m,c0,b,m)
  key = 3
ELSE IF ( mode == 'REL>CLM' ) THEN
  CALL zgemm('N','N',n,n,n,c1,crel,m,a,m,c0,w1,m)
  CALL zgemm('N','C',n,n,n,c1,w1,m,crel,m,c0,b,m)
  key = 2
ELSE IF ( mode == 'CLM>REL' ) THEN
  CALL zgemm('C','N',n,n,n,c1,crel,m,a,m,c0,w1,m)
  CALL zgemm('N','N',n,n,n,c1,w1,m,crel,m,c0,b,m)
  key = 3
ELSE IF ( mode == 'CLM>RLM' ) THEN
  CALL zgemm('N','N',n,n,n,c1,rc,m,a,m,c0,w1,m)
  CALL zgemm('N','C',n,n,n,c1,w1,m,rc,m,c0,b,m)
  key = 2
ELSE IF ( mode == 'RLM>CLM' ) THEN
  CALL zgemm('C','N',n,n,n,c1,rc,m,a,m,c0,w1,m)
  CALL zgemm('N','N',n,n,n,c1,w1,m,rc,m,c0,b,m)
  key = 2
ELSE
  WRITE (*,*) ' MODE = ',mode
  STOP 'in <ROTATE>  MODE not allowed'
END IF

IF ( ltext > 0 ) CALL cmatstr(text,ltext,b,n,m,key,key,0,1D-8,6)
!     IF ( LTEXT.GT.0 ) CALL CMATSTR(TEXT,LTEXT,B,N,M,KEY,KEY,0,1D-12,6)
END SUBROUTINE changerep
