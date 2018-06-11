    Subroutine changerep(a, mode, b, n, m, rc, crel, rrel, text, ltext)
      Use mod_datatypes, Only: dp
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
      Implicit None

! PARAMETER definitions
      Complex (Kind=dp) :: c1, c0
      Parameter (c1=(1.0E0_dp,0.0E0_dp), c0=(0.0E0_dp,0.0E0_dp))

! Dummy arguments
      Integer :: ltext, m, n
      Character (Len=7) :: mode
      Character (Len=*) :: text
      Complex (Kind=dp) :: a(m, m), b(m, m), crel(m, m), rc(m, m), rrel(m, m)

! Local variables
      Integer :: key
      Complex (Kind=dp) :: w1(m, m)


!---------------------- transform MAT from (kappa,mue) to REAL (l,ml,ms)
      If (mode=='REL>RLM') Then
        Call zgemm('N', 'N', n, n, n, c1, rrel, m, a, m, c0, w1, m)
        Call zgemm('N', 'C', n, n, n, c1, w1, m, rrel, m, c0, b, m)
        key = 2
      Else If (mode=='RLM>REL') Then
        Call zgemm('C', 'N', n, n, n, c1, rrel, m, a, m, c0, w1, m)
        Call zgemm('N', 'N', n, n, n, c1, w1, m, rrel, m, c0, b, m)
        key = 3
      Else If (mode=='REL>CLM') Then
        Call zgemm('N', 'N', n, n, n, c1, crel, m, a, m, c0, w1, m)
        Call zgemm('N', 'C', n, n, n, c1, w1, m, crel, m, c0, b, m)
        key = 2
      Else If (mode=='CLM>REL') Then
        Call zgemm('C', 'N', n, n, n, c1, crel, m, a, m, c0, w1, m)
        Call zgemm('N', 'N', n, n, n, c1, w1, m, crel, m, c0, b, m)
        key = 3
      Else If (mode=='CLM>RLM') Then
        Call zgemm('N', 'N', n, n, n, c1, rc, m, a, m, c0, w1, m)
        Call zgemm('N', 'C', n, n, n, c1, w1, m, rc, m, c0, b, m)
        key = 2
      Else If (mode=='RLM>CLM') Then
        Call zgemm('C', 'N', n, n, n, c1, rc, m, a, m, c0, w1, m)
        Call zgemm('N', 'N', n, n, n, c1, w1, m, rc, m, c0, b, m)
        key = 2
      Else
        Write (*, *) ' MODE = ', mode
        Stop 'in <ROTATE>  MODE not allowed'
      End If

      If (ltext>0) Call cmatstr(text, ltext, b, n, m, key, key, 0, 1E-8_dp, 6)
!     IF ( LTEXT.GT.0 ) CALL CMATSTR(TEXT,LTEXT,B,N,M,KEY,KEY,0,1D-12,6)
    End Subroutine
