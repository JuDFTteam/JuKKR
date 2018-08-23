module mod_changerep

contains

subroutine changerep(a, mode, b, n, m, rc, crel, rrel, text, ltext)
  ! ********************************************************************
  ! *                                                                  *
  ! *   change the representation of matrix A and store in B           *
  ! *   according to MODE:                                             *
  ! *                                                                  *
  ! *   RLM>REL   non-relat. REAL spher. harm.  >   (kappa,mue)        *
  ! *   REL>RLM   (kappa,mue)  > non-relat. REAL spher. harm.          *
  ! *   CLM>REL   non-relat. CMPLX. spher. harm.  >   (kappa,mue)      *
  ! *   REL>CLM   (kappa,mue)  > non-relat. CMPLX. spher. harm.        *
  ! *   RLM>CLM   non-relat. REAL spher. harm.  >  CMPLX. spher. harm. *
  ! *   CLM>RLM   non-relat. CMPLX. spher. harm.  >  REAL spher. harm. *
  ! *                                                                  *
  ! *   the non-relat. representations include the  spin index         *
  ! *                                                                  *
  ! *   for LTEXT > 0 the new matrix  B  is printed                    *
  ! *                                                                  *
  ! ********************************************************************
  use :: mod_datatypes, only: dp
  implicit none

  ! PARAMETER definitions
  complex (kind=dp) :: c1, c0
  parameter (c1=(1.0e0_dp,0.0e0_dp), c0=(0.0e0_dp,0.0e0_dp))

  ! Dummy arguments
  integer :: ltext, m, n
  character (len=7) :: mode
  character (len=*) :: text
  complex (kind=dp) :: a(m, m), b(m, m), crel(m, m), rc(m, m), rrel(m, m)

  ! Local variables
  integer :: key
  complex (kind=dp) :: w1(m, m)


  ! ---------------------- transform MAT from (kappa,mue) to REAL (l,ml,ms)
  if (mode=='REL>RLM') then
    call zgemm('N', 'N', n, n, n, c1, rrel, m, a, m, c0, w1, m)
    call zgemm('N', 'C', n, n, n, c1, w1, m, rrel, m, c0, b, m)
    key = 2
  else if (mode=='RLM>REL') then
    call zgemm('C', 'N', n, n, n, c1, rrel, m, a, m, c0, w1, m)
    call zgemm('N', 'N', n, n, n, c1, w1, m, rrel, m, c0, b, m)
    key = 3
  else if (mode=='REL>CLM') then
    call zgemm('N', 'N', n, n, n, c1, crel, m, a, m, c0, w1, m)
    call zgemm('N', 'C', n, n, n, c1, w1, m, crel, m, c0, b, m)
    key = 2
  else if (mode=='CLM>REL') then
    call zgemm('C', 'N', n, n, n, c1, crel, m, a, m, c0, w1, m)
    call zgemm('N', 'N', n, n, n, c1, w1, m, crel, m, c0, b, m)
    key = 3
  else if (mode=='CLM>RLM') then
    call zgemm('N', 'N', n, n, n, c1, rc, m, a, m, c0, w1, m)
    call zgemm('N', 'C', n, n, n, c1, w1, m, rc, m, c0, b, m)
    key = 2
  else if (mode=='RLM>CLM') then
    call zgemm('C', 'N', n, n, n, c1, rc, m, a, m, c0, w1, m)
    call zgemm('N', 'N', n, n, n, c1, w1, m, rc, m, c0, b, m)
    key = 2
  else
    write (*, *) ' MODE = ', mode
    stop 'in <ROTATE>  MODE not allowed'
  end if

  if (ltext>0) call cmatstr(text, ltext, b, n, m, key, key, 0, 1e-8_dp, 6)
  ! IF ( LTEXT.GT.0 ) CALL CMATSTR(TEXT,LTEXT,B,N,M,KEY,KEY,0,1D-12,6)
end subroutine changerep

end module mod_changerep
