!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_changerep
  
  private
  public :: changerep

contains

  !-------------------------------------------------------------------------------
  !> Summary: Change representation of matrix between (kappa,mue), real and complex spherical harmonics
  !> Author: 
  !> Category: KKRhost, special-functions
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Change the representation of matrix A and store in B          
  !> according to MODE:                                            
  !>                                                               
  !> RLM>REL   non-relat. REAL spher. harm.  >   (kappa,mue)       
  !> REL>RLM   (kappa,mue)  > non-relat. REAL spher. harm.         
  !> CLM>REL   non-relat. CMPLX. spher. harm.  >   (kappa,mue)     
  !> REL>CLM   (kappa,mue)  > non-relat. CMPLX. spher. harm.       
  !> RLM>CLM   non-relat. REAL spher. harm.  >  CMPLX. spher. harm.
  !> CLM>RLM   non-relat. CMPLX. spher. harm.  >  REAL spher. harm.
  !>                                                               
  !> the non-relat. representations include the  spin index        
  !>                                                               
  !> for LTEXT > 0 the new matrix  B  is printed 
  !-------------------------------------------------------------------------------
  subroutine changerep(a, mode, b, n, m, rc, crel, rrel, text, ltext)
    use :: mod_datatypes, only: dp
    use :: mod_cmatstr, only: cmatstr
    use :: mod_constants, only: cone, czero
    implicit none

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
      call zgemm('N', 'N', n, n, n, cone, rrel, m, a, m, czero, w1, m)
      call zgemm('N', 'C', n, n, n, cone, w1, m, rrel, m, czero, b, m)
      key = 2
    else if (mode=='RLM>REL') then
      call zgemm('C', 'N', n, n, n, cone, rrel, m, a, m, czero, w1, m)
      call zgemm('N', 'N', n, n, n, cone, w1, m, rrel, m, czero, b, m)
      key = 3
    else if (mode=='REL>CLM') then
      call zgemm('N', 'N', n, n, n, cone, crel, m, a, m, czero, w1, m)
      call zgemm('N', 'C', n, n, n, cone, w1, m, crel, m, czero, b, m)
      key = 2
    else if (mode=='CLM>REL') then
      call zgemm('C', 'N', n, n, n, cone, crel, m, a, m, czero, w1, m)
      call zgemm('N', 'N', n, n, n, cone, w1, m, crel, m, czero, b, m)
      key = 3
    else if (mode=='CLM>RLM') then
      call zgemm('N', 'N', n, n, n, cone, rc, m, a, m, czero, w1, m)
      call zgemm('N', 'C', n, n, n, cone, w1, m, rc, m, czero, b, m)
      key = 2
    else if (mode=='RLM>CLM') then
      call zgemm('C', 'N', n, n, n, cone, rc, m, a, m, czero, w1, m)
      call zgemm('N', 'N', n, n, n, cone, w1, m, rc, m, czero, b, m)
      key = 2
    else
      write (*, *) ' MODE = ', mode
      stop 'in <ROTATE>  MODE not allowed'
    end if

    if (ltext>0) call cmatstr(text, ltext, b, n, m, key, key, 0, 1e-8_dp, 6)
    ! IF ( LTEXT.GT.0 ) CALL CMATSTR(TEXT,LTEXT,B,N,M,KEY,KEY,0,1D-12,6)
  end subroutine changerep

end module mod_changerep
