module mod_getdmat
  
  private
  public :: getdmat

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculate projection matrices DMATT and DTILT
  !> Author: Hubert Ebert
  !> Date: 01/11/2000
  !> Category: KKRhost, k-points, coherent-potential-approx
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Calculate projection matrices   DMATT  and  DTILT
  !> for preselected atom type  IT  on preselected site  IQ   
  !>
  !>  DM = m(t)-m(c)
  !>
  !>  D~(t) = ( 1 + ( m(t) - m(c) ) * TAU )**(-1)
  !>  D(t)  = ( 1 + TAU * ( m(t) - m(c) ) )**(-1)
  !>
  !> on entry ALL matrices have to refer to the SAME frame
  !> i.e. for KMROT <> 0 prior to calling <GETDMAT> one has to
  !> - rotate TAUQ, MSSQ to the local frame        OR
  !> - rotate MSST to the global frame
  !-------------------------------------------------------------------------------
  subroutine getdmat(tauq, dmatt, dtilt, dm, n, mssq, msst, m)

    use :: mod_datatypes, only: dp
    use :: mod_constants, only: czero, cone
    implicit none

    ! Dummy arguments
    integer :: m, n
    complex (kind=dp) :: dm(m, m), dmatt(m, m), dtilt(m, m), mssq(m, m), msst(m, m), tauq(m, m)

    ! Local variables
    integer :: i, j, info, ipiv(m)
    complex (kind=dp) :: maux(m, m)

    do j = 1, n
      do i = 1, n
        dm(i, j) = msst(i, j) - mssq(i, j)
      end do
    end do

    ! -------------------------------------------
    ! ( m(t) - m(c) ) * TAU
    ! -------------------------------------------
    call zgemm('N', 'N', n, n, n, cone, dm, m, tauq, m, czero, dtilt, m)
    call zgemm('N', 'N', n, n, n, cone, tauq, m, dm, m, czero, dmatt, m)

    ! -------------------------------------------
    ! 1 + ( m(t) - m(c) ) * TAU
    ! -------------------------------------------
    do i = 1, n
      dtilt(i, i) = cone + dtilt(i, i)
      dmatt(i, i) = cone + dmatt(i, i)
    end do

    ! -------------------------------------------
    ! D~(t) = ( 1 + ( m(t) - m(c) ) * TAU )**(-1)
    ! -------------------------------------------
    call zgetrf(n, n, dtilt, m, ipiv, info)
    call zgetri(n, dtilt, m, ipiv, maux, m*m, info)

    call zgetrf(n, n, dmatt, m, ipiv, info)
    call zgetri(n, dmatt, m, ipiv, maux, m*m, info)

  end subroutine getdmat

end module mod_getdmat
