subroutine getdmat(tauq, dmatt, dtilt, dm, n, mssq, msst, m)
!   ********************************************************************
!   *                                                                  *
!   *   calculate projection matrices   DMATT  and  DTILT              *
!   *   for preselected atom type  IT  on preselected site  IQ         *
!   *                                                                  *
!   *    DM = m(t)-m(c)                                                *
!   *                                                                  *
!   *    D~(t) = ( 1 + ( m(t) - m(c) ) * TAU )**(-1)                   *
!   *    D(t)  = ( 1 + TAU * ( m(t) - m(c) ) )**(-1)                   *
!   *                                                                  *
!   *   on entry ALL matrices have to refer to the SAME frame          *
!   *   i.e. for KMROT <> 0 prior to calling <GETDMAT> one has to      *
!   *   - rotate TAUQ, MSSQ to the local frame        OR               *
!   *   - rotate MSST to the global frame                              *
!   *                                                                  *
!   * 01/11/2000 HE                                                    *
!   ********************************************************************

  implicit none

! PARAMETER definitions
  complex *16 :: c0, c1
  parameter (c0=(0.0d0,0.0d0), c1=(1.0d0,0.0d0))

! Dummy arguments
  integer :: m, n
  complex *16 :: dm(m, m), dmatt(m, m), dtilt(m, m), mssq(m, m), msst(m, m), &
    tauq(m, m)

! Local variables
  integer :: i, j, info, ipiv(m)
  complex *16 :: maux(m, m)

  do j = 1, n
    do i = 1, n
      dm(i, j) = msst(i, j) - mssq(i, j)
    end do
  end do

!     -------------------------------------------
!                   ( m(t) - m(c) ) * TAU
!     -------------------------------------------
  call zgemm('N', 'N', n, n, n, c1, dm, m, tauq, m, c0, dtilt, m)
  call zgemm('N', 'N', n, n, n, c1, tauq, m, dm, m, c0, dmatt, m)

!     -------------------------------------------
!               1 + ( m(t) - m(c) ) * TAU
!     -------------------------------------------
  do i = 1, n
    dtilt(i, i) = c1 + dtilt(i, i)
    dmatt(i, i) = c1 + dmatt(i, i)
  end do

!     -------------------------------------------
!     D~(t) = ( 1 + ( m(t) - m(c) ) * TAU )**(-1)
!     -------------------------------------------
  call zgetrf(n, n, dtilt, m, ipiv, info)
  call zgetri(n, dtilt, m, ipiv, maux, m*m, info)

  call zgetrf(n, n, dmatt, m, ipiv, info)
  call zgetri(n, dmatt, m, ipiv, maux, m*m, info)

end subroutine
