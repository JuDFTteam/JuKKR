    Subroutine getdmat(tauq, dmatt, dtilt, dm, n, mssq, msst, m)
      Use mod_datatypes, Only: dp
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

      Implicit None

! PARAMETER definitions
      Complex (Kind=dp) :: c0, c1
      Parameter (c0=(0.0E0_dp,0.0E0_dp), c1=(1.0E0_dp,0.0E0_dp))

! Dummy arguments
      Integer :: m, n
      Complex (Kind=dp) :: dm(m, m), dmatt(m, m), dtilt(m, m), mssq(m, m), &
        msst(m, m), tauq(m, m)

! Local variables
      Integer :: i, j, info, ipiv(m)
      Complex (Kind=dp) :: maux(m, m)

      Do j = 1, n
        Do i = 1, n
          dm(i, j) = msst(i, j) - mssq(i, j)
        End Do
      End Do

!     -------------------------------------------
!                   ( m(t) - m(c) ) * TAU
!     -------------------------------------------
      Call zgemm('N', 'N', n, n, n, c1, dm, m, tauq, m, c0, dtilt, m)
      Call zgemm('N', 'N', n, n, n, c1, tauq, m, dm, m, c0, dmatt, m)

!     -------------------------------------------
!               1 + ( m(t) - m(c) ) * TAU
!     -------------------------------------------
      Do i = 1, n
        dtilt(i, i) = c1 + dtilt(i, i)
        dmatt(i, i) = c1 + dmatt(i, i)
      End Do

!     -------------------------------------------
!     D~(t) = ( 1 + ( m(t) - m(c) ) * TAU )**(-1)
!     -------------------------------------------
      Call zgetrf(n, n, dtilt, m, ipiv, info)
      Call zgetri(n, dtilt, m, ipiv, maux, m*m, info)

      Call zgetrf(n, n, dmatt, m, ipiv, info)
      Call zgetri(n, dmatt, m, ipiv, maux, m*m, info)

    End Subroutine
