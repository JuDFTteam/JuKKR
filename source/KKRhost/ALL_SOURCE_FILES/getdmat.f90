SUBROUTINE getdmat(tauq,dmatt,dtilt,dm,n,mssq,msst,m)
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

IMPLICIT NONE

! PARAMETER definitions
COMPLEX*16 C0,C1
PARAMETER (C0=(0.0D0,0.0D0),C1=(1.0D0,0.0D0))

! Dummy arguments
INTEGER M,N
COMPLEX*16 DM(M,M),DMATT(M,M),DTILT(M,M),MSSQ(M,M),MSST(M,M), &
           TAUQ(M,M)

! Local variables
INTEGER I,J,INFO,IPIV(M)
COMPLEX*16 MAUX(M,M)

DO j=1,n
  DO i=1,n
    dm(i,j) = msst(i,j) - mssq(i,j)
  END DO
END DO

!     -------------------------------------------
!                   ( m(t) - m(c) ) * TAU
!     -------------------------------------------
CALL zgemm('N','N',n,n,n,c1,dm,m,tauq,m,c0,dtilt,m)
CALL zgemm('N','N',n,n,n,c1,tauq,m,dm,m,c0,dmatt,m)

!     -------------------------------------------
!               1 + ( m(t) - m(c) ) * TAU
!     -------------------------------------------
DO i = 1,n
  dtilt(i,i) = c1 + dtilt(i,i)
  dmatt(i,i) = c1 + dmatt(i,i)
END DO

!     -------------------------------------------
!     D~(t) = ( 1 + ( m(t) - m(c) ) * TAU )**(-1)
!     -------------------------------------------
CALL zgetrf(n,n,dtilt,m,ipiv,info)
CALL zgetri(n,dtilt,m,ipiv,maux,m*m,info)

CALL zgetrf(n,n,dmatt,m,ipiv,info)
CALL zgetri(n,dmatt,m,ipiv,maux,m*m,info)

END SUBROUTINE getdmat
