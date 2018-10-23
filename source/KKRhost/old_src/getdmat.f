C*==getdmat.f    processed by SPAG 6.05Rc at 20:01 on 11 Nov 2000
      SUBROUTINE GETDMAT(TAUQ,DMATT,DTILT,DM,N,MSSQ,MSST,M)
C   ********************************************************************
C   *                                                                  *
C   *   calculate projection matrices   DMATT  and  DTILT              *
C   *   for preselected atom type  IT  on preselected site  IQ         *
C   *                                                                  *
C   *    DM = m(t)-m(c)                                                *
C   *                                                                  *
C   *    D~(t) = ( 1 + ( m(t) - m(c) ) * TAU )**(-1)                   *
C   *    D(t)  = ( 1 + TAU * ( m(t) - m(c) ) )**(-1)                   *
C   *                                                                  *
C   *   on entry ALL matrices have to refer to the SAME frame          *
C   *   i.e. for KMROT <> 0 prior to calling <GETDMAT> one has to      *
C   *   - rotate TAUQ, MSSQ to the local frame        OR               *
C   *   - rotate MSST to the global frame                              *
C   *                                                                  *
C   * 01/11/2000 HE                                                    *
C   ********************************************************************
C
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
      INTEGER M,N
      COMPLEX*16 DM(M,M),DMATT(M,M),DTILT(M,M),MSSQ(M,M),MSST(M,M),
     &           TAUQ(M,M)
C
C Local variables
C
      INTEGER I,J,INFO,IPIV(M)
      COMPLEX*16 MAUX(M,M)
C
C*** End of declarations rewritten by SPAG
C
      DO J=1,N
         DO I=1,N
            DM(I,J) = MSST(I,J) - MSSQ(I,J)
         END DO
      END DO
C
C     -------------------------------------------
C                   ( m(t) - m(c) ) * TAU
C     -------------------------------------------
      CALL ZGEMM('N','N',N,N,N,C1,DM,M,TAUQ,M,C0,DTILT,M)
      CALL ZGEMM('N','N',N,N,N,C1,TAUQ,M,DM,M,C0,DMATT,M)
C
C     -------------------------------------------
C               1 + ( m(t) - m(c) ) * TAU
C     -------------------------------------------
      DO I = 1,N
         DTILT(I,I) = C1 + DTILT(I,I)
         DMATT(I,I) = C1 + DMATT(I,I)
      END DO
C
C     -------------------------------------------
C     D~(t) = ( 1 + ( m(t) - m(c) ) * TAU )**(-1)
C     -------------------------------------------
      CALL ZGETRF(N,N,DTILT,M,IPIV,INFO)
      CALL ZGETRI(N,DTILT,M,IPIV,MAUX,M*M,INFO)
C
      CALL ZGETRF(N,N,DMATT,M,IPIV,INFO)
      CALL ZGETRI(N,DMATT,M,IPIV,MAUX,M*M,INFO)
C
      END
