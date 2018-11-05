      MODULE mod_calctref
!-------------------------------------------------------------------------------
!> Summary: Calculate t-matrices for the reference system 
!> Author:
!> Category: KKRimp, single-site, reference-system 
!>           
!-------------------------------------------------------------------------------
      CONTAINS
!-------------------------------------------------------------------------------
!> Summary: Calculate t-matrices for the reference system 
!> Author:
!> Category: KKRimp, single-site, reference-system 
!>           
!-------------------------------------------------------------------------------
      SUBROUTINE CALCTREF(ERYD,VREF,RMTREF,LMAX,LMTMAT,TREFLL,
     &                    LMAXDP1,LMMAXD)
      use mod_beshan
C   ********************************************************************
C   *                                                                  *
C   *  calculates analitically the single-site scattering matrix for a *
C   *  constant potential Vo at energy E                               *
C   *                                                                  *
C   *  input: potential radius RMTREF (R)                              *
C   *         energy ERYD (E)                                          *
C   *         potential constant value VREF (Vo)                       *
C   * output: single-site matrix TREFLL                                *
C   *                                                                  *
C   *   tmat(l) =   aj(l+1,aR)j(l,bR) - bj(l,aR)j(l+1,bR)              *
C   *             - ------------------------------------               *
C   *               j(l,bR)h(l+1,aR)a - bh(l,aR)j(l+1,bR)              *
C   *                                                                  *
C   *   a = sqrt(E),  b = sqrt(E-Vo)                                   *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C     ..
C     .. Scalar arguments
      DOUBLE COMPLEX ERYD
      INTEGER LMAX,LMAXDP1,LMMAXD,LMTMAT
      DOUBLE PRECISION RMTREF,VREF
C     ..
C     .. Array arguments ..
      DOUBLE COMPLEX TREFLL(LMMAXD,LMMAXD)
C     ..
C     .. Local scalars
      INTEGER J1,L1,LM1
      DOUBLE COMPLEX A1,B1,TMATANAL
      DOUBLE COMPLEX,parameter :: CI=(0D0,1D0)
      DOUBLE COMPLEX CIOVE
      DOUBLE COMPLEX,parameter :: CZERO=(0D0,0D0)

C     ..
C     .. Local arrays
      DOUBLE COMPLEX BESSJW1(0:LMAXDP1),BESSJW2(0:LMAXDP1),
     &               BESSYW1(0:LMAXDP1),BESSYW2(0:LMAXDP1),
     &               HANKWS1(0:LMAXDP1),HANKWS2(0:LMAXDP1)
C     ..
C     .. Intrinsic functions
      INTRINSIC SQRT
C     ..
C     .. External subroutines
!       EXTERNAL BESHAN,CINIT
C     .. 
C     .. Data statement
!       DATA CI /  /
C     ..
      LMTMAT = 0
      TREFLL=CZERO
C
      A1 = SQRT(ERYD)*RMTREF
      B1 = SQRT(ERYD-VREF)*RMTREF
      CIOVE = CI/SQRT(ERYD)
C
      CALL BESHAN(HANKWS1,BESSJW1,BESSYW1,A1,LMAXDP1)
      CALL BESHAN(HANKWS2,BESSJW2,BESSYW2,B1,LMAXDP1)
C
      DO L1 = 0,LMAX
         A1 = SQRT(ERYD)*BESSJW1(L1+1)*BESSJW2(L1) - SQRT(ERYD-VREF)
     &        *BESSJW1(L1)*BESSJW2(L1+1)
C
         B1 = SQRT(ERYD)*HANKWS1(L1+1)*BESSJW2(L1) - SQRT(ERYD-VREF)
     &        *HANKWS1(L1)*BESSJW2(L1+1)
C
         TMATANAL = -CIOVE*A1/B1
         DO J1 = -L1,L1
            LM1 = L1*(L1+1) + J1 + 1
            TREFLL(LM1,LM1) = TMATANAL
         END DO
      END DO
      LMTMAT = LM1
C
      END SUBROUTINE
      END MODULE mod_calctref
