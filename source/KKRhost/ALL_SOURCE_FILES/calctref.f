      SUBROUTINE CALCTREF13(ERYD,VREF,RMTREF,LMAX,LMTMAT,TREFLL,DTREFLL,
     &                    ALPHAREF,DALPHAREF,LMAXDP1,LMMAXD)
!   ********************************************************************
!   *                                                                  *
!   *  calculates analytically the single-site scattering matrix for a *
!   *  constant potential Vo at energy E                               *
!   *                                                                  *
!   *  input: potential radius RMTREF (R)                              *
!   *         energy ERYD (E)                                          *
!   *         potential constant value VREF (Vo)                       *
!   * output: single-site matrix TREFLL                                *
!   *         energy derivative DRTREFLL                               *
!   *                                                                  *
!   *               aj(l+1,aR)j(l,bR) - bj(l,aR)j(l+1,bR)              *
!   *   tmat(l) = - ------------------------------------               *
!   *               j(l,bR)h(l+1,aR)a - bh(l,aR)j(l+1,bR)              *
!   *                                                                  *
!   *   a = sqrt(E),  b = sqrt(E-Vo)                                   *
!   *                                                                  *
!     Derivative of t: the analytical formula for the derivative of
!     spherical Bessel functions is used:
!     (imported from KKRnano by Phivos Mavropoulos 10.10.2013)
!
!     d                     l+1
!     --  j (x) = j   (x) - --- j (x)   
!     dx   l       l-1       x   l
!
!     d
!     --  j (x) = - j (x)
!     dx   0         1
!
!     which for x = sqrt(E0)*r leads to
!
!      d          r*r             (l+1)
!     --- j (x) = --- ( j   (x) - ----- j (x) )
!     dE0  l      2 x    l-1        x    l
!
!      d            r*r
!     --- j (x) = - --- j (x)
!     dE0  0        2 x  1
!
!   ********************************************************************

      IMPLICIT NONE
!     ..
!     .. Scalar arguments
!     Input:
      DOUBLE COMPLEX ERYD                 ! energy
      INTEGER LMAX,LMAXDP1,LMMAXD
      DOUBLE PRECISION VREF,RMTREF        ! repulsive potential and its radius

!     Output:
      INTEGER LMTMAT
!     ..
!     .. Array arguments ..
!     Output:
      DOUBLE COMPLEX TREFLL(LMMAXD,LMMAXD)    ! t-matrix
      DOUBLE COMPLEX DTREFLL(LMMAXD,LMMAXD)   ! energy derivative of t-matrix             ! LLY Lloyd
      DOUBLE COMPLEX ALPHAREF(0:LMAXDP1-1),DALPHAREF(0:LMAXDP1-1) ! alpha matrix and derivative   ! LLY Lloyd
!     ..
!     .. Local scalars
      INTEGER J1,L1,LM1
      DOUBLE COMPLEX A1,B1,DA1,DB1,TMATANAL,DTMATANAL
      DOUBLE COMPLEX ROOTE1,ROOTE2
      DOUBLE COMPLEX CI,CIOVE
!     ..
!     .. Local arrays
      DOUBLE COMPLEX BESSJW1(0:LMAXDP1),BESSJW2(0:LMAXDP1), ! Bessel & Hankel
     &               BESSYW1(0:LMAXDP1),BESSYW2(0:LMAXDP1),
     &               HANKWS1(0:LMAXDP1),HANKWS2(0:LMAXDP1),
     &               DBESSJW1(0:LMAXDP1),DBESSJW2(0:LMAXDP1), ! Derivatives
     &               DHANKWS1(0:LMAXDP1)
!     ..
!     .. Intrinsic functions
      INTRINSIC SQRT
!     ..
!     .. External subroutines
      EXTERNAL BESHAN,CINIT
!     .. 
!     .. Data statement
      DATA CI / (0D0,1D0) /
!     ..
      LMTMAT = 0
      CALL CINIT(LMMAXD*LMMAXD,TREFLL)
      CALL CINIT(LMMAXD*LMMAXD,DTREFLL)

      ROOTE1 = SQRT(ERYD)
      ROOTE2 = SQRT(ERYD-VREF)

      A1 = ROOTE1*RMTREF
      B1 = ROOTE2*RMTREF
      CIOVE = CI/ROOTE1

!     Bessel functions and their derivatives:
      CALL BESHAN(HANKWS1,BESSJW1,BESSYW1,A1,LMAXDP1)
      CALL BESHAN(HANKWS2,BESSJW2,BESSYW2,B1,LMAXDP1)

      DBESSJW1(0) = - BESSJW1(1)/A1
      DBESSJW2(0) = - BESSJW2(1)/B1
      DHANKWS1(0) = - HANKWS1(1)/A1

      DO L1 = 1,LMAX + 1
         DBESSJW1(L1) = (BESSJW1(L1-1) - (L1+1)*BESSJW1(L1)/A1)/A1
         DBESSJW2(L1) = (BESSJW2(L1-1) - (L1+1)*BESSJW2(L1)/B1)/B1
         DHANKWS1(L1) = (HANKWS1(L1-1) - (L1+1)*HANKWS1(L1)/A1)/A1
      END DO

      DO L1 = 0,LMAX + 1
         DBESSJW1(L1) = 0.5D0*DBESSJW1(L1)*RMTREF**2 
         DBESSJW2(L1) = 0.5D0*DBESSJW2(L1)*RMTREF**2 
         DHANKWS1(L1) = 0.5D0*DHANKWS1(L1)*RMTREF**2 
      END DO



      DO L1 = 0,LMAX
         A1 = ROOTE1      * BESSJW1(L1+1) * BESSJW2(L1) - 
     &        ROOTE2 * BESSJW1(L1)   * BESSJW2(L1+1)

         B1 = ROOTE1      * HANKWS1(L1+1) * BESSJW2(L1) - 
     &        ROOTE2 * HANKWS1(L1)   * BESSJW2(L1+1)


         DA1 = 0.5D0/ROOTE1 * BESSJW1(L1+1) * BESSJW2(L1) -
     &         0.5D0/ROOTE2 * BESSJW1(L1)   * BESSJW2(L1+1) +
     &         ROOTE1 * DBESSJW1(L1+1) * BESSJW2(L1) -
     &         ROOTE2 * DBESSJW1(L1)   * BESSJW2(L1+1) +
     &         ROOTE1 * BESSJW1(L1+1)  * DBESSJW2(L1) -
     &         ROOTE2 * BESSJW1(L1)    * DBESSJW2(L1+1)

         DB1 = 0.5D0/ROOTE1 * HANKWS1(L1+1) * BESSJW2(L1) -
     &         0.5D0/ROOTE2 * HANKWS1(L1)   * BESSJW2(L1+1) +
     &         ROOTE1 * DHANKWS1(L1+1) * BESSJW2(L1) -
     &         ROOTE2 * DHANKWS1(L1)   * BESSJW2(L1+1) +
     &         ROOTE1 * HANKWS1(L1+1)  * DBESSJW2(L1) -
     &         ROOTE2 * HANKWS1(L1)    * DBESSJW2(L1+1)


         TMATANAL = -CIOVE*A1/B1
         DTMATANAL = CI * 0.5D0 / ROOTE1**3 * A1/B1 - 
     &               CI  / ROOTE1 * ( DA1/B1 - A1*DB1/B1**2 )

         ALPHAREF(L1) = -(ROOTE2/ROOTE1)**L1/(RMTREF**2 * ROOTE1 * B1) ! Following R. Zeller
         DALPHAREF(L1)=( L1/2.D0/ROOTE2 - (L1+1)/2.D0/ROOTE1 - DB1/B1 )* 
     &                   ALPHAREF(L1)

         DO J1 = -L1,L1
            LM1 = L1*(L1+1) + J1 + 1
            TREFLL(LM1,LM1) = TMATANAL
            DTREFLL(LM1,LM1) = DTMATANAL
         END DO
      END DO
      LMTMAT = LM1

      END

