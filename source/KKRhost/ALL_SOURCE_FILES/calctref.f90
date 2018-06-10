SUBROUTINE calctref13(eryd,vref,rmtref,lmax,lmtmat,trefll,dtrefll,  &
        alpharef,dalpharef,lmaxdp1,lmmaxd)
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
!..
!.. Scalar arguments
!Input:
DOUBLE COMPLEX ERYD                 ! energy
INTEGER LMAX,LMAXDP1,LMMAXD
DOUBLE PRECISION VREF,RMTREF        ! repulsive potential and its radius

!Output:
INTEGER LMTMAT
!..
!.. Array arguments ..
!Output:
DOUBLE COMPLEX TREFLL(LMMAXD,LMMAXD)    ! t-matrix
DOUBLE COMPLEX DTREFLL(LMMAXD,LMMAXD)   ! energy derivative of t-matrix             ! LLY Lloyd
DOUBLE COMPLEX ALPHAREF(0:LMAXDP1-1),DALPHAREF(0:LMAXDP1-1) ! alpha matrix and derivative   ! LLY Lloyd
!..
!.. Local scalars
INTEGER J1,L1,LM1
DOUBLE COMPLEX A1,B1,DA1,DB1,TMATANAL,DTMATANAL
DOUBLE COMPLEX ROOTE1,ROOTE2
DOUBLE COMPLEX CI,CIOVE
!..
!.. Local arrays
DOUBLE COMPLEX BESSJW1(0:LMAXDP1),BESSJW2(0:LMAXDP1), &! Bessel & Hankel
               BESSYW1(0:LMAXDP1),BESSYW2(0:LMAXDP1), &
               HANKWS1(0:LMAXDP1),HANKWS2(0:LMAXDP1), &
               DBESSJW1(0:LMAXDP1),DBESSJW2(0:LMAXDP1), & ! Derivatives
               DHANKWS1(0:LMAXDP1)
!..
!.. Intrinsic functions
INTRINSIC SQRT
!..
!.. External subroutines
EXTERNAL BESHAN,CINIT
!.. 
!.. Data statement
DATA CI / (0D0,1D0) /
!..
lmtmat = 0
CALL cinit(lmmaxd*lmmaxd,trefll)
CALL cinit(lmmaxd*lmmaxd,dtrefll)

roote1 = SQRT(eryd)
roote2 = SQRT(eryd-vref)

a1 = roote1*rmtref
b1 = roote2*rmtref
ciove = ci/roote1

!     Bessel functions and their derivatives:
CALL beshan(hankws1,bessjw1,bessyw1,a1,lmaxdp1)
CALL beshan(hankws2,bessjw2,bessyw2,b1,lmaxdp1)

dbessjw1(0) = - bessjw1(1)/a1
dbessjw2(0) = - bessjw2(1)/b1
dhankws1(0) = - hankws1(1)/a1

DO l1 = 1,lmax + 1
  dbessjw1(l1) = (bessjw1(l1-1) - (l1+1)*bessjw1(l1)/a1)/a1
  dbessjw2(l1) = (bessjw2(l1-1) - (l1+1)*bessjw2(l1)/b1)/b1
  dhankws1(l1) = (hankws1(l1-1) - (l1+1)*hankws1(l1)/a1)/a1
END DO

DO l1 = 0,lmax + 1
  dbessjw1(l1) = 0.5D0*dbessjw1(l1)*rmtref**2
  dbessjw2(l1) = 0.5D0*dbessjw2(l1)*rmtref**2
  dhankws1(l1) = 0.5D0*dhankws1(l1)*rmtref**2
END DO



DO l1 = 0,lmax
  a1 = roote1      * bessjw1(l1+1) * bessjw2(l1) -  &
      roote2 * bessjw1(l1)   * bessjw2(l1+1)
  
  b1 = roote1      * hankws1(l1+1) * bessjw2(l1) -  &
      roote2 * hankws1(l1)   * bessjw2(l1+1)
  
  
  da1 = 0.5D0/roote1 * bessjw1(l1+1) * bessjw2(l1) -  &
      0.5D0/roote2 * bessjw1(l1)   * bessjw2(l1+1) +  &
      roote1 * dbessjw1(l1+1) * bessjw2(l1) -  &
      roote2 * dbessjw1(l1)   * bessjw2(l1+1) +  &
      roote1 * bessjw1(l1+1)  * dbessjw2(l1) -  &
      roote2 * bessjw1(l1)    * dbessjw2(l1+1)
  
  db1 = 0.5D0/roote1 * hankws1(l1+1) * bessjw2(l1) -  &
      0.5D0/roote2 * hankws1(l1)   * bessjw2(l1+1) +  &
      roote1 * dhankws1(l1+1) * bessjw2(l1) -  &
      roote2 * dhankws1(l1)   * bessjw2(l1+1) +  &
      roote1 * hankws1(l1+1)  * dbessjw2(l1) -  &
      roote2 * hankws1(l1)    * dbessjw2(l1+1)
  
  
  tmatanal = -ciove*a1/b1
  dtmatanal = ci * 0.5D0 / roote1**3 * a1/b1 -  &
      ci  / roote1 * ( da1/b1 - a1*db1/b1**2 )
  
  alpharef(l1) = -(roote2/roote1)**l1/(rmtref**2 * roote1 * b1) ! Following R. Zeller
  dalpharef(l1)=( l1/2.d0/roote2 - (l1+1)/2.d0/roote1 - db1/b1 )* alpharef(l1)
  
  DO j1 = -l1,l1
    lm1 = l1*(l1+1) + j1 + 1
    trefll(lm1,lm1) = tmatanal
    dtrefll(lm1,lm1) = dtmatanal
  END DO
END DO
lmtmat = lm1

END SUBROUTINE calctref13

