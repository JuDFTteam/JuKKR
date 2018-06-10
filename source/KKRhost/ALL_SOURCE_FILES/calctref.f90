subroutine calctref13(eryd, vref, rmtref, lmax, lmtmat, trefll, dtrefll, &
  alpharef, dalpharef, lmaxdp1, lmmaxd)
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

  implicit none
!..
!.. Scalar arguments
!Input:
  double complex :: eryd ! energy
  integer :: lmax, lmaxdp1, lmmaxd
  double precision :: vref, rmtref ! repulsive potential and its radius

!Output:
  integer :: lmtmat
!..
!.. Array arguments ..
!Output:
  double complex :: trefll(lmmaxd, lmmaxd) ! t-matrix
  double complex :: dtrefll(lmmaxd, lmmaxd) ! energy derivative of t-matrix             ! LLY Lloyd
  double complex :: alpharef(0:lmaxdp1-1), dalpharef(0:lmaxdp1-1) ! alpha matrix and derivative   ! LLY Lloyd
!..
!.. Local scalars
  integer :: j1, l1, lm1
  double complex :: a1, b1, da1, db1, tmatanal, dtmatanal
  double complex :: roote1, roote2
  double complex :: ci, ciove
!..
!.. Local arrays
  double complex :: bessjw1(0:lmaxdp1), bessjw2(0:lmaxdp1), & ! Bessel & Hankel
    bessyw1(0:lmaxdp1), bessyw2(0:lmaxdp1), hankws1(0:lmaxdp1), &
    hankws2(0:lmaxdp1), dbessjw1(0:lmaxdp1), dbessjw2(0:lmaxdp1), & ! Derivatives
    dhankws1(0:lmaxdp1)
!..
!.. Intrinsic functions
  intrinsic :: sqrt
!..
!.. External subroutines
  external :: beshan, cinit
!.. 
!.. Data statement
  data ci/(0d0, 1d0)/
!..
  lmtmat = 0
  call cinit(lmmaxd*lmmaxd, trefll)
  call cinit(lmmaxd*lmmaxd, dtrefll)

  roote1 = sqrt(eryd)
  roote2 = sqrt(eryd-vref)

  a1 = roote1*rmtref
  b1 = roote2*rmtref
  ciove = ci/roote1

!     Bessel functions and their derivatives:
  call beshan(hankws1, bessjw1, bessyw1, a1, lmaxdp1)
  call beshan(hankws2, bessjw2, bessyw2, b1, lmaxdp1)

  dbessjw1(0) = -bessjw1(1)/a1
  dbessjw2(0) = -bessjw2(1)/b1
  dhankws1(0) = -hankws1(1)/a1

  do l1 = 1, lmax + 1
    dbessjw1(l1) = (bessjw1(l1-1)-(l1+1)*bessjw1(l1)/a1)/a1
    dbessjw2(l1) = (bessjw2(l1-1)-(l1+1)*bessjw2(l1)/b1)/b1
    dhankws1(l1) = (hankws1(l1-1)-(l1+1)*hankws1(l1)/a1)/a1
  end do

  do l1 = 0, lmax + 1
    dbessjw1(l1) = 0.5d0*dbessjw1(l1)*rmtref**2
    dbessjw2(l1) = 0.5d0*dbessjw2(l1)*rmtref**2
    dhankws1(l1) = 0.5d0*dhankws1(l1)*rmtref**2
  end do



  do l1 = 0, lmax
    a1 = roote1*bessjw1(l1+1)*bessjw2(l1) - roote2*bessjw1(l1)*bessjw2(l1+1)

    b1 = roote1*hankws1(l1+1)*bessjw2(l1) - roote2*hankws1(l1)*bessjw2(l1+1)


    da1 = 0.5d0/roote1*bessjw1(l1+1)*bessjw2(l1) - 0.5d0/roote2*bessjw1(l1)* &
      bessjw2(l1+1) + roote1*dbessjw1(l1+1)*bessjw2(l1) - &
      roote2*dbessjw1(l1)*bessjw2(l1+1) + roote1*bessjw1(l1+1)*dbessjw2(l1) - &
      roote2*bessjw1(l1)*dbessjw2(l1+1)

    db1 = 0.5d0/roote1*hankws1(l1+1)*bessjw2(l1) - 0.5d0/roote2*hankws1(l1)* &
      bessjw2(l1+1) + roote1*dhankws1(l1+1)*bessjw2(l1) - &
      roote2*dhankws1(l1)*bessjw2(l1+1) + roote1*hankws1(l1+1)*dbessjw2(l1) - &
      roote2*hankws1(l1)*dbessjw2(l1+1)


    tmatanal = -ciove*a1/b1
    dtmatanal = ci*0.5d0/roote1**3*a1/b1 - ci/roote1*(da1/b1-a1*db1/b1**2)

    alpharef(l1) = -(roote2/roote1)**l1/(rmtref**2*roote1*b1) ! Following R. Zeller
    dalpharef(l1) = (l1/2.d0/roote2-(l1+1)/2.d0/roote1-db1/b1)*alpharef(l1)

    do j1 = -l1, l1
      lm1 = l1*(l1+1) + j1 + 1
      trefll(lm1, lm1) = tmatanal
      dtrefll(lm1, lm1) = dtmatanal
    end do
  end do
  lmtmat = lm1

end subroutine

