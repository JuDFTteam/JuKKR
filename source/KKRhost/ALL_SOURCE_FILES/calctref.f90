    Subroutine calctref13(eryd, vref, rmtref, lmax, lmtmat, trefll, dtrefll, &
      alpharef, dalpharef, lmaxdp1, lmmaxd)
      Use mod_datatypes, Only: dp
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

      Implicit None
!..
!.. Scalar arguments
!Input:
      Complex (Kind=dp) :: eryd ! energy
      Integer :: lmax, lmaxdp1, lmmaxd
      Real (Kind=dp) :: vref, rmtref ! repulsive potential and its radius

!Output:
      Integer :: lmtmat
!..
!.. Array arguments ..
!Output:
      Complex (Kind=dp) :: trefll(lmmaxd, lmmaxd) ! t-matrix
      Complex (Kind=dp) :: dtrefll(lmmaxd, lmmaxd) ! energy derivative of t-matrix             ! LLY Lloyd
      Complex (Kind=dp) :: alpharef(0:lmaxdp1-1), dalpharef(0:lmaxdp1-1) ! alpha matrix and derivative   ! LLY Lloyd
!..
!.. Local scalars
      Integer :: j1, l1, lm1
      Complex (Kind=dp) :: a1, b1, da1, db1, tmatanal, dtmatanal
      Complex (Kind=dp) :: roote1, roote2
      Complex (Kind=dp) :: ci, ciove
!..
!.. Local arrays
! Bessel & Hankel
! Derivatives
      Complex (Kind=dp) :: bessjw1(0:lmaxdp1), bessjw2(0:lmaxdp1), &
        bessyw1(0:lmaxdp1), bessyw2(0:lmaxdp1), hankws1(0:lmaxdp1), &
        hankws2(0:lmaxdp1), dbessjw1(0:lmaxdp1), dbessjw2(0:lmaxdp1), &
        dhankws1(0:lmaxdp1)
!..
!.. Intrinsic functions
      Intrinsic :: sqrt
!..
!.. External subroutines
      External :: beshan, cinit
!.. 
!.. Data statement
      Data ci/(0E0_dp, 1E0_dp)/
!..
      lmtmat = 0
      Call cinit(lmmaxd*lmmaxd, trefll)
      Call cinit(lmmaxd*lmmaxd, dtrefll)

      roote1 = sqrt(eryd)
      roote2 = sqrt(eryd-vref)

      a1 = roote1*rmtref
      b1 = roote2*rmtref
      ciove = ci/roote1

!     Bessel functions and their derivatives:
      Call beshan(hankws1, bessjw1, bessyw1, a1, lmaxdp1)
      Call beshan(hankws2, bessjw2, bessyw2, b1, lmaxdp1)

      dbessjw1(0) = -bessjw1(1)/a1
      dbessjw2(0) = -bessjw2(1)/b1
      dhankws1(0) = -hankws1(1)/a1

      Do l1 = 1, lmax + 1
        dbessjw1(l1) = (bessjw1(l1-1)-(l1+1)*bessjw1(l1)/a1)/a1
        dbessjw2(l1) = (bessjw2(l1-1)-(l1+1)*bessjw2(l1)/b1)/b1
        dhankws1(l1) = (hankws1(l1-1)-(l1+1)*hankws1(l1)/a1)/a1
      End Do

      Do l1 = 0, lmax + 1
        dbessjw1(l1) = 0.5E0_dp*dbessjw1(l1)*rmtref**2
        dbessjw2(l1) = 0.5E0_dp*dbessjw2(l1)*rmtref**2
        dhankws1(l1) = 0.5E0_dp*dhankws1(l1)*rmtref**2
      End Do



      Do l1 = 0, lmax
        a1 = roote1*bessjw1(l1+1)*bessjw2(l1) - roote2*bessjw1(l1)*bessjw2(l1+ &
          1)

        b1 = roote1*hankws1(l1+1)*bessjw2(l1) - roote2*hankws1(l1)*bessjw2(l1+ &
          1)


        da1 = 0.5E0_dp/roote1*bessjw1(l1+1)*bessjw2(l1) - &
          0.5E0_dp/roote2*bessjw1(l1)*bessjw2(l1+1) + &
          roote1*dbessjw1(l1+1)*bessjw2(l1) - roote2*dbessjw1(l1)*bessjw2(l1+1 &
          ) + roote1*bessjw1(l1+1)*dbessjw2(l1) - roote2*bessjw1(l1)*dbessjw2( &
          l1+1)

        db1 = 0.5E0_dp/roote1*hankws1(l1+1)*bessjw2(l1) - &
          0.5E0_dp/roote2*hankws1(l1)*bessjw2(l1+1) + &
          roote1*dhankws1(l1+1)*bessjw2(l1) - roote2*dhankws1(l1)*bessjw2(l1+1 &
          ) + roote1*hankws1(l1+1)*dbessjw2(l1) - roote2*hankws1(l1)*dbessjw2( &
          l1+1)


        tmatanal = -ciove*a1/b1
        dtmatanal = ci*0.5E0_dp/roote1**3*a1/b1 - ci/roote1*(da1/b1-a1*db1/b1 &
          **2)

        alpharef(l1) = -(roote2/roote1)**l1/(rmtref**2*roote1*b1) ! Following R. Zeller
        dalpharef(l1) = (l1/2.E0_dp/roote2-(l1+1)/2.E0_dp/roote1-db1/b1)* &
          alpharef(l1)

        Do j1 = -l1, l1
          lm1 = l1*(l1+1) + j1 + 1
          trefll(lm1, lm1) = tmatanal
          dtrefll(lm1, lm1) = dtmatanal
        End Do
      End Do
      lmtmat = lm1

    End Subroutine
