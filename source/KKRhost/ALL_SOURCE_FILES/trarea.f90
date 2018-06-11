    Subroutine trarea(a, b, lmax)
      Use mod_datatypes, Only: dp
! from complex to real  (differenciated spherical harmonics)

!.. Parameters ..
      Real (Kind=dp) :: rtwo
      Complex (Kind=dp) :: ci
      Parameter (rtwo=1.414213562373E0_dp, ci=(0.E0_dp,1.E0_dp))
!..
!.. Scalar Arguments ..
      Integer :: lmax
!..
!.. Array Arguments ..
      Complex (Kind=dp) :: a(*)
      Real (Kind=dp) :: b(*)
!..
!.. Local Scalars ..
      Real (Kind=dp) :: sgm
      Integer :: i, l, m
!..
!.. Intrinsic Functions ..
      Intrinsic :: conjg, real

!    calculate real the spherical harmonics derivetived
      i = 0
      Do l = 0, lmax
        i = i + l + 1
        b(i) = real(a(i))
        sgm = -1.E0_dp
        Do m = 1, l
          b(i-m) = real(ci*(a(i-m)-conjg(a(i-m))))/rtwo
          b(i+m) = sgm*real((a(i+m)+conjg(a(i+m))))/rtwo
          sgm = -sgm
        End Do
        i = i + l
      End Do
      Return
    End Subroutine
