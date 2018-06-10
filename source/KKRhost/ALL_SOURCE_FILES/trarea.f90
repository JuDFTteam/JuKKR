subroutine trarea(a, b, lmax)
! from complex to real  (differenciated spherical harmonics)

!.. Parameters ..
  double precision :: rtwo
  double complex :: ci
  parameter (rtwo=1.414213562373d0, ci=(0.d0,1.d0))
!..
!.. Scalar Arguments ..
  integer :: lmax
!..
!.. Array Arguments ..
  double complex :: a(*)
  double precision :: b(*)
!..
!.. Local Scalars ..
  double precision :: sgm
  integer :: i, l, m
!..
!.. Intrinsic Functions ..
  intrinsic :: conjg, dble

!    calculate real the spherical harmonics derivetived
  i = 0
  do l = 0, lmax
    i = i + l + 1
    b(i) = dble(a(i))
    sgm = -1.d0
    do m = 1, l
      b(i-m) = dble(ci*(a(i-m)-conjg(a(i-m))))/rtwo
      b(i+m) = sgm*dble((a(i+m)+conjg(a(i+m))))/rtwo
      sgm = -sgm
    end do
    i = i + l
  end do
  return
end subroutine
