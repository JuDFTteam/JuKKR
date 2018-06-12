subroutine trarea(a, b, lmax)
  use :: mod_datatypes, only: dp
  ! from complex to real  (differenciated spherical harmonics)

  ! .. Parameters ..
  real (kind=dp) :: rtwo
  complex (kind=dp) :: ci
  parameter (rtwo=1.414213562373e0_dp, ci=(0.e0_dp,1.e0_dp))
  ! ..
  ! .. Scalar Arguments ..
  integer :: lmax
  ! ..
  ! .. Array Arguments ..
  complex (kind=dp) :: a(*)
  real (kind=dp) :: b(*)
  ! ..
  ! .. Local Scalars ..
  real (kind=dp) :: sgm
  integer :: i, l, m
  ! ..
  ! .. Intrinsic Functions ..
  intrinsic :: conjg, real

  ! calculate real the spherical harmonics derivetived
  i = 0
  do l = 0, lmax
    i = i + l + 1
    b(i) = real(a(i))
    sgm = -1.e0_dp
    do m = 1, l
      b(i-m) = real(ci*(a(i-m)-conjg(a(i-m))))/rtwo
      b(i+m) = sgm*real((a(i+m)+conjg(a(i+m))))/rtwo
      sgm = -sgm
    end do
    i = i + l
  end do
  return
end subroutine trarea
