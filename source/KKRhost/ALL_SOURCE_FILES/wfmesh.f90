subroutine wfmesh(e, ek, cvlight, nsra, z, r, s, rs, irm, irmd, lmaxd)
  use :: mod_datatypes, only: dp
  real (kind=dp), parameter :: eps=1.0D-12
  ! .. Scalar Arguments ..
  complex (kind=dp) :: e, ek
  real (kind=dp) :: cvlight, z
  integer :: irm, irmd, lmaxd, nsra
  ! ..
  ! .. Intrinsic Functions ..
  intrinsic :: real, sqrt
  ! ..
  ! .. Array Arguments ..
  real (kind=dp) :: r(irmd), rs(irmd, 0:lmaxd), s(0:lmaxd)
  ! ..
  ! .. Local Scalars ..
  real (kind=dp) :: s1
  integer :: ir, l
  ! ..
  if (nsra==1) ek = sqrt(e)
  if (nsra==2) ek = sqrt(e+e*e/(cvlight*cvlight))
  do l = 0, lmaxd

    if (nsra==2) then
      s1 = sqrt(real(l*l+l+1,kind=dp)-4.0e0_dp*z*z/(cvlight*cvlight))
      if (abs(z)<eps) s1 = real(l, kind=dp)
    else
      s1 = real(l, kind=dp)
    end if
    s(l) = s1
    rs(1, l) = 0.0e0_dp
    do ir = 2, irm
      rs(ir, l) = r(ir)**s1
    end do
    do ir = irm + 1, irmd
      rs(ir, l) = 0.0e0_dp
    end do

  end do
  return
end subroutine wfmesh
