subroutine wfmesh(e, ek, cvlight, nsra, z, r, s, rs, irm, irmd, lmaxd)
!.. Scalar Arguments ..
  double complex :: e, ek
  double precision :: cvlight, z
  integer :: irm, irmd, lmaxd, nsra
!..
!.. Intrinsic Functions ..
  intrinsic :: dble, sqrt
!..
!.. Array Arguments ..
  double precision :: r(irmd), rs(irmd, 0:lmaxd), s(0:lmaxd)
!..
!.. Local Scalars ..
  double precision :: s1
  integer :: ir, l
!..
  if (nsra==1) ek = sqrt(e)
  if (nsra==2) ek = sqrt(e+e*e/(cvlight*cvlight))
  do l = 0, lmaxd

    if (nsra==2) then
      s1 = sqrt(dble(l*l+l+1)-4.0d0*z*z/(cvlight*cvlight))
      if (z==0.0d0) s1 = dble(l)
    else
      s1 = dble(l)
    end if
    s(l) = s1
    rs(1, l) = 0.0d0
    do ir = 2, irm
      rs(ir, l) = r(ir)**s1
    end do
    do ir = irm + 1, irmd
      rs(ir, l) = 0.0d0
    end do

  end do
  return
end subroutine
