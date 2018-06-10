function cjlz(l, z)
!   ********************************************************************
!   *                                                                  *
!   *   SPHERICAL BESSEL-FUNCTION  J(L,Z)  FOR COMPLEX ARGUMENT  Z     *
!   *                  see:  e.g. MERZBACHER EQ. (10.22)               *
!   *                                                                  *
!   ********************************************************************

  implicit none

! PARAMETER definitions
  complex *16 :: c1
  parameter (c1=(1.0d0,0.0d0))
  integer :: lp2max
  parameter (lp2max=25)


! Dummy arguments
  integer :: l
  complex *16 :: z
  complex *16 :: cjlz

! Local variables
  real *8 :: dfac
  complex *16 :: dt, s(lp2max), t, zsq
  integer :: i, k, llp1

  zsq = z*z
  llp1 = l + l + 1

  if (abs(zsq/dble(llp1))<=10.d0) then

    dfac = 1.0d0
    do k = 3, llp1, 2
      dfac = dfac*dble(k)
    end do

    dt = c1
    t = c1
    do i = 2, 400, 2
      dt = -dt*zsq/dble(i*(i+llp1))
      t = t + dt
      if (abs(dt)<1.0d-10) go to 100
    end do

100 continue
    cjlz = t*z**l/dfac

  else
    if (l>23) stop '<cjlz>: l too large'

    s(2) = sin(z)/z
    if (l<=0) then
      cjlz = s(2)*z**l
      return
    end if

    s(1) = cos(z)
    do i = 3, l + 2
      s(i) = (s(i-1)*(2*i-5)-s(i-2))/zsq
    end do
    cjlz = s(l+2)*z**l

  end if
end function
