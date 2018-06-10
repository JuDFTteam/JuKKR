function cnlz(l, z)
!   ********************************************************************
!   *                                                                  *
!   *     von NEUMANN  - FUNCTION  N(L,Z)  FOR COMPLEX ARGUMENT  Z     *
!   *                  see:  e.g. MERZBACHER EQ. (10.34)               *
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
  complex *16 :: cnlz

! Local variables
  complex *16 :: cjlz

  real *8 :: dfac
  complex *16 :: dt, s(lp2max), t, zsq
  integer :: i, k, llp1

  if (l<0) then
    cnlz = cjlz(l+1, z)
    return
  end if

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
      dt = -dt*zsq/dble(i*(i-llp1))
      t = t + dt
      if (abs(dt)<1.0d-10) go to 100
    end do

100 continue
    cnlz = -t*z**(-l-1)*dfac/dble(llp1)

  else
    if (l>23) stop '<cnlz>: l too large'
    s(2) = cos(z)
    if (l<=0) then
      cnlz = -s(2)*z**(-l-1)
      return
    end if

    s(1) = -sin(z)/z
    do i = 3, l + 2
      s(i) = s(i-1)*(2*i-5) - zsq*s(i-2)
    end do
    cnlz = -s(l+2)*z**(-l-1)

  end if

end function
