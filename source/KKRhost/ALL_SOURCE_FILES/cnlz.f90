module mod_cnlz
  use :: mod_datatypes, only: dp
  private :: dp

contains

  function cnlz(l, z)
    ! ********************************************************************
    ! *                                                                  *
    ! *     von NEUMANN  - FUNCTION  N(L,Z)  FOR COMPLEX ARGUMENT  Z     *
    ! *                  see:  e.g. MERZBACHER EQ. (10.34)               *
    ! *                                                                  *
    ! ********************************************************************
    use :: mod_cjlz
    implicit none

    ! PARAMETER definitions
    complex (kind=dp) :: c1
    parameter (c1=(1.0e0_dp,0.0e0_dp))
    integer :: lp2max
    parameter (lp2max=25)

    ! Dummy arguments
    integer :: l
    complex (kind=dp) :: z
    complex (kind=dp) :: cnlz

    real (kind=dp) :: dfac
    complex (kind=dp) :: dt, s(lp2max), t, zsq
    integer :: i, k, llp1

    if (l<0) then
      cnlz = cjlz(l+1, z)
      return
    end if

    zsq = z*z
    llp1 = l + l + 1

    if (abs(zsq/real(llp1,kind=dp))<=10.e0_dp) then

      dfac = 1.0e0_dp
      do k = 3, llp1, 2
        dfac = dfac*real(k, kind=dp)
      end do

      dt = c1
      t = c1

      do i = 2, 400, 2
        dt = -dt*zsq/real(i*(i-llp1), kind=dp)
        t = t + dt
        if (abs(dt)<1.0e-10_dp) go to 100
      end do

100   continue
      cnlz = -t*z**(-l-1)*dfac/real(llp1, kind=dp)

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

  end function cnlz

end module mod_cnlz
