!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_cjlz
  
  private
  public :: cjlz

contains

  !-------------------------------------------------------------------------------
  !> Summary: Computes spherial Bessel function
  !> Author: 
  !> Category: KKRhost, dirac, special-functions
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> SPHERICAL BESSEL-FUNCTION  J(L,Z)  FOR COMPLEX ARGUMENT  Z
  !> see:  e.g. MERZBACHER EQ. (10.22)
  !-------------------------------------------------------------------------------
  function cjlz(l, z)

    use :: mod_datatypes, only: dp
    use :: mod_constants, only: cone
    implicit none

    ! PARAMETER definitions
    integer :: lp2max
    parameter (lp2max=25)


    ! Dummy arguments
    integer :: l
    complex (kind=dp) :: z
    complex (kind=dp) :: cjlz

    ! Local variables
    real (kind=dp) :: dfac
    complex (kind=dp) :: dt, s(lp2max), t, zsq
    integer :: i, k, llp1

    zsq = z*z
    llp1 = l + l + 1

    if (abs(zsq/real(llp1,kind=dp))<=10.e0_dp) then

      dfac = 1.0e0_dp
      do k = 3, llp1, 2
        dfac = dfac*real(k, kind=dp)
      end do

      dt = cone
      t = cone
      do i = 2, 400, 2
        dt = -dt*zsq/real(i*(i+llp1), kind=dp)
        t = t + dt
        if (abs(dt)<1.0e-10_dp) go to 100
      end do

100   continue
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
  end function cjlz

end module mod_cjlz
