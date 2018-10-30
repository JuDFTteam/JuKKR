!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_fplaner
  
  private
  public :: fplaner

contains

  !-------------------------------------------------------------------------------
  !> Summary: Derivative of real-space contribution to Ewald sum
  !> Author: 
  !> Category: KKRhost, geometry
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> This sub calculates the derivatives of the real
  !> space contribution to the ewald sum .
  !>
  !>              l
  !>             d     erfc(lamda*sqrt(d*d+z*z))
  !>      lim    --   ------------------------
  !>      z->0     l        sqrt(d*d+z*z)
  !>             dz
  !>
  !> Up to l = 4 (l=1,3,5,7 etc vanish)
  !-------------------------------------------------------------------------------
  subroutine fplaner(alpha, g, r)

    use :: mod_datatypes, only: dp
    use :: mod_constants, only: pi
    implicit none

    real (kind=dp), parameter :: sqpi = sqrt(pi)
    real (kind=dp) :: alpha, g(0:4), r
    integer :: l
    real (kind=dp) :: lamda, er, ex, pref

    do l = 0, 4
      g(l) = 0.e0_dp
    end do

    er = erfc(alpha)
    ex = exp(-alpha*alpha)
    lamda = alpha/r

    g(0) = er/r

    pref = sqrt(5.e0_dp/pi)/4.e0_dp
    g(2) = -pref*(er/r/r/r+ex*2.e0_dp*lamda/r/r/sqpi)

    pref = 3.e0_dp*sqrt(9.e0_dp/pi)/16.e0_dp/9.e0_dp
    g(4) = pref*(9.e0_dp*er+ex*(12.e0_dp*alpha**3+18.e0_dp*alpha)/sqpi)/r/r/r/r/r

  end subroutine fplaner

end module mod_fplaner
