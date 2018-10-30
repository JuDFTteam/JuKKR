!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Sum up the first `N` elements of the `real (kind=dp)` array `V(*)` with a stepwidth of `IV`
!> Author: 
!> Sum up the first `N` elements of the `real (kind=dp)` array `V(*)` with a stepwidth of `IV`
!------------------------------------------------------------------------------------
module mod_ssum
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Sum up the first `N` elements of the `real (kind=dp)` array `V(*)` with a stepwidth of `IV`
  !> Author: 
  !> Category: numerical-tools, KKRhost
  !> Deprecated: False 
  !> Sum up the first `N` elements of the `real (kind=dp)` array `V(*)` with a stepwidth of `IV`
  !-------------------------------------------------------------------------------
  function ssum(n, v, iv)

    implicit none
    real (kind=dp) :: ssum
    ! .. Scalar Arguments ..
    integer, intent(in) :: n  !! Number of elements to sum over
    integer, intent(in) :: iv !! Stepwidth
    real (kind=dp), dimension(*), intent(in) :: v !! Input vector
    ! ..
    ! .. Local Scalars ..
    real (kind=dp) :: vsum
    integer :: i, ibot, itop
    ! ..
    if (iv>=0) then
      ibot = 1
      itop = 1 + (n-1)*iv

    else
      ibot = 1 - (n-1)*iv
      itop = 1
    end if

    vsum = 0.0e0_dp
    do i = ibot, itop, iv
      vsum = vsum + v(i)
    end do
    ssum = vsum
  end function ssum

end module mod_ssum
