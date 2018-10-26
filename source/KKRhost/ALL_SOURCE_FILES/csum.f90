!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_csum
  
  use :: mod_datatypes, only: dp
  private
  public :: csum

contains

  !-------------------------------------------------------------------------------
  !> Summary: Sums complex array elements with given length and stepwidth
  !> Author: 
  !> Date: 19.10.95
  !> Category: KKRhost, undefined
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Sum up the first N elements of the complex
  !> array V(*) with a stepwidth of IV
  !-------------------------------------------------------------------------------
  complex (kind=dp) function csum(n, v, iv)
    use :: mod_constants, only: czero
    implicit none
    ! .. Scalar Arguments ..
    integer, intent(in) :: iv, n
    ! ..
    ! .. Array Arguments ..
    complex (kind=dp), intent(in) :: v(*)
    ! ..
    ! .. Local Scalars ..
    complex (kind=dp) :: vsum
    integer :: i, ibot, itop
    ! ..
    if (iv>=0) then
      ibot = 1
      itop = 1 + (n-1)*iv

    else
      ibot = 1 - (n-1)*iv
      itop = 1
    end if

    vsum = czero
    do i = ibot, itop, iv
      vsum = vsum + v(i)
    end do
    csum = vsum
    return
  end function csum

end module mod_csum
