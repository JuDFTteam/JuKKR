!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_cdjlzdz
  
  private
  public :: cdjlzdz

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculates derivative of Bessel function
  !> Author: 
  !> Category: KKRhost, special-functions, dirac
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Calculates derivative
  !> d j(L,Z) / dz
  !> analytically
  !-------------------------------------------------------------------------------
  function cdjlzdz(l, z, mode)
    use :: mod_datatypes, only: dp
    use :: mod_cjlz, only: cjlz
    implicit none

    ! Dummy arguments
    integer :: l, mode
    complex (kind=dp) :: z
    complex (kind=dp) :: cdjlzdz

    if (mode==1) then

      if (l==0) then

        cdjlzdz = l*cjlz(l, z)/z - cjlz(l+1, z)
      else
        cdjlzdz = (l*cjlz(l-1,z)-(l+1)*cjlz(l+1,z))/real(2*l+1, kind=dp)
        return
      end if
    else if (mode==2) then

      if (l==0) then
        cdjlzdz = l*cjlz(l, z)/z - cjlz(l+1, z)
      else
        cdjlzdz = cjlz(l-1, z) - (l+1)*cjlz(l, z)/z
        return
      end if
    else
      cdjlzdz = l*cjlz(l, z)/z - cjlz(l+1, z)
    end if
  end function cdjlzdz

end module mod_cdjlzdz
