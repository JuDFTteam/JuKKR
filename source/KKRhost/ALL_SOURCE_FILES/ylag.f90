!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: lagrangian interpolation
!> Author: People who wrote it
!> Lagrangian interpolation xi is interpolated entry into x-array n is the order 
!> of lagrangran interpolation y is array from which ylag is obtained by interpolation
!> ind is the min-i for x(i).gt.xi if ind=0,x-array will be searched imax is 
!> max index of x-and y-arrays
!------------------------------------------------------------------------------------
!> @note 07/12/94  HE  arg. IEX removed
!> @endnote
!------------------------------------------------------------------------------------
module mod_ylag
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: lagrangian interpolation
  !> Author: 
  !> Category: numerical-tools, KKRhost 
  !> Deprecated: False 
  !> Lagrangian interpolation xi is interpolated entry into x-array n is the order 
  !> of lagrangran interpolation y is array from which ylag is obtained by interpolation
  !> ind is the min-i for x(i).gt.xi if ind=0,x-array will be searched imax is 
  !> max index of x-and y-arrays
  !-------------------------------------------------------------------------------
  !> @note 07/12/94  HE  arg. IEX removed
  !> @endnote
  !-------------------------------------------------------------------------------
  function ylag(xi, x, y, ind1, n1, imax)

    implicit none

    real (kind=dp), parameter :: eps = 1.0e-12_dp

    ! Dummy arguments
    integer, intent(in) :: n1
    integer, intent(in) :: imax
    integer, intent(in) :: ind1
    real (kind=dp), intent(in) :: xi
    real (kind=dp), dimension(imax), intent(in) :: x, y
    real (kind=dp) :: ylag

    ! Local variables
    real (kind=dp) :: d, p, s, xd
    integer :: i, ind, inl, inu, j, n
    save :: d, i, ind, inl, inu, j, n, p, s, xd

    ind = ind1
    n = n1
    if (n>imax) n = imax
    if (ind>0) go to 110
    do j = 1, imax
      if (abs(xi-x(j))<eps) go to 150
      if (xi<x(j)) go to 100
    end do
    go to 120
100 continue
    ind = j
110 continue
    if (ind>1) then
    end if
    inl = ind - (n+1)/2
    if (inl<=0) inl = 1
    inu = inl + n - 1
    if (inu<=imax) go to 130
120 continue
    inl = imax - n + 1
    inu = imax
130 continue
    s = 0.0e0_dp
    p = 1.0e0_dp
    do j = inl, inu
      p = p*(xi-x(j))
      d = 1.0e0_dp
      do i = inl, inu
        if (i/=j) then
          xd = x(j)
        else
          xd = xi
        end if
        d = d*(xd-x(i))
      end do
      s = s + y(j)/d
    end do
    ylag = s*p
140 continue
    return
150 continue
    ylag = y(j)
    go to 140
  end function ylag

end module mod_ylag
