!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: This routine returns an array `y2(1:n)` of length `n` which contains the second derivatives of the interpolating function at the tabulated points `xi`.
!> Author: 
!> Given arrays `x(1:n)` and `y(1:n)` containing a tabulated function, i.e., 
!> \(y(i) = f(xi)\), with \(x1<x2<\cdots<xN\) , and given values `yp1` and `ypn`
!> for the 1st derivative of the interpolating function at points
!> 1 and `n`, respectively, this routine returns an array `y2(1:n)` of
!> length `n` which contains the second derivatives of the interpolating
!> function at the tabulated points `xi`.
!> If `yp1` and/or `ypn` are equal to `1.e30` or larger, the routine is
!> signaled to set the corresponding boundary condition for a natural
!> spline, with zero second derivative on that boundary.
!> Parameter: `NMAX` is the largest anticipated value of `n`.
!------------------------------------------------------------------------------------
!> @note Taken from "Numerical Recipes in Fortran 77", W.H.Press et al.
!> @endnote
!------------------------------------------------------------------------------------
module mod_spline
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: This routine returns an array `y2(1:n)` of length `n` which contains the second derivatives of the interpolating function at the tabulated points `xi`.
  !> Author: 
  !> Category: numerical-tools, KKRhost
  !> Deprecated: False
  !> Given arrays `x(1:n)` and `y(1:n)` containing a tabulated function, i.e., 
  !> \(y(i) = f(xi)\), with \(x1<x2<\cdots<xN\) , and given values `yp1` and `ypn`
  !> for the 1st derivative of the interpolating function at points
  !> 1 and `n`, respectively, this routine returns an array `y2(1:n)` of
  !> length `n` which contains the second derivatives of the interpolating
  !> function at the tabulated points `xi`.
  !> If `yp1` and/or `ypn` are equal to `1.e30` or larger, the routine is
  !> signaled to set the corresponding boundary condition for a natural
  !> spline, with zero second derivative on that boundary.
  !> Parameter: `NMAX` is the largest anticipated value of `n`.
  !-------------------------------------------------------------------------------
  !> @note Taken from "Numerical Recipes in Fortran 77", W.H.Press et al.
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine spline_real(nmax, x, y, n, yp1, ypn, y2)
    implicit none
    integer :: n, nmax
    real (kind=dp) :: yp1, ypn
    real (kind=dp), dimension(nmax) :: x, y, y2
    integer :: i, k
    complex (kind=dp) :: p, qn, sig, un
    complex (kind=dp), dimension(nmax) :: u

    if (n>nmax) stop 'SPLINE: n > NMAX.'
    if (abs(yp1)>0.99e30_dp) then
      ! The lower boundary condition is set either to be "natural"
      y2(1) = 0.0_dp
      u(1) = 0.0_dp
    else
      ! or else to have a specified first derivative.
      y2(1) = -0.50_dp
      u(1) = (3.0_dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    end if

    do i = 2, n - 1
      ! This is the decomposition loop of the tridiagonal algorithm. y2 and u
      ! are used for temporary storage of the decomposed factors.
      sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
      p = sig*y2(i-1) + 2.0_dp
      y2(i) = (sig-1.0_dp)/p
      u(i) = (6.0_dp*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    end do

    if (abs(ypn)>0.99e30_dp) then
      ! The upper boundary condition is set either to be "natural"
      qn = 0.0_dp
      un = 0.0_dp
    else
      ! or else to have a specified 1rst derivative.
      qn = 0.5_dp
      un = (3.0_dp/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    end if
    y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0_dp)
    do k = n - 1, 1, -1
      ! This is the backsubstitution loop of the tridiagonal algorithm.
      y2(k) = y2(k)*y2(k+1) + u(k)
    end do

    return
  end subroutine spline_real


  !-------------------------------------------------------------------------------
  !> Summary: This routine returns an array `y2(1:n)` of length `n` which contains the second derivatives of the interpolating function at the tabulated points `xi`.
  !> Author: 
  !> Category: numerical-tools, KKRhost
  !> Deprecated: False
  !> Complex version of `spline_real`
  !-------------------------------------------------------------------------------
  subroutine spline_complex(nmax, x, y, n, yp1, ypn, y2)
    implicit none
    integer :: n, nmax
    complex (kind=dp) :: yp1, ypn
    complex (kind=dp), dimension(nmax) :: y, y2
    real (kind=dp), dimension(nmax) :: x
    integer :: i, k
    real (kind=dp) :: p, qn, sig, un
    real (kind=dp), dimension(nmax) :: u

    if (n>nmax) stop 'SPLINE: n > NMAX.'
    if (abs(yp1)>0.99e30_dp) then
      ! The lower boundary condition is set either to be "natural"
      y2(1) = 0.0_dp
      u(1) = 0.0_dp
    else
      ! or else to have a specified first derivative.
      y2(1) = -0.5_dp
      u(1) = (3.0_dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    end if

    do i = 2, n - 1
      ! This is the decomposition loop of the tridiagonal algorithm. y2 and u
      ! are used for temporary storage of the decomposed factors.
      sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
      p = sig*y2(i-1) + 2.0_dp
      y2(i) = (sig-1.0_dp)/p
      u(i) = (6.0_dp*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    end do

    if (abs(ypn)>0.99e30_dp) then
      ! The upper boundary condition is set either to be "natural"
      qn = 0.0_dp
      un = 0.0_dp
    else
      ! or else to have a specified 1rst derivative.
      qn = 0.5_dp
      un = (3.0_dp/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    end if
    y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0_dp)
    do k = n - 1, 1, -1
      ! This is the backsubstitution loop of the tridiagonal algorithm.
      y2(k) = y2(k)*y2(k+1) + u(k)
    end do

    return
  end subroutine spline_complex

end module mod_spline
