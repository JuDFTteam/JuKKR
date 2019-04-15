!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: This subroutine calculates real spherical harmonics with the normalization : <y|y> =1 
!> Author: M. Weinert
!> This subroutine calculates real spherical harmonics with the normalization : <y|y> =1
!> returns also r = length of vector v
!> generate the complex spherical harmonics for the vector v
!> using a stable upward recursion in l. (see notes by m. weinert.)
!------------------------------------------------------------------------------------
!> @note Converted to real spherical harmonics B. Drittler 1987
!> @endnote
!------------------------------------------------------------------------------------
module mod_ymy
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: This subroutine calculates real spherical harmonics with the normalization : <y|y> =1 
  !> Author: M. Weinert
  !> Category: special-functions, KKRhost
  !> Deprecated: False 
  !> This subroutine calculates real spherical harmonics with the normalization : <y|y> =1
  !> returns also r = length of vector v
  !> generate the complex spherical harmonics for the vector v
  !> using a stable upward recursion in l. (see notes by m. weinert.)
  !-------------------------------------------------------------------------------
  !> @note Converted to real spherical harmonics B. Drittler 1987
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine ymy(v1, v2, v3, r, ylm, lmax)

    use :: mod_rcstop, only: rcstop
    use :: mod_constants, only : pi
    implicit none
    ! .. Parameters ..
    real (kind=dp) :: szero
    parameter (szero=1.0e-20_dp)
    ! ..
    integer, intent (in) :: lmax !! Maximum l component in wave function expansion
    real (kind=dp), intent (in) :: v1
    real (kind=dp), intent (in) :: v2
    real (kind=dp), intent (in) :: v3
    real (kind=dp), intent (out) :: r
    real (kind=dp), dimension((2*lmax+1)**2), intent (out) :: ylm !! real spherical harmonic to a given l,m 
    ! ..
    ! .. Local Scalars ..
    real (kind=dp) :: a, cd, cph, cth, fac, fpi, rtwo, sgm, sph, sth, t, xy, xyz
    integer :: i, l, m
    ! ..
    ! .. Local Arrays ..
    real (kind=dp), dimension(0:lmax) :: c
    real (kind=dp), dimension(0:lmax) :: s
    real (kind=dp), dimension(0:lmax,0:lmax) :: p
    ! ..
    fpi = 4.e0_dp*pi
    rtwo = sqrt(2.0e0_dp)

    ! --->    calculate sin and cos of theta and phi

    xy = v1**2 + v2**2
    xyz = xy + v3**2

    r = sqrt(xyz)
    if (xyz<=0.0e0_dp) then
      write (*, *) xyz
      call rcstop('ylm=0   ')

    else

      if (xy>szero*xyz) then
        xy = sqrt(xy)
        xyz = sqrt(xyz)
        cth = v3/xyz
        sth = xy/xyz
        cph = v1/xy
        sph = v2/xy

      else

        sth = 0.0e0_dp
        cth = 1.0e0_dp
        if (v3<0) cth = -1.0e0_dp
        cph = 1.0e0_dp
        sph = 0.0e0_dp
      end if
      ! --->    generate associated legendre functions for m.ge.0
      ! loop over m values
      fac = 1.0e0_dp
      do m = 0, lmax - 1
        fac = -(2*m-1)*fac
        p(m, m) = fac
        p(m+1, m) = (2*m+1)*cth*fac
        ! --->    recurse upward in l
        do l = m + 2, lmax
          p(l, m) = ((2*l-1)*cth*p(l-1,m)-(l+m-1)*p(l-2,m))/(l-m)
        end do
        fac = fac*sth
      end do
      p(lmax, lmax) = -(2*lmax-1)*fac
      ! --->    determine sin and cos of phi
      s(0) = 0.0e0_dp
      s(1) = sph
      c(0) = 1.0e0_dp
      c(1) = cph
      do m = 2, lmax
        s(m) = 2*cph*s(m-1) - s(m-2)
        c(m) = 2*cph*c(m-1) - c(m-2)
      end do
      ! --->    multiply in the normalization factors
      i = 0
      do l = 0, lmax
        i = i + l + 1
        a = sqrt((2*l+1)/fpi)
        cd = 1
        ylm(i) = a*p(l, 0)
        sgm = -rtwo
        do m = 1, l
          t = (l+1-m)*(l+m)
          cd = cd/t
          t = a*sqrt(cd)
          ylm(i+m) = sgm*t*p(l, m)*c(m)
          ylm(i-m) = sgm*t*p(l, m)*s(m)
          sgm = -sgm
        end do
        i = i + l
      end do

    end if

    return

  end subroutine ymy

end module mod_ymy
