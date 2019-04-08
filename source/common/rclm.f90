!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Transform matrix to spherical harmonics representation.
!> Author: Ph. Mavropoulos
!> Transform matrix to spherical harmonics representation.
!>
!> * (KEY=1) From real to complex spherical harmonics basis
!> * (KEY=1) From complex to real spherical harmonics basis
!>
!> Transformation matrix is \(A\): \(Y_{real} = A  Y_{cmplx} \)
!> \begin{equation}
!> A= 
!> \begin{bmatrix} I\frac{i}{\sqrt{2}} & 0 & -J\frac{i}{\sqrt{2}} \\ 0 & 1 & 0 \\ J\frac{1}{\sqrt{2}} & 0 & I\frac{1}{\sqrt{2}}\end{bmatrix}
!>\end{equation}
!> (of dimension \(2l+1\)), where \(I\) is the \(l\times l\) unit matrix,
!> \(J\) is the \(l\times l\) antidiagonal unit matrix 
!> \begin{equation}
!> J =
!> \begin{bmatrix} 0 & 0 & 1 \\ 0 & 1 & 0 \\ 1 & 0 & 0 \end{bmatrix}
!> \( A_{mm'} = \int Y_{lm} Y^{c *}_{lm'} d\Omega\)
!> \end{equation}
!> Transformation rule: Complex --> Real (VC --> VR)
!> \begin{equation}
!> VR_{mm'} = \sum_{m1,m2} A_{m,m1} VC_{m1,m2} A^*_{m'm2} 
!> \end{equation}
!> with
!> \begin{equation}
!> VC_{m1,m2} = \int Y^{c *}_{m2} V Y^{c}_{m1} 
!> \end{equation}
!> LDIM corresponds to the dimension of the array `VMAT` as `VMAT(2*LDIM+1,2*LDIM+1)`.
!> LL is the angular momentum l, corresponding to the part of `VMAT` used -- `VMAT(2*LL+1,2*LL+1)` 
!> Result returns in again in `VMAT` (original `VMAT` is destroyed) 
!------------------------------------------------------------------------------------
!> @warning Original `VMAT` is destroyed.
!> @endwarning
!------------------------------------------------------------------------------------
module mod_rclm
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Transform matrix to spherical harmonics representation.
  !> Author: Ph. Mavropoulos
  !> Category: special-functions, lda+u, numerical-tools, KKRhost
  !> Deprecated: False 
  !> Transform matrix to spherical harmonics representation.
  !>
  !> * (KEY=1) From real to complex spherical harmonics basis
  !> * (KEY=1) From complex to real spherical harmonics basis
  !>
  !> Transformation matrix is \(A\): \(Y_{real} = A  Y_{cmplx} \)
  !> \begin{equation}
  !> A= 
  !> \begin{bmatrix} I\frac{i}{\sqrt{2}} & 0 & -J\frac{i}{\sqrt{2}} \\ 0 & 1 & 0 \\ J\frac{1}{\sqrt{2}} & 0 & I\frac{1}{\sqrt{2}}\end{bmatrix}
  !>\end{equation}
  !> (of dimension \(2l+1\)), where \(I\) is the \(l\times l\) unit matrix,
  !> \(J\) is the \(l\times l\) antidiagonal unit matrix 
  !> \begin{equation}
  !> J =
  !> \begin{bmatrix} 0 & 0 & 1 \\ 0 & 1 & 0 \\ 1 & 0 & 0 \end{bmatrix}
  !> \end{equation}
  !> \( A_{mm'} = \int Y_{lm} Y^{c *}_{lm'} d\Omega\)
  !> Transformation rule: Complex --> Real (VC --> VR)
  !> \begin{equation}
  !> VR_{mm'} = \sum_{m1,m2} A_{m,m1} VC_{m1,m2} A^*_{m'm2} 
  !> \end{equation}
  !> with
  !> \begin{equation}
  !> VC_{m1,m2} = \int Y^{c *}_{m2} V Y^{c}_{m1} 
  !> \end{equation}
  !> LDIM corresponds to the dimension of the array `VMAT` as `VMAT(2*LDIM+1,2*LDIM+1)`.
  !> LL is the angular momentum l, corresponding to the part of `VMAT` used -- `VMAT(2*LL+1,2*LL+1)` 
  !> Result returns in again in `VMAT` (original `VMAT` is destroyed) 
  !-------------------------------------------------------------------------------
  !> @warning Original `VMAT` is destroyed.
  !> @endwarning
  !-------------------------------------------------------------------------------
  subroutine rclm(key, ll, ldim, vmat)

    use :: mod_cinit, only: cinit
    use :: mod_constants, only: cone,ci
    implicit none
    real (kind=dp), parameter :: eps = 1.0e-12_dp
    ! ..
    integer :: ldim, key, ll
    complex (kind=dp), dimension(2*ldim+1, 2*ldim+1) :: vmat
    ! .. Locals
    complex (kind=dp), dimension(2*ldim+1, 2*ldim+1) :: vtmp
    complex (kind=dp), dimension(2*ldim+1, 2*ldim+1) :: aa, aac
    real (kind=dp) :: ovsqrtwo
    complex (kind=dp) :: oneovrt, ciovrt, cimovrt, blj
    complex (kind=dp) :: a11, a13, a31, a33
    integer :: mdim, icall
    integer :: m1, m2, m3, m4, mm, midl

    data icall/0/
    save :: icall, oneovrt, ciovrt, cimovrt
    ! ======================================================================
    mdim = 2*ldim + 1
    if (ll>ldim) stop 'ERROR IN RCLM: LL>LDIM'
    ! ----------------------------------------------------------------------
    icall = icall + 1
    if (icall==1) then
      ovsqrtwo = 1e0_dp/sqrt(2e0_dp)
      oneovrt = ovsqrtwo*cone
      ciovrt = ovsqrtwo*ci
      cimovrt = -ciovrt
    end if
    ! ----------------------------------------------------------------------

    if (key==2) then               ! matrix elements
      a11 = ciovrt
      a13 = cimovrt
      a31 = oneovrt
      a33 = oneovrt
    else if (key==1) then          ! adjoint matrix elements
      a11 = conjg(ciovrt)
      a13 = oneovrt
      a31 = conjg(cimovrt)
      a33 = oneovrt
    else
      stop 'ERROR IN RCLM: KEY NOT EQUAL 1 OR 2'
    end if

    ! -> Construct transformation matrix AA

    call cinit(mdim*mdim, aa)
    mm = 2*ll + 1
    midl = ll + 1
    do m1 = 1, ll
      aa(m1, m1) = a11
      aa(m1, mm+1-m1) = a13
    end do

    aa(midl, midl) = cone
    do m1 = midl + 1, mm
      aa(m1, m1) = a33
      aa(m1, mm+1-m1) = a31
    end do

    ! -> Construct transformation matrix AA+

    do m2 = 1, mm
      do m1 = 1, mm
        aac(m1, m2) = conjg(aa(m2,m1))
      end do
    end do

    ! -> Multiply from left

    call cinit(mdim*mdim, vtmp)
    do m2 = 1, mm
      do m4 = 1, mm
        blj = aac(m4, m2)
        if (abs(blj)>eps) then
          do m3 = 1, mm
            vtmp(m3, m2) = vtmp(m3, m2) + vmat(m3, m4)*blj
          end do
        end if
      end do
    end do

    call cinit(mdim*mdim, vmat)
    do m2 = 1, mm
      do m3 = 1, mm
        blj = vtmp(m3, m2)
        if (abs(blj)>eps) then
          do m1 = 1, mm
            vmat(m1, m2) = vmat(m1, m2) + aa(m1, m3)*blj
          end do
        end if
      end do
    end do
  end subroutine rclm

end module mod_rclm
