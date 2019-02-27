!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Wrapper to setup the rotation matrices to transform from the local to the global frame of references
!> Author: 
!> Wrapper to setup the rotation matrices to transform from the local to the global frame of references
!------------------------------------------------------------------------------------
module mod_rotatespinframe

  use :: mod_datatypes, only: dp

  character (len=25), parameter :: spinmode= 'kkr'

  private
  public :: spinmode, rotatematrix, rotatevector

contains

  !-------------------------------------------------------------------------------
  !> Summary: Rotates a matrix in the local frame pointing in the direction of \(\phi\) and \(\theta\) to the global frame
  !> Author: 
  !> Category: numerical-tools, KKRhost
  !> Deprecated: False 
  !> Rotates a matrix in the local frame pointing in the direction of phi and theta to the global frame
  !> \begin{equation}
  !> U=
  !> \begin{bmatrix} \cos{\frac{\theta}{2}} \exp{-\frac{i\phi}{2}} & -\sin{\frac{\theta}{2}}\exp{-i\frac{i\phi}{2}} \\ \sin{\frac{\theta}{2}} \exp{ \frac{i\phi}{2}} & \cos{\frac{\theta}{2}} \exp{ \frac{i\phi}{2}}\end{bmatrix}
  !> \end{equation}
  !> `Udegga = transpose(complex conjug ( U ) )`
  !> 
  !> @note
  !>   * mode=0: 'loc->glob'
  !>   * mode=1: 'glob->loc'
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine rotatematrix(mat, theta, phi, lmmax, mode)
    implicit none
    ! interface
    integer, intent (in) :: lmmax
    integer, intent (in) :: mode
    real (kind=dp), intent (in) :: phi
    real (kind=dp), intent (in) :: theta
    complex (kind=dp), dimension (2*lmmax, 2*lmmax), intent (inout) :: mat
    ! local
    complex (kind=dp), dimension (2*lmmax, 2*lmmax) :: umat
    complex (kind=dp), dimension (2*lmmax, 2*lmmax) :: udeggamat
    complex (kind=dp), dimension (2*lmmax, 2*lmmax) :: mattemp

    ! ***********************************************************************
    ! create the rotation matrix:
    ! | cos(theta/2) exp(-i/2 phi)   -sin(theta/2) exp(-i/2 phi) |
    ! U= |                                                          |
    ! | sin(theta/2) exp( i/2 phi)    cos(theta/2) exp( i/2 phi) |

    ! Udegga = transpose(complex conjug ( U ) )
    ! ***********************************************************************


    call create_umatrix(theta, phi, lmmax, umat, udeggamat)
    ! ***********************************************************************
    ! calculate matrix in the global frame:

    ! t_glob = U * t_loc * Udegga
    ! ***********************************************************************


    if (mode==0) then              ! 'loc->glob'
      call zgemm('N', 'N', 2*lmmax, 2*lmmax, 2*lmmax, (1e0_dp,0e0_dp), mat, 2*lmmax, udeggamat, 2*lmmax, (0e0_dp,0e0_dp), mattemp, 2*lmmax)
      call zgemm('N', 'N', 2*lmmax, 2*lmmax, 2*lmmax, (1e0_dp,0e0_dp), umat, 2*lmmax, mattemp, 2*lmmax, (0e0_dp,0e0_dp), mat, 2*lmmax)
    else if (mode==1) then         ! 'glob->loc'
      call zgemm('N', 'N', 2*lmmax, 2*lmmax, 2*lmmax, (1e0_dp,0e0_dp), mat, 2*lmmax, umat, 2*lmmax, (0e0_dp,0e0_dp), mattemp, 2*lmmax)
      call zgemm('N', 'N', 2*lmmax, 2*lmmax, 2*lmmax, (1e0_dp,0e0_dp), udeggamat, 2*lmmax, mattemp, 2*lmmax, (0e0_dp,0e0_dp), mat, 2*lmmax)
    else
      stop '[rotatematrix] mode not known'
    end if

  end subroutine rotatematrix

  !-------------------------------------------------------------------------------
  !> Summary: Does the rotation from the old local to the new local spin frame reference
  !> Author: 
  !> Category: numerical-tools, KKRhost
  !> Deprecated: False 
  !> Does the rotation from the old local to the new local spin frame for densities and charges
  !> \begin{equation}
  !> \rho_{loc}(ir,lm)= W1 \rho_{glob}(ir,lm) W2
  !> \end{equation}
  !> where \(\rho\) and \(W\) are matricies in spin space 
  !-------------------------------------------------------------------------------
  subroutine rotatevector(rho2nsc,rho2ns,nrmax,lmpotd,theta,phi,theta_old,phi_old,  &
    nrmaxd)

    implicit none
    ! interface
    integer :: nrmaxd, lmpotd, nrmax
    complex (kind=dp) :: rho2nsc(nrmaxd, lmpotd, 4)
    real (kind=dp) :: rho2ns(nrmaxd, lmpotd, 4)
    real (kind=dp) :: theta, phi
    real (kind=dp) :: theta_old, phi_old
    ! local
    integer :: ir, ilm
    complex (kind=dp) :: w1(2, 2), w2(2, 2)
    complex (kind=dp) :: w1_11w2_11, w1_11w2_21 !w1_11w2_22, w1_11w2_12
    complex (kind=dp) :: w1_22w2_22, w1_22w2_12 !w1_22w2_11, w1_22w2_21 
    complex (kind=dp) :: w1_12w2_11, w1_12w2_21 !w1_12w2_22, w1_12w2_12 
    complex (kind=dp) :: w1_21w2_22, w1_21w2_12 !w1_21w2_11, w1_21w2_21

    call create_wmatrix(theta, phi, theta_old, phi_old, 1, w1, w2)

    ! here some are unused and thus commented out
    w1_11w2_11 = w1(1, 1)*w2(1, 1)
    !w1_11w2_22 = w1(1, 1)*w2(2, 2)
    !w1_11w2_12 = w1(1, 1)*w2(1, 2)
    w1_11w2_21 = w1(1, 1)*w2(2, 1)
    !w1_22w2_11 = w1(2, 2)*w2(1, 1)
    w1_22w2_22 = w1(2, 2)*w2(2, 2)
    w1_22w2_12 = w1(2, 2)*w2(1, 2)
    !w1_22w2_21 = w1(2, 2)*w2(2, 1)
    w1_12w2_11 = w1(1, 2)*w2(1, 1)
    !w1_12w2_22 = w1(1, 2)*w2(2, 2)
    !w1_12w2_12 = w1(1, 2)*w2(1, 2)
    w1_12w2_21 = w1(1, 2)*w2(2, 1)
    !w1_21w2_11 = w1(2, 1)*w2(1, 1)
    w1_21w2_22 = w1(2, 1)*w2(2, 2)
    w1_21w2_12 = w1(2, 1)*w2(1, 2)
    !w1_21w2_21 = w1(2, 1)*w2(2, 1)

    do ir = 1, nrmax
      do ilm = 1, lmpotd
        rho2ns(ir, ilm, 1) = aimag(+(rho2nsc(ir,ilm,1)*w1_11w2_11)+(rho2nsc(ir,ilm,3)*w1_11w2_21)+(rho2nsc(ir,ilm,4)*w1_12w2_11)+(rho2nsc(ir,ilm,2)*w1_12w2_21))

        rho2ns(ir, ilm, 2) = aimag(+(rho2nsc(ir,ilm,1)*w1_21w2_12)-(rho2nsc(ir,ilm,3)*w1_21w2_22)-(rho2nsc(ir,ilm,4)*w1_22w2_12)+(rho2nsc(ir,ilm,2)*w1_22w2_22))
      end do
    end do

  end subroutine rotatevector

  !-------------------------------------------------------------------------------
  !> Summary: Create the rotation matrix \(W\):
  !> Author: 
  !> Category: numerical-tools, KKRhost
  !> Deprecated: False 
  !> Create the rotation matrix \(W\)
  !> \begin{equation}
  !> W1= U_{degga_{locnew}}  U_{locold}
  !> \end{equation}
  !> \begin{equation}
  !> W2= U_{degga_{locold}}  U_{locnew}
  !> \end{equation}
  !> The rotation matrix is created such that it rotates an operator
  !> which is in a local frame (locold) to another local frame (locnew)
  !> This is done by first transforming the old local frame to the
  !> global frame using the U matrix and then transforming the global
  !> frame to the new local frame
  !> \begin{equation}
  !> A_{locnew} = W1  A_{locold}  W2
  !> \end{equation}
  !> `Udegga = transpose(complex conjug ( U ) )`
  !-------------------------------------------------------------------------------
  subroutine create_wmatrix(theta, phi, theta_old, phi_old, lmmax, wmat1, wmat2)
    implicit none
    ! interface
    real (kind=dp), intent (in) :: phi
    real (kind=dp), intent (in) :: theta
    real (kind=dp), intent (in) :: phi_old
    real (kind=dp), intent (in) :: theta_old
    integer, intent (in) :: lmmax
    complex (kind=dp), intent (out) :: wmat1(2*lmmax, 2*lmmax)
    complex (kind=dp), intent (out) :: wmat2(2*lmmax, 2*lmmax)
    ! local
    complex (kind=dp) :: umat1(2*lmmax, 2*lmmax)
    complex (kind=dp) :: udeggamat1(2*lmmax, 2*lmmax)
    complex (kind=dp) :: umat2(2*lmmax, 2*lmmax)
    complex (kind=dp) :: udeggamat2(2*lmmax, 2*lmmax)

    call create_umatrix(theta_old, phi_old, lmmax, umat1, udeggamat1)

    call create_umatrix(theta, phi, lmmax, umat2, udeggamat2)

    call zgemm('N', 'N', 2*lmmax, 2*lmmax, 2*lmmax, (1e0_dp,0e0_dp), udeggamat2, 2*lmmax, umat1, 2*lmmax, (0e0_dp,0e0_dp), wmat1, 2*lmmax)
    call zgemm('N', 'N', 2*lmmax, 2*lmmax, 2*lmmax, (1e0_dp,0e0_dp), udeggamat1, 2*lmmax, umat2, 2*lmmax, (0e0_dp,0e0_dp), wmat2, 2*lmmax)


  end subroutine create_wmatrix

  !-------------------------------------------------------------------------------
  !> Summary: Create the rotation matrix \(U\)
  !> Author: 
  !> Category: numerical-tools, KKRhost
  !> Deprecated: False 
  !> Create the rotation matrix \(U\)
  !> \begin{equation}
  !> U=
  !> \begin{bmatrix} \cos{\frac{\theta}{2}} \exp{-\frac{i\phi}{2}} & -\sin{\frac{\theta}{2}}\exp{-i\frac{i\phi}{2}} \\ \sin{\frac{\theta}{2}} \exp{ \frac{i\phi}{2}} & \cos{\frac{\theta}{2}} \exp{ \frac{i\phi}{2}}\end{bmatrix}
  !> \end{equation}
  !> `Udegga = transpose(complex conjug ( U ) )`
  !-------------------------------------------------------------------------------
  subroutine create_umatrix(theta, phi, lmmax, umat, udeggamat)

    use :: mod_constants, only: ci, czero

    implicit none
    ! interface
    real (kind=dp), intent (in) :: phi
    real (kind=dp), intent (in) :: theta
    integer, intent (in) :: lmmax
    complex (kind=dp), intent (out) :: umat(2*lmmax, 2*lmmax)
    complex (kind=dp), intent (out) :: udeggamat(2*lmmax, 2*lmmax)
    ! local
    complex (kind=dp) :: umat11, umat12, umat21, umat22
    complex (kind=dp) :: udeggamat11, udeggamat12, udeggamat21, udeggamat22
    integer :: ival

    if (spinmode=='regular') then
      umat11 = cos(theta/2.0e0_dp)*exp(-ci/2.0e0_dp*phi)
      umat12 = -sin(theta/2.0e0_dp)*exp(-ci/2.0e0_dp*phi)
      umat21 = sin(theta/2.0e0_dp)*exp(ci/2.0e0_dp*phi)
      umat22 = cos(theta/2.0e0_dp)*exp(ci/2.0e0_dp*phi)
    else if (spinmode=='kkr') then
      umat11 = cos(theta/2.0e0_dp)*exp(ci/2.0e0_dp*phi)
      umat12 = sin(theta/2.0e0_dp)*exp(ci/2.0e0_dp*phi)
      umat21 = -sin(theta/2.0e0_dp)*exp(-ci/2.0e0_dp*phi)
      umat22 = cos(theta/2.0e0_dp)*exp(-ci/2.0e0_dp*phi)
    else
      stop '[create_Umatrix] mode not known'
    end if

    umat = czero
    do ival = 1, lmmax
      umat(ival, ival) = umat11
      umat(ival, lmmax+ival) = umat12
      umat(lmmax+ival, ival) = umat21
      umat(lmmax+ival, lmmax+ival) = umat22
    end do

    if (spinmode=='regular') then
      udeggamat11 = cos(theta/2.0e0_dp)*exp(ci/2.0e0_dp*phi)
      udeggamat12 = sin(theta/2.0e0_dp)*exp(-ci/2.0e0_dp*phi)
      udeggamat21 = -sin(theta/2.0e0_dp)*exp(ci/2.0e0_dp*phi)
      udeggamat22 = cos(theta/2.0e0_dp)*exp(-ci/2.0e0_dp*phi)
    else if (spinmode=='kkr') then
      udeggamat11 = cos(theta/2.0e0_dp)*exp(-ci/2.0e0_dp*phi)
      udeggamat12 = -sin(theta/2.0e0_dp)*exp(ci/2.0e0_dp*phi)
      udeggamat21 = sin(theta/2.0e0_dp)*exp(-ci/2.0e0_dp*phi)
      udeggamat22 = cos(theta/2.0e0_dp)*exp(ci/2.0e0_dp*phi)
    else
      stop '[create_Umatrix] mode not known'
    end if

    udeggamat = czero
    do ival = 1, lmmax
      udeggamat(ival, ival) = udeggamat11
      udeggamat(ival, lmmax+ival) = udeggamat12
      udeggamat(lmmax+ival, ival) = udeggamat21
      udeggamat(lmmax+ival, lmmax+ival) = udeggamat22
    end do

  end subroutine create_umatrix

end module mod_rotatespinframe
