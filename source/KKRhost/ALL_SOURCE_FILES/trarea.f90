!------------------------------------------------------------------------------------
!> Summary: From complex to real (differenciated spherical harmonics)
!> Author: 
!> From complex to real (differenciated spherical harmonics)
!------------------------------------------------------------------------------------
module mod_trarea
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: From complex to real (differenciated spherical harmonics)
  !> Author: 
  !> Category: special-functions, numerical-tools, KKRhost
  !> Deprecated: False 
  !> From complex to real (differenciated spherical harmonics)
  !-------------------------------------------------------------------------------
  subroutine trarea(a, b, lmax)

    use :: mod_constants, only: ci

    ! .. Parameters ..
    real (kind=dp) :: rtwo
    parameter (rtwo=1.414213562373e0_dp)
    ! .. Input variables
    integer, intent(in) :: lmax !! Maximum l component in wave function expansion
    complex (kind=dp), dimension(*), intent(in) :: a !! Input complex spherical harmonics
    ! .. Output variables
    real (kind=dp), dimension(*), intent(out) :: b  !! Output real derivative of an spherical harmonic
    ! ..
    ! .. Local Scalars ..
    real (kind=dp) :: sgm
    integer :: i, l, m
    ! ..
    ! .. Intrinsic Functions ..
    intrinsic :: conjg, real

    ! calculate real the spherical harmonics derivetived
    i = 0
    do l = 0, lmax
      i = i + l + 1
      b(i) = real(a(i))
      sgm = -1.e0_dp
      do m = 1, l
        b(i-m) = real(ci*(a(i-m)-conjg(a(i-m))))/rtwo
        b(i+m) = sgm*real((a(i+m)+conjg(a(i+m))))/rtwo
        sgm = -sgm
      end do
      i = i + l
    end do
    return
  end subroutine trarea

end module mod_trarea
