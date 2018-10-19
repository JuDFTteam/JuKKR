!------------------------------------------------------------------------------------
!> Summary: Creation of the radial mesh for the old-solver 
!> Author: 
!> Creation of the radial mesh for the old-solver 
!------------------------------------------------------------------------------------
module mod_wfmesh
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Creation of the radial mesh for the old-solver 
  !> Author: 
  !> Category: radial-mesh, KKRhost
  !> Deprecated: False 
  !> Creation of the radial mesh for the old-solver 
  !-------------------------------------------------------------------------------
  subroutine wfmesh(e,ek,cvlight,nsra,z,r,s,rs,irm,irmd,lmaxd)

    real (kind=dp), parameter :: eps = 1.0e-12_dp
    ! .. Scalar Arguments ..
    integer, intent(in) :: irm  !! Current radial point
    integer, intent(in) :: irmd !! Maximum number of radial points
    integer, intent(in) :: nsra
    integer, intent(in) :: lmaxd !! Maximum l component in wave function expansion
    real (kind=dp), intent(in) :: z       !! Nuclear charge
    real (kind=dp), intent(in) :: cvlight !! Speed of light divided by the fine structure constant
    complex (kind=dp), intent(in) :: e
    real (kind=dp), dimension(irmd), intent(in) :: r !! Set of real space vectors (in a.u.)
    ! .. In/Out variables
    complex (kind=dp), intent(inout) :: ek
    real (kind=dp), dimension(0:lmaxd), intent(inout) :: s
    real (kind=dp), dimension(irmd, 0:lmaxd), intent(inout) :: rs
    ! .. Intrinsic Functions ..
    intrinsic :: real, sqrt
    ! .. Local Scalars ..
    real (kind=dp) :: s1
    integer :: ir, l
    ! ..
    if (nsra==1) ek = sqrt(e)
    if (nsra==2) ek = sqrt(e+e*e/(cvlight*cvlight))
    do l = 0, lmaxd

      if (nsra==2) then
        s1 = sqrt(real(l*l+l+1,kind=dp)-4.0e0_dp*z*z/(cvlight*cvlight))
        if (abs(z)<eps) s1 = real(l, kind=dp)
      else
        s1 = real(l, kind=dp)
      end if
      s(l) = s1
      rs(1, l) = 0.0e0_dp
      do ir = 2, irm
        rs(ir, l) = r(ir)**s1
      end do
      do ir = irm + 1, irmd
        rs(ir, l) = 0.0e0_dp
      end do

    end do
    return
  end subroutine wfmesh

end module mod_wfmesh
