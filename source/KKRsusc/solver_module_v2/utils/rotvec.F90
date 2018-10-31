  subroutine rotvec(u0,u1,rotmat)
! Rotation matrix rotating vector u0 to vector u1
! Rodrigues rotation formula

  use global, only: iodb, i4b, r8b

  implicit none

  real(kind=r8b), intent(in)  :: u0(3)
  real(kind=r8b), intent(in)  :: u1(3)
  real(kind=r8b), intent(out) :: rotmat(3,3)
! ---------------------------------------
  real(kind=r8b),   parameter :: tol = 1.d-6, rtol = 1.d-8
  real(kind=r8b)    :: axis(3), u0len, u1len, u01len, cosa, sina, cross(3,3), ui(3), fac
  integer(kind=i4b) :: i, j

! Lengths of the vectors (if not unit vectors)
  u0len = sqrt(dot_product(u0,u0))
  u1len = sqrt(dot_product(u1,u1))
  if (u0len < tol .or. u1len < tol) stop 'rotvec: check u0 and u1'
  cosa = dot_product(u0,u1)/(u0len*u1len)
! Rotation axis: axis = u0 x u1 / | u0 x u1 |
  axis(1) = u0(2)*u1(3) - u0(3)*u1(2)
  axis(2) = u0(3)*u1(1) - u0(1)*u1(3)
  axis(3) = u0(1)*u1(2) - u0(2)*u1(1)
! chop small numbers
  where (abs(axis) < tol) axis = 0.d0
  u01len = sqrt(dot_product(axis,axis))
! ----------------------------------------------------------------------
  if (u01len > rtol) then
!   noncollinear: no problem
    sina = u01len/(u0len*u1len)
    axis = axis/u01len
    fac = 1.d0 - cosa
  else
!   collinear: handle zeros
    if (cosa > 0.d0) then
!     cosa = +1 => no rotation
      cosa = +1.d0; sina = 0.d0
      axis(:) = 0.d0
      fac = 0.d0
      write(iodb,'("rotvec: identity")')
    else
!     cosa = -1 => inversion
      cosa = -1.d0; sina = 0.d0
      fac = 2.d0
!     choose a rotation axis
      i = minloc(abs(u0),dim=1)
      ui(:) = 0.d0
      ui(i) = 1.d0
      axis(1) = u0(2)*ui(3) - u0(3)*ui(2)
      axis(2) = u0(3)*ui(1) - u0(1)*ui(3)
      axis(3) = u0(1)*ui(2) - u0(2)*ui(1)
      u01len = sqrt(dot_product(axis,axis))
      axis(:) = axis(:)/u01len
      write(iodb,'("rotvec: inversion")')
    end if
  end if
! ----------------------------------------------------------------------
  write(iodb,'("cosa, sina, rot axis=",5es16.8)') cosa, sina, axis
! cross product matrix
  cross(:,:) = 0.d0
  cross(2,1) =  axis(3); cross(1,2) = -axis(3)
  cross(3,1) = -axis(2); cross(1,3) =  axis(2)
  cross(3,2) =  axis(1); cross(2,3) = -axis(1)
! rotation matrix
  rotmat(:,:) = sina*cross(:,:) + fac*matmul(cross,cross)
  do i=1,3
    rotmat(i,i) = rotmat(i,i) + 1.d0
  end do
! All done!
  end subroutine rotvec

