  subroutine spin_rotation(u0,u1,pauli,spinrot)
! Constructs the spin rotation from u0 to u1
  use global, only: i4b, r8b, c8b, iodb

  implicit none

! Initial direction
  real(kind=r8b),    intent(in)  :: u0(3)
! Final direction
  real(kind=r8b),    intent(in)  :: u1(3)
! Pauli matrices: x, y, z, 0
  complex(kind=c8b), intent(in)  :: pauli(2,2,4)
! Spin rotation matrix
  complex(kind=c8b), intent(out) :: spinrot(2,2)
! ----------------------------------------------------------------------
  real(kind=r8b),    parameter :: tol = 1.d-6
  complex(kind=c8b), parameter :: iu = (0.d0,1.d0)
  real(kind=r8b)    :: axis(3), u0len, u1len, cosa, alpha, sina

! Lengths of the vectors (if not unit vectors)
  u0len = sqrt(dot_product(u0,u0))
  u1len = sqrt(dot_product(u1,u1))
! Angle between them
  alpha = acos(dot_product(u0,u1)/(u0len*u1len))
! Trig
  cosa = cos(0.5d0*alpha)
  sina = sin(0.5d0*alpha)
! Rotation axis
  axis(1) = u1(2)*u0(3) - u1(3)*u0(2)
  axis(2) = u1(3)*u0(1) - u1(1)*u0(3)
  axis(3) = u1(1)*u0(2) - u1(2)*u0(1)
  u0len = sqrt(dot_product(axis,axis))
! Collinear?
  if (u0len < tol) then
    if (abs(cosa) < tol) then
      spinrot = pauli(:,:,1)  ! sigma_x <=> flip spin
    else
      spinrot = pauli(:,:,4)  ! sigma_0 <=> unit matrix
    end if
    return
  end if
  axis = axis/u0len
!  write(iodb,'("spinrotation: angle, axis=",4f16.8)') 45.d0*alpha/atan(1.d0), axis(1:3)
! Spin rotation matrix
  spinrot = cosa*pauli(:,:,4) + iu*sina*(axis(1)*pauli(:,:,1) + axis(2)*pauli(:,:,2) + axis(3)*pauli(:,:,3))
! All done!
  end subroutine spin_rotation
