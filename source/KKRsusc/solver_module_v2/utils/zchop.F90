  subroutine zchop(z,tol)
! set to zero real and imaginary part of complex numbers below tol
  use global, only: r8b, c8b

  implicit none

  complex(kind=c8b), intent(inout) :: z
  real(kind=r8b),    intent(in)    :: tol
! ----------------------------------------------------------------------
  real(kind=r8b) :: re, im

  re = real(z); im = aimag(z)
  if (abs(re) < tol) re = 0.d0
  if (abs(im) < tol) im = 0.d0
  z = cmplx(re,im)
! All done!
  end subroutine zchop
