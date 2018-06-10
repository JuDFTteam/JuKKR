subroutine simp3(f, fint, istart, iend, drdi)
!-----------------------------------------------------------------------
!     this subroutine does an integration from istart to iend of
!     the real function f with an extended 3-point-simpson :

!                          r(istart)

!                       fint = { f(r') dr'

!                           r(iend)

!-----------------------------------------------------------------------
!.. Scalar Arguments ..
  double precision :: fint
  integer :: iend, istart
!..
!.. Array Arguments ..
  double precision :: drdi(*), f(*)
!..
!.. Local Scalars ..
  double precision :: a1, a2
  integer :: i, ist
!..
!.. Intrinsic Functions ..
  intrinsic :: mod
!..
  a1 = 4.0d0/3.0d0
  a2 = 2.0d0/3.0d0

!---> initialize fint

  if (mod(iend-istart,2)==0) then
    fint = f(istart)*drdi(istart)/3.0d0
    ist = istart + 1

  else
    fint = (f(istart+3)*drdi(istart+3)-5.0d0*f(istart+2)*drdi(istart+2)+ &
      19.0d0*f(istart+1)*drdi(istart+1)+9.0d0*f(istart)*drdi(istart))/24.0d0 + &
      f(istart+1)*drdi(istart+1)/3.0d0
    ist = istart + 2
  end if

!---> calculate with an extended 3-point-simpson

  do i = ist, iend - 1, 2
    fint = fint + a1*f(i)*drdi(i) + a2*f(i+1)*drdi(i+1)
  end do
  fint = fint - f(iend)*drdi(iend)/3.0d0

end subroutine
