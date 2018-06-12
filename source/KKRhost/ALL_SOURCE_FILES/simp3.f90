subroutine simp3(f, fint, istart, iend, drdi)
  use :: mod_datatypes, only: dp
  ! -----------------------------------------------------------------------
  ! this subroutine does an integration from istart to iend of
  ! the real function f with an extended 3-point-simpson :

  ! r(istart)

  ! fint = { f(r') dr'

  ! r(iend)

  ! -----------------------------------------------------------------------
  ! .. Scalar Arguments ..
  real (kind=dp) :: fint
  integer :: iend, istart
  ! ..
  ! .. Array Arguments ..
  real (kind=dp) :: drdi(*), f(*)
  ! ..
  ! .. Local Scalars ..
  real (kind=dp) :: a1, a2
  integer :: i, ist
  ! ..
  ! .. Intrinsic Functions ..
  intrinsic :: mod
  ! ..
  a1 = 4.0e0_dp/3.0e0_dp
  a2 = 2.0e0_dp/3.0e0_dp

  ! ---> initialize fint

  if (mod(iend-istart,2)==0) then
    fint = f(istart)*drdi(istart)/3.0e0_dp
    ist = istart + 1

  else
    fint = (f(istart+3)*drdi(istart+3)-5.0e0_dp*f(istart+2)*drdi(istart+2)+ &
      19.0e0_dp*f(istart+1)*drdi(istart+1)+9.0e0_dp*f(istart)*drdi(istart))/ &
      24.0e0_dp + f(istart+1)*drdi(istart+1)/3.0e0_dp
    ist = istart + 2
  end if

  ! ---> calculate with an extended 3-point-simpson

  do i = ist, iend - 1, 2
    fint = fint + a1*f(i)*drdi(i) + a2*f(i+1)*drdi(i+1)
  end do
  fint = fint - f(iend)*drdi(iend)/3.0e0_dp

end subroutine simp3
