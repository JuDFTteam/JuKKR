module mod_rinit

contains

subroutine rinit(n, a)
  use :: mod_datatypes, only: dp
  ! **********************************************************************
  ! * Setting the first N values of a real (kind=dp) array A to zero   *
  ! **********************************************************************
  ! ..
  ! .. Arguments ..
  integer :: n
  real (kind=dp) :: a(*)
  ! ..
  ! .. Locals ..
  integer :: i, m, mp1
  real (kind=dp) :: dzero
  ! ..
  data dzero/0.0e0_dp/
  ! ..
  ! ..
  m = mod(n, 5)
  if (m/=0) then
    do i = 1, m
      a(i) = dzero
    end do
    if (n<5) return
  end if
  mp1 = m + 1
  do i = mp1, n, 5
    a(i) = dzero
    a(i+1) = dzero
    a(i+2) = dzero
    a(i+3) = dzero
    a(i+4) = dzero
  end do

end subroutine rinit

end module mod_rinit
