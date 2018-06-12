function ssum(n, v, iv)
  use :: mod_datatypes, only: dp
  ! **********************************************************************
  ! sum up the first N elements of the real (kind=dp)
  ! array V(*) with a stepwidth of IV
  ! ----------------------------------------------------------------------
  implicit none
  real (kind=dp) :: ssum
  ! .. Scalar Arguments ..
  integer :: iv, n
  ! ..
  ! .. Array Arguments ..
  real (kind=dp) :: v(*)
  ! ..
  ! .. Local Scalars ..
  real (kind=dp) :: vsum
  integer :: i, ibot, itop
  ! ..
  if (iv>=0) then
    ibot = 1
    itop = 1 + (n-1)*iv

  else
    ibot = 1 - (n-1)*iv
    itop = 1
  end if

  vsum = 0.0e0_dp
  do i = ibot, itop, iv
    vsum = vsum + v(i)
  end do
  ssum = vsum
end function ssum
