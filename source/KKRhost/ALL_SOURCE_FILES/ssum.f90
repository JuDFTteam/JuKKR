double precision function ssum(n, v, iv)
! **********************************************************************
!        sum up the first N elements of the double precision
!        array V(*) with a stepwidth of IV
! ----------------------------------------------------------------------
  implicit none
!.. Scalar Arguments ..
  integer :: iv, n
!..
!.. Array Arguments ..
  double precision :: v(*)
!..
!.. Local Scalars ..
  double precision :: vsum
  integer :: i, ibot, itop
!..
  if (iv>=0) then
    ibot = 1
    itop = 1 + (n-1)*iv

  else
    ibot = 1 - (n-1)*iv
    itop = 1
  end if

  vsum = 0.0d0
  do i = ibot, itop, iv
    vsum = vsum + v(i)
  end do
  ssum = vsum
end function
