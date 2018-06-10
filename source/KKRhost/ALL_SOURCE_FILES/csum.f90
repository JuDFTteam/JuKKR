! 19.10.95 *************************************************************
complex *16 function csum(n, v, iv)
! **********************************************************************
!        sum up the first N elements of the double complex
!        array V(*) with a stepwidth of IV
! ----------------------------------------------------------------------
!.. Scalar Arguments ..
  integer :: iv, n
!..
!.. Array Arguments ..
  double complex :: v(*)
!..
!.. Local Scalars ..
  double complex :: vsum
  integer :: i, ibot, itop
!..
  if (iv>=0) then
    ibot = 1
    itop = 1 + (n-1)*iv

  else
    ibot = 1 - (n-1)*iv
    itop = 1
  end if

  vsum = (0d0, 0d0)
  do i = ibot, itop, iv
    vsum = vsum + v(i)
  end do
  csum = vsum
  return
end function
