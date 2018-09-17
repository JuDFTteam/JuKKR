module mod_csum
  use :: mod_datatypes, only: dp
  private :: dp

contains

  ! 19.10.95 *************************************************************
  function csum(n, v, iv)
    complex (kind=dp) :: csum
    ! **********************************************************************
    ! sum up the first N elements of the complex (kind=dp)
    ! array V(*) with a stepwidth of IV
    ! ----------------------------------------------------------------------
    ! .. Scalar Arguments ..
    integer :: iv, n
    ! ..
    ! .. Array Arguments ..
    complex (kind=dp) :: v(*)
    ! ..
    ! .. Local Scalars ..
    complex (kind=dp) :: vsum
    integer :: i, ibot, itop
    ! ..
    if (iv>=0) then
      ibot = 1
      itop = 1 + (n-1)*iv

    else
      ibot = 1 - (n-1)*iv
      itop = 1
    end if

    vsum = (0e0_dp, 0e0_dp)
    do i = ibot, itop, iv
      vsum = vsum + v(i)
    end do
    csum = vsum
    return
  end function csum

end module mod_csum
