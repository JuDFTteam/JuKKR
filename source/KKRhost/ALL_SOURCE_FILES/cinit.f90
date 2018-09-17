module mod_cinit
  use :: mod_datatypes, only: dp
  private :: dp

contains

  subroutine cinit(n, a)
    ! **********************************************************************
    ! * Setting the first N values of a complex (kind=dp) array A to zero     *
    ! **********************************************************************
    ! ..
    ! .. Arguments ..
    integer :: n
    complex (kind=dp) :: a(*)
    ! ..
    ! .. Locals
    integer :: i, m, mp1
    complex (kind=dp) :: czero
    ! ..
    data czero/(0.0e0_dp, 0.0e0_dp)/
    ! ..
    m = mod(n, 5)
    if (m/=0) then
      do i = 1, m
        a(i) = czero
      end do
      if (n<5) return
    end if
    mp1 = m + 1
    do i = mp1, n, 5
      a(i) = czero
      a(i+1) = czero
      a(i+2) = czero
      a(i+3) = czero
      a(i+4) = czero
    end do

  end subroutine cinit

end module mod_cinit
