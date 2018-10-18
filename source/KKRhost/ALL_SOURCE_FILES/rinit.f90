!------------------------------------------------------------------------------------
!> Summary: Setting the first `N` values of a `real (kind=dp)` array `A` to zero
!> Author: 
!> Setting the first `N` values of a `real (kind=dp)` array `A` to zero
!------------------------------------------------------------------------------------
!> @note Maybe it can be replaced by calls susch as `A(1:N)=0.0d0`
!> @endnote
!------------------------------------------------------------------------------------
module mod_rinit
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Setting the first `N` values of a `real (kind=dp)` array `A` to zero
  !> Author: 
  !> Category: numerical-tools, KKRhost
  !> Deprecated: False
  !> Setting the first `N` values of a `real (kind=dp)` array `A` to zero
  !-------------------------------------------------------------------------------
  !> @note Maybe it can be replaced by calls susch as `A(1:N)=0.0d0`
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine rinit(n, a)

    implicit none
    ! .. Arguments ..
    integer, intent(in) :: n !! Number of entries to set to zero
    real (kind=dp), dimension(*), intent(inout) :: a !! Array which entries will be set to zero
    ! ..
    ! .. Locals ..
    integer :: i, m, mp1
    real (kind=dp) :: dzero
    ! ..
    data dzero/0.0e0_dp/
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
