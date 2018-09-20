module mod_cinit
  
  private
  public :: cinit

contains

  !-------------------------------------------------------------------------------
  !> Summary: Initialize array elements to complex zero
  !> Author: 
  !> Category: KKRhost, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Setting the first N values of a complex (kind=dp) array A to zero
  !> @note
  !> This can easily be removed and replace by a single line
  !>  `a(:n) = czero`
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine cinit(n, a)
    use :: mod_datatypes, only: dp
    use :: mod_constants, only: czero
    implicit none
    ! ..
    ! .. Arguments ..
    integer :: n
    complex (kind=dp) :: a(*)
    ! ..
    ! .. Locals
    integer :: i, m, mp1
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
