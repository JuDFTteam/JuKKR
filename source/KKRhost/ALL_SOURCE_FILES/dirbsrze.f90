module mod_dirbsrze

  private
  public :: dirbsrze

contains

  !-------------------------------------------------------------------------------
  !> Summary: Radial function support for Burlisch-Stoer method
  !> Author: 
  !> Category: KKRhost, dirac, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Diagonal rational function extrapolation to support the
  !> Burlisch-Stoer method
  !>
  !> See: numerical recipes chapter 15.4
  !-------------------------------------------------------------------------------
  subroutine dirbsrze(iest, xest, yest, yz, dy, nv, nuse)

    use :: mod_datatypes, only: dp
    implicit none

    ! PARAMETER definitions
    integer :: ncfmax, iseqmax, nusemax
    parameter (ncfmax=8, iseqmax=30, nusemax=7)

    ! Dummy arguments
    integer :: iest, nuse, nv
    real (kind=dp) :: xest
    complex (kind=dp) :: dy(nv), yest(nv), yz(nv)

    ! Local variables
    complex (kind=dp) :: b, b1, c, d(ncfmax, nusemax), ddy, v, yy
    real (kind=dp) :: fx(nusemax), x(iseqmax)
    integer :: j, k, m1
    save :: b, b1, c, d, ddy, fx, j, k, m1, v, x, yy

    x(iest) = xest
    if (iest==1) then
      do j = 1, nv
        yz(j) = yest(j)
        d(j, 1) = yest(j)
        dy(j) = yest(j)
      end do
    else
      m1 = min(iest, nuse)
      do k = 1, m1 - 1
        fx(k+1) = x(iest-k)/xest
      end do
      do j = 1, nv
        yy = yest(j)
        v = d(j, 1)
        c = yy
        d(j, 1) = yy
        do k = 2, m1
          b1 = fx(k)*v
          b = b1 - c
          if (b/=0._dp) then
            b = (c-v)/b
            ddy = c*b
            c = b1*b
          else
            ddy = v
          end if
          v = d(j, k)
          d(j, k) = ddy
          yy = yy + ddy
        end do
        dy(j) = ddy
        yz(j) = yy
      end do
    end if
  end subroutine dirbsrze

end module mod_dirbsrze
