module mod_csimpk

  private
  public :: csimpk

contains

  !-------------------------------------------------------------------------------
  !> Summary: Integration of complex function with extended 3-point Simpson
  !> Author: 
  !> Category: KKRhost, radial-grid, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> This subroutine does an integration up to rcut of an
  !> complex function cf with an extended 3-point-simpson:
  !>
  !> rcut
  !> cfint = { cf(r') dr'
  !> 0
  !>
  !> Modified for functions with kinks - at each kink the
  !> integration is restarted.
  !>
  !> @warning
  !> Input `cf` is destroyed!
  !> @endwarning
  !-------------------------------------------------------------------------------
  subroutine csimpk(cf, cfint, ipan, ircut, drdi)

    use :: mod_datatypes, only: dp
    use :: mod_csum, only: csum
    implicit none

    ! .. Scalar Arguments ..
    complex (kind=dp) :: cfint
    integer :: ipan
    ! ..
    ! .. Array Arguments ..
    complex (kind=dp) :: cf(*)
    real (kind=dp) :: drdi(*)
    integer :: ircut(0:ipan)
    ! ..
    ! .. Local Scalars ..
    real (kind=dp) :: a1, a2
    integer :: i, ien, ip, ist, n
    ! ..
    a1 = 4.0e0_dp/3.0e0_dp
    a2 = 2.0e0_dp/3.0e0_dp
    cfint = 0.0e0_dp

    do ip = 1, ipan

      ! ---> loop over kinks

      ist = ircut(ip-1) + 1
      ien = ircut(ip)

      do i = ist, ien
        cf(i) = cf(i)*drdi(i)
      end do

      if (mod(ien-ist,2)==0) then
        cfint = cfint + (cf(ist)-cf(ien))/3.0e0_dp
        ist = ist + 1
        n = (ien-ist+1)/2
      else
        ! ---> four point lagrange integration for the first step
        cfint = cfint + (9.0e0_dp*cf(ist)+19.0e0_dp*cf(ist+1)-5.0e0_dp*cf(ist+2)+cf(ist+3))/24.0e0_dp + (cf(ist+1)-cf(ien))/3.0e0_dp
        ist = ist + 2
        n = (ien-ist+1)/2
      end if

      ! ---> calculate with an extended 3-point-simpson

      cfint = cfint + a1*csum(n, cf(ist), 2) + a2*csum(n, cf(ist+1), 2)
    end do

  end subroutine csimpk

end module mod_csimpk
