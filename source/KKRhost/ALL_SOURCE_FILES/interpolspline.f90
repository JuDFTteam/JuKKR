module mod_interpolspline
  
  private
  public :: interpolspline

contains

  !-------------------------------------------------------------------------------
  !> Summary: Spline interpolation of potentail to new radial mesh
  !> Author: 
  !> Category: KKRhost, potential, radial-grid
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !-------------------------------------------------------------------------------
  subroutine interpolspline(rmesh, rmeshnew, vpot, vpotnew, nrmax, nrmaxnew)

    use :: mod_datatypes, only: dp
    use :: mod_spline_real, only: spline_real
    use :: mod_splint_real, only: splint_real
    implicit none

    ! interface
    integer :: nrmax
    integer :: nrmaxnew
    real (kind=dp) :: rmesh(nrmax)
    real (kind=dp) :: rmeshnew(nrmaxnew)
    real (kind=dp) :: vpot(nrmax)
    real (kind=dp) :: vpotnew(nrmaxnew)
    ! local
    real (kind=dp) :: maxa
    real (kind=dp) :: spline(nrmax)
    real (kind=dp) :: parsum, parsumderiv, r0
    integer :: ir

    maxa = 1.e35_dp
    call spline_real(nrmax, rmesh, vpot, nrmax, maxa, maxa, spline)
    ! CALL SPLINE(IRMDJJ,R,VM2Z,NR,maxa,maxa,VM2ZB)

    do ir = 1, nrmaxnew
      r0 = rmeshnew(ir)
      call splint_real(rmesh, vpot, spline, nrmax, r0, parsum, parsumderiv)
      vpotnew(ir) = parsum
    end do
  end subroutine interpolspline

end module mod_interpolspline
