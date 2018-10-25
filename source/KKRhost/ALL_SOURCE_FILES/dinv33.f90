module mod_dinv33
  
  private
  public :: dinv33

contains

  !-------------------------------------------------------------------------------
  !> Summary: Invert 3x3 matrix
  !> Author: 
  !> Category: KKRhost, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Inverts 3X3 matrix
  !> 
  !> Inputs:
  !>   matrix:input matrix
  !>   iopt  :if 0, usual inverse
  !>             1, transpose of inverse
  !>             2, 2*pi*inverse
  !>             3, 2*pi*transpose of inverse
  !>
  !> Outputs:
  !>   invers:as modified according to iopt
  !>   det   :determinant, (or det/2*pi if iopt=2,3)
  !>
  !> Remarks:
  !>  To generate reciprocal lattice vectors, call dinv33(plat,3,plat)
  !-------------------------------------------------------------------------------
  subroutine dinv33(matrix, iopt, invers, det)

    use :: mod_datatypes, only: dp
    use :: mod_cross, only: cross
    use :: mod_ddot1, only: ddot1
    use :: mod_dscal1, only: dscal1
    use :: mod_dswap1, only: dswap1
    use :: mod_constants, only: pi
    implicit none

    integer :: iopt
    real (kind=dp) :: matrix(3, 3), invers(3, 3), det
    ! Local parameters:
    integer :: i, j
    real (kind=dp) :: twopi
    parameter (twopi=2.e0_dp*pi)

    call cross(matrix(1,2), matrix(1,3), invers(1,1))
    call cross(matrix(1,3), matrix(1,1), invers(1,2))
    call cross(matrix(1,1), matrix(1,2), invers(1,3))
    det = ddot1(3, matrix, 1, invers, 1)
    if (iopt>=2) det = det/twopi
    if (mod(iopt,2)==0) then
      do i = 1, 3
        do j = i + 1, 3
          call dswap1(1, invers(i,j), 1, invers(j,i), 1)
        end do
      end do
    end if
    call dscal1(9, 1.e0_dp/det, invers, 1)
  end subroutine dinv33

end module mod_dinv33
