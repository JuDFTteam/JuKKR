!------------------------------------------------------------------------------------
 !> Summary: Checks if a set of vectors are lattice vectors
 !> Author: 
!> Checks if a set of vectors are lattice vectors
!------------------------------------------------------------------------------------
module mod_latvec
  use :: mod_datatypes, only: dp
  private :: dp

contains
    !-------------------------------------------------------------------------------
    !> Summary: Checks if a set of vectors are lattice vectors
    !> Author: 
    !> Category: geometry, k-points, KKRhost 
    !> Deprecated: False 
    !> Checks if a set of vectors are lattice vectors
    !-------------------------------------------------------------------------------
  logical function latvec(n, qlat, vec)
    ! - Checks if a set of vectors are lattice vectors
    ! ----------------------------------------------------------------------
    ! i Inputs:
    ! i   n     :number of vectors
    ! i   qlat  :primitive translation vectors in reciprocal space
    ! i   vec   :double-precision vector
    ! o Outputs:
    ! o   latvec:.true. if all vectors are lattice vectors
    ! r Remarks:
    ! ----------------------------------------------------------------------
    implicit none
    ! Passed parameters:
    integer, intent(in) :: n !! Number of vectors
    real (kind=dp), dimension(3,*), intent(in) :: vec   !! Double precision vector
    real (kind=dp), dimension(3,*), intent(in) :: qlat  !! Primitive translation vectors in reciprocal space
    ! Local parameters:
    integer :: i, m
    real (kind=dp) :: tol, vdiff
    parameter (tol=1.e-3_dp)
    ! Common block:
    ! Intrinsic functions:
    intrinsic :: abs, anint

    latvec = .false.
    do i = 1, n
      do m = 1, 3
        vdiff = vec(1, i)*qlat(1, m) + vec(2, i)*qlat(2, m) + vec(3, i)*qlat(3, m)
        vdiff = abs(vdiff-anint(vdiff))
        if (vdiff>tol) return
      end do
    end do
    latvec = .true.
  end function latvec

end module mod_latvec
