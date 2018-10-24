module mod_cross

  private
  public :: cross

contains

  !-------------------------------------------------------------------------------
  !> Summary: Computes cross product of two vectors
  !> Author: 
  !> Category: KKRhost, undefined
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Cross product cr12 = X1 cross X2
  !> 
  !> Inputs:
  !>   x1    :first vector to multiply
  !>   x2    :second vector to multiply
  !>
  !> Outputs:
  !>   cr12  :cross product
  !>
  !> @note
  !> the `crospr` subroutine is a duplication of this and should be removed
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine cross(x1, x2, cr12)

    use :: mod_datatypes, only: dp
    implicit none
    ! Passed parameters:
    real (kind=dp) :: x1(3)
    real (kind=dp) :: x2(3)
    real (kind=dp) :: cr12(3)

    cr12(1) = x1(2)*x2(3) - x1(3)*x2(2)
    cr12(2) = x1(3)*x2(1) - x1(1)*x2(3)
    cr12(3) = x1(1)*x2(2) - x1(2)*x2(1)

  end subroutine cross

end module mod_cross
