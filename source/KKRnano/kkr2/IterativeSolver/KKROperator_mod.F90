module KKROperator_mod
  use OperatorT_mod

  use SparseMatrixDescription_mod
  use ClusterInfo_mod
  use MultScatData_mod

  implicit none

  type, extends(OperatorT) :: KKROperator
    private
    type(MultScatData), pointer :: ms => null()
    contains
    procedure :: apply  => apply_KKROperator
  end type

  contains

  !----------------------------------------------------------------------------
  !> Applies Operator on mat_X and returns result in mat_AX
  subroutine apply_KKROperator(self, mat_X, mat_AX)
    use vbrmv_mat_mod, only: multiply_vbr
    class(KKROperator) :: self
    double complex, intent(in)  :: mat_X(:,:)
    double complex, intent(out) :: mat_AX(:,:)

    ! perform sparse VBR matrix * dense matrix
    call multiply_vbr(self%ms%GLLH, mat_X, mat_AX, self%ms%sparse)
  end subroutine

  !---------------------------------------------------------------------------
  subroutine associate_ms_workspace(self, ms)
    class(KKROperator) :: self
    type(MultScatData), target :: ms

    if (.not. associated(self%ms)) self%ms => ms

  end subroutine

  !----------------------------------------------------------------------------
  ! Return a reference to the multiple scattering workspace.
  function get_ms_workspace(self) result(ms)
    class(KKROperator) :: self
    type(MultScatData), pointer :: ms

    ms => self%ms
  end function

end module
