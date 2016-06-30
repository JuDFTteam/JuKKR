!> Module to apply KKR coefficient matrix.
!>
!> The setup of the KKR coefficient matrix is rather complicated.
!> Therefore the needed data is stored in a 'MultScatData' struct.
!>
!> *) One has to get a reference (pointer) to the MultScatData workspace by
!>    using 'get_ms_workspace' and set up the workspace properly (routine kkrmat01)
!> *) Then one can apply the KKR coefficient matrix on any dense matrix using
!>    'apply'

module KKROperator_mod
  use MultScatData_mod, only: MultScatData
  use OperatorT_mod, only: OperatorT
  implicit none
  private
  public :: KKROperator, create, destroy, multiply

  !> Represents the operator/matrix (1 - \Delta T G_ref).
  type, extends(OperatorT) :: KKROperator
    type(MultScatData) :: ms
    contains
      procedure :: apply => multiply_KKROperator
  endtype

  interface create
    module procedure create_KKROperator
  endinterface
  
  interface destroy
    module procedure destroy_KKROperator
  endinterface
  
  interface multiply
    module procedure multiply_KKROperator
  endinterface

  contains

  subroutine create_KKROperator(self)
    class(KKROperator) :: self
  endsubroutine ! create

  subroutine destroy_KKROperator(self)
    use MultScatData_mod, only: destroy
    class(KKROperator) :: self

    call destroy(self%ms)
  endsubroutine ! destroy

  !----------------------------------------------------------------------------
  !> Applies Operator on mat_X and returns result in mat_AX.
  subroutine multiply_KKROperator(self, mat_X, mat_AX)
    class(KKROperator) :: self
    double complex, intent(in)  :: mat_X(:,:)
    double complex, intent(out) :: mat_AX(:,:)

    ! perform sparse VBR matrix * dense matrix
    call multiply_vbr(self%ms%GLLH, mat_X, mat_AX, self%ms%sparse)
  endsubroutine ! apply

  subroutine multiply_vbr(A, x, Ax, sparse)
    use vbrmv_mat_mod, only: vbrmv_mat
    use SparseMatrixDescription_mod, only: SparseMatrixDescription
    double complex, intent(in)  :: A(:), x(:,:)
    double complex, intent(out) :: Ax(:,:)
    type(SparseMatrixDescription), intent(in) :: sparse

    call vbrmv_mat(sparse%blk_nrows, sparse%ia, sparse%ja, sparse%ka, &
                   A, sparse%kvstr, x, Ax, &
                   sparse%max_blockdim, sparse%max_blocks_per_row)

  endsubroutine ! multiply_vbr

endmodule ! KKROperator_mod
