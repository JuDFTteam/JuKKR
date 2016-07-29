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
  use SparseMatrixDescription_mod, only: SparseMatrixDescription
  use ClusterInfo_mod, only: ClusterInfo
  implicit none
  private
  public :: KKROperator, create, destroy, multiply

  !> Represents the operator/matrix (1 - \Delta T G_ref).
  type :: KKROperator
    integer :: lmmaxd
    integer :: naez
    type(SparseMatrixDescription) :: sparse
    double complex, allocatable :: GLLh(:)
    double complex, allocatable :: dGLLh(:)
    double complex, allocatable :: mat_B(:,:) ! ToDo: make it a sparse operator since it is mostly zero or an implicit action of subtracting mat_B
    double complex, allocatable :: mat_X(:,:)
    integer(kind=2), allocatable :: atom_indices(:) !< a copy of the atom indices
    type(ClusterInfo), pointer :: cluster_info
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

  subroutine create_KKROperator(self, cluster_info, lmmaxd, atom_indices)
    use TEST_lcutoff_mod, only: lmax_array
    use fillKKRMatrix_mod, only: getKKRMatrixStructure
    use SparseMatrixDescription_mod, only: create, getNNZ, getNrows

    type(KKROperator), intent(inout) :: self
    type(ClusterInfo), target, intent(in) :: cluster_info
    integer, intent(in) :: lmmaxd
    integer(kind=2), intent(in) :: atom_indices(:)

    integer :: sum_cluster, nCols, nRows

    self%cluster_info => cluster_info

    sum_cluster = sum(cluster_info%numn0_trc)
    self%naez = size(cluster_info%indn0_trc, 2)
!   naclsd = size(cluster_info%indn0_trc, 1) ! not used
    self%lmmaxd = lmmaxd

    allocate(self%atom_indices, source=atom_indices) ! local truncation zone indices of the source atoms

    call create(self%sparse, self%naez, sum_cluster)

    call getKKRMatrixStructure(lmax_array, cluster_info%numn0_trc, cluster_info%indn0_trc, self%sparse)

    nRows = getNrows(self%sparse)!, self%naez)
    nCols = lmmaxd*size(atom_indices)
    
    allocate(self%mat_B(nRows,nCols))
    allocate(self%mat_X(nRows,nCols))
    allocate(self%GLLh(getNNZ(self%sparse))) ! allocate memory for sparse matrix
    allocate(self%dGLLh(getNNZ(self%sparse))) ! allocate memory for derivative
    
  endsubroutine ! create


  subroutine destroy_KKROperator(self)
    use SparseMatrixDescription_mod, only: destroy
    type(KKROperator), intent(inout) :: self

    integer :: ist ! ignore status
    
    deallocate(self%GLLh, stat=ist)
    deallocate(self%dGLLh, stat=ist)
    deallocate(self%mat_X, stat=ist)
    deallocate(self%mat_B, stat=ist)

    call destroy(self%sparse)

    deallocate(self%atom_indices, stat=ist)
    nullify(self%cluster_info)
  endsubroutine ! destroy

  !----------------------------------------------------------------------------
  !> Applies Operator on mat_X and returns result in mat_AX.
  subroutine multiply_KKROperator(self, mat_X, mat_AX)
    type(KKROperator) :: self
    double complex, intent(in)  :: mat_X(:,:)
    double complex, intent(out) :: mat_AX(:,:)

    ! perform sparse VBR matrix * dense matrix
    call multiply_vbr(self%GLLH, mat_X, mat_AX, self%sparse)
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
