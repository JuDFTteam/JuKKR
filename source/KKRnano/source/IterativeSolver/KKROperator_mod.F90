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
  use bsrmm_mod, only: bsrMultPlan, destroy
  implicit none
  private
  public :: KKROperator, create, destroy, multiply
  public :: make_multiplication_plan

  !> Represents the operator/matrix (1 - \Delta T G_ref).
  type :: KKROperator
    type(SparseMatrixDescription) :: bsr_A, bsr_X, bsr_B ! ToDo: introduce B with less entries than X
    double complex, allocatable :: mat_A(:,:,:,:) !< dim(fastBlockDim,slowBlockDim,bsr_A%nnzb,0:Lly)
    double complex, allocatable :: mat_B(:,:,:)   !< dim(fastBlockDim,slowBlockDim,bsr_B%nnzb)
    double complex, allocatable :: mat_X(:,:,:)   !< dim(fastBlockDim,slowBlockDim,bsr_X%nnzb)
    integer(kind=2), allocatable :: atom_indices(:) !< local truncation zone indices of the source atoms
    type(ClusterInfo), pointer :: cluster
    type(bsrMultPlan) :: plan ! this plan allows to achieve performance when multiplying mat_A to mat_X
    integer, allocatable :: B_subset_of_X(:) !< dim(B%nnzb) this plan allows to subtract mat_B from mat_X
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

#ifndef NDEBUG
  integer(kind=2), allocatable, protected, public :: local_atom_indices(:) ! see above
#endif

  contains

  subroutine create_KKROperator(self, cluster, lmsd, atom_indices, Lly)
    use bsrmm_mod, only: bsr_times_bsr ! planning
    use Truncation_mod, only: lmax_a_array
    use fillKKRMatrix_mod, only: getKKRMatrixStructure, getKKRSolutionStructure, getRightHandSideStructure
    use SparseMatrixDescription_mod, only: SparseMatrixDescription, subset

    type(KKROperator), intent(inout) :: self
    type(ClusterInfo), target, intent(in) :: cluster
    integer, intent(in) :: lmsd, Lly
    integer(kind=2), intent(in) :: atom_indices(:) !< local truncation zone indices of the source atoms

    integer :: nCols, nRows, nLloyd, ist, nRHS_group

    self%cluster => cluster

    allocate(self%atom_indices(size(atom_indices))) ! local truncation zone indices of the source atoms
    self%atom_indices = atom_indices ! copy

#ifndef NDEBUG
    deallocate(local_atom_indices, stat=ist) ! ignore status
    allocate(local_atom_indices(size(atom_indices)), stat=ist) ! see above
    local_atom_indices(:) = atom_indices(:) ! make a copy that we can use for DEBUG purposes by use KKROperator_mod, only: local_atom_indices
#endif

    ! create block sparse structure of matrix A
    call getKKRMatrixStructure(self%bsr_A, cluster%numn0, cluster%indn0)
    
    ! create block sparse structure of solution X
    call getKKRSolutionStructure(self%bsr_X, lmax_a_array)

    nRHS_group = 1; if (self%bsr_X%nCols == 1) nRHS_group = size(atom_indices) ! rectangluar shaped -- does not conform with a correct parallelization of truncation for num_local_atoms > 1
    
    ! create a very sparse right hand side descriptor
    call getRightHandSideStructure(self%bsr_B, atom_indices, nRHS_group)
    allocate(self%B_subset_of_X(self%bsr_B%nnzb))
    ist = subset(set=self%bsr_X, sub=self%bsr_B, list=self%B_subset_of_X)
    if (ist /= 0) stop 'KKROperator: creation of subset list failed!'

#ifdef  TRANSPOSE_TO_ROW_MAJOR
#define NCOLS nRows
#define NROWS nCols
    write(*,*) "warning: transposition of KKROperator-quantities"
#endif

    nRows = lmsd ! here we can introduce memory alignment
    nCols = lmsd*nRHS_group

    allocate(self%mat_B(NROWS,NCOLS,self%bsr_B%nnzb))
    allocate(self%mat_X(NROWS,NCOLS,self%bsr_X%nnzb))

    nCols = lmsd
    nLloyd = min(max(0, Lly), 1) ! for the energy derivative needed in Lloyd''s formula

    allocate(self%mat_A(NROWS,NCOLS,self%bsr_A%nnzb,0:nLloyd)) ! allocate memory for the KKR operator

    call make_multiplication_plan(self)
    
  endsubroutine ! create
    
  subroutine make_multiplication_plan(self)
    use bsrmm_mod, only: bsr_times_bsr ! planning
    type(KKROperator), intent(inout) :: self

    ! plan the operation
    call bsr_times_bsr(self%plan, self%bsr_A%RowStart, self%bsr_A%ColIndex, &
    shape(self%mat_A), self%bsr_X%RowStart, self%bsr_X%ColIndex, shape(self%mat_X))

  endsubroutine ! make_plan


  elemental subroutine destroy_KKROperator(self)
    use SparseMatrixDescription_mod, only: destroy
    type(KKROperator), intent(inout) :: self

    integer :: ist ! ignore status
    
    deallocate(self%mat_A, stat=ist)
    deallocate(self%mat_X, stat=ist)
    deallocate(self%mat_B, stat=ist)

    call destroy(self%bsr_A)
    call destroy(self%bsr_X)
    call destroy(self%bsr_B)

    call destroy(self%plan)
    deallocate(self%B_subset_of_X, stat=ist)
    
    deallocate(self%atom_indices, stat=ist)
    nullify(self%cluster)
  endsubroutine ! destroy

  !----------------------------------------------------------------------------
  !> Applies Operator on mat_X and returns result in mat_AX.
  subroutine multiply_KKROperator(self, mat_X, mat_AX, nFlop)
    use bsrmm_mod, only: bsr_times_bsr

    type(KKROperator), intent(in) :: self
    double complex, intent(in)  :: mat_X(:,:,:)
    double complex, intent(out) :: mat_AX(:,:,:)
    integer(kind=8), intent(inout) :: nFlop

    if (self%plan%mtasks > 0) then

      ! perform BSR matrix * BSR matrix: bsr_times_bsr(plan, Y, A, X) according to plan
      call bsr_times_bsr(self%plan, mat_AX, self%mat_A(:,:,:,0), mat_X)
      nFlop = nFlop + self%plan%nFlop

    else
    
      ! perform BSR matrix * BSR matrix: bsr_times_bsr(Y, ia, ja, A, ix, jx, X, nFlop) spontaneously
      call bsr_times_bsr(mat_AX, self%bsr_A%RowStart, self%bsr_A%ColIndex, self%mat_A(:,:,:,0), self%bsr_X%RowStart, self%bsr_X%ColIndex, mat_X, nFlop)
      
    endif
    
  endsubroutine ! apply
  
endmodule ! KKROperator_mod
