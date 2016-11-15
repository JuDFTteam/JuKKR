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
  use bsrmm_mod, only: bsrMultPlan
  implicit none
  private
  public :: KKROperator, create, destroy, multiply

  !> Represents the operator/matrix (1 - \Delta T G_ref).
  type :: KKROperator
    integer :: lmmaxd
    type(SparseMatrixDescription) :: bsr_A, bsr_X!, bsr_B ToDo: introduce B with less entries than X
    double complex, allocatable :: mat_A(:,:,:,:) !< dim(fastBlockDim,slowBlockDim,bsr_A%nnzb,0:Lly)
    double complex, allocatable :: mat_B(:,:,:)   !< dim(fastBlockDim,slowBlockDim,bsr_B%nnzb)
    double complex, allocatable :: mat_X(:,:,:)   !< dim(fastBlockDim,slowBlockDim,bsr_X%nnzb)
    integer(kind=2), allocatable :: atom_indices(:) !< local truncation zone indices of the source atoms
    type(ClusterInfo), pointer :: cluster
    type(bsrMultPlan) :: plan
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

  subroutine create_KKROperator(self, cluster, lmmaxd, atom_indices, Lly)
    use bsrmm_mod, only: bsr_times_bsr ! planning
    use TEST_lcutoff_mod, only: lmax_a_array
    use fillKKRMatrix_mod, only: getKKRMatrixStructure, getKKRSolutionStructure

    type(KKROperator), intent(inout) :: self
    type(ClusterInfo), target, intent(in) :: cluster
    integer, intent(in) :: lmmaxd, Lly
    integer(kind=2), intent(in) :: atom_indices(:) !< local truncation zone indices of the source atoms

    integer :: nCols, nRows, nBlocks, nLloyd, ist

    self%cluster => cluster
    self%lmmaxd = lmmaxd

#ifndef __GFORTRAN__
    allocate(self%atom_indices, source=atom_indices) ! local truncation zone indices of the source atoms
#else
    allocate(self%atom_indices(size(atom_indices))) ! local truncation zone indices of the source atoms
    self%atom_indices = atom_indices ! copy
#endif

#ifndef NDEBUG
    deallocate(local_atom_indices, stat=ist) ! ignore status
    allocate(local_atom_indices(size(atom_indices)), stat=ist) ! see above
    local_atom_indices(:) = atom_indices(:) ! make a copy that we can use for DEBUG purposes by use KKROperator_mod, only: local_atom_indices
#endif


    ! create block sparse structure of matrix A
    call getKKRMatrixStructure(cluster%numn0, cluster%indn0, self%bsr_A)
    
    ! create block sparse structure of solution X
    call getKKRSolutionStructure(lmax_a_array, self%bsr_X)

    nRows = lmmaxd ! here we can introduce memory alignment
    if (self%bsr_X%nCols == 1) then
      nCols = lmmaxd*size(atom_indices) ! rectangluar shaped -- does not conform with a correct parallelization of truncation for num_local_atoms > 1
    else
      nCols = lmmaxd
    endif
    nBlocks = self%bsr_X%nnzb

    allocate(self%mat_B(nRows,nCols,nBlocks))
    allocate(self%mat_X(nRows,nCols,nBlocks))

    nRows = lmmaxd ! here we can introduce memory alignment
    nCols = lmmaxd
    nBlocks = self%bsr_A%nnzb
    nLloyd = min(max(0, Lly), 1) ! for the energy derivative needed in Lloyd''s formula

    allocate(self%mat_A(nRows,nCols,nBlocks,0:nLloyd)) ! allocate memory for the KKR operator

    ! plan the operation
    call bsr_times_bsr(self%plan, self%bsr_A%RowStart, self%bsr_A%ColIndex, shape(self%mat_A), self%bsr_X%RowStart, self%bsr_X%ColIndex, shape(self%mat_X))
    
  endsubroutine ! create


  subroutine destroy_KKROperator(self)
    use SparseMatrixDescription_mod, only: destroy
    type(KKROperator), intent(inout) :: self

    integer :: ist ! ignore status
    
    deallocate(self%mat_A, stat=ist)
    deallocate(self%mat_X, stat=ist)
    deallocate(self%mat_B, stat=ist)

    call destroy(self%bsr_A)
    call destroy(self%bsr_X)
!   call destroy(self%bsr_B)

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
