!> Workspace for solution of the Dyson equation

module MultScatData_mod
  use SparseMatrixDescription_mod, only: SparseMatrixDescription
  use ClusterInfo_mod, only: ClusterInfo
  implicit none
  private
  public :: MultScatData, create, destroy

  type MultScatData
    integer :: lmmaxd
    integer :: naez
    type(SparseMatrixDescription) :: sparse
    double complex, allocatable :: GLLh(:)
    double complex, allocatable :: mat_B(:,:)
    double complex, allocatable :: mat_X(:,:)
    integer, allocatable :: atom_indices(:) !< a copy of the atom indices
    type(ClusterInfo), pointer :: cluster_info
  endtype
  
  interface create
    module procedure createMultScatData
  endinterface
  
  interface destroy
    module procedure destroyMultScatData
  endinterface
  
  contains

  !------------------------------------------------------------------------------
  !> Create workspace for multiple scattering calculation.
  subroutine createMultScatData(ms, cluster_info, lmmaxd, atom_indices)
    use TEST_lcutoff_mod, only: lm_array
    use fillKKRMatrix_mod, only: getKKRMatrixStructure
    use SparseMatrixDescription_mod, only: create, getNNZ, getNrows

    type(MultScatData), intent(inout) :: ms
    type(ClusterInfo), target, intent(in) :: cluster_info
    integer, intent(in) :: lmmaxd
    integer, intent(in) :: atom_indices(:)

    integer :: sum_cluster, naez, naclsd, nCols, nRows

    ms%cluster_info => cluster_info

    sum_cluster = sum(cluster_info%numn0_trc)
    naez   = size(cluster_info%indn0_trc, 1)
    naclsd = size(cluster_info%indn0_trc, 2) ! not used
    ms%lmmaxd = lmmaxd
    ms%naez = naez

    allocate(ms%atom_indices, source=atom_indices)

    call create(ms%sparse, naez, sum_cluster)

    call getKKRMatrixStructure(lm_array, cluster_info%numn0_trc, cluster_info%indn0_trc, ms%sparse)

    nRows = getNrows(ms%sparse, naez)
    nCols = lmmaxd*size(atom_indices)
    
    allocate(ms%mat_B(nRows,nCols))
    allocate(ms%mat_X(nRows,nCols))
    allocate(ms%GLLh(getNNZ(ms%sparse))) ! allocate memory for sparse matrix
    
  endsubroutine ! create

  !------------------------------------------------------------------------------
  !> Destroys workspace for multiple scattering calculation.
  elemental subroutine destroyMultScatData(ms)
    use SparseMatrixDescription_mod, only: destroy
    type(MultScatData), intent(inout) :: ms

    integer :: ist ! ignore status
    
    deallocate(ms%GLLh, stat=ist)
    deallocate(ms%mat_X, stat=ist)
    deallocate(ms%mat_B, stat=ist)

    call destroy(ms%sparse)

    deallocate(ms%atom_indices, stat=ist)
    nullify(ms%cluster_info)
  endsubroutine ! destroy

endmodule ! MultScatData_mod
