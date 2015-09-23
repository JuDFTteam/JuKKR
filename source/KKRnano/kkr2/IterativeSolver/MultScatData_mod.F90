!> Workspace for solution of the Dyson equation

module MultScatData_mod
  use SparseMatrixDescription_mod, only: SparseMatrixDescription
  use ClusterInfo_mod, only: ClusterInfo
  implicit none
  private
  public :: MultScatData, create, destroy
  public :: createMultScatData, destroyMultScatData ! deprecated

  type MultScatData
    type(SparseMatrixDescription) :: sparse

    double complex, allocatable :: mat_B(:,:)
    double complex, allocatable :: mat_X(:,:)
    double complex, allocatable :: EIKRP(:)
    double complex, allocatable :: EIKRM(:)
    double complex, allocatable :: GLLH(:)

    integer, allocatable :: atom_indices(:)
    integer :: lmmaxd
    integer :: naez
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
    use TEST_lcutoff_mod, only: lmarray
    use fillKKRMatrix_mod, only: getKKRMatrixStructure
    use SparseMatrixDescription_mod, only: createSparseMatrixDescription, getNNZ

    type(MultScatData), intent(inout) :: ms
    type(ClusterInfo), target, intent(in) :: cluster_info
    integer, intent(in) :: lmmaxd
    integer, intent(in) :: atom_indices(:)

    integer :: sum_cluster, naez, naclsd

    ms%cluster_info => cluster_info

    sum_cluster = sum(cluster_info%numn0_trc)
    naez = size(cluster_info%indn0_trc, 1)
    naclsd = size(cluster_info%indn0_trc, 2)
    ms%lmmaxd = lmmaxd
    ms%naez = naez

    allocate(ms%atom_indices, source=atom_indices)

    call createSparseMatrixDescription(ms%sparse, naez, sum_cluster)

    call getKKRMatrixStructure(lmarray, cluster_info%numn0_trc, cluster_info%indn0_trc, ms%sparse)

    allocate(ms%mat_B(ms%sparse%kvstr(naez+1)-1,LMMAXD*size(atom_indices)))
    allocate(ms%mat_X(ms%sparse%kvstr(naez+1)-1,LMMAXD*size(atom_indices)))

    allocate(ms%GLLH(getNNZ(ms%sparse))) ! allocate memory for sparse matrix

    allocate(ms%eikrm(naclsd))
    allocate(ms%eikrp(naclsd))
  endsubroutine ! create

  !------------------------------------------------------------------------------
  !> Destroys workspace for multiple scattering calculation.
  subroutine destroyMultScatData(ms)
    use SparseMatrixDescription_mod, only: destroySparseMatrixDescription
    type(MultScatData), intent(inout) :: ms

    integer :: ist ! ignore status
    
    deallocate(ms%eikrp, stat=ist)
    deallocate(ms%eikrm, stat=ist)
    deallocate(ms%GLLH, stat=ist)
    deallocate(ms%mat_X, stat=ist)
    deallocate(ms%mat_B, stat=ist)

    call destroySparseMatrixDescription(ms%sparse)

    deallocate(ms%atom_indices, stat=ist)
  endsubroutine ! destroy

endmodule MultScatData_mod
