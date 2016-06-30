!> Workspace for solution of the Dyson equation

module MultScatData_mod
  use SparseMatrixDescription_mod, only: SparseMatrixDescription, getNNZ, &
  create, destroy
  use ClusterInfo_mod, only: ClusterInfo
  implicit none
  private
  public :: MultScatData, create, destroy

  type MultScatData
    type(SparseMatrixDescription) :: sparse
    double complex, allocatable :: GLLh(:)
    double complex, allocatable :: dGLLh(:)
    double complex, allocatable :: mat_B(:,:) ! ToDo: make it a sparse operator since it is mostly zero or an implicit action of subtracting mat_B
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

CONTAINS

!------------------------------------------------------------------------------
!> Create workspace for multiple scattering calculation.
subroutine createMultScatData(ms, cluster_info, lmmaxd, atom_indices)
  use TEST_lcutoff_mod, only: lmax_array
  use fillKKRMatrix_mod, only: getKKRMatrixStructure
  
  implicit none
  type (MultScatData), intent(inout) :: ms
  type (ClusterInfo), target  :: cluster_info
  integer, intent(in) :: lmmaxd
  integer, intent(in) :: atom_indices(:)

  integer :: sum_cluster
  integer :: naez
  integer :: naclsd

  ms%cluster_info => cluster_info

  sum_cluster = sum(cluster_info%numn0_trc)
  naez = size(cluster_info%indn0_trc, 1)
  naclsd = size(cluster_info%indn0_trc, 2)
  ms%lmmaxd = lmmaxd
  ms%naez = naez

  allocate(ms%atom_indices, source=atom_indices)

  call create(ms%sparse, naez, sum_cluster)

  call getKKRMatrixStructure(lmax_array, cluster_info%numn0_trc, &
                             cluster_info%indn0_trc, ms%sparse)

  allocate(ms%mat_B(ms%sparse%kvstr(naez+1)-1,LMMAXD * size(atom_indices)))
  allocate(ms%mat_X(ms%sparse%kvstr(naez+1)-1,LMMAXD * size(atom_indices)))

  ! allocate memory for sparse matrix
  allocate(ms%GLLH(getNNZ(ms%sparse)))
  allocate(ms%DGLLH(getNNZ(ms%sparse)))

  allocate(ms%eikrm(naclsd))
  allocate(ms%eikrp(naclsd))
end subroutine

!------------------------------------------------------------------------------
!> Destroys workspace for multiple scattering calculation.
subroutine destroyMultScatData(ms)
  implicit none
  type (MultScatData), intent(inout) :: ms

  deallocate(ms%eikrp)
  deallocate(ms%eikrm)
  deallocate(ms%GLLH)
  deallocate(ms%DGLLH)
  deallocate(ms%mat_X)
  deallocate(ms%mat_B)

  call destroy(ms%sparse)

  deallocate(ms%atom_indices)
end subroutine

end module MultScatData_mod
