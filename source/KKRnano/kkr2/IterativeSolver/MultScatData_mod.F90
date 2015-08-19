!> Workspace for solution of the Dyson equation

module MultScatData_mod
  use SparseMatrixDescription_mod
  use ClusterInfo_mod
  implicit none

  type MultScatData
    type (SparseMatrixDescription) :: sparse

    double complex, dimension(:,:), allocatable :: mat_B
    double complex, dimension(:,:), allocatable :: mat_X
    double complex, allocatable :: EIKRP(:)
    double complex, allocatable :: EIKRM(:)
    double complex, allocatable, dimension(:) :: GLLH
    double complex, allocatable, dimension(:) :: DGLLH

    integer, allocatable :: atom_indices(:)
    integer :: lmmaxd
    integer :: naez
    type (ClusterInfo), pointer :: cluster_info

  end type

CONTAINS

!------------------------------------------------------------------------------
!> Create workspace for multiple scattering calculation.
subroutine createMultScatData(ms, cluster_info, lmmaxd, atom_indices)
  use TEST_lcutoff_mod, only: lmarray
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

  call createSparseMatrixDescription(ms%sparse, naez, sum_cluster)

  call getKKRMatrixStructure(lmarray, cluster_info%numn0_trc, &
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

  call destroySparseMatrixDescription(ms%sparse)

  deallocate(ms%atom_indices)
end subroutine

end module MultScatData_mod
