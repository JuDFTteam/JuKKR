!------------------------------------------------------------------------------
!> Storage for tetrahedra data.
module tetrahedra_common
  use shape_constants_mod
  implicit none

  save

  integer, dimension(:), allocatable :: ISIGNU  !< ??? sign for rotation sense ???
  integer, dimension(:), allocatable :: NTT     !< number of tetrahedra for polygon
  real(kind=DP), dimension(:), allocatable  ::RD       !< distances pyramid footpoint to edge
  real(kind=DP), dimension(:), allocatable  ::R0       !< foot points of perpendicular of polygons
  real(kind=DP), dimension(:), allocatable  ::FA       !< tetrahedron angle, phi-angle corresponding to 1st vertex  FA < FD < FB
  real(kind=DP), dimension(:), allocatable  ::FB       !< tetrahedron angle, phi-angle corresponding to 2nd vertex
  real(kind=DP), dimension(:), allocatable  ::FD       !< tetrahedron angle, phi-angle corresponding to foot point between 1st and 2nd vertex

  contains

  subroutine createtetra(nfaced, nvtotd)
    implicit none
    integer :: nfaced, nvtotd
    allocate(ISIGNU(NVTOTD))
    allocate(NTT(NFACED))
    allocate(RD(NVTOTD))
    allocate(R0(NFACED))
    allocate(FA(NVTOTD))
    allocate(FB(NVTOTD))
    allocate(FD(NVTOTD))
  end subroutine

  subroutine destroytetra()
    implicit none
    deallocate(ISIGNU)
    deallocate(NTT)
    deallocate(RD)
    deallocate(R0)
    deallocate(FA)
    deallocate(FB)
    deallocate(FD)

  end subroutine

end module
