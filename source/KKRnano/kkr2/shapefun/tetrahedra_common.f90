!------------------------------------------------------------------------------
!> Storage for tetrahedra data.
module tetrahedra_common
  use shape_constants_mod, only: DP
  implicit none
  public

  integer, allocatable :: ISIGNU(:)   !< ??? sign for rotation sense ???
  integer, allocatable :: NTT(:)      !< number of tetrahedra for polygon
  real(kind=DP), allocatable :: RD(:) !< distances pyramid footpoint to edge
  real(kind=DP), allocatable :: R0(:) !< foot points of perpendicular of polygons
  real(kind=DP), allocatable :: FA(:) !< tetrahedron angle, phi-angle corresponding to 1st vertex  FA < FD < FB
  real(kind=DP), allocatable :: FB(:) !< tetrahedron angle, phi-angle corresponding to 2nd vertex
  real(kind=DP), allocatable :: FD(:) !< tetrahedron angle, phi-angle corresponding to foot point between 1st and 2nd vertex

  contains

  subroutine createtetra(nfaced, nvtotd)
    integer, intent(in) :: nfaced, nvtotd
    allocate(ISIGNU(NVTOTD))
    allocate(NTT(NFACED))
    allocate(RD(NVTOTD))
    allocate(R0(NFACED))
    allocate(FA(NVTOTD))
    allocate(FB(NVTOTD))
    allocate(FD(NVTOTD))
  end subroutine

  subroutine destroytetra()
    deallocate(ISIGNU)
    deallocate(NTT)
    deallocate(RD)
    deallocate(R0)
    deallocate(FA)
    deallocate(FB)
    deallocate(FD)
  end subroutine

end module
