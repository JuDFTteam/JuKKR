!------------------------------------------------------------------------------
!> Storage for tetrahedra data.
module tetrahedra_common
  use shape_constants_mod, only: dp
  implicit none
  public

  integer, allocatable :: ntt(:)      !< number of tetrahedra for polygon
  real(kind=dp), allocatable :: r0(:) !< foot points of perpendicular of polygons
  
  real(kind=dp), allocatable :: rd(:) !< distances pyramid footpoint to edge
  real(kind=dp), allocatable :: fa(:) !< tetrahedron angle, phi-angle corresponding to 1st vertex  fa < fd < fb
  real(kind=dp), allocatable :: fb(:) !< tetrahedron angle, phi-angle corresponding to 2nd vertex
  real(kind=dp), allocatable :: fd(:) !< tetrahedron angle, phi-angle corresponding to foot point between 1st and 2nd vertex
  integer, allocatable :: isignu(:)   !< ??? sign for rotation sense ???

  contains

  integer function createtetra(nfaced, nvtotd) result(ist)
    integer, intent(in) :: nfaced, nvtotd
    allocate(ntt(nfaced), r0(nfaced), &
    rd(nvtotd), fa(nvtotd), fb(nvtotd), fd(nvtotd), isignu(nvtotd), stat=ist)
  endfunction ! create

  integer function destroytetra() result(ist)
    deallocate(ntt, r0, rd, fa, fb, fd, isignu, stat=ist)
  endfunction ! destroy

endmodule ! tetrahedra_common
