!------------------------------------------------------------------------------
!> Storage for polygon face data and additional information on vertices
module PolygonFaces_mod
  implicit none
  public

!   integer, allocatable :: ntt(:)      !< number of tetrahedra for polygon
!   double precision, allocatable :: r0(:) !< foot points of perpendicular of polygons
!   double precision, allocatable :: alpha(:)    !< Euler angles alpha to rotate faces perpendicular to z-axis
!   double precision, allocatable :: beta(:)     !< Euler angles beta
!   double precision, allocatable :: gamma(:)    !< Euler angles gamma
  
  type PolygonFace
    integer :: ntt
    double precision :: r0
    double precision :: euler(1:3) ! former angles alpha beta gamma
  endtype
  
  type(PolygonFace), allocatable :: face(:)
  
  double precision, allocatable :: rd(:) !< distances pyramid footpoint to edge
  double precision, allocatable :: fa(:) !< tetrahedron angle, phi-angle corresponding to 1st vertex  fa < fd < fb
  double precision, allocatable :: fb(:) !< tetrahedron angle, phi-angle corresponding to 2nd vertex
  double precision, allocatable :: fd(:) !< tetrahedron angle, phi-angle corresponding to foot point between 1st and 2nd vertex
  integer, allocatable :: isignu(:)   !< ??? sign for rotation sense ???
  
  contains

  integer function createtetra(nfaced, nvtotd) result(ist)
    integer, intent(in) :: nfaced, nvtotd
    allocate(&!ntt(nfaced), r0(nfaced), alpha(nfaced), beta(nfaced), gamma(nfaced), &
      face(nfaced), &
      rd(nvtotd), fa(nvtotd), fb(nvtotd), fd(nvtotd), isignu(nvtotd), stat=ist)
  endfunction ! create

  integer function destroytetra() result(ist)
    deallocate(&!ntt, r0, alpha, beta, gamma, &
      face, &
      rd, fa, fb, fd, isignu, stat=ist)
  endfunction ! destroy

endmodule ! PolygonFaces_mod
