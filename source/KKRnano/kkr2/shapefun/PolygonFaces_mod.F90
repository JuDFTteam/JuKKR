!------------------------------------------------------------------------------
!> Storage for polygon face data and additional information on vertices
module PolygonFaces_mod
  implicit none
  private
  public :: TetrahedronAngles, PolygonFace, destroy

  type TetrahedronAngles
    double precision :: rd !< distances pyramid footpoint to edge
    double precision :: fa !< tetrahedron angle, phi-angle corresponding to 1st vertex  fa < fd < fb
    double precision :: fb !< tetrahedron angle, phi-angle corresponding to 2nd vertex
    double precision :: fd !< tetrahedron angle, phi-angle corresponding to foot point between 1st and 2nd vertex
    integer(kind=1) :: isignu  !< sign for rotation sense
  endtype

  type PolygonFace
    integer :: ntt              !< number of tetrahedra for polygon
    double precision :: r0      !< foot points of perpendicular of polygons
    double precision :: Euler(1:3) !< Euler angles alpha beta gamma to rotate faces perpendicular to z-axis
    type(TetrahedronAngles), allocatable :: ta(:)
  endtype

  interface destroy
    module procedure destroyPolygonFace, destroyTetrahedronAngles
  endinterface  
  
  contains
  
  elemental subroutine destroyPolygonFace(self)
    type(PolygonFace), intent(inout) :: self
    integer :: ist
    call destroy(self%ta) 
    deallocate(self%ta, stat=ist)
  endsubroutine ! destroy

  elemental subroutine destroyTetrahedronAngles(self)
    type(TetrahedronAngles), intent(in) :: self
    ! nothing to be deallocated   
  endsubroutine ! destroy

endmodule ! PolygonFaces_mod
