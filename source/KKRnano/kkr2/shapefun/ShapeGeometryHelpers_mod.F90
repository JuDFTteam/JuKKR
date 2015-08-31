#include "DebugHelpers/test_macros.h"
!>    Auxillary module needed for shape function calculation.

module ShapeGeometryHelpers_mod
  implicit none
  private

  public :: polchk, perp, nrm2, operator(.dot.)
  
  interface operator(.dot.)
    module procedure inner_product_v3
  endinterface
  
  interface nrm2
    module procedure norm2_squared_v3
  endinterface
  

  contains
!---------------------------------------------------------------------
!>    this subroutine reads the coordinates of the vertices of each
!>    (polygon)  face of  a convex polyhedron and  checks  if these
!>    vertices arranged  consecutively define a  polygon. then  the
!>    subroutine  determines  the  vertices  and  the  edges of the
!>    polyhedron and checks if  the  number  of  vertices  plus the
!>    number of faces equals the number of edges plus 2.
!>
!>    the angular sum of the polygons is checked. it has to be (n-2)*pi
!>    n is the number of vertices.
!
!     ----------------------------------------------------------------
! ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  = 
  subroutine polchk(nface, nvertices, vert, tolvdist)
    use shape_constants_mod, only: pi
    integer, intent(in) :: nface
    integer, intent(in) :: nvertices(:)
    double precision, intent(in) :: vert(:,:,:) ! vert(3,nvertd,nfaced)
    double precision, intent(in) :: tolvdist

    integer :: ivert, ivertp,ivertm,ivrt,iedge,nvrt,nedge, iface,nvert
    double precision :: arg, down, up, fisum
    double precision :: vrt0(3), vrtp(3), vrtm(3)
    double precision, allocatable :: v1(:,:), v2(:,:), vrt(:,:)
    integer :: nfaced, nedged, nvertd, nvrtd
    logical :: new
    
    
    nfaced = size(nvertices)
    CHECKASSERT(size(vert, 3) == nfaced)
    nvertd = size(vert, 2)
    nvrtd  = nfaced*nvertd
    nedged = nvrtd+nfaced-2

    allocate(v1(3,nedged), v2(3,nedged), vrt(3,nvrtd))
    ! todo: change mode from storing copies of vertices to storing indicies ivert and iface only!  

    nvrt = 0
    nedge = 0
    do iface = 1, nface
      nvert = nvertices(iface)
      fisum = (nvert - 2)*pi
!
!------> treatment of vertices
!
      do ivert = 1, nvert
        vrt0 = vert(1:3,ivert,iface)
        new = .true.
        do ivrt = 1, nvrt
          if (nrm2(vrt0 - vrt(1:3,ivrt)) < tolvdist) new = .false. ! 0:drop this vertex
        enddo ! ivrt
        
        if (new) then
          nvrt = nvrt+1
          if (nvrt > nvrtd) stop 'increase nvrtd'
          vrt(1:3,nvrt) = vert(1:3,ivert,iface)
        endif ! new
        
        ivertp = ivert+1; if (ivert == nvert) ivertp = 1        !!! alternative: ivertp = modulo(ivert+1-1, nvert)+1
        vrtp = vert(1:3,ivertp,iface)
        ivertm = ivert-1; if (ivert == 1) ivertm = nvert        !!! alternative: ivertm = modulo(ivert-1-1, nvert)+1
        vrtm = vert(1:3,ivertm,iface) ! check if the  consecutive vertices define a polygon

        down = sqrt(nrm2(vrtp - vrt0)*nrm2(vrtm - vrt0))
        up = (vrtp - vrt0) .dot. (vrtm - vrt0)
        
        if (down < tolvdist) then
          write(6,*) 'down',down
          stop 'identical consecutive vertices'
        endif
        
        arg = up/down
        if (abs(arg) >= 1.d0) arg = sign(1.d0,arg)
        ! arg = min(max(-1.d0, arg), 1.d0) ! new formulation
        fisum = fisum - acos(arg)
!
!------> treatment of edges
!
        new = .true. ! 1:save all different edges
        do iedge = 1, nedge
          if (nrm2(vrt0 - v1(1:3,iedge)) < tolvdist) then
            if(nrm2(vrtp - v2(1:3,iedge)) < tolvdist) new = .false. ! 0:do not save
          elseif(nrm2(vrt0 - v2(1:3,iedge)) < tolvdist) then
            if(nrm2(vrtp - v1(1:3,iedge)) < tolvdist) new = .false. ! 0:do not save
          endif
        enddo ! iedge
        
        if (new) then
          nedge = nedge+1
          if (nedge > nedged) stop 'insufficient nedged'
          v1(1:3,nedge) = vert(1:3,ivert ,iface)
          v2(1:3,nedge) = vert(1:3,ivertp,iface)
        endif ! new
        
      enddo ! ivert
      
      if (fisum > 1.d-6) then
        write(6,*) 'fisum  = ',fisum
        write(6,*) 'iface  = ',iface
        write(6,*) 'nedge  = ',nedge
        write(6,*) 'nvert  = ',nvert
        write(6,*) 'nvrt   = ',nvrt
        do ivert = 1, nvrt
          write(6,'(3F16.9)') vert(1:3,ivert,iface)
        enddo ! ivert
        stop 'not consecutive vertices of a polygon'
      endif ! fisum
      
    enddo ! iface
    if ((nvrt+nface) /= (nedge+2)) then
      write(6,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
      write(6,*) '    serious warning from shape      '
      write(6,*) '   >>  stop illegal polyhedron      '
      write(6,*) 'nvrt = ',nvrt,' ; nface = ',nface,' ; nedge = ',nedge
    endif
    
  endsubroutine polchk

!-----------------------------------------------------------------------
!>    given  two  distinct  points   r1 , r2, this  routine calculates
!>    the coordinates  of the foot of  the  perpendicular from a point
!>    r0 to the line joining r1   and  r2. the logical variable inside
!>    gives the additional information whether the foot of the perpen-
!>    dicular lies within the segment or not.
!-----------------------------------------------------------------------
  subroutine perp(r0, r1, r2, tolvdist, rd, inside)
    double precision, intent(in) :: r0(3), r1(3), r2(3)
    double precision, intent(in) :: tolvdist
    double precision, intent(out) :: rd(3)
    logical,          intent(out) :: inside
!---------------------------------------------------------------------
    double precision :: dxyz(3), dabc(3), s, d
    
    dxyz = r2 - r1    
    s = r0 .dot. dxyz
    d = nrm2(dxyz)
    
    if (d < tolvdist**2) then
      write(6,fmt="(///33x,'from perp:   identical points'/33x,2('(',3e14.6,')',3x))") r1(1:3), r2(1:3)
      stop
    endif
    
    dabc(1) = s*dxyz(1) + dxyz(2)*(r1(1)*r2(2) - r1(2)*r2(1)) + dxyz(3)*(r1(1)*r2(3) - r1(3)*r2(1))
    dabc(2) = s*dxyz(2) + dxyz(3)*(r1(2)*r2(3) - r1(3)*r2(2)) + dxyz(1)*(r1(2)*r2(1) - r1(1)*r2(2))
    dabc(3) = s*dxyz(3) + dxyz(1)*(r1(3)*r2(1) - r1(1)*r2(3)) + dxyz(2)*(r1(3)*r2(2) - r1(2)*r2(3))
    rd(1:3) = dabc/d

    inside = (d - max(nrm2(rd - r1), nrm2(rd - r2)) > tolvdist)
  endsubroutine perp
  
  double precision function inner_product_v3(v, w) result(vxw)
    double precision, intent(in) :: v(3), w(3)
    vxw = v(1)*w(1) + v(2)*w(2) + v(3)*w(3)
  endfunction
  
  double precision function norm2_squared_v3(v) result(vxv)
    double precision, intent(in) :: v(3)
    vxv = v(1)*v(1) + v(2)*v(2) + v(3)*v(3)
  endfunction
  

endmodule ShapeGeometryHelpers_mod
