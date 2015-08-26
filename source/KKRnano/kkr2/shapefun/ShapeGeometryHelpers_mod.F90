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
  subroutine polchk(nface,nvertices,xvert,yvert,zvert,tolvdist)
    use shape_constants_mod, only: pi
!     integer, parameter :: nedged = nvrtd+nfaced-2
!     integer   nvertices(nfaced)
!     double precision    xvert(nvertd,nfaced),yvert(nvertd,nfaced),zvert(nvertd,nfaced)

    integer, intent(in) :: nface
    integer, intent(in) :: nvertices(:)
    double precision, intent(in) :: xvert(:,:),yvert(:,:),zvert(:,:) ! todo: unite into an array of shape vert(3,nvertd,nfaced)
    double precision, intent(in) :: tolvdist

    integer :: ivert,inew,ivertp,ivertm,ivrt,iedge,nvrt,nedge, iface,nvert
    double precision :: arg,a1,a2,down,up,fisum,t
    double precision :: vrt0(3), vrtp(3), vrtm(3)
    double precision, allocatable :: v1(:,:),v2(:,:),v(:,:),vrt(:,:) ! v1(3,nedged),v2(3,nedged),v(3,nvertd),vrt(3,nvrtd)
    integer :: nfaced, nedged, nvertd, nvrtd

    nfaced = size(nvertices)
    nvertd = size(xvert,1)
    nvrtd  = nfaced*nvertd
    nedged = nvrtd+nfaced-2

    allocate(v1(3,nedged), v2(3,nedged), v(3,nvertd), vrt(3,nvrtd))

    nvrt = 0
    nedge = 0
    do iface = 1, nface
      nvert = nvertices(iface)
      fisum = (nvert-2)*pi
      do ivert = 1, nvert
        v(1:3,ivert) = [xvert(ivert,iface), yvert(ivert,iface), zvert(ivert,iface)]
      enddo ! ivert
!
!------> t r e a t m e n t   o f   v e r t i c e s
!
      do ivert = 1, nvert
        vrt0 = v(1:3,ivert)
        inew = 1 ! 1:save all different vertices
        do ivrt = 1, nvrt
          if (nrm2(vrt0 - vrt(1:3,ivrt)) < tolvdist) inew = 0 ! 0:drop this vertex
        enddo ! ivrt
        
        if (inew == 1) then
          nvrt = nvrt+1
          if(nvrt > nvrtd) stop 'increase nvrtd'
          vrt(1:3,nvrt) = v(1:3,ivert)
        endif ! inew
        
        ivertp = ivert+1; if (ivert == nvert) ivertp = 1        !!! alternative: ivertp = modulo(ivert+1-1, nvert)+1
        vrtp = v(1:3,ivertp)
        ivertm = ivert-1; if (ivert == 1) ivertm = nvert        !!! alternative: ivertm = modulo(ivert-1-1, nvert)+1
        vrtm = v(1:3,ivertm) ! check if the  consecutive vertices define a polygon

        a1 = nrm2(vrtp - vrt0)
        a2 = nrm2(vrtm - vrt0)
        up = (vrtp - vrt0) .dot. (vrtm - vrt0)
        down = sqrt(a1*a2)
        
        if(down >= tolvdist)   then
          arg = up/down
          if (abs(arg) >= 1.d0) arg = sign(1.d0,arg)
          fisum = fisum - acos(arg)
        else
          write(6,*) 'down',down
          stop 'identical consecutive vertices'
        endif
!
!------> t r e a t m e n t   o f   e d g e s
!
        inew = 1 ! 1:save all different edges
        do iedge = 1, nedge
          if (nrm2(vrt0 - v1(1:3,iedge)) < tolvdist) then
            if (nrm2(vrtp - v2(1:3,iedge)) < tolvdist) inew = 0 ! 0:do not save
          else
            if (nrm2(vrt0 - v2(1:3,iedge)) < tolvdist) then
              if (nrm2(vrtp - v1(1:3,iedge)) < tolvdist) inew = 0 ! 0:do not save
            endif
          endif
        enddo ! iedge
        
        if (inew == 1) then
          nedge = nedge+1
          if (nedge > nedged) stop 'insufficient nedged'
          v1(1:3,nedge) = v(1:3,ivert)
          v2(1:3,nedge) = v(1:3,ivertp)
        endif ! inew
        
      enddo ! ivert
      if (fisum > 1d-6) then
        write(6,*) 'fisum  = ',fisum
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
  subroutine perp(r0,r1,r2,rd,tolvdist,inside)
    logical, intent(out) :: inside ! result
    double precision, intent(out) :: rd(3)
    double precision, intent(in) :: r0(3), r1(3), r2(3)
    double precision, intent(in) :: tolvdist
!---------------------------------------------------------------------
    double precision :: dxyz(3), dabc(3), s, d
    
    dxyz = r2 - r1    
    s = r0 .dot. dxyz
    d = nrm2(dxyz)
    
    if (d < tolvdist*tolvdist) then
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
