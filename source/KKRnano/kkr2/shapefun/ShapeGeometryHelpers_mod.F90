!>    Auxillary module needed for shape function calculation.

module ShapeGeometryHelpers_mod
  implicit none
  private

  public :: polchk, perp

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
!     integer, parameter :: nedged = nvrtd+nfaced-2
!     integer   nvertices(nfaced)
!     double precision    xvert(nvertd,nfaced),yvert(nvertd,nfaced),zvert(nvertd,nfaced)

    integer, intent(in) :: nface
    integer, intent(in) :: nvertices(:)
    double precision, intent(in) :: xvert(:,:),yvert(:,:),zvert(:,:)
    double precision, intent(in) :: tolvdist

    integer :: ivert,inew,ivertp,ivertm,ivrt,iedge,nvrt,nedge
    integer :: iface,nvert
    double precision :: arg,a1,a2,down,up,fisum,t
    double precision :: vrtx,vrty,vrtz,vrtpx,vrtpy,vrtpz,vrtmx,vrtmy,vrtmz
    double precision :: pi314
    double precision, allocatable :: v1(:,:),v2(:,:),v(:,:),vrt(:,:) ! v1(3,nedged),v2(3,nedged),v(3,nvertd),vrt(3,nvrtd)
    integer :: nfaced, nedged, nvertd, nvrtd

    nfaced  =  size(nvertices)
    nvertd  =  size(xvert,1)
    nvrtd   =  nfaced*nvertd
    nedged  =  nvrtd+nfaced-2

    allocate(v1(3,nedged),v2(3,nedged),v(3,nvertd),vrt(3,nvrtd))

    pi314  =  4.d0*atan(1.d0)

    nvrt = 0
    nedge = 0
    do iface = 1,nface
      nvert = nvertices(iface)
      fisum = (nvert-2)*pi314
      do ivert = 1, nvert
        v(1,ivert) = xvert(ivert,iface)
        v(2,ivert) = yvert(ivert,iface)
        v(3,ivert) = zvert(ivert,iface)
      enddo ! ivert
!
!------> t r e a t m e n t   o f   v e r t i c e s
!
      do ivert = 1, nvert
        vrtx = v(1,ivert)
        vrty = v(2,ivert)
        vrtz = v(3,ivert)
        inew = 1 ! 1:save all different vertices
        do ivrt = 1,nvrt
          t = (vrtx-vrt(1,ivrt))**2 + (vrty-vrt(2,ivrt))**2 + (vrtz-vrt(3,ivrt))**2
          if(t < tolvdist) inew = 0
        enddo ! ivrt
        
        if (inew  ==  1) then
          nvrt = nvrt+1
          if(nvrt > nvrtd) stop 'increase nvrtd'
          vrt(1,nvrt) = v(1,ivert)
          vrt(2,nvrt) = v(2,ivert)
          vrt(3,nvrt) = v(3,ivert)
        endif
        ivertp = ivert+1
        if (ivert == nvert) ivertp = 1
        vrtpx = v(1,ivertp)
        vrtpy = v(2,ivertp)
        vrtpz = v(3,ivertp)
        ivertm = ivert-1
        if (ivert == 1) ivertm = nvert
        vrtmx = v(1,ivertm)
        vrtmy = v(2,ivertm)               ! check if the  consecutive
        vrtmz = v(3,ivertm)               ! vertices define a polygon
        a1 = sqrt((vrtpx-vrtx)**2 + (vrtpy-vrty)**2 + (vrtpz-vrtz)**2)
        a2 = sqrt((vrtmx-vrtx)**2 + (vrtmy-vrty)**2 + (vrtmz-vrtz)**2)
        down = a1*a2
        up = (vrtpx-vrtx)*(vrtmx-vrtx) + (vrtpy-vrty)*(vrtmy-vrty) + (vrtpz-vrtz)*(vrtmz-vrtz)
        ! write(6,*)
        ! write(6,*) vrtx,vrty,vrtz
        ! write(6,*) vrtmx,vrtmy,vrtmz
        ! write(6,*) vrtpx,vrtpy,vrtpz
        ! write(6,*) 'fisum ',ivert,a1,a2,up,arg,acos(arg)
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
          t = (vrtx-v1(1,iedge))**2 + (vrty-v1(2,iedge))**2 + (vrtz-v1(3,iedge))**2
          if (t < tolvdist) then
            t = (vrtpx-v2(1,iedge))**2 + (vrtpy-v2(2,iedge))**2 + (vrtpz-v2(3,iedge))**2
            if(t < tolvdist) inew = 0
          else
            t = (vrtx-v2(1,iedge))**2 + (vrty-v2(2,iedge))**2 + (vrtz-v2(3,iedge))**2
            if(t < tolvdist) then
              t = (vrtpx-v1(1,iedge))**2 + (vrtpy-v1(2,iedge))**2 + (vrtpz-v1(3,iedge))**2
              if(t < tolvdist) inew = 0
            endif
          endif
        enddo ! iedge
        if (inew == 1) then
          nedge = nedge+1
          if(nedge > nedged) stop 'insufficient nedged'
          v1(1,nedge) = v(1,ivert )
          v1(2,nedge) = v(2,ivert )
          v1(3,nedge) = v(3,ivert )
          v2(1,nedge) = v(1,ivertp)
          v2(2,nedge) = v(2,ivertp)
          v2(3,nedge) = v(3,ivertp)
        endif ! inew
      enddo ! ivert
      if(fisum > 1d-6) then
        write(6,*) 'fisum  = ',fisum
        stop 'not consecutive vertices of a polygon'
      end if
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
    double precision :: dx,dy,dz,s,d,da,db,dc,d1,d2,co
    
    dx = r2(1) - r1(1)
    dy = r2(2) - r1(2)
    dz = r2(3) - r1(3)
    s = r0(1)*dx + r0(2)*dy + r0(3)*dz
    d = dx*dx + dy*dy + dz*dz
    if(d < tolvdist*tolvdist) then
      write(6,fmt="(///33x,'from perp:   identical points'/33x,2('(',3e14.6,')',3x))") r1(1:3), r2(1:3)
      stop
    endif      
    da = s*dx + dy*(r1(1)*r2(2) - r1(2)*r2(1)) + dz*(r1(1)*r2(3) - r1(3)*r2(1))
    db = s*dy + dz*(r1(2)*r2(3) - r1(3)*r2(2)) + dx*(r1(2)*r2(1) - r1(1)*r2(2))
    dc = s*dz + dx*(r1(3)*r2(1) - r1(1)*r2(3)) + dy*(r1(3)*r2(2) - r1(2)*r2(3))
    rd(1) = da/d
    rd(2) = db/d
    rd(3) = dc/d
    d1 = (rd(1) - r1(1))**2 + (rd(2) - r1(2))**2 + (rd(3) - r1(3))**2
    d2 = (rd(1) - r2(1))**2 + (rd(2) - r2(2))**2 + (rd(3) - r2(3))**2

    co = d - max(d1, d2)
    inside = (co > tolvdist)
  endsubroutine perp

endmodule ShapeGeometryHelpers_mod
