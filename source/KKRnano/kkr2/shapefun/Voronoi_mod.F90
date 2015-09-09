
module Voronoi_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  public :: Voronoi_construction

  interface normal_plane
    module procedure normal_plane0f, normal_planef
  endinterface
  
  contains
  
  
  !***********************************************************************
  ! Given a cluster of atomic positions at RVEC(3,NVEC), this subroutine
  ! returns information about the Voronoi cell around the origin. It is
  ! supposed, of course, that the origin corresponds to an atomic position
  ! which is not included in the RVEC(3,N). The information returned is:
  !
  ! RMT: Muffin-tin radius (radius of inscribed sphere centered at the 
  !      origin).
  !
  ! ROUT: Radius of circumscribed sphere centered at the origin.
  !
  ! VOLUME: Volume of the Voronoi cell. 
  !
  ! NFACE: Number of faces of the cell.
  !
  ! A3(NFACE),B3(NFACE),C3(NFACE),D3(NFACE): Coefficients defining the
  ! faces via A3*x+B3*y+C3*z=D3. The arrays are filled in the first
  ! NFACE positions with the information, the rest can be garbage.
  !
  ! NVERT(NFACE): Number of vertices corresponting to each face.
  !
  ! vert(1,NVERTMAX,NFACE),vert(2,NVERTMAX,NFACE),vert(3,NVERTMAX,NFACE): 
  ! Coordinates of theese vertices.
  !
  ! The Voronoi construction performed here allows for different than
  ! 50%/50% bisections. For this, each atomic site, positioned say at
  ! vector r(i), is assigned a weight w(i). Then a point r in space
  ! belongs to the cell i if, for all j, 
  ! |r-r(i)|**2 - w(i) < |r-r(j)|**2 - w(j).    (1)
  ! The points for which the inequality becomes an equality is the 
  ! bisector, and it can be easily shown that it is a plane perpendicular
  ! to the vector r(j)-r(i).
  ! One can parametrize the segment connecting r(i) and r(j) using a
  ! parameter t, 0 < t < 1. For t=0 we are at r(i), for t=1 at r(j), and
  ! in-between we have a linear dependence of the distance from r(i) with
  ! t. I.e., the position vector is r = r(i)*(1-t) + r(j)*t.
  ! The point where the bisector cuts this segment will then be at
  ! t=(1/2)*( (w(i)-w(j))/dist**2 + 1 )         (2)
  ! where dist is the distance of the two points. As a special case we see
  ! that, for w(i)=w(j), the bisector will be in the middle, else the
  ! atomic site with the bigger weight gains more space.
  !
  ! The above procedure is guaranteed to divide space into tesselating 
  ! (:=space-filling) convex polyhedra (one of their sides could be at 
  ! infinity in special cases). The convexity of the polyhedra is 
  ! guaranteed by the constuction, i.e. as mutual cut of half-spaces, and
  ! the property of tesselation is guaranteed by the fact that every point
  ! in space is assigned to some atom.
  !
  ! However, this procedure might place the polyhedron in such a way that
  ! its own atomic site is not contained in it. This can happen if there
  ! is a big enough difference in the weights, whence, from eq.(2) above,
  ! t can become <0 or >1 (then the bisector has passed one of the two
  ! points). This cannot be allowed in a KKR calculation, and therefore
  ! each time it is checked that this is not the case. If such a case
  ! occurs, a different set of weights must be chosen.
  !
  !
  ! Uses subroutines NORMALPLANE, POLYHEDRON, and function DISTPLANE.
  subroutine Voronoi_construction(nvec,rvec,nvertmax,nfaced,weight0,weight,tolvdist,tolarea, &
     rmt,rout,volume, nface, planes, nvert, vert, output)
    use Sorting_mod, only: dsort

    integer, intent(in) :: nvec, nvertmax, nfaced
    double precision, intent(in) :: rvec(:,:) ! cluster positions without the origin among them
    double precision, intent(in) :: weight0, weight(:)  ! Weight of central atom, and of all others (dimensioned as RVEC). 
    double precision, intent(in) :: tolvdist ! Max. tolerance for distance of two vertices
    double precision, intent(in) :: tolarea  ! Max. tolerance for area of polygon face
    logical, intent(in) :: output ! test output

    integer, intent(out) :: nface, nvert(:)
    double precision, intent(out) :: volume, rmt, rout
    double precision, intent(out) :: planes(0:,:)
    double precision, intent(out) :: vert(:,:,:) ! (3,nvertmax,nfaced) ! vertices
    
    integer :: ivec, iface, ivert, i, nplane
    double precision :: tetrvol,rsq,tau, v1(3), v2(3), v3(3)
    double precision :: facearea(nfaced),trianglearea   
    double precision :: temp
    integer :: isort(nfaced) ! Index for sorting
    double precision :: rsort(nfaced)              ! Aux. function for sorting
    double precision, allocatable :: ptmp(:,:), vtmp(:,:,:) ! For sorting
    integer, allocatable :: ntmp(:)                        ! For sorting

    !---------------------------------------------------------------
    ! Check that the origin is not included in RVEC.
    do ivec = 1, nvec
      if (sum(abs(rvec(1:3,ivec))) < 1.d-12) die_here("vector"+ivec+"is zero!")
    enddo ! ivec
    !---------------------------------------------------------------
    ! Define the planes as normal to the vectors RVEC, passing from t as in eq. (2) above:
    do ivec = 1, nvec
      rsq = rvec(1,ivec)**2 + rvec(2,ivec)**2 + rvec(3,ivec)**2
      tau = 0.5d0*((weight0 - weight(ivec))/rsq + 1.d0)
      planes(0:3,ivec) = normal_plane(rvec(1:3,ivec), tau)
    enddo ! ivec
    !---------------------------------------------------------------
    ! Find the Voronoi polyhedron.
    nplane = nvec
    if (output) write(*,*) 'Entering POLYHEDRON08'
    call polyhedron08(nplane, nvertmax, nfaced, tolvdist, tolarea, planes, nface, nvert, vert, output)
    if (output) write(*,*) 'Exited POLYHEDRON08'

    !---------------------------------------------------------------
    ! Calculate the volume as sum of the volumes of all tetrahedra 
    ! connecting the origin to the faces. Use for each tetrahedron
    ! volume = det((r0-r1),(r0-r2),(r0-r3))/6, where r0 is here the
    ! origin and r1,r2,r3 the vectors of the 3 other vertices.
    ! Algorithm requires that the face is a convex polygon with the 
    ! vertices ordered (doesn't matter if they are clock- or anticlockwise).
    volume = 0.d0
    do iface = 1, nface
      v1(1:3) = vert(1:3,1,iface)
      facearea(iface) = 0.d0
      do ivert = 2, nvert(iface)-1
          v2(1:3) = vert(1:3,ivert,iface)
          v3(1:3) = vert(1:3,ivert+1,iface)
          tetrvol = v1(1)*(v2(2)*v3(3)-v3(2)*v2(3))+v2(1)*(v3(2)*v1(3)-v1(2)*v3(3))+v3(1)*(v1(2)*v2(3)-v2(2)*v1(3))
          volume = volume + abs(tetrvol)
          trianglearea = 0.5d0 * sqrt( &
                ( v1(1)*v2(2) + v2(1)*v3(2) + v3(1)*v1(2) - v2(2)*v3(1) - v3(2)*v1(1) - v1(2)*v2(1))**2 &
              + ( v1(2)*v2(3) + v2(2)*v3(3) + v3(2)*v1(3) - v2(3)*v3(2) - v3(3)*v1(2) - v1(3)*v2(2))**2 &
              + ( v1(3)*v2(1) + v2(3)*v3(1) + v3(3)*v1(1) - v2(1)*v3(3) - v3(1)*v1(3) - v1(1)*v2(3))**2 )
          facearea(iface) = facearea(iface)+ trianglearea
      enddo ! ivert
    enddo ! iface
    volume = volume/6.d0

    if (output) then
      write(*,*) ' Polyhedron properties '
      write(*,*) ' Number of faces : ',nface
    endif ! output

    if (output) then
      do iface=1,nface
          write(*,fmt="(' Face ',i4,'   has ',i4,' vertices ','; Area= ',e12.4)") iface, nvert(iface), facearea(iface)
          do ivert = 1, nvert(iface)
            write(*,fmt="(i5,4e16.8)") ivert, vert(1:3,ivert,iface)
          enddo ! ivert
          write(*,fmt="(' Face coefficients:',4e16.8)") planes(1:3,iface), planes(0,iface)
      enddo ! iface
      write(*,*) 'The Volume is : ',volume
    endif ! output

    !---------------------------------------------------------------
    ! Find RMT:
    iface = 1
    rmt = dist_plane(planes(1:3,iface), planes(0,iface))
    do iface = 2, nface
      temp = dist_plane(planes(1:3,iface), planes(0,iface))
      if (temp < rmt) rmt = temp
    enddo ! iface
    
    ! Fint ROUT, the largest radius of all vertices
    rout = 0.d0
    do iface = 1,nface
      do ivert = 1, nvert(iface)
          temp = sum(vert(1:3,ivert,iface)**2)
          if (temp > rout) rout = temp
      enddo ! ivert
    enddo ! iface
    rout = sqrt(rout)

    if (output) write(*,fmt="('Voronoi subroutine: RMT=',e16.8,'; ROUT=',e16.8,'; RATIO=',f12.2,' %')") rmt,rout,rmt*100/rout

    !---------------------------------------------------------------
    ! Sort faces:
    do iface = 1, nface
      rsort(iface) = 1.d8*dist_plane(planes(1:3,iface), planes(0,iface)) + dot_product([0.1d0, 1.d2, 1.d5], planes(1:3,iface))
    enddo ! iface
    
    call dsort(rsort, isort, nface)
    ! Rearrange using a temporary arrays ptmp, vtmp, ntmp
    allocate(ptmp(0:3,nface), vtmp(3,nvertmax,nface), ntmp(nface))
    ptmp(0:3,:) = planes(0:3,1:nface)
    vtmp(:,:,:) = vert(:,:,1:nface)
    ntmp(:)     = nvert(1:nface)
    do iface = 1, nface
      i = isort(iface)
      planes(:,i) = ptmp(:,  iface) ! (0:3, ...)
      vert(:,:,i) = vtmp(:,:,iface) ! (1:3,:, ...)
      nvert   (i) = ntmp    (iface)
    enddo ! iface
    
    deallocate(ptmp, vtmp, ntmp)
    !---------------------------------------------------------------

  endsubroutine Voronoi_construction
  

  subroutine polyhedron08(nplane, nvertmax, nfaced, tolvdist,tolarea, planes, nface, nvert, vert, output)
    ! given a set of planes, defined by a3*x+b3*y+c3*z=d3 and defining
    ! a convex part of space (the minimal one containing the origin, 
    ! usually a ws-polyhedron), this subroutine returns the actual faces of
    ! the polyhedron, discarding the planes that do not contain faces. also,
    ! the coordinates of the verticess of the faces xvert,yvert,zvert and their 
    ! number nvert per face are returned. the coefficients of the actual
    ! faces are returned in the same arrays a3,b3,c3, and d3.

    logical, intent(in) :: output
    integer, intent(in) :: nplane                ! initial number of planes.
    integer, intent(in) :: nvertmax              ! max. number of vertices per plane.
    integer, intent(in) :: nfaced                ! max. number of faces.
    double precision, intent(in) :: tolvdist              ! max. tolerance for distance of two vertices
    double precision, intent(in) :: tolarea               ! max. tolerance for area of polygon face
    double precision, intent(inout) :: planes(0:,:)  ! coefs. defining the planes, dimensioned >= nplane.
    double precision, intent(inout) :: vert(:,:,:) ! (3,nvertmax,nfaced) ! cartesian coords. of vertices for each plane (2nd index is for planes).
    integer, intent(out) :: nvert(:)  ! number of vertices found for each face
    integer, intent(out) :: nface     ! number of faces found (with nvert>0).

    !---------------------------------------------------------------
    ! find all faces and vertices of the polyhedron. on output, the vertices of each face are sorted (clockwise or anticlockwise).
    if (output) write(*,*) 'entering vertex3d'
    call vertex3d(nplane, planes, nvertmax, tolvdist, nface, nvert, vert, output)
    if (output) write(*,*) 'vertex3d found',nface,' faces with >3 vertices.'
    !---------------------------------------------------------------
    ! analyze the faces and vertices of the polyhedron.
    ! use criteria for rejecting faces that are too small or vertices that are too close to each other.
    ! on output, number of faces and vertices may be reduced after some rejections have taken place.
    if (output) write(*,*) 'entering analyzevert3d'
    call analyzevert3d(nvertmax, nfaced, tolvdist, tolarea, nplane, nface, nvert, vert, planes, output)
    if (output) write(*,*) 'analyzevert3d accepted',nface,' faces.'

  endsubroutine polyhedron08


    !***********************************************************************
  subroutine vertex3d(nplane, planes, nvertmax, tolvdist, nface, nvert, vert, output)
    ! given a set of planes, defined by a3*x+b3*y+c3*z=d3 and defining
    ! a convex part of space (the minimal one containing the origin, 
    ! usually a ws-polyhedron), this subroutine returns the vertices
    ! of this polyhedron in cartesian coordinates. for the planes that
    ! are not faces of the polyhedron, a value nvert(iplane)=0 is returned.
    ! the total no. of faces found is returned as nface.

    integer, intent(in) :: nplane ! number of planes.
    integer, intent(in) :: nvertmax              ! max. number of vertices per plane.
    double precision, intent(in) :: planes(0:,:)  ! coefs. defining the planes, dimensioned >= nplane.
    double precision, intent(in) :: tolvdist               ! min. distance between vertices
    integer, intent(out) :: nvert(:)  ! number of vertices found for each face
    integer, intent(out) :: nface     ! number of faces found (with nvert>0).
    double precision, intent(inout) :: vert(:,:,: ) ! (3,nvertmax,*) ! cartesian coords. of vertices for each plane ! (2nd index is for planes).
    logical, intent(in) :: output

    integer ip1,ip2,ip3,ipl,kpl ! plane indices
    integer ivert                    ! vertex index
    double precision cut(1:3) ! cut point of three planes.
    double precision det(0:3) ! determinants of 3x3 system for xcut,ycut,zcut.
    double precision tolvdist2, distance2   ! a distance criterium of two points in space
    ! the following are for sorting the vertices of each face:
    double precision v1(3),v2(3),v3(3)   ! auxiliary vectors...
    double precision sinfiv1v2,cosfiv1v2 ! ...and their inner and outer products
    double precision fi(nvertmax)        ! ...and also their relative angles.
    double precision uv(3),vl,sn         ! unit vector, length, sign of sin(fi)
    logical laccept ! determining whether a cut point is inside the polyhedron.
    !---------------------------------------------------------------
    ! check & initialize
    if (nplane < 4) die_here("nplane was only"+nplane)

    tolvdist2 = tolvdist**2
    
    nvert(1:nplane) = 0
    !===============================================================
    ! start loop over all planes that can be cut:
    plane1: do ip1 = 1, nplane
      ! start loop over all other planes:
      plane2: do ip2 = 1,nplane
      if (ip2 == ip1) cycle plane2
        ! start loop over all other-other (!) planes. do from ip2+1 to 
        ! nplane so that no pair is considered twice.
        plane3: do ip3 = ip2+1, nplane        ! nikos  ip2+1,nplane
          if (ip3 == ip1) cycle plane3
            !     if (ip3 == ip2) goto 100 ! added by nikos
            ! solve the 3x3 system to find the cut point.
            !-----------------------------------
    !-----------------------------------
#define a3(I) planes(1,I)
#define b3(I) planes(2,I)
#define c3(I) planes(3,I)
#define d3(I) planes(0,I)

    det(0) = a3(ip1)*(b3(ip2)*c3(ip3) - b3(ip3)*c3(ip2)) + a3(ip2)*(b3(ip3)*c3(ip1) - b3(ip1)*c3(ip3)) + a3(ip3)*(b3(ip1)*c3(ip2) - b3(ip2)*c3(ip1))

    !---------------------------------------------------------------
    if (dabs(det(0)) > 1.d-12) then ! there is a cut point 

    det(1) = d3(ip1)*(b3(ip2)*c3(ip3) - b3(ip3)*c3(ip2)) + d3(ip2)*(b3(ip3)*c3(ip1) - b3(ip1)*c3(ip3)) + d3(ip3)*(b3(ip1)*c3(ip2) - b3(ip2)*c3(ip1))
    det(2) = a3(ip1)*(d3(ip2)*c3(ip3) - d3(ip3)*c3(ip2)) + a3(ip2)*(d3(ip3)*c3(ip1) - d3(ip1)*c3(ip3)) + a3(ip3)*(d3(ip1)*c3(ip2) - d3(ip2)*c3(ip1))
    det(3) = a3(ip1)*(b3(ip2)*d3(ip3) - b3(ip3)*d3(ip2)) + a3(ip2)*(b3(ip3)*d3(ip1) - b3(ip1)*d3(ip3)) + a3(ip3)*(b3(ip1)*d3(ip2) - b3(ip2)*d3(ip1))

    cut(1:3) = det(1:3)/det(0)

    !     write(*,333) ip1,ip2,ip3,cut
    !333  format('cutting point of planes ',3i5,':',3d15.7)
    !-----------------------------------
    ! accept this cut point as a vertex, if it belongs to the polyhedron. so,
    ! make a loop over all other (than ip1,2,3) planes:
    laccept = .true.
    do kpl = 1, nplane
      if (kpl == ip1 .or. kpl == ip2 .or. kpl == ip3) cycle
      laccept = laccept .and. halfspace(a3(kpl),b3(kpl),c3(kpl),d3(kpl),cut(1),cut(2),cut(3))
    enddo ! kpl

#undef d3
#undef c3
#undef b3
#undef a3
    !-----------------------------------
            !-----------------------------------
            if (laccept) then
            ! if the cut point found belongs to the cell, we accept it unless it has
            ! occured before for this face (ip1). such a situation is possible
            ! when 4 or more planes pass through the same point (e.g. for the vertices
            ! of the fcc ws-cell). so...
              do ivert = 1,nvert(ip1)
                distance2 = (vert(1,ivert,ip1) - cut(1))**2 + (vert(2,ivert,ip1) - cut(2))**2 + (vert(3,ivert,ip1) - cut(3))**2
                if (distance2 < tolvdist2) then
                  laccept = .false. ! vertex is too close to a previous one.
                  exit              ! jump loop, no need to continue.
                endif
              enddo ! ivert
              
            endif ! laccept

            ! now we are ready to add the point to the vertex list.
            if (laccept) then
              nvert(ip1) = nvert(ip1) + 1
              vert(1:3,nvert(ip1),ip1) = cut(1:3)
            endif

          endif ! (dabs(det) > 1.d-12)
          !---------------------------------------------------------------

        enddo plane3 ! ip3
      enddo plane2 ! ip2
    !     write(*,*) 'number of vertices for plane ',ip1,'  :',nvert(ip1)
    enddo plane1 ! ip1

    !===============================================================
    ! each plane should finally have either at least 3 vertices, if it is a
    ! face of the polyhedron, or none at all. check this:
    if (output) then
      do ipl = 1, nplane
        if (nvert(ipl) == 1 .or. nvert(ipl) == 2) warn(6, "there is a problem with the vertices for plane #"-ipl-', only'+nvert(ipl)+'vertices were found!')
      enddo ! ipl
    endif ! output

    !===============================================================
    ! for each face of the polyhedron, sort the vertices in a consecutive order
    ! as vertices of the polygon. the order is not necessarily mathematically
    ! positive.
    nface = 0
    do ipl = 1, nplane
      if (nvert(ipl) >= 3) then
        nface = nface + 1      ! count the faces
        fi(1) = -4.d0          ! just a number smaller than -pi.
      ! unit vector in the direction of first vertex:
        vl = sqrt(vert(1,1,ipl)**2 + vert(2,1,ipl)**2 + vert(3,1,ipl)**2)
        uv(1:3) = vert(1:3,1,ipl)/vl

      ! define the vector connecting the first vertex to the (now-) second:
        v1(1:3) = vert(1:3,2,ipl) - vert(1:3,1,ipl)

        do ivert = 2,nvert(ipl)
    ! define the vector connecting the first vertex to the current one:
          v2(1:3) = vert(1:3,ivert,ipl) - vert(1:3,1,ipl)
    ! find the angle fi between v1 and v2 
    ! ( always, -pi < fi < pi from the definition of datan2 )
          cosfiv1v2 = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
          v3 = crospr(v1, v2) ! cross product = |v1|*|v2|*sinfi
          sinfiv1v2 = sqrt(v3(1)*v3(1) + v3(2)*v3(2) + v3(3)*v3(3))
    ! sign of sinfi is defined with respect to unit vector uv (see above)
          sn = uv(1)*v3(1) + uv(2)*v3(2) + uv(3)*v3(3)
          if (sn < 0) sinfiv1v2 = -sinfiv1v2

          if (sinfiv1v2 == 0.d0 .and. cosfiv1v2 == 0.d0) &
            die_here("found two identical vertex points")
!             then
!       ! point falls exactly on 1st vertex...
!               fi(ivert) = -4.d0
!               die_here("found two identical vertex points")
!       ! ...while it shouldn't ! (this was checked earlier)
!             else
!               fi(ivert) = datan2(sinfiv1v2, cosfiv1v2)
!             endif
          fi(ivert) = datan2(sinfiv1v2, cosfiv1v2)
           
        enddo ! ivert = 3,nvert(ipl)

        ! store with respect to the angle found:
        call sortvertices(nvert(ipl), fi, vert(:,:,ipl)) ! sort the vertices in-place

      endif ! (nvert(ipl) >= 3)
    enddo ! ipl = 1, nplane

    endsubroutine vertex3d

    !------------------------------------------------------------------------------
    subroutine analyzevert3d(nvertmax, nfaced, tolvdist, tolarea, nplane, nface, nvert, vert, planes, output)
    ! analyze the faces and vertices of the polyhedron.
    ! use criteria for rejecting faces that are too small or vertices that are too close to each other.
    ! on output, number of faces and vertices may be reduced after some rejections have taken place.
    logical, intent(in) :: output
    integer, intent(in) :: nvertmax, nfaced
    integer, intent(in) :: nplane ! number of planes
    integer, intent(inout) :: nface ! number of faces (usually much less than nplane)
    integer, intent(inout) :: nvert(:) ! (nfaced) ! number of vertices for each face
    double precision, intent(inout) :: vert(:,:,:) ! (3,nvertmax,nfaced)
    double precision, intent(inout) :: planes(0:,:) ! coefs. defining the planes, to be reordered at end
    double precision, intent(in) :: tolvdist, tolarea  ! max. tolerance for distance of two vertices and area of face.

    double precision vdist2, tolvdist2 !, vdist ! distance between consecutive vertices
    logical lacceptvert(nvertmax),lacctot,lfoundnext,lthisisthelast
    logical lacceptface(nfaced)
!   integer newindexface(nfaced)
    integer iface,ivert,ivert2,inext,iplane
    integer ifacenewcount,ivertnewcount
    double precision dv(3),x1,x2,x3,y1,y2,y3,z1,z2,z3,x4,y4,z4,trianglearea
    double precision facearea(nfaced)
    double precision det,detsum

    tolvdist2 = tolvdist**2
    
    ! first analyze vertices.
    ! reject doubles (also vertices which fall almost on the previous vertex).
    do iplane = 1, nplane
      if (nvert(iplane) == 0) cycle ! no need checking this one, proceed to next!

      lacctot = .true.
    ! first vertices 1st with 2nd, 2nd with 3rd, etc...
      lacceptvert(1) = .true.                 ! first vertex is always accepted.
      ivert = 1
      lthisisthelast = .false.   ! this flag will go up at the last acceptable vertex.

    ! double loop: first loop is over all vertices; 
    ! but if during the loop vertices are found that have to be rejected, they are jumped over.
      do while (ivert < nvert(iplane) .and. .not. lthisisthelast)
        lfoundnext = .false.    ! this flag will become true when the next acceptable vertex is found.
        ivert2 = ivert + 1
    ! second loop is over subsequent vertices (i.e., vertices ivert2 > ivert).
    ! look for the first acceptable-as-next vertex, but do not go beyond last vertex.
    ! stop loop as soon as the acceptable next vertex is found.
        do while (.not. lfoundnext .and. ivert2 <= nvert(iplane))
          dv(1:3) = vert(1:3,ivert2,iplane) - vert(1:3,ivert,iplane)
          vdist2 = dv(1)**2 + dv(2)**2 + dv(3)**2
!             vdist = sqrt(vdist2)

          if (vdist2 >= tolvdist2) then
            lacceptvert(ivert2) = .true.   ! vertex is to be accepted
            inext = ivert2                 ! set this as the next vertex
            lfoundnext = .true.            ! and we have a winner, exit loop.
          else
            lacceptvert(ivert2) = .false.  ! remember that vertex is to be rejected later
            lacctot = .false.              ! remember that at least one vertex has to be rejected.
            ivert2 = ivert2 + 1            ! now compare to the next vertex
          endif  
        enddo ! while

        if (.not. lfoundnext) lthisisthelast = .true.  ! if there is no next acceptable vertex,
                                                       ! then this was the last one. jump out.
        ivert = inext
      enddo ! while

    ! ...and now 1st with last to close the cycle:
      ivert = 1
      ivert2 = nvert(iplane)
      dv(1:3) = vert(1:3,ivert2,iplane) - vert(1:3,ivert,iplane)
      vdist2 = dv(1)**2 + dv(2)**2 + dv(3)**2
!     vdist = sqrt(vdist2)
      if (vdist2 >= tolvdist2) then
        lacceptvert(ivert2) = .true.        ! vertex is to be accepted
      else
        lacceptvert(ivert2) = .false.       ! remember that vertex is to be rejected later
        lacctot = .false.                   ! remember that at least one vertex has to be rejected.
      endif ! vdist > tolvdist


    ! reject vertices which were found inappropriate and re-index vertices in each plane:
      if (.not. lacctot) then
        ivertnewcount = 0
        do ivert = 1, nvert(iplane)
          if (lacceptvert(ivert)) then
            ivertnewcount = ivertnewcount + 1    ! one more vertex to accept 
            if (ivertnewcount /= ivert) then     ! otherwise the correct value is already at the correct place
              vert(1:3,ivertnewcount,iplane) = vert(1:3,ivert,iplane) ! re-index vertex
            endif
          endif ! lacceptvert(ivert)
        enddo ! ivert
        nvert(iplane) = ivertnewcount
      endif ! not lacctot

    enddo ! iplane



    ! now analyze faces, reject faces with less than three vertices and faces of very small area.
    do iplane = 1, nplane
      if (nvert(iplane) >= 3) then  ! calculate area
        x1 = vert(1,1,iplane)
        y1 = vert(2,1,iplane)
        z1 = vert(3,1,iplane)
        facearea(iplane) = 0.d0
        do ivert = 2, nvert(iplane)-1
          x2 = vert(1,ivert,iplane)
          y2 = vert(2,ivert,iplane)
          z2 = vert(3,ivert,iplane)
          x3 = vert(1,ivert+1,iplane)
          y3 = vert(2,ivert+1,iplane)
          z3 = vert(3,ivert+1,iplane)
          trianglearea = 0.5d0 * dabs( (x2 - x1)*(x3 - x1) + (y2 - y1)*(y3 - y1) + (z2 - z1)*(z3 - z1) )  ! formula incorrect? e.r.
          facearea(iplane) = facearea(iplane)+ trianglearea
        enddo ! ivert

        if (output) write(*,fmt="('analyzevert3d: face',i5,' has area',e12.4)") iplane,facearea(iplane)

        if (facearea(iplane) >= tolarea) then
          lacceptface(iplane) = .true.
        else
          lacceptface(iplane) = .false.  ! reject faces with small area
          if (output) write(*,fmt="('face will be rejected ; max. area tolerance=',e12.4)") tolarea 
        endif
          
      else
        lacceptface(iplane) = .false.  ! reject planes with less than 3 vertices
        if (output) write(*,fmt="('plane',i5,' has only',i3,' vertices and is rejected')") iplane,nvert(iplane)
      endif

    enddo ! iplane


    ! re-order the faces so that the accepted ones are in the first nface positions (nface is recalculated);
    ! the rest of the array entries are not taken care of, and can contain garbage.
    ifacenewcount = 0
    do iplane = 1, nplane
      if (lacceptface(iplane)) then
        ifacenewcount = ifacenewcount + 1  ! one more face to accept
        if (ifacenewcount /= iplane) then   ! otherwise the correct value is already at the correct place
          nvert(ifacenewcount) = nvert(iplane) ! re-index face vertex number
          do ivert = 1, nvert(iplane)
            vert(1:3,ivert,ifacenewcount) = vert(1:3,ivert,iplane) ! re-index face vertices
          enddo ! ivert
          planes(0:3,ifacenewcount) = planes(0:3,iplane) ! re-index face equation parameters
        endif ! ifacenewcount /= iplane
      endif ! lacceptface(iplane)
    enddo ! iplane
    nface = ifacenewcount

    ! check for every face that all vertices lie on the same plane by checking linear dependence
    do iface = 1, nface
      x2 = vert(1,2,iface) - vert(1,1,iface)
      y2 = vert(1,2,iface) - vert(2,1,iface)
      z2 = vert(1,2,iface) - vert(3,1,iface)
      x3 = vert(1,3,iface) - vert(1,1,iface)
      y3 = vert(1,3,iface) - vert(2,1,iface)
      z3 = vert(1,3,iface) - vert(3,1,iface)
      detsum = 0.d0
      do ivert = 4, nvert(iface)
        x4 = vert(1,ivert,iface) - vert(1,1,iface)
        y4 = vert(1,ivert,iface) - vert(2,1,iface)
        z4 = vert(1,ivert,iface) - vert(3,1,iface)
        det = x2*(y3*z4 - y4*z3) + y2*(z3*x4 - z4*x3) + z2*(x3*y4 - x4*y3)
        detsum = detsum + abs(det)
        if (abs(det) > 1.d-16 .and. output) write(*,fmt="('error from analyzevert3d: vertices not on single plane. iface=',i5,' ivert=',i5,' determinant=',e12.4)") iface,ivert,det
      enddo ! ivert
      if (output) write(*,fmt="('analyzevert3d: checking that vertices lie on plane. iface=',i5,' ; determinants sum to be zero=',e12.4)") iface,detsum
    enddo ! iface

  endsubroutine analyzevert3d

  !***********************************************************************
  logical function halfspace(a,b,c,d,x,y,z)
  ! given a plane a*x+b*y+c*z=d, and a point (x,y,z) in space, this 
  ! function takes the value true if (x,y,z) lies in the half-space 
  ! defined by the plane and the origin (0,0,0) (including the plane 
  ! itself). else, the value false is returned.
  !
  ! the criterion used is that the inner product of the vector (x,y,z) 
  ! with the vector d connecting the origin to the plane vertically be 
  ! less than or equal to d**2:  (d_x,d_y,d_z)*(x,y,z) =< d**2.
    double precision, intent(in) :: a,b,c,d,x,y,z

    if (dabs(a)+dabs(b)+dabs(c) < 1.d-80) die_here('halfspace: a,b,c too small.')

    halfspace = (d*(a*x + b*y + c*z) <= d*d)
    ! (re-checked 31may2008 fm)

  endfunction ! halfspace
  
  logical function half_space(plane, vec)
    double precision, intent(in) :: plane(0:3), vec(1:3)
    half_space = halfspace(plane(1), plane(2), plane(3), plane(0), vec(1), vec(2), vec(3))
  endfunction ! halfspace

  !***********************************************************************
  subroutine normalplane(x1,y1,z1,x2,y2,z2,tau,a,b,c,d)
    ! given two points in space, r1=(x1,y1,z1) and r2=(x2,y2,z2), this
    ! subroutine returns the coefficients defining a plane through the 
    ! equation a*x+b*y+c*z=d, which is normal to the vector r2-r1 and passes
    ! through the point (1.-tau)*r1 + tau*r2 (tau thus being a parameter
    ! defining how close the plane is to each of the two points).
    double precision, intent(in) :: x1,y1,z1,x2,y2,z2,tau
    double precision, intent(out) :: a,b,c,d
    
    double precision :: onemtau
    ! the plane is defined as 
    ! (a,b,c)*(x-x1,y-y1,z-z1)=const=
    !                         =(distance from r1 to (1.-tau)*r1 + tau*r2)**2
    ! so a,b,c are the coords. of a vector connecting the point r1 to
    ! the point (1.-tau)*r1 + tau*r2.
    onemtau = 1.d0 - tau

    a = onemtau*x1 + tau*x2
    b = onemtau*y1 + tau*y2
    c = onemtau*z1 + tau*z2
    d = a*(a + x1) + b*(b + y1) + c*(c + z1)

  endsubroutine ! normal_plane
  
  !***********************************************************************
  subroutine normalplane0(x1,y1,z1,tau, a,b,c,d)
    ! given a point in space, r1=(x1,y1,z1), this
    ! subroutine returns the coefficients defining a plane through the 
    ! equation a*x+b*y+c*z=d, which is normal to the vector r1 and passes
    ! through the point tau*r1 (tau thus being a parameter
    ! defining how close the plane is to the point).
    double precision, intent(in) :: x1,y1,z1,tau
    double precision, intent(out) :: a,b,c,d
    ! the plane is defined by
    ! (a,b,c)*(x,y,z) = d = (tau * r1)**2
    ! so a,b,c are the coords. of the vector tau * r1.
    ! if tau=0 (plane passes through the origin), then d=0.

    if (tau /= 0.d0) then
      a = tau*x1
      b = tau*y1
      c = tau*z1
      d = a*a + b*b + c*c
    else
      a = x1
      b = y1
      c = z1
      d = 0.d0
    endif

  endsubroutine ! normal_plane
  
  function normal_plane0f(v1, tau) result(plane)
    double precision :: plane(0:3)
    double precision, intent(in) :: v1(3), tau
    
    call normalplane0(v1(1), v1(2), v1(3), tau, plane(1), plane(2), plane(3), plane(0))

  endfunction ! normal_plane

  function normal_planef(v1, v2, tau) result(plane)
    double precision :: plane(0:3)
    double precision, intent(in) :: v1(3), v2(3), tau
    
    call normalplane(v1(1), v1(2), v1(3), v2(1), v2(2), v2(3), tau, plane(1), plane(2), plane(3), plane(0))

  endfunction ! normal_plane
  

  subroutine sortvertices(n, s, xyz)
    ! sorts the array s(n) in ascending order using straight insertion. 
    ! the arrays z(n), y(n), and z(n) follow.
    ! on output, arrays s, x, y, and z return sorted.
    integer, intent(in) :: n
    double precision, intent(inout) :: s(:), xyz(:,:)
    double precision :: tmp(0:3)
    integer :: i, j

    outer: do j = 2, n
      tmp = [s(j), xyz(1,j), xyz(2,j), xyz(3,j)]
      do i = j-1, 1, -1
          if (s(i) <= tmp(0)) then
            s(i+1) = tmp(0)
            xyz(1:3,i+1) = tmp(1:3)
            cycle outer
          endif
          s(i+1) = s(i)
          xyz(1:3,i+1) = xyz(1:3,i)
      enddo ! i
      s(1) = tmp(0)
      xyz(1:3,1) = tmp(1:3)
    enddo outer ! j

  endsubroutine ! sort vertices

  function crospr(x, y) result(z)
    ! crosp computes the cross product of x and y returning it into z.
    double precision :: z(3) ! result
    double precision, intent(in) :: x(3), y(3)
    z(1) = x(2)*y(3) - x(3)*y(2)
    z(2) = x(3)*y(1) - x(1)*y(3)
    z(3) = x(1)*y(2) - x(2)*y(1)
  endfunction crospr

  double precision function distplane(a,b,c,d)
    ! returns the distance of a plane a*x+b*y+c*z=d to the origin.
    double precision, intent(in) :: a,b,c,d
    double precision :: abcsq

    abcsq = a*a + b*b + c*c

    assert(abcsq >= 1.d-100)

    distplane = dabs(d)/dsqrt(abcsq)  
    ! or: distplane = dsqrt(d*d/abcsq)  

  endfunction distplane

  double precision function dist_plane(abc, d)
    ! returns the distance of a plane a*x+b*y+c*z=d to the origin.
    double precision, intent(in) :: abc(3), d
    double precision :: abcsq

    abcsq = abc(1)**2 + abc(2)**2 + abc(3)**2

    assert(abcsq >= 1.d-100)

    dist_plane = dsqrt(d*d/abcsq)  

  endfunction dist_plane

endmodule ! Voronoi_mod
