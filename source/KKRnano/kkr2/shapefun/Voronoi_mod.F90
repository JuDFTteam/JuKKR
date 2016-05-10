
module Voronoi_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  public :: Voronoi_construction
  public :: jellstart12

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
  subroutine Voronoi_construction(nvec, rvec, nvertmax, nfaced, weight0, weights, tolvdist, tolarea, &
     rmt, rout, volume, nface, planes, nvert, vert, atom_id, output)
    use Sorting_mod, only: dsort

    integer, intent(in) :: nvec, nvertmax, nfaced
    double precision, intent(in) :: rvec(:,:) ! cluster positions without the origin among them
    double precision, intent(in) :: weight0, weights(:)  ! Weight of central atom, and of all others (dimensioned as RVEC). 
    double precision, intent(in) :: tolvdist ! Max. tolerance for distance of two vertices
    double precision, intent(in) :: tolarea  ! Max. tolerance for area of polygon face
    integer, intent(in) :: atom_id
    logical, intent(in) :: output ! test output

    integer, intent(out) :: nface, nvert(:)
    double precision, intent(out) :: volume, rmt, rout
    double precision, intent(out) :: planes(0:,:)
    double precision, intent(out) :: vert(:,:,:) ! (3,nvertmax,nfaced) ! vertices
    
    integer :: ivec, iface, ivert, i, nplane
    double precision :: tetrvol, rsq, tau, v1(3), v2(3), v3(3), rout2, dist
    double precision :: facearea(nfaced), trianglearea   
    integer :: isort(nfaced) ! Index for sorting
    double precision :: rsort(nfaced)              ! Aux. function for sorting
    double precision, allocatable :: ptmp(:,:), vtmp(:,:,:) ! For sorting
    integer, allocatable :: ntmp(:)                        ! For sorting


    ! Check that the origin is not included in RVEC.
    do ivec = 1, nvec
      if (sum(abs(rvec(1:3,ivec))) < 1.d-12) die_here("vector #"-ivec+"is zero!")
    enddo ! ivec

    ! Define the planes as normal to the vectors RVEC, passing from t as in eq. (2) above:
    do ivec = 1, nvec
      rsq = sum(rvec(1:3,ivec)**2)
      tau = 0.5d0*((weight0 - weights(ivec))/rsq + 1.d0)
      planes(0:3,ivec) = normal_plane(rvec(1:3,ivec), tau)
    enddo ! ivec

    ! Find the Voronoi polyhedron.
    nplane = nvec
    if (output) write(*,'(a,i0)') ' Entering POLYHEDRON08 for atom #',atom_id
    call polyhedron08(nplane, nvertmax, nfaced, tolvdist, tolarea, planes, nface, nvert, vert, output)
    if (output) write(*,'(a,i0)') ' Exited POLYHEDRON08 for atom #',atom_id


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
          facearea(iface) = facearea(iface) + trianglearea
      enddo ! ivert
    enddo ! iface
    volume = volume/6.d0

    if (output) then
      write(*,*) ' Polyhedron properties '
      write(*,*) ' Number of faces : ',nface
      do iface = 1, nface
        write(*,fmt="(' Face ',i4,'   has ',i4,' vertices ','; Area= ',e12.4)") iface, nvert(iface), facearea(iface)
        do ivert = 1, nvert(iface)
          write(*,fmt="(i5,4e16.8)") ivert, vert(1:3,ivert,iface)
        enddo ! ivert
        write(*,fmt="(' Face coefficients:',4e16.8)") planes(1:3,iface), planes(0,iface)
      enddo ! iface
      write(*,*) 'The Volume is : ',volume
    endif ! output
    
    ! Fint rout, the largest radius of all vertices
    rout2 = 0.d0
    do iface = 1, nface
      do ivert = 1, nvert(iface)
        rout2 = max(rout2, sum(vert(1:3,ivert,iface)**2))
      enddo ! ivert
    enddo ! iface
    rout = sqrt(rout2)

    ! Find rmt and preprare sorting of faces
    rmt = dist_plane(planes(0:3,1)) ! init with the distance of the 1st plane
    do iface = 1, nface
      dist = dist_plane(planes(0:3,iface))
      rsort(iface) = 1.d9*dist + dot_product([1.d0, 1.d3, 1.d6], planes(1:3,iface)) ! first criterion is the distance, 2nd is z-coords, 3rd y, 4th x
      rmt = min(rmt, dist)
    enddo ! iface

    if (output) &
      write(*, fmt="('Voronoi subroutine: RMT=',e16.8,'; ROUT=',e16.8,'; RATIO=',f12.2,' %, atom #',i0)") rmt,rout,rmt*100/rout,atom_id
    
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

  endsubroutine ! Voronoi_construction
  

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

  endsubroutine ! polyhedron08


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

    integer :: p1, p2, p3, ip ! plane indices
    integer :: ivert          ! vertex index
    integer :: nwarn
    double precision :: cut(1:3) ! cut point of three planes.
    double precision :: det(0:3) ! determinants of 3x3 system
    double precision :: tolvdist2, distance2   ! a distance criterium of two points in space
    ! the following are for sorting the vertices of each face:
    double precision :: v1(3), v2(3), v3(3)   ! auxiliary vectors...
    double precision :: sinfiv1v2, cosfiv1v2 ! ...and their inner and outer products
    double precision :: phi(nvertmax)        ! ...and also their relative angles.
    double precision :: uv(3), vl, sn         ! unit vector, length, sign of sin(phi)
    logical :: laccept ! determining whether a cut point is inside the polyhedron.

    if (nplane < 4) die_here("nplane was only"+nplane)

    nwarn = 0
    tolvdist2 = tolvdist**2
    
    nvert(1:nplane) = 0
    !===============================================================
    ! start loop over all planes that can be cut:
    plane1: do p1 = 1, nplane
      ! start loop over all other planes:
      plane2: do p2 = 1, nplane
      if (p2 == p1) cycle plane2
        ! start loop over all other-other (!) planes. do from p2+1 to 
        ! nplane so that no pair is considered twice.
        plane3: do p3 = p2+1, nplane        ! nikos  p2+1,nplane
          if (p3 == p1) cycle plane3
          !     if (p3 == p2) goto 100 ! added by nikos
          ! solve the 3x3 system to find the cut point.
          !-----------------------------------
#define a(I) planes(1,I)
#define b(I) planes(2,I)
#define c(I) planes(3,I)
#define d(I) planes(0,I)

          det(0) = a(p1)*(b(p2)*c(p3) - b(p3)*c(p2)) + a(p2)*(b(p3)*c(p1) - b(p1)*c(p3)) + a(p3)*(b(p1)*c(p2) - b(p2)*c(p1))

          if (dabs(det(0)) <= 1.d-12) cycle plane3 ! there is no cut point 

          det(1) = d(p1)*(b(p2)*c(p3) - b(p3)*c(p2)) + d(p2)*(b(p3)*c(p1) - b(p1)*c(p3)) + d(p3)*(b(p1)*c(p2) - b(p2)*c(p1))
          det(2) = a(p1)*(d(p2)*c(p3) - d(p3)*c(p2)) + a(p2)*(d(p3)*c(p1) - d(p1)*c(p3)) + a(p3)*(d(p1)*c(p2) - d(p2)*c(p1))
          det(3) = a(p1)*(b(p2)*d(p3) - b(p3)*d(p2)) + a(p2)*(b(p3)*d(p1) - b(p1)*d(p3)) + a(p3)*(b(p1)*d(p2) - b(p2)*d(p1))
          
#undef d
#undef c
#undef b
#undef a

          cut(1:3) = det(1:3)/det(0)

          !     write(*,333) p1,p2,p3,cut
          !333  format('cutting point of planes ',3i5,':',3d15.7)
          !-----------------------------------
          ! accept this cut point as a vertex, if it belongs to the polyhedron. so,
          ! make a loop over all other (than p1,2,3) planes:
          laccept = .true.
          do ip = 1, nplane
            if (ip == p1 .or. ip == p2 .or. ip == p3) cycle
      !     laccept = laccept .and. halfspace(a(ip),b(ip),c(ip),d(ip),cut(1),cut(2),cut(3))
            laccept = laccept .and. half_space(planes(0:3,ip), cut(1:3))
          enddo ! ip

          !-----------------------------------
          if (laccept) then
            ! if the cut point found belongs to the cell, we accept it unless it has
            ! occured before for this face (p1). such a situation is possible
            ! when 4 or more planes pass through the same point (e.g. for the vertices
            ! of the fcc ws-cell). so...
            do ivert = 1, nvert(p1)
              distance2 = sum((vert(1:3,ivert,p1) - cut(1:3))**2)
              if (distance2 < tolvdist2) then
                laccept = .false. ! vertex is too close to a previous one.
                exit              ! jump loop, no need to continue.
              endif
            enddo ! ivert
            
          endif ! laccept

          ! now we are ready to add the point to the vertex list.
          if (laccept) then
            nvert(p1) = nvert(p1) + 1
            vert(1:3,nvert(p1),p1) = cut(1:3)
          endif

        enddo plane3 ! p3
      enddo plane2 ! p2
    !     write(*,*) 'number of vertices for plane ',p1,'  :',nvert(p1)
    enddo plane1 ! p1

    !===============================================================
    ! each plane should finally have either at least 3 vertices, if it is a
    ! face of the polyhedron, or none at all. check this:
    do ip = 1, nplane
      if (nvert(ip) == 1 .or. nvert(ip) == 2) then
        if (output) then ! todo: this should happen in the very verbose case
          warn(6, "there is a problem with the vertices for plane #"-ip-', only'+nvert(ip)+'vertices were found!')
        else
          nwarn = nwarn+1
        endif
      endif
    enddo ! ip

    !===============================================================
    ! for each face of the polyhedron, sort the vertices in a consecutive order
    ! as vertices of the polygon. the order is not necessarily mathematically
    ! positive.
    nface = 0
    do ip = 1, nplane
      if (nvert(ip) >= 3) then
        nface = nface + 1      ! count the faces
        phi(1) = -4.d0          ! just a number smaller than -pi.
      ! unit vector in the direction of first vertex:
        vl = sqrt(sum(vert(1:3,1,ip)**2))
        uv(1:3) = vert(1:3,1,ip)/vl

      ! define the vector connecting the first vertex to the (now-) second:
        v1(1:3) = vert(1:3,2,ip) - vert(1:3,1,ip)

        do ivert = 2, nvert(ip)
    ! define the vector connecting the first vertex to the current one:
          v2(1:3) = vert(1:3,ivert,ip) - vert(1:3,1,ip)
    ! find the angle phi between v1 and v2
    ! ( always, -pi < phi < pi from the definition of datan2 )
          cosfiv1v2 = dot_product(v1(1:3), v2(1:3))
          v3 = crospr(v1, v2) ! cross product = |v1|*|v2|*sinfi
          sinfiv1v2 = sqrt(sum(v3(1:3)**2))
    ! sign of sinfi is defined with respect to unit vector uv (see above)
          sn = dot_product(uv(1:3), v3(1:3))
          if (sn < 0) sinfiv1v2 = -sinfiv1v2

          if (sinfiv1v2 == 0.d0 .and. cosfiv1v2 == 0.d0) die_here("found two identical vertex points")
!             then
!       ! point falls exactly on 1st vertex...
!               phi(ivert) = -4.d0
!               die_here("found two identical vertex points")
!       ! ...while it shouldn't ! (this was checked earlier)
!             else
!               phi(ivert) = datan2(sinfiv1v2, cosfiv1v2)
!             endif
          phi(ivert) = datan2(sinfiv1v2, cosfiv1v2)
           
        enddo ! ivert = 3,nvert(ip)

        ! store with respect to the angle found:
        call sortvertices(nvert(ip), phi(1:nvert(ip)), vert(:,1:nvert(ip),ip)) ! sort the vertices in-place

      endif ! (nvert(ip) >= 3)
    enddo ! ip = 1, nplane

    if (output .and. nwarn > 0) & !  ! todo: this should happen in the medium verbose case
      warn(6, "there is a problem with the vertices for"+nwarn+"planes: only 1 or 2 vertices were found!")
    
    endsubroutine ! vertex3d

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

    double precision :: vdist2, tolvdist2 ! distance between consecutive vertices
    logical :: lacceptvert(nvertmax),lacctot,lfoundnext,lthisisthelast
    logical :: lacceptface(nfaced)
    integer :: iface, ivert, ivert2, inext, iplane, ifacenewcount, ivertnewcount
    double precision :: dv(3), v1(3), v2(3), v3(3), v4(3), trianglearea
    double precision :: facearea(nfaced)
    double precision :: det, detsum

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
        v1(1:3) = vert(1:3,1,iplane)
        facearea(iplane) = 0.d0
        do ivert = 2, nvert(iplane)-1
          v2(1:3) = vert(1:3,ivert,iplane)
          v3(1:3) = vert(1:3,ivert+1,iplane)
          trianglearea = 0.5d0*abs((v2(1) - v1(1))*(v3(1) - v1(1)) + (v2(2) - v1(2))*(v3(2) - v1(2)) + (v2(3) - v1(3))*(v3(3) - v1(3)))  ! formula incorrect? e.r.
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
      v2(1:3) = [vert(1,2,iface) - vert(1,1,iface), vert(1,2,iface) - vert(2,1,iface), vert(1,2,iface) - vert(3,1,iface)]
      v3(1:3) = [vert(1,3,iface) - vert(1,1,iface), vert(1,3,iface) - vert(2,1,iface), vert(1,3,iface) - vert(3,1,iface)]
      detsum = 0.d0
      do ivert = 4, nvert(iface)
        v4(1:3) = vert(1,ivert,iface) - vert(1:3,1,iface) ! yesss, the left argument is a scalar here
        det = v2(1)*(v3(2)*v4(3) - v4(2)*v3(3)) + v2(2)*(v3(3)*v4(1) - v4(3)*v3(1)) + v2(3)*(v3(1)*v4(2) - v4(1)*v3(2))
        detsum = detsum + abs(det)
        if (abs(det) > 1.d-16 .and. output) write(*, fmt="('error from analyzevert3d: vertices not on single plane. iface=',i5,' ivert=',i5,' determinant=',e12.4)") iface,ivert,det
      enddo ! ivert
      if (output) write(*,fmt="('analyzevert3d: checking that vertices lie on plane. iface=',i5,' ; determinants sum to be zero=',e12.4)") iface,detsum
    enddo ! iface

  endsubroutine ! analyzevert3d

!   !***********************************************************************
!   logical function halfspace(a,b,c,d,x,y,z)
!   ! given a plane a*x+b*y+c*z=d, and a point (x,y,z) in space, this 
!   ! function takes the value true if (x,y,z) lies in the half-space 
!   ! defined by the plane and the origin (0,0,0) (including the plane 
!   ! itself). else, the value false is returned.
!   !
!   ! the criterion used is that the inner product of the vector (x,y,z) 
!   ! with the vector d connecting the origin to the plane vertically be 
!   ! less than or equal to d**2:  (d_x,d_y,d_z)*(x,y,z) =< d**2.
!     double precision, intent(in) :: a,b,c,d,x,y,z
! 
!     if (dabs(a)+dabs(b)+dabs(c) < 1.d-80) die_here('halfspace: a,b,c too small.')
! 
!     halfspace = (d*(a*x + b*y + c*z) <= d*d) !!! re-checked 31may2008 FM
! 
!   endfunction ! halfspace
  
  logical function half_space(plane, vec)
    double precision, intent(in) :: plane(0:3), vec(1:3)
!   half_space = halfspace(plane(1), plane(2), plane(3), plane(0), vec(1), vec(2), vec(3))
    if (sum(abs(plane(1:3))) < 1.d-80) die_here('halfspace: a,b,c too small.')

    half_space = (plane(0)*dot_product(plane(1:3), vec(1:3)) <= plane(0)**2) !!! re-checked 31may2008 FM
  endfunction ! halfspace
  
  
  function normal_plane0f(v1, tau) result(plane)
    double precision :: plane(0:3) ! former [d,a,b,c]
    double precision, intent(in) :: v1(3), tau
!     ! given a point in space, r1=(v1(1),v1(2),v1(3)), this
!     ! subroutine returns the coefficients defining a plane through the 
!     ! equation a*x+b*y+c*z=d, which is normal to the vector r1 and passes
!     ! through the point tau*r1 (tau thus being a parameter
!     ! defining how close the plane is to the point).
    
!     ! so a,b,c are the coords. of the vector tau * r1.
!     ! if tau=0 (plane passes through the origin), then d=0.

    if (tau /= 0.d0) then
      plane(1:3) = tau*v1(1:3)
      plane(0) = sum(plane(1:3)**2)
    else
      plane(1:3) = v1(1:3)
      plane(0) = 0.d0
    endif

  endfunction ! normal_plane
  
  function normal_planef(v1, v2, tau) result(plane)
    double precision :: plane(0:3) ! former [d,a,b,c]
    double precision, intent(in) :: v1(3), v2(3), tau
    ! given two points in space, r1=(v1(1),v1(2),v1(3)) and r2=(v2(1),v2(2),v2(3)), this
    ! subroutine returns the coefficients defining a plane through the 
    ! equation a*x+b*y+c*z=d, which is normal to the vector r2-r1 and passes
    ! through the point (1.-tau)*r1 + tau*r2 (tau thus being a parameter
    ! defining how close the plane is to each of the two points).
    
    ! the plane is defined as 
    ! (a,b,c)*(x-v1(1),y-v1(2),z-v1(3))=const=
    !                         =(distance from r1 to (1.-tau)*r1 + tau*r2)**2
    ! so a,b,c are the coords. of a vector connecting the point r1 to
    ! the point (1.-tau)*r1 + tau*r2.

    plane(1:3) = (1.d0 - tau)*v1(1:3) + tau*v2(1:3)
    plane(0) = dot_product(plane(1:3), plane(1:3) + v1(1:3))

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
  endfunction ! crospr

 
  double precision function dist_plane(dabc)
    ! returns the distance of a plane a*x+b*y+c*z=d to the origin.
    double precision, intent(in) :: dabc(0:3)
    double precision :: abcsq

    abcsq = sum(dabc(1:3)**2)
    assert(abcsq >= 1.d-100)
    dist_plane = sqrt(dabc(0)**2/abcsq)  

  endfunction ! dist_plane


!---------- Routines for creation of Jellium potentials -----------

SUBROUTINE jellstart12(nspin,ins,natoms,z,idshape,  &
        rwscl,rmtcl,meshn,xrn,drn,  &
        irws,irns,  &
        alatnew,qbound,dims,atom_index)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2016-01-21  Time: 16:35:08
 
! ******************************************************
! * This subroutine reads a jellium potential from the database
! * file. and interpolates to the new mesh
! ******************************************************

use DimParams_mod, only: DimParams

IMPLICIT NONE

INTEGER, INTENT(IN)                      :: nspin
INTEGER, INTENT(IN)                      :: ins
INTEGER, INTENT(IN)                      :: natoms
REAL*8, INTENT(IN)                       :: z(:)
INTEGER, INTENT(IN)                      :: idshape(*)
REAL*8, INTENT(IN)                       :: rwscl(*)
REAL*8, INTENT(IN)                       :: rmtcl(*)
INTEGER, INTENT(IN)                      :: meshn(:)
REAL*8, INTENT(IN)                       :: xrn(:,:)
REAL*8, INTENT(IN)                       :: drn(:,:)
INTEGER, INTENT(IN)                      :: irws(:)
INTEGER, INTENT(IN)                      :: irns(:)
REAL*8, INTENT(IN)                       :: alatnew
REAL*8, INTENT(IN OUT)                   :: qbound
type(DimParams), intent(in)              :: dims
INTEGER, INTENT(IN)                      :: atom_index


! Parameters that are no longer taken from 'inc.geometry' but depend on dims
INTEGER :: npotd,lmpotd,irmind,inslpd,lmxspd,irmdjj

!     ..
!     .. Scalar Arguments ..
REAL*8           efermi
INTEGER :: kshape
!     ..
!     .. Array Arguments ..

INTEGER :: lcore(20),ncore
!     ..
!     .. Local Scalars ..
REAL*8           ea,s1,z1,  &
    vbc(2),ecore1(20),maxa,aout,bout,rmtout,  &
    parsum,parsumderiv,r0,rmaxout,rmtnew
REAL*8           rws0,br,za,zvali,einf,ar,amsh
INTEGER :: i,ir,iri,  &
    ispin,  &
    ipot,id,lm,lm1,irnsout,irmtout,irwsout,  &
    nr,iat,ncore1,lcore1(20),irc,nz,  &
    irs1,nsec,nzvali,nc,ii,i1,i2
LOGICAL, allocatable :: potlm(:)
!     ..
!     .. Local Arrays ..
REAL*8, allocatable  ::    u(:),drdi(:),ecore(:),  &
                           rmesh(:),vins(:,:),  &
                           vm2z(:),vinsout(:,:),  &
                           vm2zout(:),vm2zb(:),rout(:), vinsb(:,:),drdiout(:),  &
                           work(:,:),ra(:)
CHARACTER (LEN=40) :: baner
CHARACTER (LEN=4) :: aaaa,tran

CHARACTER (LEN=4) :: elem_file(0:113)
CHARACTER (LEN=26) :: atompot
CHARACTER (LEN=2) :: txtc(20)
character(len=17) :: filename

DATA elem_file/'Vac0',  &
    'H_01','He02','Li03','Be04','B_05','C_06','N_07','O_08', 'F_09','Ne10',  &
    'Na11','Mg12','Al13','Si14','P_15','S_16','Cl17','Ar18',  &
    'K_19','Ca20','Sc21','Ti22', 'V_23','Cr24','Mn25','Fe26','Co27','Ni28',  &
    'Cu29','Zn30', 'Ga31','Ge32','As33','Se34','Br35','Kr36','Rb37','Sr38',  &
    'Y_39','Zr40', 'Nb41','Mo42','Tc43','Ru44','Rh45','Pd46','Ag47','Cd48',  &
    'In49','Sn50', 'Sb51','Te52','I_53','Xe54','Cs55','Ba56','La57','Ce58',  &
    'Pr59','Nd60', 'Pm61','Sm62','Eu63','Gd64','Tb65','Dy66','Ho67','Er68',  &
    'Tm69','Yb70', 'Lu71','Hf72','Ta73','W_74','Re75','Os76','Ir77','Pt78',  &
    'Au79','Hg80', 'Tl81','Pb82','Bi83','Po84','At85','Rn68','Fr87','Ra88',  &
    'Ac89','Th90', 'Pa91','U_92','Np93','Pu94','Am95','Cm96','Bk97','Cf98',  &
    'Es99','Fm__', 'Md__','No__','Lr__','Rf__','Db__','Sg__','Bh__','Hs__',  &
    'Mt__','Uun_', 'Uuu_','Uub_','NoE_'/
!     --------------------------------------------------------------

! set parameters depending on 'dims' and 'num_local_atoms'
npotd=dims%nspind*natoms
lmpotd= (dims%lpot+1)**2
irmind=dims%irmd-dims%irnsd
inslpd= (dims%irnsd+1)*lmpotd
lmxspd= (2*dims%lpot+1)**2
irmdjj=1501
kshape=2  ! always full-pot calculations
WRITE(*,*) 'dims%irmd', dims%irmd

! allocate arrays
allocate(u(dims%irmd))
allocate(drdi(irmdjj))
allocate(ecore(20))
allocate(rmesh(irmdjj))
allocate(vins(irmind:dims%irmd,lmpotd))
allocate(vm2z(irmdjj))
allocate(vinsout(irmind:dims%irmd,lmpotd))
allocate(vm2zout(dims%irmd))
allocate(vm2zb(irmdjj))
allocate(rout(dims%irmd))
allocate(vinsb(dims%irmd,lmpotd))
allocate(drdiout(dims%irmd))
allocate(work(dims%irmd,lmpotd))
allocate(ra(dims%irmd))
allocate(potlm(lmpotd))


WRITE(6,*) ' ****  READING  POTENTIAL  **** '

!OPEN(19,STATUS='UNKNOWN',FILE='output.pot')
write(filename, fmt="(a,i7.7)") "potential.",atom_index
open(19, file=filename, form="formatted", action='write')

DO i2=1,lmpotd
  DO i1=irmind,dims%irmd
    vins(i1,i2) = 0.d0
  END DO
END DO
DO i1=1,irmdjj
  vm2z(i1) = 0.d0
END DO

DO iat = 1,natoms
  DO ispin=1,nspin
    DO lm=1,lmpotd
      potlm(lm) =.false.
    END DO
    ipot =  nspin* (iat-1) + ispin
    
! Find out what atom is needed
    
    nz = z(iat)
    IF (((nz >= 24.AND.nz <= 28).OR.(nz >= 57.AND.nz <= 70))  &
          .AND.ispin == 2) THEN
      atompot = 'ElementDataBase/'//elem_file(nz)//'.pots2'
    ELSE
      atompot = 'ElementDataBase/'//elem_file(nz)//'.pot  '
    END IF
    WRITE(6,*) 'Using database ....: ',atompot
    OPEN(21,STATUS='OLD',FILE=atompot,ERR=1010)
!           IRWS1 =  NR
    
! --------------------------------------------------------------------
    
    efermi = .409241D+00
    vbc(1)    = .500D0
    vbc(2)    = .500D0
!           read potential from jellium
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    READ(21,141) baner,aaaa,aaaa
    READ(21,142) rws0,s1,irs1,br
    READ(21,142) za,zvali,nsec,einf
!     Calculate number of core states
    nz = za             ! make integer
    nzvali = zvali      ! make integer
    nc = nz - nzvali
    ncore = 0
    IF (nc == 2 ) ncore = 1 ! 1s
    IF (nc == 4 ) ncore = 2 ! 1s2s
    IF (nc == 10) ncore = 3 ! 1s2s2p
    IF (nc == 12) ncore = 4 ! 1s2s2p3s
    IF (nc == 18) ncore = 5 ! 1s2s2p3s3p
    IF (nc == 28) ncore = 6 ! 1s2s2p3s3p3d
    IF (nc == 30) ncore = 7 ! 1s2s2p3s3p3d4s
    IF (nc == 36) ncore = 8 ! 1s2s2p3s3p3d4s4p
    IF (nc == 46) ncore = 9 ! 1s2s2p3s3p3d4s4p4d
    IF (nc == 48) ncore = 10 ! 1s2s2p3s3p3d4s4p4d4s
    IF (nc == 54) ncore = 11 ! 1s2s2p3s3p3d4s4p4d4s4p
    IF (nc == 68) ncore = 12 ! 1s2s2p3s3p3d4s4p4d4s4p4f
    IF (nc == 78) ncore = 13 ! 1s2s2p3s3p3d4s4p4d4s4p4f5d
    IF (nc == 80) ncore = 14 ! 1s2s2p3s3p3d4s4p4d4s4p4f5d6s
    IF (nc == 86) ncore = 15 ! 1s2s2p3s3p3d4s4p4d4s4p4f5d6s4p
    WRITE(6,*) '*************************************'
    WRITE(6,*) '   Potential Interpolation Program   '
    WRITE(6,*) '   Using the Jellium Database v1.0   '
    WRITE(6,*) '*************************************'
    WRITE(6,163) efermi
    WRITE(6,161) za
    WRITE (6,162) ncore
    READ(21,133)(lcore(nc),txtc(nc),ecore(nc),nc=1,ncore)
    WRITE(6,*) ' ** Position of the Core States ** '
    DO i=1,ncore
      WRITE(6,135) lcore(i),txtc(i),ecore(i)
      IF (txtc(i) == 's ') lcore(i) = 0
      IF (txtc(i) == 'p ') lcore(i) = 1
      IF (txtc(i) == 'd ') lcore(i) = 2
      IF (txtc(i) == 'f ') lcore(i) = 3
    END DO
    ncore1 = ncore
    DO i=1,ncore1
      lcore1(i) = lcore(i)
      ecore1(i) = ecore(i)
    END DO
    WRITE(6,134) einf
    WRITE(6,*) '**********************************************'
    
    READ(21,131)(vm2z(ii),ii=1,irs1)
    READ(21,132)tran
    CLOSE (21)
!     make mesh r0
    ar = LOG(s1/br+1.d0)/FLOAT(irs1-1)
    ea=EXP(ar)
    amsh=1.d0
    rmesh(1)=0.d0
    drdi(1)=br*ar
    DO  i=2,irs1
      amsh=amsh*ea
      rmesh(i)=br*(amsh-1.d0)
      drdi(i)=drdi(1)*amsh
    END DO
    WRITE(6,*) 'Jellium Potential Read In ',irs1,' points'
    nr = irs1
    z1 = za
    
    131        FORMAT(4D15.8)
    132        FORMAT(a4)
    133        FORMAT(i3,a2,d15.8)
    135        FORMAT(i3,a2,f15.6,' Ry')
    134        FORMAT('All other states are above :',f8.4,' Ry in Energy')
    141        FORMAT(3X,a40,3X,a4,3X,a4)
    142        FORMAT(7X,f8.4,7X,f8.4,7X,i5,7X,f8.4)
    161        FORMAT('Potential Atomic Number :',f7.2)
    162        FORMAT('Number of Core   States :',i4)
    163        FORMAT('Jellium Fermi Energy :',f10.5,' Ry')
! --------------------------------------------------------------------
    
!     The input mesh has been constructed. Now construct the output mesh.
    
    id = idshape(iat)
    rout(1) = 0.d0
    aout = 0.025d0
    rmaxout = rwscl(id)
    rmtout  = rmtcl(id)
    irwsout = irws(iat)
    irmtout = irws(iat) - meshn(id)
    irnsout = irns(iat)  ! 22.1.12 Changed from IRNS(ID) to IRNS(IAT)
    
    
    IF (ins == 0) THEN
      bout = rmaxout / (EXP(aout*REAL(irwsout-1))-1.0D0)
      DO ir=2,irwsout
        ea = EXP(aout*REAL(ir-1))
        rout(ir) = bout* (ea-1.0D0)
        drdiout(ir) = aout*bout*ea
      END DO
      DO i=1,irwsout
        IF (rout(i) < rmtout) irmtout = i
      END DO
      IF (MOD(irmtout,2) == 0) irmtout = irmtout+1
      rmtnew = rout(irmtout)
      rmaxout = rout(irwsout)
    ELSE
      bout = rmtout /  (EXP(aout*REAL(irmtout-1))-1.0D0)
      DO ir=2,irmtout
        ea = EXP(aout*REAL(ir-1))
        rout(ir) = bout* (ea-1.0D0)
        drdiout(ir) = aout*bout*ea
      END DO
      DO iri=1,meshn(id)
        ir = iri + irmtout
        rout(ir) = alatnew*xrn(iri,id)   ! scaling is always 1.0d0
        drdiout(ir) = alatnew*drn(iri,id)
      END DO
      rmtnew = rout(irmtout)
      rmaxout = rout(irwsout)
    END IF
    
!  Ok now interpolate
    
    
    maxa = 1.d35
    CALL spline(irmdjj,rmesh,vm2z,nr,maxa,maxa,vm2zb)
    
! OK with spline
    
    vm2zout(1) = vm2z(1)
    DO ir = 2,irwsout
      r0 = rout(ir)
      CALL splint(rmesh,vm2z,vm2zb,nr,r0,parsum,parsumderiv)
      vm2zout(ir) = parsum
    END DO
    
    
    
    
    IF (ins > 0) THEN
      irc = irwsout - irnsout
      DO lm1=1,lmpotd
        DO ir = irc,irwsout
          vinsout(ir,lm1) = 0.d0
        END DO
      END DO
    END IF
    
    CALL ritesone12(19,ispin,z1,alatnew,rmtout,rmtnew,rmaxout,  &
        rout,drdiout,vm2zout,irwsout,aout,bout,ins,irnsout,  &
        vinsout,qbound,irwsout,kshape,efermi,vbc,  &
        ecore1,lcore1,ncore1,elem_file(nz),nspin,dims)
    
! Next atom or next spin
    
  END DO
END DO

RETURN


!1000  WRITE(6,*) 'Error read file......... ',i13
WRITE(6,*) 'Error occured on atom... ',iat
STOP
1010  WRITE(6,*) ' Error in JELLSTART '
WRITE(6,*) ' Potential.............',elem_file(nz)
WRITE(6,*) ' Does not exist in the database'
STOP
8000 FORMAT (a40)
8010 FORMAT (3X,a4,26X,f8.3)
8011 FORMAT ('#  ',a4,'POTENTIAL             Z = ',f8.3)
8012 FORMAT ('#  ',a4,'POTENTIAL SPIN UP     Z=  ',f8.3)
8013 FORMAT ('#  ',a4,'POTENTIAL SPIN DOWN   Z=  ',f8.3)
8020  FORMAT(a40)
8030 FORMAT (4F12.8)
8040 FORMAT (1X,3I6)
8050 FORMAT (2D15.8)
8060 FORMAT (3F12.8)
8070 FORMAT (2I5)
9051 FORMAT (1P,4D20.12)
9000 FORMAT (7A4,6X,'  exc:',a24,3X,a10)
9010 FORMAT (3F12.8)
9020 FORMAT (f10.5,/,f10.5,2F15.10)
9030 FORMAT (i3,/,2D15.8,/,2I2)
9140 FORMAT (i5,1P,d20.11)
9040 FORMAT (f10.5,/,f10.5,2F15.10)
9050 FORMAT (i3,/,2D15.8,/,2I2)
9060 FORMAT (1P,2D15.6,1P,d15.8)
9160 FORMAT (10I5)
9061 FORMAT (1P,5D15.8)
9070 FORMAT (i5,1P,d20.11)
! 9080 FORMAT (10x,20a4)
9080 FORMAT (' < ',20A4)
9081 FORMAT (' <#',20A4)
9090 FORMAT (10I5)
9100 FORMAT (1P,4D20.13)
END SUBROUTINE jellstart12



SUBROUTINE ritesone12(ifile,is,z,alat,rmt,rmtnew,rws,  &
        r,drdi,vm2z,irws,a,b,ins,irns,  &
        vins,qbound,irc,kshape,efermi,vbc,ecore,  &
        lcore,ncore,elem_name,nspin,dims)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2016-01-21  Time: 16:59:42
 
! ************************************************************************
!      this subroutine stores in 'ifile' the necessary results
!      (potentials e.t.c.) to start self-consistency iterations

!      modified for the full potential case - if ins .gt. 0 there
!       is written a different potential card
!       if the sum of absolute values of an lm component of vins (non
!       spher. potential) is less than the given rms error qbound this
!       component will not be stored .

!        (see to subroutine start , where most of the arrays are
!         described)

!                            modified by b. drittler  aug. 1988
!-----------------------------------------------------------------------
!     .. Parameters ..


use DimParams_mod, only: DimParams

IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ifile
INTEGER, INTENT(IN)                      :: is
REAL*8, INTENT(IN)                       :: z
REAL*8, INTENT(IN)                       :: alat
REAL*8, INTENT(IN)                       :: rmt
REAL*8, INTENT(IN)                       :: rmtnew
REAL*8, INTENT(IN)                       :: rws
REAL*8, INTENT(IN)                       :: r(:)
REAL*8, INTENT(IN)                       :: drdi(:)
REAL*8, INTENT(IN)                       :: vm2z(:)
INTEGER, INTENT(IN)                      :: irws
REAL*8, INTENT(IN)                       :: a
REAL*8, INTENT(IN)                       :: b
INTEGER, INTENT(IN)                      :: ins
INTEGER, INTENT(IN)                      :: irns
REAL*8, INTENT(IN)                       :: qbound
INTEGER, INTENT(IN)                      :: irc
INTEGER, INTENT(IN)                      :: kshape
REAL*8, INTENT(IN)                       :: efermi
REAL*8, INTENT(IN)                       :: vbc(2)
REAL*8, INTENT(IN)                       :: ecore(20)
INTEGER, INTENT(IN)                      :: lcore(20)
INTEGER, INTENT(IN)                      :: ncore
CHARACTER (LEN=4), INTENT(IN)            :: elem_name
INTEGER, INTENT(IN)                      :: nspin
type(DimParams), intent(in)              :: dims
REAL*8, INTENT(IN)                       :: vins((dims%irmd-dims%irnsd):dims%irmd,(dims%lpot+1)**2)

!     ..
!     .. Local Scalars ..
REAL*8           a1,b1,rmax,rmt1,rmtnw1,rv,sum,z1
INTEGER :: icore,inew,ir,irmin,irns1,isave,j,lm,lmnr,lmpot, &
ncore1,nr,lpot
!     ..
!     .. Local Arrays ..
REAL*8, allocatable :: dradi(:),ecore1(:),ra(:),vm2za(:)
INTEGER :: lcore1(20)
!     ..
!     .. Intrinsic Functions ..
INTRINSIC SQRT
!     ..
! -------------------------------------------------------------------

allocate(dradi(dims%irmd))
allocate(ecore1(20))
allocate(ra(dims%irmd))
allocate(vm2za(dims%irmd))

lpot=dims%lpot
isave = 1
inew  = 1



lmpot = (lpot+1)* (lpot+1)

rmt1 = rmt
rmtnw1 = rmtnew
z1 = z
rmax = rws
IF (kshape == 0) THEN
  nr = irws
  
ELSE
  nr = irc
END IF

irns1 = irns
irmin = nr - irns1
a1 = a
b1 = b
ncore1 = ncore

DO  j = 1,nr
  ra(j) = r(j)
  dradi(j) = drdi(j)
  
!--->       store only lm=1 component of the potential
  
  vm2za(j) = vm2z(j)
END DO

IF (ncore1 >= 1) THEN
  
  DO  j = 1,ncore1
    lcore1(j) = lcore(j)
    ecore1(j) = ecore(j)
  END DO
END IF


IF (nspin == 1) THEN
  WRITE (ifile,FMT=8999) elem_name,''
ELSE
  IF (is == 1) THEN
    WRITE (ifile,FMT=8995) elem_name,''
  ELSE
    WRITE (ifile,FMT=8996) elem_name,''
  END IF
END IF
!     WRITE (IFILE,FMT=9000) (ITITLE(I),I=1,7),TXC(KXC+1)
WRITE (ifile,FMT=9010) rmt1,alat,rmtnw1
WRITE (ifile,FMT=9020) z1,rmax,efermi,vbc(is)
IF (nr <= 999) THEN
  WRITE (ifile,FMT=9030) nr,a1,b1,ncore1,inew
ELSE
  WRITE (ifile,FMT=9031) nr,a1,b1,ncore1,inew
END IF
IF (ncore1 >= 1) WRITE (ifile,FMT=9040) (lcore1(icore),  &
    ecore1(icore),icore=1,ncore1)


IF (ins == 0 ) THEN
  
!--->       store only the spherically averaged potential
!           (in mt or as - case)
!           this is done always for the host
  
  IF (inew == 0) THEN
    WRITE (ifile,FMT=9050) (ra(ir),dradi(ir),vm2za(ir),ir=1,nr)
  ELSE
    WRITE (ifile,FMT=9051) (vm2za(ir),ir=1,nr)
  END IF
  
ELSE
  
!--->     store the full potential , but the non spherical contribution
!         only from irns1 up to irws1 ;
!         remember that the lm = 1 contribution is multiplied
!         by a factor 1/sqrt(4 pi)
  
  WRITE (ifile,FMT=9060) nr,irns1,lmpot,isave
  WRITE (ifile,FMT=9070) (vm2za(ir),ir=1,nr)
  IF (lpot > 0) THEN
    lmnr = 1
    DO  lm = 2,lmpot
      sum = 0.0D0
      DO  ir = irmin,nr
        rv = vins(ir,lm)*ra(ir)
        sum = sum + rv*rv*dradi(ir)
      END DO
      
      IF (SQRT(sum) > qbound) THEN
        lmnr = lmnr + 1
        WRITE (ifile,FMT=9060) lm
        WRITE (ifile,FMT=9070) (vins(ir,lm),ir=irmin,nr)
      END IF
      
    END DO
    
!--->         write a one to mark the end
    
    IF (lmnr < lmpot) WRITE (ifile,FMT=9060) isave
  END IF
  
END IF
!write(6,*) ' Potential finished go on'
50   CONTINUE
60 CONTINUE
8995 FORMAT (a4,' POTENTIAL SPIN DOWN',10X,'  exc:',a24)
8996 FORMAT (a4,' POTENTIAL SPIN UP  ',10X,'  exc:',a24)
8999 FORMAT (a4,' POTENTIAL ',19X,'  exc:',a24)
9000 FORMAT (7A4,6X,'  exc:',a24,3X,a10)
9010 FORMAT (3F12.8)
9020 FORMAT (f10.5,/,f10.5,2F15.10)
9030 FORMAT (i3,/,2D15.8,/,2I2)
9031 FORMAT (i4,/,2D15.8,/,2I2)
9040 FORMAT (i5,1P,d20.11)
9050 FORMAT (1P,2D15.6,1P,d15.8)
9051 FORMAT (1P,4D20.12)
9060 FORMAT (10I5)
9070 FORMAT (1P,4D20.13)
END SUBROUTINE ritesone12


!***********************************************************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2016-01-12  Time: 14:48:44

SUBROUTINE spline(nmax,x,y,n,yp1,ypn,y2)

IMPLICIT NONE

INTEGER, INTENT(IN)                      :: nmax
REAL*8, INTENT(IN)                       :: x(nmax)
REAL*8, INTENT(IN)                       :: y(nmax)
INTEGER, INTENT(IN)                      :: n
REAL*8, INTENT(IN)                       :: yp1
REAL*8, INTENT(IN)                       :: ypn
REAL*8, INTENT(OUT)                      :: y2(nmax)


! Given arrays x(1:n) and  y(1:n) containing a tabulated function,
! i.e., y i = f(xi), with x1<x2<...<xN , and given values yp1 and ypn
! for the 1rst derivative of the interpolating function at points
! 1 and n, respectively, this routine returns an array y2(1:n) of
! length n which contains the second derivatives of the interpolating
! function at the tabulated points xi.
! If yp1 and/or ypn are equal to 1.e30 or larger, the routine is
! signaled to set the corresponding boundary condition for a natural
! spline, with zero second derivative on that boundary.
! Parameter: NMAX is the largest anticipated value of n.
! Taken from "Numerical Recipes in Fortran 77", W.H.Press et al.
INTEGER :: i,k
REAL*8          p,qn,sig,un,u(nmax)

IF (n > nmax) STOP 'SPLINE: n > NMAX.'
IF (yp1 > 0.99D30) THEN
! The lower boundary condition is set either to be "natural"
  y2(1) = 0.d0
  u(1) = 0.d0
ELSE
! or else to have a specified first derivative.
  y2(1) = -0.5D0
  u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
END IF

DO i = 2,n-1
! This is the decomposition loop of the tridiagonal algorithm. y2 and u
! are used for temporary storage of the decomposed factors.
  sig = (x(i)-x(i-1)) / (x(i+1)-x(i-1))
  p = sig * y2(i-1) + 2.d0
  y2(i) = (sig-1.d0)/p
  u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))  &
      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1)) / p
END DO

IF (ypn > 0.99D30) THEN
! The upper boundary condition is set either to be "natural"
  qn = 0.d0
  un = 0.d0
ELSE
! or else to have a specified 1rst derivative.
  qn = 0.5D0
  un = (3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
END IF
y2(n) = (un-qn*u(n-1)) / (qn*y2(n-1)+1.d0)
DO k = n-1,1,-1
! This is the backsubstitution loop of the tridiagonal algorithm.
  y2(k)=y2(k)*y2(k+1)+u(k)
END DO

RETURN
END SUBROUTINE spline


!***********************************************************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2016-01-12  Time: 14:48:43

SUBROUTINE splint(xa,ya,y2a,n,x,y,yderiv)

IMPLICIT NONE

REAL*8, INTENT(IN)                       :: xa(*)
REAL*8, INTENT(IN)                       :: ya(*)
REAL*8, INTENT(IN OUT)                   :: y2a(*)
INTEGER, INTENT(IN)                      :: n
REAL*8, INTENT(IN)                       :: x
REAL*8, INTENT(OUT)                      :: y
REAL*8, INTENT(OUT)                      :: yderiv


! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a
! function (with the xai's in order), and given the array y2a(1:n), which
! is the output from spline above, and given a value of x, this routine
! returns a cubic-spline interpolated value y and the derivative yderiv.
! Taken from "Numerical Recipes in Fortran 77", W.H.Press et al.
INTEGER :: k,khi,klo
REAL*8         a,b,h
! We will  nd the right place in the table by means of bisection.
! This is optimal if sequential calls to this routine are at random
! values of x. If sequential calls are in order, and closely
! spaced, one would do better to store previous values of
! klo and khi and test if they remain appropriate on the
! next call.
klo=1
khi=n
1    IF (khi-klo > 1) THEN
  k=(khi+klo)/2
  IF(xa(k) > x)THEN
    khi=k
  ELSE
    klo=k
  END IF
  GO TO 1
END IF
! klo and khi now bracket the input value of x.
h=xa(khi)-xa(klo)
! The xa's must be distinct.
IF (h == 0.d0) PAUSE 'bad xa input in splint'
! Cubic spline polynomial is now evaluated.
a = (xa(khi)-x)/h
b = (x-xa(klo))/h
y = a*ya(klo) + b*ya(khi) +  &
    ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi)) * (h**2)/6.d0
yderiv = (ya(khi)-ya(klo))/h -  &
    ((3.d0*a*a-1.d0)*y2a(klo) - (3.d0*b*b-1.d0)*y2a(khi))*h/6.d0

RETURN
END SUBROUTINE splint


endmodule ! Voronoi_mod
