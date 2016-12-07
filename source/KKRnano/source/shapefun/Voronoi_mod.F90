
module Voronoi_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  public :: Voronoi_construction
  
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
     rmt, rout, volume, nface, planes, nvert, vert, atom_id)
    use Sorting_mod, only: dsort

    integer, intent(in) :: nvec, nvertmax, nfaced
    double precision, intent(in) :: rvec(:,:) !> cluster positions without the origin among them
    double precision, intent(in) :: weight0, weights(:) !> Weight of central atom, and of all others (dimensioned as RVEC). 
    double precision, intent(in) :: tolvdist !> Max. tolerance for distance of two vertices
    double precision, intent(in) :: tolarea  !> Max. tolerance for area of polygon face
    integer, intent(in) :: atom_id !> global atom id passed for more helpful error messages
    integer, intent(out) :: nface !> number of faces
    integer, intent(out) :: nvert(:) !> number of vertices
    double precision, intent(out) :: volume !> volume of the Voronoi cell
    double precision, intent(out) :: rmt !> smallest distance of a plane
    double precision, intent(out) :: rout !> largest distance of a plane contributing
    double precision, intent(out) :: planes(0:,:) !> planes
    double precision, intent(out) :: vert(:,:,:) ! (3,nvertmax,nfaced) !> vertices
    
    ! locals
#ifdef DEBUGSHAPEFUNCTIONS
    logical, parameter :: output = .true. !> activates logging for debug
#else
    logical, parameter :: output = .false. !> deactivates logging for debug
#endif
    integer :: ivec, iface, ivert, i, nplane, isort(nfaced)
    double precision :: tetrvol, rsq, tau, rout2, dist, trianglearea
    double precision :: facearea(nfaced), rsort(nfaced), v1(3), v2(3), v3(3)
    double precision, allocatable :: ptmp(:,:), vtmp(:,:,:) ! For sorting
    integer, allocatable :: ntmp(:) ! For sorting

    ! Check that the origin is not included in RVEC.
    do ivec = 1, nvec
      if (sum(abs(rvec(1:3,ivec))) < 1.d-12) die_here("vector #"-ivec+"is zero! atom#"-atom_id)
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
          v2(1:3) = vert(1:3,ivert  ,iface)
          v3(1:3) = vert(1:3,ivert+1,iface)
          tetrvol = v1(1)*(v2(2)*v3(3) - v3(2)*v2(3)) + v2(1)*(v3(2)*v1(3) - v1(2)*v3(3)) + v3(1)*(v1(2)*v2(3) - v2(2)*v1(3))
          volume = volume + abs(tetrvol)
          trianglearea = 0.5d0 * sqrt( &
                ( v1(1)*v2(2) + v2(1)*v3(2) + v3(1)*v1(2) - v2(2)*v3(1) - v3(2)*v1(1) - v1(2)*v2(1) )**2 &
              + ( v1(2)*v2(3) + v2(2)*v3(3) + v3(2)*v1(3) - v2(3)*v3(2) - v3(3)*v1(2) - v1(3)*v2(2) )**2 &
              + ( v1(3)*v2(1) + v2(3)*v3(1) + v3(3)*v1(1) - v2(1)*v3(3) - v3(1)*v1(3) - v1(1)*v2(3) )**2 )
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
    rout2 = 0.d0 ! init with zero
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

    if (output) write(*, fmt="('Voronoi subroutine: RMT=',e16.8,'; ROUT=',e16.8,'; RATIO=',f12.2,' %, atom #',i0)") rmt,rout,rmt*100/rout,atom_id
    
    call dsort(rsort, isort, nface)
    ! Rearrange using a temporary arrays ptmp, vtmp, ntmp
    allocate(ptmp(0:3,nface), vtmp(3,nvertmax,nface), ntmp(nface))
    ptmp(0:3,:) = planes(0:3,1:nface)
    vtmp(:,:,:) = vert(:,:,  1:nface)
    ntmp    (:) = nvert     (1:nface)
    do iface = 1, nface
      i = isort(iface)
      planes(:,i) = ptmp(:,  iface) ! (0:3,  ...)
      vert(:,:,i) = vtmp(:,:,iface) ! (1:3,:,...)
      nvert   (i) = ntmp    (iface) ! (      ...)
    enddo ! iface    
    deallocate(ptmp, vtmp, ntmp, stat=i) ! ignore status

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
          vdist2 = dv(1)*dv(1) + dv(2)*dv(2) + dv(3)*dv(3)
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
      vdist2 = dv(1)*dv(1) + dv(2)*dv(2) + dv(3)*dv(3)
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

  logical function half_space(p, v)
    double precision, intent(in) :: p(0:3), v(1:3)
#ifndef NDEBUG
    if (sum(abs(p(1:3))) < 1.d-80) die_here('halfspace: a,b,c too small.')
#endif
    half_space = (p(0)*(p(1)*v(1) + p(2)*v(2) + p(3)*v(3)) <= p(0)*p(0))
  endfunction ! halfspace
  
  function normal_plane(v1, tau) result(p)
    double precision :: p(0:3) ! result, former [d,a,b,c]
    double precision, intent(in) :: v1(3), tau
!     ! given a point in space, r1=(v1(1),v1(2),v1(3)), this
!     ! subroutine returns the coefficients defining a plane through the 
!     ! equation a*x+b*y+c*z=d, which is normal to the vector r1 and passes
!     ! through the point tau*r1 (tau thus being a parameter
!     ! defining how close the plane is to the point).
    
!     ! so a,b,c are the coords. of the vector tau * r1.
!     ! if tau=0 (plane passes through the origin), then d=0.

    if (tau /= 0.d0) then
      p(1:3) = tau*v1(1:3)
      p(0)   = p(1)*p(1) + p(2)*p(2) + p(3)*p(3)
    else
      p(0)   = 0.d0
      p(1:3) = v1(1:3)
    endif

  endfunction ! normal_plane
  

  subroutine sortvertices(n, s, xyz)
    ! sorts the array s(n) in ascending order using straight insertion. 
    ! the arrays z(n), y(n), and z(n) follow.
    ! on output, arrays s, x, y, and z return sorted.
    integer, intent(in) :: n
    double precision, intent(inout) :: s(:), xyz(:,:)
    double precision :: sxyz(0:3) ! temp. swap
    integer :: i, j

    outer: do j = 2, n
      sxyz(0) = s(j) ; sxyz(1:3) = xyz(1:3,j)
      do i = j-1, 1, -1
        if (s(i) <= sxyz(0)) then
          s(i+1) = sxyz(0) ; xyz(1:3,i+1) = sxyz(1:3)
          cycle outer
        endif
        s(i+1) = s(i) ; xyz(1:3,i+1) = xyz(1:3,i)
      enddo ! i
      s(1) = sxyz(0) ; xyz(1:3,1) = sxyz(1:3)
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

 
  double precision function dist_plane(p)
    ! returns the distance of a plane a*x+b*y+c*z=d to the origin.
    double precision, intent(in) :: p(0:3)
    double precision :: d2

    d2 = p(1)*p(1) + p(2)*p(2) + p(3)*p(3)
    assert(d2 >= 1.d-100)
    dist_plane = abs(p(0))/sqrt(d2)

  endfunction ! dist_plane


endmodule ! Voronoi_mod
