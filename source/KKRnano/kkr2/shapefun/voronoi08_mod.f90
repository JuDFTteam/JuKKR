module voronoi08_mod
  implicit none
  private
  public :: distplane, normalplane0, polyhedron08, dsort

  interface normal_plane
    module procedure normal_plane0f, normal_planef, normalplane, normalplane0
  endinterface
  
  contains

  subroutine polyhedron08(nplane, nvertmax, nfaced, tolvdist,tolarea, a3,b3,c3,d3, nface, nvert, xvert,yvert,zvert, output)
    ! given a set of planes, defined by a3*x+b3*y+c3*z=d3 and defining
    ! a convex part of space (the minimal one containing the origin, 
    ! usually a ws-polyhedron), this subroutine returns the actual faces of
    ! the polyhedron, discarding the planes that do not contain faces. also,
    ! the coordinates of the verticess of the faces xvert,yvert,zvert and their 
    ! number nvert per face are returned. the coefficients of the actual
    ! faces are returned in the same arrays a3,b3,c3, and d3.

    logical, intent(in) :: output

    ! input:
    integer, intent(in) :: nplane                ! initial number of planes.
    integer, intent(in) :: nvertmax              ! max. number of vertices per plane.
    integer, intent(in) :: nfaced                ! max. number of faces.
    double precision, intent(in) :: tolvdist              ! max. tolerance for distance of two vertices
    double precision, intent(in) :: tolarea               ! max. tolerance for area of polygon face
    double precision, intent(inout) :: a3(*),b3(*),c3(*),d3(*)  ! coefs. defining the planes, dimensioned >= nplane.
    double precision, intent(inout) :: xvert(nvertmax,nfaced), yvert(nvertmax,nfaced), zvert(nvertmax,nfaced) ! cartesian coords. of vertices for each plane (2nd index is for planes).
    integer, intent(out) :: nvert(*)  ! number of vertices found for each face
    integer, intent(out) :: nface     ! number of faces found (with nvert>0).
    
!   integer :: nverttot ! total number of vertices

    !---------------------------------------------------------------
    ! find all faces and vertices of the polyhedron. on output, the vertices of each face are sorted (clockwise or anticlockwise).
    if (output) write(*,*) 'entering vertex3d'
    call vertex3d(nplane, a3, b3, c3, d3, nvertmax, tolvdist, nface, nvert, xvert, yvert, zvert, output)
    if (output) write(*,*) 'vertex3d found',nface,' faces with >3 vertices.'
    !---------------------------------------------------------------
    ! analyze the faces and vertices of the polyhedron.
    ! use criteria for rejecting faces that are too small or vertices that are too close to each other.
    ! on output, number of faces and vertices may be reduced after some rejections have taken place.
    if (output) write(*,*) 'entering analyzevert3d'
    call analyzevert3d(nvertmax, nfaced, tolvdist, tolarea, nplane, nface, nvert, xvert,yvert,zvert, a3,b3,c3,d3, output)
    if (output) write(*,*) 'analyzevert3d accepted',nface,' faces.'

  endsubroutine polyhedron08


!***********************************************************************
subroutine vertex3d(nplane, a3,b3,c3,d3, nvertmax, tolvdist, nface, nvert, xvert,yvert,zvert, output)
! given a set of planes, defined by a3*x+b3*y+c3*z=d3 and defining
! a convex part of space (the minimal one containing the origin, 
! usually a ws-polyhedron), this subroutine returns the vertices
! of this polyhedron in cartesian coordinates. for the planes that
! are not faces of the polyhedron, a value nvert(iplane)=0 is returned.
! the total no. of faces found is returned as nface.

integer, intent(in) :: nplane ! number of planes.
integer, intent(in) :: nvertmax              ! max. number of vertices per plane.
double precision, intent(in) :: a3(*),b3(*),c3(*),d3(*)  ! coefs. defining the planes, dimensioned >= nplane.
double precision, intent(in) :: tolvdist               ! min. distance between vertices
integer, intent(out) :: nvert(*)  ! number of vertices found for each face
integer, intent(out) :: nface     ! number of faces found (with nvert>0).
double precision, intent(inout) :: xvert(nvertmax,*), yvert(nvertmax,*), zvert(nvertmax,*) ! cartesian coords. of vertices for each plane ! (2nd index is for planes).
logical, intent(in) :: output

integer ip1,ip2,ip3,ipl,kpl ! plane indices
integer ivert                    ! vertex index
double precision          xcut,ycut,zcut ! cut point of three planes.
double precision          det,detx,dety,detz ! determinants of 3x3 system for xcut,ycut,zcut.
double precision          distance   ! a distance criterium of two points in space
! the following are for sorting the vertices of each face:
double precision          v1(3),v2(3),v3(3)   ! auxiliary vectors...
double precision          sinfiv1v2,cosfiv1v2 ! ...and their inner and outer products
double precision          fi(nvertmax)        ! ...and also their relative angles.
double precision          uv(3),vl,sn         ! unit vector, length, sign of sin(fi)
logical laccept ! determining whether a cut point is inside the polyhedron.
!---------------------------------------------------------------
! check & initialize
if (nplane < 4) then
  write(*,*) 'vert3d: error:nplane was only',nplane
  stop
endif

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
det  = a3(ip1)*(b3(ip2)*c3(ip3) - b3(ip3)*c3(ip2)) + a3(ip2)*(b3(ip3)*c3(ip1) - b3(ip1)*c3(ip3)) + a3(ip3)*(b3(ip1)*c3(ip2) - b3(ip2)*c3(ip1))

!---------------------------------------------------------------
if (dabs(det) > 1.d-12) then ! there is a cut point 

detx = d3(ip1)*(b3(ip2)*c3(ip3) - b3(ip3)*c3(ip2)) + d3(ip2)*(b3(ip3)*c3(ip1) - b3(ip1)*c3(ip3)) + d3(ip3)*(b3(ip1)*c3(ip2) - b3(ip2)*c3(ip1))
dety = a3(ip1)*(d3(ip2)*c3(ip3) - d3(ip3)*c3(ip2)) + a3(ip2)*(d3(ip3)*c3(ip1) - d3(ip1)*c3(ip3)) + a3(ip3)*(d3(ip1)*c3(ip2) - d3(ip2)*c3(ip1))
detz = a3(ip1)*(b3(ip2)*d3(ip3) - b3(ip3)*d3(ip2)) + a3(ip2)*(b3(ip3)*d3(ip1) - b3(ip1)*d3(ip3)) + a3(ip3)*(b3(ip1)*d3(ip2) - b3(ip2)*d3(ip1))

xcut = detx/det
ycut = dety/det
zcut = detz/det

!     write(6,333) ip1,ip2,ip3,xcut,ycut,zcut
!333  format('cutting point of planes ',3i5,':',3d15.7)
!-----------------------------------
! accept this cut point as a vertex, if it belongs to the polyhedron. so,
! make a loop over all other (than ip1,2,3) planes:
laccept = .true.
do kpl = 1, nplane
   if (kpl == ip1 .or. kpl == ip2 .or. kpl == ip3) cycle
   laccept = laccept .and. halfspace(a3(kpl),b3(kpl),c3(kpl),d3(kpl),xcut,ycut,zcut)
enddo ! kpl
!-----------------------------------
if (laccept) then
! if the cut point found belongs to the cell, we accept it unless it has
! occured before for this face (ip1). such a situation is possible
! when 4 or more planes pass through the same point (e.g. for the vertices
! of the fcc ws-cell). so...
   do ivert = 1,nvert(ip1)
      distance = (xvert(ivert,ip1) - xcut)**2 + (yvert(ivert,ip1) - ycut)**2 + (zvert(ivert,ip1) - zcut)**2
      distance = dsqrt(distance)
      if (distance < tolvdist ) then
         laccept = .false. ! vertex is too close to a previous one.
         exit              ! jump loop, no need to continue.
      endif
   enddo ! ivert
endif ! laccept

! now we are ready to add the point to the vertex list.
if (laccept) then
   nvert(ip1) = nvert(ip1) + 1
   xvert(nvert(ip1),ip1) = xcut
   yvert(nvert(ip1),ip1) = ycut
   zvert(nvert(ip1),ip1) = zcut
endif

endif ! (dabs(det) > 1.d-12)
!---------------------------------------------------------------

enddo plane3 ! ip3
enddo plane2 ! ip2
!     write(6,*) 'number of vertices for plane ',ip1,'  :',nvert(ip1)
enddo plane1 ! ip1

!===============================================================
! each plane should finally have either at least 3 vertices, if it is a
! face of the polyhedron, or none at all. check this:
if (output) then
  do ipl = 1, nplane
    if (nvert(ipl) == 1 .or. nvert(ipl) == 2) then
      write(*,*) 'vertex3d: error:there is a problem with the vertices.'
      write(*,*) 'for plane',ipl,' ,only ',nvert(ipl),' vertices were found.'
    endif
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
    vl = dsqrt(xvert(1,ipl)**2 + yvert(1,ipl)**2 + zvert(1,ipl)**2)
    uv(1) = xvert(1,ipl)/vl
    uv(2) = yvert(1,ipl)/vl
    uv(3) = zvert(1,ipl)/vl

  ! define the vector connecting the first vertex to the (now-) second:
    v1(1) = xvert(2,ipl) - xvert(1,ipl)
    v1(2) = yvert(2,ipl) - yvert(1,ipl)
    v1(3) = zvert(2,ipl) - zvert(1,ipl)

    do ivert = 2,nvert(ipl)
  ! define the vector connecting the first vertex to the current one:
        v2(1) = xvert(ivert,ipl) - xvert(1,ipl)
        v2(2) = yvert(ivert,ipl) - yvert(1,ipl)
        v2(3) = zvert(ivert,ipl) - zvert(1,ipl)
  ! find the angle fi between v1 and v2 
  ! ( always, -pi < fi < pi from the definition of datan2 )
        cosfiv1v2 = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
        v3 = crospr(v1, v2) ! cross product = |v1|*|v2|*sinfi
        sinfiv1v2 = dsqrt(v3(1)*v3(1) + v3(2)*v3(2) + v3(3)*v3(3))
  ! sign of sinfi is defined with respect to unit vector uv (see above)
        sn = uv(1)*v3(1) + uv(2)*v3(2) + uv(3)*v3(3)
        if (sn < 0) sinfiv1v2 = -sinfiv1v2

        if (sinfiv1v2 == 0.d0 .and. cosfiv1v2 == 0.d0) then
  ! point falls exactly on 1st vertex...
          fi(ivert) = -4.d0
          write(*,*) 'vertex3d: error: found two identical vertex points'
  ! ...while it shouldn't ! (this was checked earlier)
          stop
        else
          fi(ivert) = datan2(sinfiv1v2, cosfiv1v2)
        endif

    enddo ! ivert = 3,nvert(ipl)

  ! store with respect to the angle found:
    call sortvertices(nvert(ipl), fi, xvert(1,ipl), yvert(1,ipl), zvert(1,ipl))

  endif ! (nvert(ipl) >= 3)
enddo ! ipl = 1, nplane

endsubroutine vertex3d

!------------------------------------------------------------------------------
subroutine analyzevert3d(nvertmax, nfaced, tolvdist, tolarea, nplane, nface, nvert, xvert,yvert,zvert, a3,b3,c3,d3, output)
! analyze the faces and vertices of the polyhedron.
! use criteria for rejecting faces that are too small or vertices that are too close to each other.
! on output, number of faces and vertices may be reduced after some rejections have taken place.
logical, intent(in) :: output
integer, intent(in) :: nvertmax, nfaced
integer, intent(in) :: nplane ! number of planes
integer, intent(inout) :: nface ! number of faces (usually much less than nplane)
integer, intent(inout) :: nvert(nfaced) ! number of vertices for each face
double precision, intent(inout) :: xvert(nvertmax,nfaced), yvert(nvertmax,nfaced), zvert(nvertmax,nfaced)
double precision, intent(in) :: tolvdist, tolarea  ! max. tolerance for distance of two vertices and area of face.
double precision, intent(inout) :: a3(*),b3(*),c3(*),d3(*) ! coefs. defining the planes, to be reordered at end


double precision vdist             ! distance between consecutive vertices
logical lacceptvert(nvertmax),lacctot,lfoundnext,lthisisthelast
logical lacceptface(nfaced)
integer newindexface(nfaced)
integer iface,ivert,ivert2,inext,iplane
integer ifacenewcount,ivertnewcount
double precision dx,dy,dz,x1,x2,x3,y1,y2,y3,z1,z2,z3,trianglearea
double precision facearea(nfaced)
double precision x4,y4,z4,det,detsum

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
         dx = xvert(ivert2,iplane) -  xvert(ivert,iplane)
         dy = yvert(ivert2,iplane) -  yvert(ivert,iplane)
         dz = zvert(ivert2,iplane) -  zvert(ivert,iplane)
         vdist = dsqrt(dx*dx + dy*dy + dz*dz)

         if (vdist >= tolvdist) then
            lacceptvert(ivert2) = .true.   ! vertex is to be accepted
            inext = ivert2                 ! set this as the next vertex
            lfoundnext = .true.            ! and we have a winner, exit loop.
         else
            lacceptvert(ivert2) = .false.  ! remember that vertex is to be rejected later
            lacctot = .false.              ! remember that at least one vertex has to be rejected.
            ivert2 = ivert2 + 1            ! now compare to the next vertex
        endif  
      enddo ! while

      if (.not.lfoundnext) lthisisthelast = .true. ! if there is no next acceptable vertex,
                                                        ! then this was the last one. jump out.
      ivert = inext
      
   enddo ! while

! ...and now 1st with last to close the cycle:
   ivert = 1
   ivert2 = nvert(iplane)
   dx = xvert(ivert2,iplane) - xvert(ivert,iplane)
   dy = yvert(ivert2,iplane) - yvert(ivert,iplane)
   dz = zvert(ivert2,iplane) - zvert(ivert,iplane)
   vdist = dsqrt(dx*dx + dy*dy + dz*dz)
   if (vdist >= tolvdist) then
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
               xvert(ivertnewcount,iplane) = xvert(ivert,iplane) ! re-index vertex
               yvert(ivertnewcount,iplane) = yvert(ivert,iplane)
               zvert(ivertnewcount,iplane) = zvert(ivert,iplane)
            endif
         endif ! lacceptvert(ivert)
      enddo ! ivert
      nvert(iplane) = ivertnewcount
   endif ! not lacctot

enddo ! iplane



! now analyze faces, reject faces with less than three vertices and faces of very small area.
do iplane = 1, nplane
   if (nvert(iplane) >= 3) then  ! calculate area
      x1 = xvert(1,iplane)
      y1 = yvert(1,iplane)
      z1 = zvert(1,iplane)
      facearea(iplane) = 0.d0
      do ivert = 2, nvert(iplane)-1
         x2 = xvert(ivert,iplane)
         y2 = yvert(ivert,iplane)
         z2 = zvert(ivert,iplane)
         x3 = xvert(ivert+1,iplane)
         y3 = yvert(ivert+1,iplane)
         z3 = zvert(ivert+1,iplane)
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
            xvert(ivert,ifacenewcount) = xvert(ivert,iplane) ! re-index face vertices
            yvert(ivert,ifacenewcount) = yvert(ivert,iplane) ! re-index face vertices
            zvert(ivert,ifacenewcount) = zvert(ivert,iplane) ! re-index face vertices
         enddo ! ivert
         a3(ifacenewcount) = a3(iplane) ! re-index face equation parameters
         b3(ifacenewcount) = b3(iplane) ! re-index face equation parameters
         c3(ifacenewcount) = c3(iplane) ! re-index face equation parameters
         d3(ifacenewcount) = d3(iplane) ! re-index face equation parameters
      endif ! ifacenewcount /= iplane
   endif ! lacceptface(iplane)
enddo ! iplane
nface = ifacenewcount

! check for every face that all veritces lie on the same plane
! by checking linear dependence
do iface = 1, nface
   x2 = xvert(2,iface) - xvert(1,iface)
   y2 = xvert(2,iface) - yvert(1,iface)
   z2 = xvert(2,iface) - zvert(1,iface)
   x3 = xvert(3,iface) - xvert(1,iface)
   y3 = xvert(3,iface) - yvert(1,iface)
   z3 = xvert(3,iface) - zvert(1,iface)
   detsum = 0.d0
   do ivert = 4, nvert(iface)
      x4 = xvert(ivert,iface) - xvert(1,iface)
      y4 = xvert(ivert,iface) - yvert(1,iface)
      z4 = xvert(ivert,iface) - zvert(1,iface)
      det = x2*(y3*z4 - y4*z3) + y2*(z3*x4 - z4*x3) + z2*(x3*y4 - x4*y3)
      detsum = detsum + dabs(det)
      if (dabs(det) > 1.d-16 .and. output) write(*,fmt="('error from analyzevert3d: vertices not on single plane. iface=',i5,' ivert=',i5,' determinant=',e12.4)") iface,ivert,det
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

    if (dabs(a)+dabs(b)+dabs(c) < 1.d-80) stop 'halfspace: a,b,c too small.'

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

  endsubroutine ! normalplane
  
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

  endsubroutine ! normalplane
  
  function normal_plane0f(v1, tau) result(plane)
    double precision :: plane(0:3)
    double precision, intent(in) :: v1(3), tau
    
    call normalplane0(v1(1), v1(2), v1(3), tau, plane(1), plane(2), plane(3), plane(0))

  endfunction ! normalplane

  function normal_planef(v1, v2, tau) result(plane)
    double precision :: plane(0:3)
    double precision, intent(in) :: v1(3), v2(3), tau
    
    call normalplane(v1(1), v1(2), v1(3), v2(1), v2(2), v2(3), tau, plane(1), plane(2), plane(3), plane(0))

  endfunction ! normalplane
  

!***********************************************************************
  subroutine sortvertices(n,s,x,y,z)
  ! sorts the array s(n) in ascending order using straight insertion. 
  ! the arrays z(n), y(n), and z(n) follow.
  ! on output, arrays s, x, y, and z return sorted.
    integer, intent(in) :: n
    double precision, intent(inout) :: s(*), x(*),y(*),z(*)
    double precision :: tmp(0:3)
    integer :: i, j

    outer: do j = 2, n
      tmp = [s(j), x(j), y(j), z(j)]
      do i = j-1, 1, -1
          if (s(i) <= tmp(0)) then
            s(i+1) = tmp(0)
            x(i+1) = tmp(1)
            y(i+1) = tmp(2)
            z(i+1) = tmp(3)
            cycle outer
          endif
          s(i+1) = s(i)
          x(i+1) = x(i)
          y(i+1) = y(i)
          z(i+1) = z(i)
      enddo ! i
      s(1) = tmp(0)
      x(1) = tmp(1)
      y(1) = tmp(2)
      z(1) = tmp(3)
    enddo outer ! j

  endsubroutine ! sort vertices

! ************************************************************************
  function crospr(x, y) result(z)
! ************************************************************************
!     crosp computes the cross product of x and y returning it into z.
! ------------------------------------------------------------------------
  double precision :: z(3) ! result
  double precision, intent(in) :: x(3), y(3)
  z(1) = x(2)*y(3) - x(3)*y(2)
  z(2) = x(3)*y(1) - x(1)*y(3)
  z(3) = x(1)*y(2) - x(2)*y(1)
  endfunction crospr

!***********************************************************************
double precision function distplane(a,b,c,d)
! returns the distance of a plane a*x+b*y+c*z=d to the origin.
  double precision, intent(in) :: a,b,c,d
  double precision :: abcsq

  abcsq = a*a + b*b + c*c

  if (abcsq < 1.d-100) stop 'distplane'

  distplane = dabs(d)/dsqrt(abcsq)  
  ! or: distplane = dsqrt(d*d/abcsq)  

endfunction distplane

! ************************************************************************
subroutine dsort(w, ind, nmax, pos)
! ************************************************************************
!     p.zahn, april 96
!     w   is the original array returned unchanged
!     ind is an array that holds the new positions
!     nmax number of ellements to be sorted
!     pos the position where the first element is found
! ------------------------------------------------------------------------
  integer, intent(in) :: nmax
  double precision, intent(in) :: w(*)
  integer, intent(out) :: ind(*)
  integer, intent(out), optional :: pos
  
  double precision, parameter :: bound = 1.0d-12
  double precision :: diff
  integer :: i, ii, j, jj, k
  ! ------------------------------------------------------------------------
  do i = 1, nmax
    ind(1:nmax) = i
  enddo ! i

  j = 1
  do while (j < nmax/3)
    j = 3*j+1
  enddo ! while

  do while (j > 1)
    j = j/3
    jj = 1
    do while (jj == 1)
      jj = 0
      do k = 1, nmax-j
        diff = abs(w(ind(k)) - w(ind(k+j)))
        if (w(ind(k)) > w(ind(k+j)) .and. diff > bound) then
          ii       = ind(k)
          ind(k)   = ind(k+j)
          ind(k+j) = ii
          jj = 1
        endif
      enddo ! k=1,nmax-j
   enddo ! while (jj == 1)
  enddo ! while (j > 1)

  if (present(pos)) then
    do i = 1, nmax
      if (ind(i) == 1) pos = i
    enddo ! i
  endif

endsubroutine ! dsort

endmodule ! voronoi08_mod
