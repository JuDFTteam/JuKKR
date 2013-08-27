!***********************************************************************
subroutine voronoi08( &
     nvec,rvec,nvertmax,nfaced,weight0,weight,tolvdist,tolarea, &
     rmt,rout,volume,nface,a3,b3,c3,d3,nvert,xvert,yvert,zvert, &
     output)
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
  ! XVERT(NVERTMAX,NFACE),YVERT(NVERTMAX,NFACE),ZVERT(NVERTMAX,NFACE): 
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

  use voronoi08_mod ! import helper routines
  implicit none

  ! Input:
  integer  nvec,nvertmax,nfaced
  real*8   rvec(3,nfaced)
  real*8   weight0,weight(nfaced)  ! Weight of central atom, 
  !                          !   and of all others (dimensioned as RVEC). 
  real*8  tolvdist              ! Max. tolerance for distance of two vertices
  real*8  tolarea               ! Max. tolerance for area of polygon face
  logical output  ! test output

  ! Output:
  integer nface
  integer nvert(nfaced)
  real*8              volume
  real*8              a3(nfaced),b3(nfaced),c3(nfaced),d3(nfaced)
  real*8              xvert(nvertmax,nfaced),yvert(nvertmax,nfaced), &
       zvert(nvertmax,nfaced)
  real*8              rmt,rout
  ! Inside:
  integer ivec,iface,ivert,i,nverttot
  integer nplane
  real*8           x1,y1,z1,x2,y2,z2,x3,y3,z3,tetrvol,rsq,tau
  real*8           facearea(nfaced),trianglearea   
  real*8           temp

  integer isort(nfaced),pos,inew    ! Index for sorting
  real*8 rsort(nfaced)              ! Aux. function for sorting
  real*8 atmp(nfaced),btmp(nfaced),ctmp(nfaced),dtmp(nfaced) ! For sorting
  real*8 xvtmp(nvertmax,nfaced),yvtmp(nvertmax,nfaced)       ! For sorting
  real*8 zvtmp(nvertmax,nfaced)                              ! For sorting
  integer nvtmp(nfaced)                                      ! For sorting
  integer npanel

  !---------------------------------------------------------------
  ! Check that the origin is not included in RVEC.
  do ivec = 1,nvec
     if (dabs(rvec(1,ivec))+dabs(rvec(2,ivec))+dabs(rvec(3,ivec)) &
          .lt.1.d-12) then
        write(*,*) 'VORONOI: Vector',ivec,'is zero.'
        stop 'VORONOI'
     endif
  enddo
  !---------------------------------------------------------------
  ! Define the planes as normal to the vectors RVEC, passing from t as
  ! in eq. (2) above:
  do ivec = 1,nvec
     rsq = rvec(1,ivec) * rvec(1,ivec) &
          + rvec(2,ivec) * rvec(2,ivec) &
          + rvec(3,ivec) * rvec(3,ivec)
     tau = 0.5d0*((weight0-weight(ivec))/rsq + 1.d0)
     call normalplane0( &
          rvec(1,ivec),rvec(2,ivec),rvec(3,ivec),tau, &
          a3(ivec),b3(ivec),c3(ivec),d3(ivec))
  enddo
  !---------------------------------------------------------------
  ! Find the Voronoi polyhedron.
  nplane = nvec
  if (output) write(*,*) 'Entering POLYHEDRON08'
  call polyhedron08( &
       nplane,nvertmax,nfaced,tolvdist,tolarea, &
       a3,b3,c3,d3, &
       nface,nvert,xvert,yvert,zvert,output)
  if (output) write(*,*) 'Exited POLYHEDRON08'

  !---------------------------------------------------------------
  ! Calculate the volume as sum of the volumes of all tetrahedra 
  ! connecting the origin to the faces. Use for each tetrahedron
  ! volume = det((r0-r1),(r0-r2),(r0-r3))/6, where r0 is here the
  ! origin and r1,r2,r3 the vectors of the 3 other vertices.
  ! Algorithm requires that the face is a convex polygon with the 
  ! vertices ordered (doesn't matter if they are clock- or anticlockwise).
  volume = 0.d0
  do iface = 1,nface
     x1 = xvert(1,iface)
     y1 = yvert(1,iface)
     z1 = zvert(1,iface)
     facearea(iface) = 0.d0
     do ivert = 2,nvert(iface)-1
        x2 = xvert(ivert,iface)
        y2 = yvert(ivert,iface)
        z2 = zvert(ivert,iface)
        x3 = xvert(ivert+1,iface)
        y3 = yvert(ivert+1,iface)
        z3 = zvert(ivert+1,iface)
        tetrvol = x1*(y2*z3-y3*z2)+x2*(y3*z1-y1*z3)+x3*(y1*z2-y2*z1)
        volume = volume + dabs(tetrvol)
        trianglearea = 0.5d0 * dsqrt( &
             ( x1*y2 + x2*y3 + x3*y1 - y2*x3 - y3*x1 - y1*x2)**2 &
             + ( y1*z2 + y2*z3 + y3*z1 - z2*y3 - z3*y1 - z1*y2)**2 &
             + ( z1*x2 + z2*x3 + z3*x1 - x2*z3 - x3*z1 - x1*z2)**2 )
        ! wrong           TRIANGLEAREA = 0.5d0 * DABS(
        !     &           (X2-X1)*(X3-X1)+(Y2-Y1)*(Y3-Y1)+(Z2-Z1)*(Z3-Z1) )
        facearea(iface) = facearea(iface)+ trianglearea
     enddo
  enddo
  volume = volume/6.d0

  write(6,*) ' Polyhedron properties '
  write(6,*) ' Number of faces : ',nface

  if (output) then
     do iface=1,nface
        write(6,201) iface,nvert(iface),facearea(iface)
201     format(' Face ',i4,'   has ',i4,' vertices ','; Area= ',e12.4)   ! 
        do ivert=1,nvert(iface)
           write(*,9000) ivert,xvert(ivert,iface),yvert(ivert,iface) &
                ,zvert(ivert,iface)
        enddo
        write(*,9010) a3(iface),b3(iface),c3(iface),d3(iface)
     end do
  endif

  write(6,*) 'The Volume is : ',volume

9000 format(i5,4e16.8)
9010 format(' Face coefficients:',4e16.8)
  !---------------------------------------------------------------
  ! Find RMT:
  rmt = distplane(a3(1),b3(1),c3(1),d3(1))
  do iface = 2,nface
     temp = distplane(a3(iface),b3(iface),c3(iface),d3(iface))
     if (temp.lt.rmt) rmt = temp
  enddo
  ! Fint ROUT:
  rout = 0.d0
  do iface = 1,nface
     do ivert = 1,nvert(iface)
        temp = xvert(ivert,iface)*xvert(ivert,iface) &
             + yvert(ivert,iface)*yvert(ivert,iface) &
             + zvert(ivert,iface)*zvert(ivert,iface)
        if (temp.gt.rout) rout = temp
     enddo
  enddo
  rout = dsqrt(rout)
  write(*,9020) rmt,rout,rmt*100/rout
9020 format('Voronoi subroutine: RMT=',e16.8,'; ROUT=',e16.8, &
       '; RATIO=',f12.2,' %')
  !---------------------------------------------------------------
  ! Sort faces:
  do iface = 1,nface
     rsort(iface) = &
          distplane(a3(iface),b3(iface),c3(iface),d3(iface))
     rsort(iface) = 100000000.d0  * rsort(iface) + &
          100000.d0  * c3(iface)    + &
          100.d0  * b3(iface)    + &
          0.1d0 * a3(iface)
  enddo
  call dsort(rsort,isort,nface,pos)
  !     Rearrange using a temporary array
  atmp(:) = a3(:)
  btmp(:) = b3(:)
  ctmp(:) = c3(:)
  dtmp(:) = d3(:)
  xvtmp(:,:) = xvert(:,:)
  yvtmp(:,:) = yvert(:,:)
  zvtmp(:,:) = zvert(:,:)
  nvtmp(1:nface) = nvert(1:nface)
  do iface = 1,nface
     inew = isort(iface)
     a3(inew) = atmp(iface)
     b3(inew) = btmp(iface)
     c3(inew) = ctmp(iface)
     d3(inew) = dtmp(iface)
     xvert(:,inew) = xvtmp(:,iface)
     yvert(:,inew) = yvtmp(:,iface)
     zvert(:,inew) = zvtmp(:,iface)
     nvert(inew) = nvtmp(iface)
  enddo
  !---------------------------------------------------------------


  return
end subroutine voronoi08











