!> Module to determine radial mesh panel positions.

module ShapeCriticalPoints_mod
  implicit none
  private
  public :: criticalShapePoints
  
  interface rotate ! change management
    module procedure rotate_common, rotate_arg
  endinterface
  
  interface euler ! change management
    module procedure euler_common, euler_arg
  endinterface
  
  contains

!------------------------------------------------------------------------------
!> Finds critical points for mesh generation and performs geometrical tests.

! *) khcp ? not needed

! TODO: add output parameters to argument list

! TODO: find out whats needed in subsequent part

!    &                 TOLVDIST,   ! Max. tolerance for distance of two vertices
!    &                 TOLEULER,   ! Used in calculation of Euler angles, subr. EULER
!    &                 NMIN,       ! Min. number of points in panel
!    &                 NVERTICES8,XVERT8,YVERT8,ZVERT8,NFACE8,LMAX8,
!    &                 DLT8,KEYPAN8,NM8

!> @param[in]  AFACE coefficients of plane equations of cell faces
!> @param[in]  BFACE
!> @param[in]  CFACE
!> @param[in]  DFACE
!> @param[in]  TOLVDIST numerical tolerance for geometrical constructions
!> @param[in]  TOLEULER numerical tolerance for plane rotations
!> @param[in]  NVERTICES number of vertices of cell
!> @param[in]  XVERT x-coordinates of cell vertices
!> @param[in]  YVERT y-coordinates of cell vertices
!> @param[in]  ZVERT z-coordinates of cell vertices
!> @param[in]  NFACE number of cell faces
!> @param[in]  LMAX calculate Shape-functions up to this cutoff
!> @param[out] NPAN number of panels, depending on critical points on radial mesh
!> @param[out] CRT critical points (where panels have to be placed)
!> Note: the smallest critical point corresponds to the muffin-tin radius

!=====================================================================
subroutine criticalShapePoints(aface,bface,cface,dface, &
  tolvdist, toleuler, &
  nvertices,xvert,yvert,zvert,nface,lmax, &
  npan, crt, npand)

  use shape_constants_mod, only: verbosity, check_geometry, isumd, lmaxd1, dp, pi
  use tetrahedra_common, only: ntt, fa, fb, fd, rd, isignu
  use shapegeometryhelpers_mod, only: polchk

  integer, intent(in) :: npand
  integer, intent(in) :: nvertices(:) ! (nfaced)
  double precision, intent(in) :: aface(:), bface(:), cface(:), dface(:) ! (nfaced) --> collect to (0:3,nfaced)
  double precision, intent(in) :: xvert(:,:), yvert(:,:), zvert(:,:) ! (nvertd,nfaced) --> collect to (3,nvertd,nfaced)
  
  integer, intent(out) :: npan
  double precision, intent(out) :: crt(npand)
  
  !-----------------------------------------------------------------------
  
  integer ::     iface,ipan,isum,is0
  integer ::     iv,ivert,ivtot,l,lmax
  integer ::     nface,nvert,nvtot,ibmax
  double precision ::      sq3o3,coa
  double precision ::      a1,a2,a3,a4
  double precision ::      tolvdist,toleuler
  double precision, allocatable :: v(:,:)
  double precision :: z(3)
  integer :: nvertd

! #define k_HCP  
  
  nvertd = size(xvert,1)
  allocate(v(3,nvertd))

  !-----------------------------------------
  !  this call does some geometrical tests (n.stefanou 98)
  if (check_geometry) call polchk(nface,nvertices,xvert,yvert,zvert,tolvdist)
  
#ifdef k_HCP 
    coa = 2.d0
    sq3o3 = sqrt(3.d0)/3.d0
#endif

  ibmax = (lmax+1)**2
  isum = 0
  do l = 0, lmax
    is0 = (2*l+1)**2
    isum = isum + is0
  enddo ! l

  if (isum > isumd .or. lmax > lmaxd1)  then
    ! check needed for routine dreal
    write(6,fmt="(23x,'from main : isum=',i7,'  greater than dimensioned',i7/23x,'       or   lmax=',i7,'  greater than dimensioned',i7)") isum,isumd,lmax,lmaxd1
    stop
  endif

  ipan=0
  ivtot=0
  !.......................................................................
  !     s t o r a g e            i n    c o m m o n        b l o c k s
  !     c a l c u l a t i o n    o f    r o t a t i o n    m a t r i c e s
  !.......................................................................
  do iface = 1, nface

    z(1:3) = [aface(iface), bface(iface), cface(iface)]/dface(iface)

#ifdef k_HCP 
    z(1) = z(1)*sq3o3
    z(3) = z(3)*8.d0/(coa*3.d0)
#endif      

    do ivert = 1, nvertices(iface)
      v(1:3,ivert) = [xvert(ivert,iface), yvert(ivert,iface), zvert(ivert,iface)]

#ifdef k_HCP 
      v(1,ivert) = v(1,ivert)*sq3o3
      v(3,ivert) = v(3,ivert)*coa
#endif
    enddo ! ivert

    call crit(iface,nvertices(iface),v,z,ipan,ivtot,toleuler,tolvdist,crt,npand)

    if (verbosity > 0) then
      write(6,fmt="(/10x,i3,'-th pyramid subdivided in ',i3,' tetrahedra')") iface,ntt(iface)
    endif

  enddo ! iface ! end of loop over faces

  !.......................................................................
  !     d e f i n i t i o n    o f    t h e    s u i t a b l e    m e s h
  !.......................................................................
  nvtot = ivtot
  npan = ipan

  if (verbosity > 1) then
    write(6,fmt="(//15x,'fa/pi',5x,'fb/pi',5x,'fd/pi',6x,'rd',8x,'isignu'/)")
    do iv = 1, nvtot
      write(6,fmt="(i10,4f10.4,i10)") iv,fa(iv)/pi,fb(iv)/pi,fd(iv)/pi,rd(iv),isignu(iv)
    enddo ! iv
  endif
  
end subroutine


!------------------------------------------------------------------------------
!>    Calculates critical points where panels have to be placed.
!>    @brief
!>    Call this routine for each face. Give face index as argument IFACE.
!>
!>    It also performs the subdivision in tetrahedra, which
!>    are then stored in module tetrahedras_common
!>    (needed for shape-function integration)
!>    (should be separate routine)
!>
!>    depends on: EULER, PERP, ROTATE (only needed for this routine!!!)
!>
!>    @param[in]     IFACE face index
!>    @param[in]     NVERT number of vertices
!>    @param[in]     V     coordinates of vertices V(3, NVERTD)
!>    @param[in,out] Z     coefficients determining plane of face Z(3) - changed on output!
!>                         Z(1) * X + Z(2) * Y + Z(3) * Z = 1 (ONE!)
!>    @param[in,out] IPAN  panel counter, pass 0 for 1st call
!>    @param[in,out] IVTOT tetrahedron index: has to be 0 for 1st call
!>    @param[in]     TOLEULER tolerance for Euler angles
!>    @param[in]     TOLVDIST tolerance for distances
!>    @param[in,out] CRT   array of critical points, CRT(NPAND)
!>    @param[in]     NPAND maximal number of panels allowed
subroutine CRIT(iface,nvert,v,z,ipan,ivtot,toleuler,tolvdist,crt, npand)
  !-----------------------------------------------------------------------
  !     this routine calculates the critical points 'crt' of the shape
  !     functions due to the face: z(1)*x + z(2)*y + z(3)*z = 1
  !     the face is rotated through the appropriate euler angles to be
  !     perpendicular to the z-axis. a further subdivision of the cen-
  !     tral pyramid into elementary tetrahedra is performed. the  ne-
  !     cessary quantities for the calculation are stored in common.
  !-----------------------------------------------------------------------
  use shape_constants_mod, only: pi, verbosity, dp
  use tetrahedra_common, only: ntt, r0, rd, isignu, fd, fa, fb
  use angles_common, only: alpha, beta, gamma
  use shapegeometryhelpers_mod, only: perp, nrm2, operator(.dot.)

  integer, intent(in) :: npand, iface, nvert
  integer, intent(inout) :: ipan, ivtot
  double precision, intent(in) :: toleuler, tolvdist
  double precision, intent(in) :: v(:,:) ! (3,nvertd)
  double precision, intent(inout) :: z(3)
  double precision, intent(inout) :: crt(npand)

  integer :: i, ix, icorn, ivert, inew, ip, ivertp, ivert1, iback
  double precision :: arg, a1, a2, a3, cf1, cf2, cf3, co, crrt, dd, down, d1, d2, ff, f1, f2, omega, rdd, s, sf1, sf2, sf3, up, xj, yj, zmod2, zvmod

  integer, allocatable :: in(:) ! (nvertd)
  double precision, allocatable :: vz(:,:) ! (3,nvertd)
  double precision :: rdv(3)
  double precision, parameter :: origin(3) = 0.d0
  logical :: inside
  double precision, parameter :: tol_small = 1d-6
  double precision, parameter :: tol_large = 1d-4

  !-----------------------------------------------------------------------

  allocate(in(size(v,2)), vz(3, size(v,2)))

  ntt(iface) = 0

  if (verbosity > 0) write(6,203) iface,z(1:3)

  zmod2 = nrm2(z)

  if (zmod2 <= tol_small) goto 100  ! check if normal vector of plane is valid

  s = 2d0*pi

  z = z/zmod2

  ix = 1
  zvmod = sqrt(nrm2(v(1:3,1)-z))

  if (zvmod < tol_small) ix = 2

! call euler(z,v(:,ix),iface,toleuler) ! store euler angles in common block
  call euler(z,v(:,ix),toleuler, alpha(iface), beta(iface), gamma(iface)) ! get euler angles directly

  if (verbosity > 0) write(6,204) alpha(iface)/pi,beta(iface)/pi,gamma(iface)/pi

!   call rotate(v,vz,iface,nvert) ! pass the index iface and access the angles in the common block
  call rotate(v, vz, nvert, alpha(iface), beta(iface), gamma(iface)) ! pass the angles directly

  r0(iface) = 1.d0/sqrt(zmod2)
  icorn = 0 ! 0: is not a corner

  if (verbosity > 0) write(6,207)

  do ivert = 1, nvert
    if (abs(r0(iface)-vz(3,ivert)) > tol_small) goto 101  !check if vertices lie in same plane
    !.......................................................................
    !     d i s t a n c e s   o f   v e r t i c e s   f r o m   c e n t e r
    !.......................................................................
    crrt = sqrt(nrm2(vz(1:3,ivert)))

    !     check if a new panel has been found (at radial point crrt)
    inew = 1
    do ip = 1, ipan
      if (abs(crrt - crt(ip)) < tol_small) inew = 0
    enddo ! ip

    if (inew == 1) then
      ipan = ipan+1
      if (ipan > npand) goto 102

      crt(ipan) = crrt
    endif ! new
    ivertp=ivert+1

    if (ivert == nvert) ivertp = 1

    !.......................................................................
    !     d i s t a n c e s   o f   e d g e s   f r o m   c e n t e r
    !.......................................................................

    call perp(origin,vz(1:3,ivert),vz(1:3,ivertp),rdv,tolvdist,inside)

    rdd = sqrt(nrm2(rdv)) ! footpoint of line origin-to-edge

    if (inside) then  ! add a new panel only when footpoint is contained on edge
      inew = 1
      do ip = 1, ipan
        if (abs(rdd - crt(ip)) < tol_large) inew = 0
      enddo ! ip
      if (inew == 1) then
        ipan = ipan+1
        if (ipan > npand) goto 102
        crt(ipan) = rdd
      endif ! new
    endif ! inside
    
    a1 = sqrt(vz(1,ivert)**2 + vz(2,ivert )**2)
    a2 = sqrt(vz(1,ivertp)**2 + vz(2,ivertp)**2)
    up = vz(1,ivert)*vz(1,ivertp) + vz(2,ivert)*vz(2,ivertp)
    down = a1*a2

    if (down > tol_small) then ! true if not a corner
      arg = up/down
      if (abs(arg) >= 1d0) arg = sign(1d0, arg)
      omega = acos(arg)
      s = s - omega

      if (abs(omega - pi) > tol_small) then
        !.......................................................................
        !     s u b d i v i s i o n    i n t o    t e t r a h e d r a
        !.......................................................................
        ntt(iface) = ntt(iface)+1
        ivtot = ivtot+1

        if (verbosity > 0) then
          write(6,205) ivtot,ivert, vz(1:3,ivert )
          write(6,206)       ivertp,vz(1:3,ivertp)
        endif

        a3 = sqrt(rdv(1)*rdv(1) + rdv(2)*rdv(2))
        rd(ivtot) = rdd
        isignu(ivtot) = 1
        
        cf1 = vz(1,ivert )/a1
        cf2 = vz(1,ivertp)/a2
        sf1 = vz(2,ivert )/a1
        sf2 = vz(2,ivertp)/a2
        cf3 = rdv(1)/a3
        sf3 = rdv(2)/a3

        if (abs(sf1) < toleuler .and. abs(cf1+1d0) < toleuler) then
          f1 = pi
        else
          f1 = 2d0*atan2(sf1, cf1+1d0)
        endif

        if (abs(sf2) < toleuler .and. abs(cf2+1d0) < toleuler) then
          f2 = pi
        else
          f2 = 2d0*atan2(sf2, cf2+1d0)
        endif

        if (abs(sf3) < toleuler .and. abs(cf3+1d0) < toleuler) then
          fd(ivtot) = pi
        else
          fd(ivtot) = 2d0*atan2(sf3, cf3+1d0)
        endif

        fa(ivtot) = min(f1, f2)
        fb(ivtot) = max(f1, f2)
        if ((fb(ivtot) - fa(ivtot)) > pi) then
          ff = fa(ivtot) + 2d0*pi ! increase new fb
          fa(ivtot) = fb(ivtot) ! swap
          fb(ivtot) = ff
        endif
        if ((fa(ivtot) - fd(ivtot)) > pi) fd(ivtot) =  2d0*pi + fd(ivtot)
        if ((fd(ivtot) - fa(ivtot)) > pi) fd(ivtot) = -2d0*pi + fd(ivtot)
      endif
    else  ! down > tol_small
      icorn = 1 ! is a corner
    endif ! down > tol_small
  enddo ! ivert ! end of vertex loop
  
  !.......................................................................
  !     f o o t   o f   t h e    p e r p end i c u l a r   to    t h e
  !     f a c e   o u t s i d e   o r  i n s i d e   t h e   p o l y g o n
  !.......................................................................
  if (s < tol_small .or. icorn == 1) then
    inew = 1
    do ip = 1, ipan
      if (abs(r0(iface)-crt(ip)) < tol_large) inew=0
    enddo ! ip
    
    if (inew == 1) then
      ipan = ipan+1
      if (ipan > npand) goto 102
      crt(ipan) = r0(iface)
    endif ! new
  else
    do ivert1 = 1, nvert
      in(ivert1) = 0
      do ivert = 1, nvert
        ivertp = ivert+1

        if (ivert == nvert) ivertp=1

        if (ivert == ivert1 .or. ivertp == ivert1) goto 7
        down = vz(2,ivert1)*(vz(1,ivertp) - vz(1,ivert)) - vz(1,ivert1)*(vz(2,ivertp)-vz(2,ivert))

        if (abs(down) <= tol_small) goto 7

        up  = vz(1,ivert1)*(vz(2,ivert)*(vz(1,ivertp) + vz(1,ivert)) - vz(1,ivert)*(vz(2,ivertp) + vz(2,ivert)))
        xj = up/down
        yj = xj*vz(2,ivert1)/vz(1,ivert1)
        dd = (vz(1,ivertp) - vz(1,ivert))**2 + (vz(2,ivertp) - vz(2,ivert))**2
        d1 = (xj - vz(1,ivert ))**2 + (yj - vz(2,ivert ))**2
        d2 = (xj - vz(1,ivertp))**2 + (yj - vz(2,ivertp))**2

        co = dd - max(d1, d2)

        if (co > tol_small) then
          in(ivert1) = 1
          goto 6
        endif

7     continue
      enddo
6   continue
    enddo
    iback = ivtot-nvert

    do ivert = 1, nvert
      iback = iback+1
      ivertp = ivert+1; if (ivert == nvert) ivertp = 1
      if (in(ivert) == 0 .and. in(ivertp) == 0) isignu(iback) = -1
    enddo ! ivert

  endif
  return
100 write(6,fmt="(//13x,'fatal error from crit: the',i3,'-th face of the polyhedron passes through the center'  /13x,'(',3e14.7,' )')") iface,z(1:3)
  stop
101 write(6,fmt="(//13x,'fatal error from crit: the vertices of the',i3,'-th rotated polygon do not lie on the plane:'  ,e13.6,' *z = 1'/30(/13x, 3e13.6))") iface,r0(iface),vz(1:3,1:nvert)
  stop
102 write(6,fmt="(//13x,'error from crit: number of panels=',i5,' greater than dimensioned='  ,i5)") ipan,npand
  stop
203 format(//80('*')/3x,'face:',i3,' equation:',f10.4,'*x +',f10.4,'*y +',f10.4,'*z  =  1')
204 format(3x,'rotation angles  :',3(f10.4,4x)/)
205 format(i5,'       vz(',i2,')  =  (',3f10.4,' )')
206 format(5x,'       vz(',i2,')  =  (',3f10.4,' )')
207 format(/'tetrahedron',14x,'coordinates'/11('*'),14x,11('*')/)
end subroutine crit

!-----------------------------------------------------------------------
!>    given two distinct points (z(1),z(2),z(3)) and (xx(1),xx(2),xx(3))
!>    this routine defines  a local coordinate  system with the  z- axis
!>    passing through (z(1),z(2),z(3))  and the x- axis parallel to  the
!>    vector : (xx(1)-z(1),xx(2)-z(2),xx(3)-z(3)).
!>    the euler angles rotating this local coordinate system back to the
!>    original frame of reference are calculated  and stored  in common.
!-----------------------------------------------------------------------
subroutine euler_common(z, xx, iface, toleuler)
  use shape_constants_mod, only: pi, dp
  use angles_common, only: alpha, beta, gamma
  integer, intent(in) :: iface
  double precision, intent(in) :: xx(3)
  double precision, intent(inout) :: z(3)
  double precision, intent(in) :: toleuler   ! introduced by phivos (05.2008) to account for inaccuracies.
  ! earlier, 1.d-5 was hard-coded at the places in this subr. where tol is used
  !     data toleuler /1.d-10/

  double precision :: rx,rz,s,p,rzp,sa,ca,sg,cg
  double precision :: x(3), y(3)
  double precision, parameter :: tolerance = 1d-6

  !-----------------------------------------------------------------------
  if (iface > size(alpha)) then
    write(6,fmt="(//13x,'number of faces:',i5,' greater than dimensioned',i5)") iface,size(alpha)
    stop
  endif
  x = xx - z
  rx = sqrt(x(1)*x(1) + x(2)*x(2) + x(3)*x(3))
  rz = sqrt(z(1)*z(1) + z(2)*z(2) + z(3)*z(3))
  s = x(1)*z(1) + x(2)*z(2) + x(3)*z(3)

  if (rx < tolerance .or. rz < tolerance .or. s > tolerance) then
    write(6,fmt="(/13x,'from euler,illegal vectors:'/13x,2(' (',3e13.6,' )'))") x(1:3), z(1:3)
    stop
  endif

  p = sqrt(z(1)*z(1) + z(2)*z(2))

  x = x/rx
  z = z/rz

  alpha(iface) = 0.d0
  beta(iface) = acos(z(3))

  if (p < toleuler) then
  
    sg = -z(3)*x(2)
    cg =  z(3)*x(1)

    if (abs(sg) < toleuler .and. abs(cg+1d0) < toleuler)  then
      gamma(iface) = pi
    else
      gamma(iface) = 2d0*atan2(sg, cg+1d0)
    endif

  else

    rzp = rz/p
    ! y = z .cross. x or x .cross. z?
    y(1) = z(2)*x(3) - z(3)*x(2)
    y(2) = z(3)*x(1) - z(1)*x(3)
    y(3) = z(1)*x(2) - z(2)*x(1)
    
    sa = y(3)*rzp
    ca = x(3)*rzp

    if (abs(sa) < toleuler .and. abs(ca+1.d0) < toleuler)  then
      alpha(iface) = pi
    else
      alpha(iface) = 2d0*atan2(sa, ca+1d0)
    endif

    sg= z(2)*rzp
    cg=-z(1)*rzp

    if (abs(sg) < toleuler .and. abs(cg+1d0) < toleuler)  then
      gamma(iface) = pi
    else
      gamma(iface) = 2d0*atan2(sg, cg+1d0)
    endif

  endif
    
  z = z*rz

endsubroutine ! euler


subroutine euler_arg(z, xx, toleuler, alpha, beta, gamma)
  use shape_constants_mod, only: pi, dp
  double precision, intent(in) :: xx(3)
  double precision, intent(inout) :: z(3)
  double precision, intent(in) :: toleuler   ! introduced by phivos (05.2008) to account for inaccuracies.
  ! earlier, 1.d-5 was hard-coded at the places in this subr. where tol is used
  !     data toleuler /1.d-10/
  double precision, intent(out) :: alpha, beta, gamma

  double precision :: rx,rz,s,p,rzp,sa,ca,sg,cg
  double precision :: x(3), y(3)
  double precision, parameter :: tolerance = 1d-6

  !-----------------------------------------------------------------------
  x = xx - z
  rx = sqrt(x(1)*x(1) + x(2)*x(2) + x(3)*x(3))
  rz = sqrt(z(1)*z(1) + z(2)*z(2) + z(3)*z(3))
  s = x(1)*z(1) + x(2)*z(2) + x(3)*z(3)

  if (rx < tolerance .or. rz < tolerance .or. s > tolerance) then
    write(6,fmt="(/13x,'from euler,illegal vectors:'/13x,2(' (',3e13.6,' )'))") x(1:3), z(1:3)
    stop
  endif

  p = sqrt(z(1)*z(1) + z(2)*z(2))

  x = x/rx
  z = z/rz

  alpha = 0.d0
  beta = acos(z(3))

  if (p < toleuler) then
  
    sg = -z(3)*x(2)
    cg =  z(3)*x(1)

    if (abs(sg) < toleuler .and. abs(cg+1d0) < toleuler)  then
      gamma = pi
    else
      gamma = 2d0*atan2(sg, cg+1d0)
    endif

  else

    rzp = rz/p
    ! y = z .cross. x or x .cross. z?
    y(1) = z(2)*x(3) - z(3)*x(2)
    y(2) = z(3)*x(1) - z(1)*x(3)
    y(3) = z(1)*x(2) - z(2)*x(1)
    
    sa = y(3)*rzp
    ca = x(3)*rzp

    if (abs(sa) < toleuler .and. abs(ca+1.d0) < toleuler)  then
      alpha = pi
    else
      alpha = 2d0*atan2(sa, ca+1d0)
    endif

    sg= z(2)*rzp
    cg=-z(1)*rzp

    if (abs(sg) < toleuler .and. abs(cg+1d0) < toleuler)  then
      gamma = pi
    else
      gamma = 2d0*atan2(sg, cg+1d0)
    endif

  endif
    
  z = z*rz

endsubroutine ! euler




!-----------------------------------------------------------------------
!> rotation by euler angles given in module angles_common.
!>    this routine performs the rotation of nvert vectors through the
!>    euler angles: alpha(iface),beta(iface),gamma(iface).
!>    v (i,ivert) : input   vectors
!>    vz(i,ivert) : rotated vectors
!-----------------------------------------------------------------------
subroutine rotate_common(v, vz, iface, nvert)
  use shape_constants_mod, only: pi, dp
  use angles_common, only: alpha, beta, gamma
  integer, intent(in) :: nvert
  double precision, intent(in) :: v(3,*) ! (3,nvertd)
  double precision, intent(out) :: vz(3,*) ! (3,nvertd)
  integer, intent(in) :: iface

  integer :: j, ivert
  double precision :: sa,sb,sg,ca,cb,cg
  double precision :: a(3,3)

  ca = cos(alpha(iface))
  sa = sin(alpha(iface))
  cb = cos(beta(iface))
  sb = sin(beta(iface))
  cg = cos(gamma(iface))
  sg = sin(gamma(iface))
  
  a(1,1) = ca*cb*cg - sa*sg
  a(2,1) = sa*cb*cg + ca*sg
  a(3,1) = -sb*cg
  a(1,2) = -ca*cb*sg - sa*cg
  a(2,2) = -sa*cb*sg + ca*cg
  a(3,2) = sb*sg
  a(1,3) = ca*sb
  a(2,3) = sa*sb
  a(3,3) = cb
  
  do ivert = 1, nvert
    vz(1:3,ivert) = 0.d0
    do j = 1, 3
      vz(1:3,ivert) = vz(1:3,ivert) + a(1:3,j)*v(j,ivert)
    enddo ! j
  enddo ! ivert
endsubroutine ! rotate

subroutine rotate_arg(v, vz, nvert, alpha, beta, gamma)
  use shape_constants_mod, only: pi
  integer, intent(in) :: nvert
  double precision, intent(in) :: alpha, beta, gamma
  double precision, intent(in) :: v(3,*) ! (3,nvertd)
  double precision, intent(out) :: vz(3,*) ! (3,nvertd)

  integer :: j, ivert
  double precision :: sa,sb,sg,ca,cb,cg
  double precision :: a(3,3)

  ca = cos(alpha)
  sa = sin(alpha)
  cb = cos(beta)
  sb = sin(beta)
  cg = cos(gamma)
  sg = sin(gamma)
  
  a(1,1) = ca*cb*cg - sa*sg
  a(2,1) = sa*cb*cg + ca*sg
  a(3,1) = -sb*cg
  a(1,2) = -ca*cb*sg - sa*cg
  a(2,2) = -sa*cb*sg + ca*cg
  a(3,2) = sb*sg
  a(1,3) = ca*sb
  a(2,3) = sa*sb
  a(3,3) = cb
  
  do ivert = 1, nvert
    vz(1:3,ivert) = 0.d0
    do j = 1, 3
      vz(1:3,ivert) = vz(1:3,ivert) + a(1:3,j)*v(j,ivert)
    enddo ! j
  enddo ! ivert
endsubroutine ! rotate

endmodule ShapeCriticalPoints_mod
