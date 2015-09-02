#include "macros.h"
!> Module to determine radial mesh panel positions.

module ShapeCriticalPoints_mod
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  public :: criticalShapePoints
  
  contains

  !------------------------------------------------------------------------------
  !> Finds critical points for mesh generation and performs geometrical tests.
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
  subroutine criticalShapePoints(planes, tolvdist, toleuler, nvertices, vert, nface,lmax, faces, npan, crt)

    use shape_constants_mod, only: verbosity, check_geometry, pi !, isumd, lmaxd1
    use PolygonFaces_mod, only: PolygonFace
    use ShapeGeometryHelpers_mod, only: polchk

    integer, intent(in) :: nvertices(:) ! (nfaced)
    double precision, intent(in) :: planes(0:,:) ! (0:3,nfaced)
    double precision, intent(in) :: vert(:,:,:) ! (3,nvertd,nfaced)
    integer, intent(in) :: nface, lmax
    double precision, intent(in) :: tolvdist, toleuler
    
    type(PolygonFace), intent(out) :: faces(:)
    integer, intent(out) :: npan
    double precision, intent(out) :: crt(:)
    
    
    !-----------------------------------------------------------------------
    
!   double precision :: sq3o3, coa
    double precision :: z(3)
    integer :: iface, itt
    integer :: iv, ivtot, l
    integer :: nvertd
    integer :: ist

    !-----------------------------------------
    !  this call does some geometrical tests (n.stefanou 98)
    if (check_geometry) call polchk(nface, nvertices, vert, tolvdist)

! #define k_HCP 
#ifdef k_HCP
    ! hexagonal closed package - lattice structure
    coa = 2.d0
    sq3o3 = sqrt(3.d0)/3.d0
    die_here("case K_HCP is decatived!")
#endif

    npan = 0 ! will be modified by routine critical_points
    ivtot = 0 ! deprecated, for comparable verbose output only
    !.......................................................................
    ! calculation of rotation matrices
    !.......................................................................
    do iface = 1, nface

      z(1:3) = planes(1:3,iface)/planes(0,iface) ! plane normal

! #ifdef k_HCP 
!       z(1) = z(1)*sq3o3
!       z(3) = z(3)*8.d0/(coa*3.d0)
! #endif      
! 
!       do iv = 1, nvertices(iface)
!         v(1:3,iv) = vert(1:3,iv,iface)
! 
! #ifdef k_HCP 
!         v(1,iv) = v(1,iv)*sq3o3
!         v(3,iv) = v(3,iv)*coa
! #endif
!       enddo ! iv

      call critical_points(faces(iface), nvertices(iface), vert(:,:,iface), z, npan, ivtot, toleuler, tolvdist, crt, face_index=iface)

      if (verbosity > 0) write(6,fmt="(/10x,i3,'-th pyramid subdivided in ',i3,' tetrahedra')") iface,faces(iface)%ntt

    enddo ! iface ! end of loop over faces

!     deallocate(v, stat=ist)
    
    !.......................................................................
    !     definition of the suitable mesh
    !.......................................................................

    if (verbosity > 1) then
      write(6,fmt="(//15x,'fa/pi',5x,'fb/pi',5x,'fd/pi',6x,'rd',8x,'isignu'/)")
      iv = 0
      do iface = 1, nface
        do itt = 1, faces(iface)%ntt
          iv = iv+1
#define tet faces(iface)%ta(itt)
          write(6,fmt="(i10,4f10.4,i10)") iv, [tet%fa, tet%fb, tet%fd]/pi, tet%rd, tet%isignu
#undef tet
        enddo ! itt
      enddo ! iface
    endif ! verbosity
    
  endsubroutine criticalShapePoints

  
  double precision function get_angle(sine, cosine, tolerance) result(angle)
    use shape_constants_mod, only: pi
    double precision, intent(in) :: sine, cosine, tolerance
    if (abs(sine) < tolerance .and. abs(cosine + 1.d0) < tolerance) then
      angle = pi
    else
      angle = 2.d0*atan2(sine, cosine + 1.d0)
    endif
  endfunction get_angle
  

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
  subroutine critical_points(face, nvert, v, z, ipan, ivtot, toleuler, tolvdist, crt, face_index)
    !-----------------------------------------------------------------------
    !     this routine calculates the critical points 'crt' of the shape
    !     functions due to the face: z(1)*x + z(2)*y + z(3)*z = 1
    !     the face is rotated through the appropriate euler angles to be
    !     perpendicular to the z-axis. a further subdivision of the cen-
    !     tral pyramid into elementary tetrahedra is performed. the  ne-
    !     cessary quantities for the calculation are stored in common.
    !-----------------------------------------------------------------------
    use shape_constants_mod, only: pi, verbosity
    use PolygonFaces_mod, only: TetrahedronAngles, PolygonFace
    use shapegeometryhelpers_mod, only: perp, nrm2, operator(.dot.)

    integer, intent(in) :: nvert
    type(PolygonFace), intent(inout) :: face
    integer, intent(inout) :: ipan, ivtot ! ivtot can be removed when we do not try to have the verbose output compatible with older versions
    double precision, intent(in) :: toleuler, tolvdist
    double precision, intent(in) :: v(:,:) ! (3,nvertd)
    double precision, intent(inout) :: z(3)
    double precision, intent(inout) :: crt(:)
    integer, intent(in) :: face_index ! is only needed for verbose output since we pass the face descriptor to this routine
    
    integer :: iv, ivert, ivertp, jv, npand
    double precision :: arg, a1, a2, a3, cf1, cf2, cf3, co, crrt, dd, down, d1, d2, ff, f1, f2, f3, omega, rdd, s, sf1, sf2, sf3, up, xj, yj, zmod2

    logical(kind=1), allocatable :: in(:) ! (nvertd)
    double precision, allocatable :: vz(:,:) ! (3,nvertd)
    double precision :: rdv(3)
    logical :: inside, new, corner
    double precision, parameter :: origin(3) = 0.d0
    double precision, parameter :: tol_small = 1d-6
    double precision, parameter :: tol_large = 1d-4
    type(TetrahedronAngles) :: ta(nvert), t1
    integer :: ist
    
    !-----------------------------------------------------------------------
    npand = size(crt)
    
    allocate(in(nvert), vz(3,nvert), stat=ist)


    if (verbosity > 0) write(6,fmt="(//80('*')/3x,'face:',i3,' equation:',f10.4,'*x +',f10.4,'*y +',f10.4,'*z  =  1')") face_index,z(1:3)

    zmod2 = nrm2(z)

    if (zmod2 <= tol_small) & ! check if normal vector of plane is valid
      die_here("the"+face_index-"-th face of the polyhedron passes through the center"+z)
 
    s = 2.d0*pi

    z = z/zmod2

    iv = 1; if (nrm2(v(1:3,1) - z) < tol_small**2) iv = 2 ! if the norm of the first vector is too small, use the second one
    face%euler = euler_angles(z, v(:,iv), toleuler) ! get euler angles directly

    if (verbosity > 0) write(6,fmt="(3x,'rotation angles  :',3(f10.4,4x)/)") face%euler(1:3)/pi

    call rotate(face%euler(1:3), nvert, v(1:3,1:nvert), vz) ! pass the angles directly

    face%r0 = 1.d0/sqrt(zmod2)
    
    corner = .false. ! is not a corner

    if (verbosity > 0) write(6,fmt="(/'tetrahedron',14x,'coordinates'/11('*'),14x,11('*')/)")

    face%ntt = 0
    
    do ivert = 1, nvert
    
      if (abs(face%r0 - vz(3,ivert)) > tol_small) & ! check if vertices lie in same plane
        die_here("the vertices of the"+face_index-"-th rotated polygon do not lie on the plane, r0 ="+face%r0+" z ="+vz(3,ivert))
      !.......................................................................
      !____distances__of__vertices__from__center
      !.......................................................................
      crrt = sqrt(nrm2(vz(1:3,ivert)))

      !     check if a new panel has been found (at radial point crrt)
      new = all(abs(crrt - crt(1:ipan)) >= tol_small)

      if (new) then
        ipan = ipan + 1 ! add a new critical point
        if (ipan > npand) then
          die_here("number of panels ="+ipan+" is greater than dimensioned ="+npand)
        else
          crt(ipan) = crrt
        endif
      endif ! new
      
      ivertp = modulo(ivert, nvert) + 1

      !.......................................................................
      !____distances__of__edges__from__center
      !.......................................................................

      call perp(origin, vz(1:3,ivert), vz(1:3,ivertp), tolvdist, rdv, inside)

      rdd = sqrt(nrm2(rdv)) ! footpoint of line origin-to-edge

      if (inside) then  
      
        ! add a new panel only when footpoint is contained on edge
        new = all(abs(rdd - crt(1:ipan)) >= tol_large)
        
        if (new) then
          ipan = ipan + 1 ! add a new critical point
          if (ipan > npand) then
            die_here("number of panels ="+ipan+" is greater than dimensioned ="+npand)
          else
            crt(ipan) = rdd
          endif
        endif ! new
        
      endif ! inside
      
      a1 = sqrt(vz(1,ivert)**2 + vz(2,ivert )**2)
      a2 = sqrt(vz(1,ivertp)**2 + vz(2,ivertp)**2)
      up = vz(1,ivert)*vz(1,ivertp) + vz(2,ivert)*vz(2,ivertp)
      down = a1*a2

      if (down > tol_small) then ! true if not a corner
        arg = up/down
        if (abs(arg) >= 1.d0) arg = sign(1.d0, arg)
        omega = acos(arg)
        s = s - omega

        if (abs(omega - pi) > tol_small) then
          !.......................................................................
          !____subdivision____into____tetrahedra
          !.......................................................................
          face%ntt = face%ntt + 1 ! count up the number of accepted tetrahedra
          
          ivtot = ivtot + 1 ! can be removed when we do not try to have the verbose output compatible with older versions
          if (verbosity > 0) then
            write(6,fmt="(i5,'       vz(',i2,')  =  (',3f10.4,' )')") ivtot,ivert, vz(1:3,ivert )
            write(6,fmt="(5x,'       vz(',i2,')  =  (',3f10.4,' )')")       ivertp,vz(1:3,ivertp)
          endif

          a3 = sqrt(rdv(1)*rdv(1) + rdv(2)*rdv(2))
          
          cf1 = vz(1,ivert )/a1
          cf2 = vz(1,ivertp)/a2
          sf1 = vz(2,ivert )/a1
          sf2 = vz(2,ivertp)/a2
          cf3 = rdv(1)/a3
          sf3 = rdv(2)/a3

          f1 = get_angle(sf1, cf1, tolerance=toleuler)
          f2 = get_angle(sf2, cf2, tolerance=toleuler)
          f3 = get_angle(sf3, cf3, tolerance=toleuler)

          t1%fa = min(f1, f2)
          t1%fb = max(f1, f2)
          t1%fd = f3
          
          if ((t1%fb - t1%fa) > pi) then
            ff = t1%fa + 2.d0*pi ! increase new fb
            t1%fa = t1%fb ! swap
            t1%fb = ff
          endif
          if ((t1%fa - t1%fd) > pi) t1%fd =  2.d0*pi + t1%fd
          if ((t1%fd - t1%fa) > pi) t1%fd = -2.d0*pi + t1%fd
          
          t1%isignu = 1
          t1%rd = rdd
!         t1%rupsq = sqrt(rdd**2 - face%r0**2) ! could be introduced here
          
          ta(face%ntt) = t1 ! copy
          
        endif ! |omega - pi| > tol_small
        
      else  ! down > tol_small
      
        corner = .true. ! is a corner
        
      endif ! down > tol_small
      
    enddo ! ivert ! end of vertex loop
    
    
    deallocate(face%ta, stat=ist)
    allocate(face%ta(face%ntt), stat=ist)
    face%ta = ta(1:face%ntt) ! copy
    
    !.......................................................................
    ! foot_of_the_perpendicular_to_the_face_outside_or_inside_the_polygon
    !.......................................................................
    if (s < tol_small .or. corner) then
    
      new = all(abs(face%r0 - crt(1:ipan)) >= tol_large)
      
      if (new) then
        ipan = ipan + 1 ! add a new critical point
        if (ipan > npand) then
          die_here("number of panels ="+ipan+" is greater than dimensioned ="+npand)
        else
          crt(ipan) = face%r0 ! found a new critical point
        endif
      endif ! new
      
    else  ! s < tol_small .or. corner
    
      do jv = 1, nvert
        in(jv) = .false.
        do iv = 1, nvert
        
          ivertp = modulo(iv, nvert) + 1

          if (iv /= jv .and. ivertp /= jv) then
          
            down = vz(2,jv)*(vz(1,ivertp) - vz(1,iv)) - vz(1,jv)*(vz(2,ivertp) - vz(2,iv))

            if (abs(down) > tol_small) then

              up  = vz(1,jv)*(vz(2,iv)*(vz(1,ivertp) + vz(1,iv)) - vz(1,iv)*(vz(2,ivertp) + vz(2,iv)))
              xj = up/down
              yj = xj*vz(2,jv)/vz(1,jv)
              dd = (vz(1,ivertp) - vz(1,iv))**2 + (vz(2,ivertp) - vz(2,iv))**2
              d1 = (xj - vz(1,iv ))**2 + (yj - vz(2,iv ))**2
              d2 = (xj - vz(1,ivertp))**2 + (yj - vz(2,ivertp))**2

              co = dd - max(d1, d2)

              if (co > tol_small) then
                in(jv) = .true.
                exit ! iv loop (inner loop)
              endif

            endif ! abs(down) > tol_small
          endif
  
        enddo ! iv
      enddo ! jv
      
      do ivert = 1, nvert
        ivertp = modulo(ivert, nvert) + 1
        if (.not. in(ivert) .and. .not. in(ivertp)) then
          if (ivert > face%ntt) then
            die_here("a sign could not be stored! the face has"+face%ntt+"tetrahedra, but we try to access #"-ivert)
          else
            face%ta(ivert)%isignu = -1
          endif
        endif
      enddo ! ivert

    endif ! s < tol_small .or. corner
    
    deallocate(in, vz, stat=ist)
  endsubroutine critical_points

  !-----------------------------------------------------------------------
  !>    given two distinct points (z(1),z(2),z(3)) and (xx(1),xx(2),xx(3))
  !>    this routine defines  a local coordinate  system with the  z- axis
  !>    passing through (z(1),z(2),z(3))  and the x- axis parallel to  the
  !>    vector : (xx(1)-z(1),xx(2)-z(2),xx(3)-z(3)).
  !>    the euler angles rotating this local coordinate system back to the
  !>    original frame of reference are calculated  and stored  in common.
  !-----------------------------------------------------------------------
  function euler_angles(zz, xx, toleuler) result(alpha_beta_gamma)
    use shape_constants_mod, only: pi
!     use ShapeGeometryHelpers, only: nrm2, operator(.dot.) ! todo

    double precision, intent(in) :: xx(3)
    double precision, intent(in) :: zz(3)
    double precision, intent(in) :: toleuler ! introduced by phivos (05.2008) to account for inaccuracies.
    ! earlier, 1.d-5 was hard-coded at the places in this subr. where tol is used data toleuler /1.d-10/
    double precision :: alpha_beta_gamma(3) ! result

    double precision :: rx,rz,s,p,rzp,sa,ca,sg,cg
    double precision :: x(3), y(3), z(3)
    double precision, parameter :: tolerance = 1d-6

    !-----------------------------------------------------------------------
    z = zz
    x = xx - z
    rx = sqrt(x(1)*x(1) + x(2)*x(2) + x(3)*x(3)) ! sqrt(nrm2(x))
    rz = sqrt(z(1)*z(1) + z(2)*z(2) + z(3)*z(3)) ! sqrt(nrm2(z))
    s = x(1)*z(1) + x(2)*z(2) + x(3)*z(3) ! x .dot. z

    if (rx < tolerance .or. rz < tolerance .or. s > tolerance) &
      die_here("Euler_angles: illegal vectors x ="+x+" z ="+z)

    p = sqrt(z(1)*z(1) + z(2)*z(2))

    x = x/rx
    z = z/rz

    alpha_beta_gamma(2) = acos(z(3)) ! beta

    if (p < toleuler) then
    
      sg = -z(3)*x(2)
      cg =  z(3)*x(1)

      alpha_beta_gamma(1) = 0.d0             ! alpha
      alpha_beta_gamma(3) = get_angle(sg, cg, tolerance=toleuler) ! gamma

    else  ! p < toleuler

      rzp = rz/p
      ! y = z .cross. x or x .cross. z?
      y(1) = z(2)*x(3) - z(3)*x(2)
      y(2) = z(3)*x(1) - z(1)*x(3)
      y(3) = z(1)*x(2) - z(2)*x(1)
      
      sa = y(3)*rzp
      ca = x(3)*rzp

      alpha_beta_gamma(1) = get_angle(sa, ca, tolerance=toleuler) ! alpha

      sg =  z(2)*rzp
      cg = -z(1)*rzp

      alpha_beta_gamma(3) = get_angle(sg, cg, tolerance=toleuler) ! gamma

    endif ! p < toleuler

  endfunction ! euler_angles


  !-----------------------------------------------------------------------
  !> rotation by euler angles given in module angles_common.
  !>    this routine performs the rotation of nvert vectors through the
  !>    euler angles: input alpha, beta, gamma
  !>    v (i,ivert) : input   vectors
  !>    vz(i,ivert) : rotated vectors
  !-----------------------------------------------------------------------
  subroutine rotate(abg, nvert, v, vz)
    use shape_constants_mod, only: pi
    
    double precision, intent(in) :: abg(3) ! former alpha, beta, gamma
    integer, intent(in) :: nvert
    double precision, intent(in) :: v(3,*) ! (3,nvertd)
    double precision, intent(out) :: vz(3,*) ! (3,nvertd)

    integer :: j, ivert
    double precision :: sn(3), cs(3), a(3,3)

    cs = cos(abg(1:3))
    sn = sin(abg(1:3))
    
    a(1,1) = cs(1)*cs(2)*cs(3) - sn(1)*sn(3)
    a(2,1) = sn(1)*cs(2)*cs(3) + cs(1)*sn(3)
    a(3,1) = -sn(2)*cs(3)
    a(1,2) = -cs(1)*cs(2)*sn(3) - sn(1)*cs(3)
    a(2,2) = -sn(1)*cs(2)*sn(3) + cs(1)*cs(3)
    a(3,2) = sn(2)*sn(3)
    a(1,3) = cs(1)*sn(2)
    a(2,3) = sn(1)*sn(2)
    a(3,3) = cs(2)
    
    do ivert = 1, nvert
      vz(1:3,ivert) = 0.d0
      do j = 1, 3
        vz(1:3,ivert) = vz(1:3,ivert) + a(1:3,j)*v(j,ivert)
      enddo ! j
    enddo ! ivert
    
  endsubroutine ! rotate

endmodule ShapeCriticalPoints_mod
