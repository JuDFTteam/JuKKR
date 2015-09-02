#include "DebugHelpers/test_macros.h"

!> Auxillary module for shape function calculation.

module ShapeIntegration_mod
  implicit none
  private
  public :: shapeIntegration

  contains

!------------------------------------------------------------------------------
!> performs the shape-function integration.
!> @brief
!> needs information from module/common block replacement "PolygonFaces_mod"
!> @param[in]  lmax calculate shape function up to lmax
!> @param[in]  nface number of cell faces
!> @param[out] meshn number of radial mesh points
!> @param[in]  xrn radial mesh points
!> @param[in]  dlt step size for gauss-legendre integration
!> @param[out] thetas_s on output it contains the non-zero shape-functions
!!                      dimension meshnd x ibmaxd
!> @param[out] lmifun_s array with lm-indeices for which the shapefunction is non-zero
!> @param[out] nfun number of non-zero shape-functions
!> @param[in]  meshnd
!> @param[in]  ibmaxd

subroutine shapeIntegration(lmax, face, meshn, xrn, dlt, thetas_s, lmifun_s, nfun, meshnd, ibmaxd)

  use shape_constants_mod, only: pi !, lmaxd1!, icd, iced!, isumd
  use PolygonFaces_mod, only: PolygonFace, TetrahedronAngles
  use shapeIntegrationhelpers_mod, only: pintg, ccoef, d_real

  integer, intent(in) :: lmax, meshn
  double precision, intent(in) :: xrn(meshnd)
  double precision, intent(in) :: dlt
  integer, intent(in) :: meshnd
  integer, intent(in) :: ibmaxd
  type(PolygonFace), intent(in) :: face(:)
  integer, intent(out) :: lmifun_s(ibmaxd)
  double precision, intent(out) :: thetas_s(meshnd,ibmaxd)
  integer, intent(out) :: nfun
  

  double precision :: cl_table(0:lmax,0:lmax,0:lmax) !, cl(icd) 
  double precision :: c_table(0:lmax,0:lmax) !, c(iced) 

  double precision :: rap, rdown, arg1, arg2, fk, fl, fpisq, rupsq

  integer :: iface, itet, m, isu, isu0, l, mp, k, lp, ilm, mlm, ir, ist, isumd, nface 
!   integer :: ib, ic0, ice0, i, ic, ice, imax, mo, ip, ipmax
  
  double precision, allocatable :: dmatl(:,:) ! (isumd,nface)
  type(TetrahedronAngles) :: t1
  
  double precision :: s(-lmax:lmax,0:lmax), s1(-lmax:lmax,0:lmax) ! this storage format uses only (lmax+1)**2 of (lmax+1)*(2*lmax+1) elements so roughly 55%
  double precision :: smm(-lmax:lmax)

  integer :: imo(-lmax:lmax) ! index translation from m=-l,-l+1,...,l-1,l to 0,1,-1,...,l,-l since the dmatl make use of this ordering
  
  ! local automatic arrays
  logical(kind=1) :: nontrivial(ibmaxd) ! the trivial shape function is zero from ilm>1 and sqrt(4*pi) for ilm==1
  double precision :: b(ibmaxd)
  integer, parameter :: idmatl_RECOMPUTE = 0 ! recompute every time, calls d_real meshn*nface times
  integer, parameter :: idmatl_MEMORIZE  = 1 ! needs some memory but calls d_real only nface times
  integer, parameter :: idmatl = idmatl_MEMORIZE ! must be in {0, 1}
  
  CHECKASSERT(idmatl == idmatl**2) ! check that idmatl is in {0, 1}
  
  fpisq = sqrt(4.d0*pi)

  mlm = (lmax+1)**2

  thetas_s = 0.d0

  imo = 0
  do m = 1, lmax
    imo( m) = 2*m-1
    imo(-m) = 2*m
  enddo ! m
  
  !.......................................................................
  !     expansion coefficients
  !.......................................................................
  call ccoef(lmax, cl_table, c_table) ! sets the coefficients: 
  
  nface = size(face, 1) ! number of faces
  
  isumd = (lmax*(3+4*(lmax+1)*(lmax+2)))/3+1 !< number of compressed rotation matrix elements = sum_l=0..lmax (2*l+1)^2
  allocate(dmatl(isumd,0:nface*idmatl), stat=ist)
  
  if (idmatl == idmatl_MEMORIZE) then
    do iface = 1, nface
      ! calculate transformation matrices for spherical harmonics (once only)
      call d_real(lmax, face(iface)%euler, dmatl(:,iface))
    enddo ! iface
  endif ! memorize
  
  nontrivial(1:mlm) = .false. ! init

  !===================== split ??? =======================================

  !.......................................................................
  !     loop over radial mesh points
  !.......................................................................
  meshloop: do ir = 1, meshn
  
    b(1:mlm) = 0.d0 ! init shape function temporary at this radius
    
    !.......................................................................
    !     loop over pyramids
    !.......................................................................
py: do iface = 1, nface

      if (xrn(ir) <= face(iface)%r0) then
        ! the radius is smaller than the distance of the base points from the origin.
        ! so, no contribution to the shapefunctions is expected (except for l=0,m=0 which is treated analytically).
        cycle py ! jump to the head of the loop over pyramids (loop counter iface)
        
      endif ! xrn(ir) <= r0
      
      arg1 = face(iface)%r0/xrn(ir)
      s = 0.d0 ! init s
      
      !.......................................................................
      !     loop over tetrahedra
      !.......................................................................
      do itet = 1, face(iface)%ntt

        t1 = face(iface)%ta(itet) ! get a work copy of the TetrahedronAngles

        if (xrn(ir) <= t1%rd) then
        
          call pintg(t1%fa, t1%fb, dlt, s1, lmax, t1%isignu, arg1, t1%fd, 0)
          s = s + s1
          
        else  ! xrn(ir) <= t1%rd
        
          ! rupsq could be precomputed since it depends only on properties of the tetrahedron and the face
          rupsq = sqrt(t1%rd**2 - face(iface)%r0**2) 
          
          rdown = sqrt(xrn(ir)**2 - face(iface)%r0**2)
          rap  = rupsq/rdown
          arg2 = rupsq/face(iface)%r0
          fk = min(max(t1%fa, t1%fd - acos(rap)), t1%fb)
          fl = min(max(t1%fa, t1%fd + acos(rap)), t1%fb)

          call pintg(t1%fa, fk, dlt, s1, lmax, t1%isignu, arg1, t1%fd, 0)
          s = s + s1
          call pintg(fk,    fl, dlt, s1, lmax, t1%isignu, arg2, t1%fd, 1)
          s = s + s1
          call pintg(fl, t1%fb, dlt, s1, lmax, t1%isignu, arg1, t1%fd, 0)
          s = s + s1
          
        endif ! xrn(ir) <= t1%rd

      enddo ! itet ! tetraeder loop

      !.......................................................................
      !  integral expansion back-rotation
      !.......................................................................

      if (idmatl /= idmatl_MEMORIZE) then
        ! calculate transformation matrices for spherical harmonics (for each ir again!)
        call d_real(lmax, face(iface)%euler, dmatl(:,0))
      endif ! recompute every time

      do l = 0, lmax
      
        ! transform to summed coefficients
        
        ! cl_table and c_table are symmetric w.r.t. m <--> -m, so no negative indices
        
        smm = 0.d0
        do m = -l, l
!         smm(m) = 0.d0
          do k = l, (l+abs(m)+1)/2, -1
            lp = 2*k-l-abs(m)
            smm(m) = smm(m) + cl_table(abs(m),lp,l)*s(m,lp)
          enddo ! k
          smm(m) = smm(m)*c_table(abs(m),l) ! scale
        enddo ! m

        isu0 = ((l-1)*(3+4*l*(l+1)))/3+1
        isu = isu0
        
        ! rotate back in the l-subspace and add to the shape functions
        ! still uses the m-ordering 0,1,-1,...,l,-l for m and mp
!         imax = 1
!         do m = 0, l
!           do i = 1, imax
!             mo = (3-2*i)*m ! mapping i=1 --> mo=m and i=2 --> mo=-m (reconstruction of signed m quantum number)
!             ilm = l*l+l+mo+1
!             ipmax = 1
!             do mp = 0, l
!               do ip = 1, ipmax
!                 isu = isu+1
!                 b(ilm) = b(ilm) + sm(ip,mp)*dmatl(isu,idmatl*iface)
!               enddo ! ip
!               ipmax = 2
!             enddo ! mp
!           enddo ! i
!           imax = 2
!         enddo ! m

        do m = -l, l
          ilm = l*l+l+m+1
          do mp = -l, l
            isu = isu0 + 1 + (2*l+1)*imo(m) + imo(mp) ! index for dmatl(isu) since it is ordered m=0,1,-1,...,l,-l
!           isu = isu0 + (2*l+1)*(m+l) + (mp+l) + 1 ! index for dmatl(isu) if it were ordered m=-l...l
            b(ilm) = b(ilm) + smm(mp)*dmatl(isu,idmatl*iface)
          enddo ! mp
        enddo ! m
        
      enddo ! l

    enddo py
    !.......................................................................
    !     defines and saves shapefunctions
    !.......................................................................
    b(1:mlm) = -b(1:mlm)/fpisq ! scale
    b(1) = fpisq + b(1) ! add constant sqrt(4*pi) in the l=0,m=0 channel
    
    nontrivial = nontrivial .or. (abs(b(:)) > 1d-6)
    
    thetas_s(ir,1:mlm) = b(1:mlm)
    
  enddo meshloop
  
  deallocate(dmatl, stat=ist)

  !now rearrange thetas_s array that it contains only non-zero shapefunctions, this is done "in-place"

  lmifun_s = 0 ! lm-index of this shape function

  nfun = 0 ! number of non-zero shape function
  do ilm = 1, mlm ! this loop must run ordered in serial
    if (nontrivial(ilm)) then

      nfun = nfun + 1
      lmifun_s(nfun) = ilm

      if (nfun < ilm) thetas_s(1:meshn,nfun) = thetas_s(1:meshn,ilm) ! compress the storage by overwriting functions that are zero or close to zero

    else
      thetas_s(1:meshn,ilm) = 0.d0 ! clear to exact zeros
    endif ! nontrivial
  enddo ! ilm
  
endsubroutine shapeIntegration

endmodule ShapeIntegration_mod
