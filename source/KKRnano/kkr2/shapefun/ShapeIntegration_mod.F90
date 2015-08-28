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

subroutine shapeIntegration(lmax, nface, meshn, xrn, dlt, thetas_s, lmifun_s, nfun, meshnd, ibmaxd)

  use shape_constants_mod, only: pi, lmaxd1, isumd, icd, iced
! use PolygonFaces_mod, only: rd, fa, fb, fd, isignu ! vertex properties read-only
! use PolygonFaces_mod, only: ntt, r0, alpha, beta, gamma ! face properties read-only ! deprecated
  use PolygonFaces_mod, only: face ! face properties
  use PolygonFaces_mod, only: TetrahedronAngles
  
  use shapeIntegrationhelpers_mod, only: pintg, ccoef, d_real

  integer, intent(in) :: lmax, nface, meshn
  real*8, intent(in) :: xrn(meshnd)
  real*8, intent(in) :: dlt
  real*8, intent(out) :: thetas_s(meshnd,ibmaxd)
  integer, intent(out) :: lmifun_s(ibmaxd)
  integer, intent(in) :: meshnd
  integer, intent(in) :: ibmaxd
  integer, intent(out) :: nfun

  real*8 :: dmatl(isumd)
  real*8 :: cl(icd)
  real*8 :: c(iced), c_table(-lmaxd1:lmaxd1,0:lmaxd1)

  real*8 :: rap
  real*8 :: rdown

  real*8 :: arg1
  real*8 :: arg2

  integer :: ivtot, iface, itet, i, m, ic, ice, isu, isu0, l, mp, k0, k, is, imax, mo, ip, ipmax, ilm, mlm, ir
!   integer :: ib, ic0, ice0
  
  real*8, allocatable :: rupsq_precomputed(:) ! dim nvtotd
  real*8, allocatable :: dmatl_precomputed(:,:)
  type(TetrahedronAngles) :: t1
  
  real*8 :: s(-lmaxd1:lmaxd1,0:lmaxd1), s1(-lmaxd1:lmaxd1,0:lmaxd1) ! this storage format usese only (lmaxd1+1)**2 of (lmaxd1+1)*(2*lmaxd1+1) elements so roughly 55%
  real*8 :: sm(2,0:lmaxd1) ! this could be reshaped into sm(-lmaxd1:lmaxd1)
!   real*8 :: smm(-lmaxd1:lmaxd1)
  real*8 :: fk, fl, fpisq, rupsq

  ! local automatic arrays
  logical(kind=1) :: nonzero(ibmaxd)
  real*8 :: b(ibmaxd)
  logical, parameter :: precompute_dmtl = .true. ! needs some memory but reduces the number of calls to d_real from meshn*nface to nface
  logical, parameter :: precompute_rupsq = .false.
  
!   if (precompute_rupsq) allocate(rupsq_precomputed(size(rd)))
  if (precompute_dmtl) allocate(dmatl_precomputed(isumd,nface))
  
  fpisq = sqrt(4.d0*pi)

  mlm = (lmax+1)**2

  thetas_s = 0.d0

  !.......................................................................
  !     expansion coefficients
  !.......................................................................
  call ccoef(lmax, cl,c ) ! sets the coefficients: 
  !  cl in m-order (0,(m,-m, m=1..l), l=0,lmax) ??
  ! and c in order ((m, m=l...0), l=0,lmax)
  
  ! copy c into a useful format
  c_table = 0.d0
  ice = 0
  do l = 0, lmax
    do mp = l, 0, -1
      ice = ice+1
      c_table(-mp,l) = c(ice) ! scaling factor
      c_table( mp,l) = c(ice) ! scaling factor
    enddo ! mp
  enddo ! l
  
  ! todo: do the equivalent for cl
  
  
  ivtot = 0
  do iface = 1, nface
  
    do itet = 1, face(iface)%ntt ! ntt(iface)
      ivtot = ivtot+1
!     rupsq_precomputed(ivtot) = sqrt(rd(ivtot)**2 - face(iface)%r0**2) ! r0(iface)**2) ! obviously |rd| >= |r0| 
    enddo ! itet
    
    if (precompute_dmtl) then
      ! calculate transformation matrices for spherical harmonics
!       call d_real(lmax,alpha(iface),beta(iface),gamma(iface),dmatl_precomputed(:,iface),isumd,lmaxd1)
!       call d_real(lmax, [alpha(iface), beta(iface), gamma(iface)], dmatl_precomputed(:,iface), isumd, lmaxd1)
      call d_real(lmax, face(iface)%euler, dmatl_precomputed(:,iface), isumd, lmaxd1)
    endif ! precompute_dmtl
    
  enddo ! iface
  nonzero(1:mlm) = .false. ! init 

  !===================== split ??? =======================================

  !.......................................................................
  !     loop over radial mesh points
  !.......................................................................
  meshloop: do ir = 1, meshn
    b(1:mlm) = 0.d0
!     ivtot = 0
    
    !.......................................................................
    !     loop over pyramids
    !.......................................................................
py: do iface = 1, nface

      if (xrn(ir) <= face(iface)%r0) then ! r0(iface)) then
        ! the radius is smaller than the distance of the base points from the origin.
        ! so, no contribution to the shapefunctions is expected (except for l=0,m=0 which is treated analytically.
!         ivtot = ivtot + face(iface)%ntt ! ntt(iface) ! fast forward total vertex index
        cycle py ! jump to the head of the loop over pyramids (loop counter iface)
        
      endif ! r <= r0
      
      arg1 = face(iface)%r0/xrn(ir) ! r0(iface)/xrn(ir)
      s = 0.d0 ! init s
      
      !.......................................................................
      !     loop over tetrahedra
      !.......................................................................
      do itet = 1, face(iface)%ntt ! ntt(iface)
!         ivtot = ivtot+1 ! forward total vertex index
!         
!         if (xrn(ir) <= rd(ivtot)) then
!         
!           call pintg(fa(ivtot), fb(ivtot), dlt, s1, lmax, isignu(ivtot), arg1, fd(ivtot), 0)
!           s = s + s1
!           
!         else  ! r <= rd(ivtot)
!         
!           rdown = sqrt(xrn(ir)**2 - face(iface)%r0**2) ! r0(iface)**2)
!           if (precompute_rupsq) then
!             rupsq = rupsq_precomputed(ivtot)
!           else
!             rupsq = sqrt(rd(ivtot)**2 - face(iface)%r0**2)
!           endif
!           rap  = rupsq/rdown
!           arg2 = rupsq/face(iface)%r0 ! /r0(iface)
!           fk = min(fb(ivtot), max(fa(ivtot), fd(ivtot) - acos(rap)))
!           fl = min(fb(ivtot), max(fa(ivtot), fd(ivtot) + acos(rap)))
!           
!           call pintg(fa(ivtot), fk,        dlt, s1, lmax, isignu(ivtot), arg1, fd(ivtot), 0)
!           s = s + s1
!           call pintg(fk,        fl,        dlt, s1, lmax, isignu(ivtot), arg2, fd(ivtot), 1)
!           s = s + s1
!           call pintg(fl,        fb(ivtot), dlt, s1, lmax, isignu(ivtot), arg1, fd(ivtot), 0)
!           s = s + s1
!           
!         endif ! r <= rd(ivtot)

        t1 = face(iface)%ta(itet) ! copy, todo: text-pointer

        if (xrn(ir) <= t1%rd) then
        
          call pintg(t1%fa, t1%fb, dlt, s1, lmax, t1%isignu, arg1, t1%fd, 0)
          s = s + s1
          
        else  ! r <= t1%rd
        
          rdown = sqrt(xrn(ir)**2 - face(iface)%r0**2) ! r0(iface)**2)
          rupsq = sqrt(t1%rd**2 - face(iface)%r0**2)
          rap  = rupsq/rdown
          arg2 = rupsq/face(iface)%r0
          fk = min(t1%fb, max(t1%fa, t1%fd - acos(rap)))
          fl = min(t1%fb, max(t1%fa, t1%fd + acos(rap)))
          
          call pintg(t1%fa, fk, dlt, s1, lmax, t1%isignu, arg1, t1%fd, 0)
          s = s + s1
          call pintg(fk,    fl, dlt, s1, lmax, t1%isignu, arg2, t1%fd, 1)
          s = s + s1
          call pintg(fl, t1%fb, dlt, s1, lmax, t1%isignu, arg1, t1%fd, 0)
          s = s + s1
          
        endif ! r <= t1%rd

      enddo ! itet ! tetraeder loop

      !.......................................................................
      !  integral expansion back-rotation
      !.......................................................................

      if (precompute_dmtl) then
        dmatl = dmatl_precomputed(:,iface)
      else  ! precompute
        ! calculate transformation matrices for spherical harmonics (for each ir again!)
        ! call d_real(lmax, [alpha(iface), beta(iface), gamma(iface)], dmatl, isumd, lmaxd1)
        call d_real(lmax, face(iface)%euler, dmatl, isumd, lmaxd1)
      endif ! precompute
      
!     ib = 0
      ic = 0
!     ice = 0
!     isu = 0

      do l = 0, lmax
      
!         CHECKASSERT(ib == l**2) ! works
!       ib = ib+l+1
        isu0 = ((l-1)*(3+4*l*(l+1)))/3+1
        isu = isu0
        
!       ic0 = ic ! store offset
        
!       ice0 = ice ! store offset
!       CHECKASSERT(ice0 == (l*(l+1))/2) ! works
!         ice0 = (l*(l+1))/2
!         ice = ice0
        
!         smm = 0.d0
!         do m = -l, l
!           smm(m) = 0.d0
!           do k = l, (l+abs(m)+1)/2, -1
!             is = 2*k-l-abs(m)
! !           ic = ic+1
! !           smm(m) = smm(m) + cl(ic)*s(m,is)
!           enddo ! k
!           smm(m) = smm(m)*c(ice0+1)
!         enddo ! m
        
        ! transform to summed coeffcients
        
        do mp = l, 0, -1
          sm(1:2,mp) = 0.d0
          k0 = (l+mp+1)/2
          do k = l, k0, -1
            is = 2*k-l-mp
            ic = ic+1
            sm(1:2,mp) = sm(1:2,mp) + cl(ic)*[s( mp,is), s(-mp,is)]
          enddo ! k
!           ice = ice+1
!           sm(1:2,mp) = sm(1:2,mp)*c(ice) ! scale
          sm(1,mp) = sm(1,mp)*c_table( mp,l) ! scale
          sm(2,mp) = sm(2,mp)*c_table(-mp,l) ! scale
        enddo ! mp
        sm(2,0) = 0.d0 ! corrrect since there is no m=-0
        
        
        ! rotate back in the l-subspace and add to the shape functions
        imax = 1
        do m = 0, l
          do i = 1, imax
            mo = (3-2*i)*m ! mapping i=1 --> mo=m and i=2 --> mo=-m (reconstruction of signed m quantum number)
            ilm = l*l+l+mo+1
            ipmax = 1
            do mp = 0, l
              do ip = 1, ipmax
                isu = isu+1
                b(ilm) = b(ilm) + sm(ip,mp)*dmatl(isu)
              enddo ! ip
              ipmax = 2
            enddo ! mp
          enddo ! i
          imax = 2
        enddo ! m
        
      enddo ! l

    enddo py
    !.......................................................................
    !     defines and saves shapefunctions
    !.......................................................................
    b(1:mlm) = -b(1:mlm)/fpisq ! scale
    b(1) = fpisq + b(1) ! add constant sqrt(4*pi) in the l=0,m=0 channel
    
    nonzero = nonzero .or. (abs(b(:)) > 1d-6)
    thetas_s(ir,1:mlm) = b(1:mlm)
    
  enddo meshloop
  
  if (precompute_dmtl) deallocate(dmatl_precomputed)
  if (precompute_rupsq) deallocate(rupsq_precomputed)

  !now rearrange thetas_s array that it contains only non-zero shapefunctions
  !this is done "in-place"

  lmifun_s = 0 ! lm-index of this shape function

  nfun = 0 ! number of non-zero shape function
  do ilm = 1, mlm ! this loop must run ordered in serial
    if (nonzero(ilm)) then

      nfun = nfun + 1
      lmifun_s(nfun) = ilm

      if (nfun < ilm) thetas_s(1:meshn,nfun) = thetas_s(1:meshn,ilm) ! compress the storage by overwriting functions that are zero or close to zero

    else
      thetas_s(1:meshn,ilm) = 0.d0 ! clear to exact zeros
    endif ! nonzero
  enddo ! ilm
  
endsubroutine shapeIntegration

endmodule ShapeIntegration_mod
