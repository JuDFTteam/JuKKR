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
!> needs information from module/common block replacement "tetrahedra_common"
!> needs information from module/common block replacement "angles_common"
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
  use tetrahedra_common, only: rd, fa, fb, fd, isignu ! vertex properties
  use tetrahedra_common, only: ntt, r0 ! face properties
  use angles_common, only: alpha, beta, gamma ! face properties
  use shapeintegrationhelpers_mod, only: pintg, ccoef, d_real

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
  real*8 :: c(iced)

  real*8 :: r
  real*8 :: rap
  real*8 :: rdown

  real*8 :: arg1
  real*8 :: arg2

  integer :: ivtot, iface, itet, i, m, ic, ice, ib, isu, l, mp, k0, k, is, imax, mo, ip, ipmax, jlm
  integer :: ilm, mlm, ir

  real*8, allocatable :: rupsq_precomputed(:) ! dim nvtotd
  real*8, allocatable :: dmatl_precomputed(:,:)

  real*8 :: s(-lmaxd1:lmaxd1,0:lmaxd1), s1(-lmaxd1:lmaxd1,0:lmaxd1) ! this storage format usese only (lmaxd1+1)**2 of (lmaxd1+1)*(2*lmaxd1+1) elements so roughly 55%
  real*8 :: sm(2,0:lmaxd1) ! this could be reshaped into sm(-lmaxd1:lmaxd1)
  real*8 :: fk
  real*8 :: fl
  real*8 :: fpisq

  ! local automatic arrays
! integer :: lofm(ibmaxd), mofm(ibmaxd)
  logical(kind=1) :: nonzero(ibmaxd)
  real*8 :: b(ibmaxd)
  logical, parameter :: precompute_dmtl = .true. ! needs some memory but reduces the number of calls to d_real from meshn*nface to nface
  
  ! nonzero index array, value is 1 for lm for which shapefunction is non-zero
  ! otherwise the value is 0

  
  allocate(rupsq_precomputed(size(rd))) ! dim: nvtotd, sum(ntt) should suffice, or?
 
  if (precompute_dmtl) allocate(dmatl_precomputed(isumd,nface))
  
!   ilm = 0
!   do l = 0, lmax
!     do m = -l, l
!       ilm = ilm+1
!       mofm(ilm) = m
!       lofm(ilm) = l
!     enddo ! m
!   enddo ! l
  
  fpisq = sqrt(4.d0*pi)

  mlm = (lmax+1)**2

  thetas_s = 0.d0

  !.......................................................................
  !     e x p a n s i o n    c o e f f i c i e n t s
  !.......................................................................
  call ccoef(lmax,cl,c) ! sets the coefficients in m-order 0,1,-1,...,l,-l
  ivtot = 0
  do iface = 1, nface
  
    do itet = 1, ntt(iface)
      ivtot = ivtot+1
      rupsq_precomputed(ivtot) = sqrt(rd(ivtot)**2 - r0(iface)**2) ! obviously |rd| >= |r0| 
    enddo ! itet
    
    if (precompute_dmtl) & ! calculate transformation matrices for spherical harmonics
      call d_real(lmax,alpha(iface),beta(iface),gamma(iface),dmatl_precomputed(:,iface),isumd,lmaxd1)
    
  enddo ! iface
  nonzero(1:mlm) = .false. ! init 

  !===================== split ??? =======================================

  !.......................................................................
  !     l o o p    o v e r    r a d i a l    m e s h    p o i n t s
  !.......................................................................
  meshloop: do ir = 1, meshn
    r = xrn(ir) ! radius
    b(1:mlm) = 0.d0
    ivtot = 0
    
    !.......................................................................
    !     l o o p    o v e r    p y r a m i d s
    !.......................................................................
py: do iface = 1, nface

      if (r <= r0(iface)) then
        ! the radius is smaller than the distance of the base points from the origin.
        ! so, no contribution to the shapefunctions is expected (except for l=0,m=0 which is treated analytically.
        ivtot = ivtot + ntt(iface) ! fast forward total vertex index
        cycle py ! jump to the head of the loop over pyramids (loop counter iface)
        
      endif ! r <= r0
      
      arg1 = r0(iface)/r
      s = 0.d0 ! init s
      
      !.......................................................................
      !     l o o p     o v e r     t e t r a h e d r a
      !.......................................................................
      do itet = 1, ntt(iface)
        ivtot = ivtot+1 ! forward total vertex index
        
        if (r <= rd(ivtot)) then
        
          call pintg(fa(ivtot),fb(ivtot),dlt,s1,lmax,isignu(ivtot), arg1,fd(ivtot),0)
          s = s + s1
          
        else  ! r <= rd(ivtot)
        
          rdown = sqrt(r**2 - r0(iface)**2)
          rap  = rupsq_precomputed(ivtot)/rdown
          arg2 = rupsq_precomputed(ivtot)/r0(iface)
          fk = min(fb(ivtot), max(fa(ivtot), fd(ivtot) - acos(rap)))
          fl = min(fb(ivtot), max(fa(ivtot), fd(ivtot) + acos(rap)))
          
          call pintg(fa(ivtot),fk,dlt,s1,lmax,isignu(ivtot), arg1,fd(ivtot),0)
          s = s + s1
          call pintg(fk       ,fl,dlt,s1,lmax,isignu(ivtot), arg2,fd(ivtot),1)
          s = s + s1
          call pintg(fl,fb(ivtot),dlt,s1,lmax,isignu(ivtot), arg1,fd(ivtot),0)
          s = s + s1
          
        endif ! r <= rd(ivtot)

      enddo ! itet ! tetraeder loop

      !.......................................................................
      !  i n t e g r a l   e x p a n s i o n        b a c k - r o t a t i o n
      !.......................................................................

      if (precompute_dmtl) then
        dmatl = dmatl_precomputed(:,iface)
      else  ! precompute
        ! calculate transformation matrices for spherical harmonics
        call d_real(lmax,alpha(iface),beta(iface),gamma(iface),dmatl,isumd,lmaxd1)
      endif ! precompute
      
      ib = 0
      ic = 0
      ice = 0

      isu = 0
      do l = 0, lmax
      
!         CHECKASSERT(ib == l**2)
        ib = ib+l+1
        
        do mp = l, 1, -1
          sm(1:2,mp) = 0.d0
          ice = ice+1
          k0 = (l+mp+1)/2
          do k = l, k0, -1
            is = 2*k-l-mp
            ic = ic+1
            sm(1:2,mp) = sm(1:2,mp) + cl(ic)*[s( mp,is), s(-mp,is)]
          enddo ! k
          sm(1:2,mp) = sm(1:2,mp)*c(ice) ! scale
        enddo ! mp
        
        do mp = 0, 0
          sm(1:1,mp) = 0.d0
          ice = ice+1
          k0 = (l+mp+1)/2
          do k = l, k0, -1
            is = 2*k-l-mp
            ic = ic+1
            sm(1:1,mp) = sm(1:1,mp) + cl(ic)*[s( mp,is)]
          enddo ! k
          sm(1:1,mp) = sm(1:1,mp)*c(ice) ! scale
        enddo ! mp
        
        imax = 1
        do m = 0, l
          do i = 1, imax
            mo = (3-2*i)*m ! mapping i=1 --> mo=m and i=2 --> mo=-m
            ilm = ib+mo ! assume that ib == l*l+l+1 ==> ilm = l*l+l+mo+1
!             CHECKASSERT(ilm == l*l+l+mo+1)
!             lofm(ilm) = l
!             mofm(ilm) = mo
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
        
        ib = ib+l
      enddo ! l

    enddo py
    !.......................................................................
    !     d e f i n e s   a n d    s a v e s   s h a p e    f u n c t i o n s
    !.......................................................................
    b(1:mlm) = -b(1:mlm)/fpisq ! scale
    b(1) = fpisq + b(1) ! add constant sqrt(4*pi) in the l=0,m=0 channel
    
    nonzero = nonzero .or. (abs(b(:)) > 1d-6)
    thetas_s(ir,1:mlm) = b(1:mlm)
    
  enddo meshloop
  
  if (precompute_dmtl) deallocate(dmatl_precomputed)

  !now rearrange thetas_s array that it contains only non-zero shapefunctions
  !this is done "in-place"

  lmifun_s = 0 ! lm-index of this shape function

  jlm = 0
  do ilm = 1, mlm
    if (nonzero(ilm)) then

      jlm = jlm + 1
      lmifun_s(jlm) = ilm ! lofm(ilm)*lofm(ilm) + lofm(ilm) + mofm(ilm) + 1

      if (jlm < ilm) then
        thetas_s(1:meshn,jlm) = thetas_s(1:meshn,ilm)
      endif ! jlm < ilm

    else
      thetas_s(1:meshn,ilm) = 0.d0
    endif ! nonzero
  enddo ! ilm

  nfun = count(nonzero(1:mlm)) ! count non-zero shape functions

endsubroutine shapeIntegration

endmodule ShapeIntegration_mod
