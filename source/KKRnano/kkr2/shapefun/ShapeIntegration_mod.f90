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

  integer :: ivtot, iface, ntet, itet, i, m, ic, ice, ib, isu, l, mp, k0, k, is, imax, mo, ip, ipmax, jbm
  integer :: ibm, ibmax, ir

  real*8, allocatable :: rupsq(:) ! dim nvtotd

  real*8 :: s(-lmaxd1:lmaxd1,0:lmaxd1), s1(-lmaxd1:lmaxd1,0:lmaxd1) ! this storage format usese only (lmaxd1+1)**2 of (lmaxd1+1)*(2*lmaxd1+1) elements so roughly 55%
  real*8 :: sm(2,0:lmaxd1) ! this could be reshaped into sm(-lmaxd1:lmaxd1)
  real*8 :: fk
  real*8 :: fl
  real*8 :: fpisq

  ! local automatic arrays
  integer :: lofm(ibmaxd)
  integer :: mofm(ibmaxd)
  logical(kind=1) :: nonzero(ibmaxd)
  real*8 :: b(ibmaxd)
  
  ! nonzero index array, value is 1 for lm for which shapefunction is non-zero
  ! otherwise the value is 0

  allocate(rupsq(size(rd))) ! dim: nvtotd

  fpisq = sqrt(4.d0*pi)

  ibmax = (lmax+1)**2

  thetas_s = 0.d0

  !.......................................................................
  !     e x p a n s i o n    c o e f f i c i e n t s
  !.......................................................................
  call ccoef(lmax,cl,c)
  ivtot = 0
  do iface = 1, nface
    ntet = ntt(iface)
    do itet = 1, ntet
      ivtot = ivtot+1
      rupsq(ivtot) = sqrt(rd(ivtot)**2 - r0(iface)**2) ! obviously |rd| >= |r0| 
    enddo ! itet
  enddo ! iface
  nonzero(1:ibmax) = .false. ! init 

  !===================== split ??? =======================================

  !.......................................................................
  !     l o o p    o v e r    r a d i a l    m e s h    p o i n t s
  !.......................................................................
  meshloop: do ir = 1, meshn
    r = xrn(ir) ! radius
    b(1:ibmax) = 0.d0
    ivtot = 0
    
    !.......................................................................
    !     l o o p    o v e r    p y r a m i d s
    !.......................................................................
py: do iface = 1, nface
      ntet = ntt(iface)

      if (r <= r0(iface)) then
        ivtot = ivtot+ntet
        cycle py
        
      endif ! r <= r0
      
      arg1 = r0(iface)/r
      rdown = sqrt(r**2 - r0(iface)**2)
!       do i = 0, lmax
!         s(0,i) = 0.d0
!       enddo ! i
!       do m = 1, lmax
!         do i = 0, lmax-m
!           s(-m,i) = 0.d0
!           s( m,i) = 0.d0
!         enddo ! i
!       enddo ! m
      s = 0.d0 ! init s
      
      !.......................................................................
      !     l o o p     o v e r     t e t r a h e d r a
      !.......................................................................
      do itet = 1, ntet
        ivtot = ivtot+1
        
        if (r <= rd(ivtot)) then
        
          call pintg(fa(ivtot),fb(ivtot),dlt,s1,lmax,isignu(ivtot), arg1,fd(ivtot),0)
!           do i = 0, lmax
!             s(0,i) = s(0,i) + s1(0,i)
!           enddo ! i
!           do m = 1, lmax
!             do i = 0, lmax-m
!               s(-m,i) = s(-m,i) + s1(-m,i)
!               s( m,i) = s( m,i) + s1( m,i)
!             enddo ! i
!           enddo ! m
          s = s + s1
          
        else  ! r <= rd(ivtot)
        
          rap  = rupsq(ivtot)/rdown
          arg2 = rupsq(ivtot)/r0(iface)
          fk = fd(ivtot) - acos(rap)
          fl = fd(ivtot) + acos(rap)

          fk = max(fa(ivtot), fk)
          fl = max(fa(ivtot), fl)
          fk = min(fb(ivtot), fk)
          fl = min(fb(ivtot), fl)
          call pintg(fa(ivtot),fk,dlt,s1,lmax,isignu(ivtot), arg1,fd(ivtot),0)
          s = s + s1
          call pintg(fk       ,fl,dlt,s1,lmax,isignu(ivtot), arg2,fd(ivtot),1)
          s = s + s1
          call pintg(fl,fb(ivtot),dlt,s1,lmax,isignu(ivtot), arg1,fd(ivtot),0)
          s = s + s1
          
!           do i = 0, lmax
!             s(0,i) = s(0,i) + s1(0,i) + s2(0,i) + s3(0,i)
!           enddo ! i
!           do m = 1, lmax
!             do i = 0, lmax-m
!               s(-m,i) = s(-m,i) + s1(-m,i) + s2(-m,i) + s3(-m,i)
!               s( m,i) = s( m,i) + s1( m,i) + s2( m,i) + s3( m,i)
!             enddo ! i
!           enddo ! m
          
        endif ! r <= rd(ivtot)

      enddo ! itet ! tetraeder loop

      !.......................................................................
      !  i n t e g r a l   e x p a n s i o n        b a c k - r o t a t i o n
      !.......................................................................

      ib = 0
      ic = 0
      ice = 0
      ! calculate transformation matrices for spherical harmonics
      call d_real(lmax,alpha(iface),beta(iface),gamma(iface),dmatl,isumd,lmaxd1)

      isu = 0
      do l = 0, lmax
      
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
        
        sm(1,0) = 0.d0
        ice = ice+1
        k0 = (l+1)/2
        do k = l, k0, -1
          is = 2*k-l
          ic = ic+1
          sm(1,0) = sm(1,0) + cl(ic)*s(0,is)
        enddo ! k
        sm(1,0) = sm(1,0)*c(ice) ! scale
        
        imax = 1
        do m = 0, l
          do i = 1, imax
            mo = (3-2*i)*m
            ibm = ib+mo
            lofm(ibm) = l
            mofm(ibm) = mo
            ipmax = 1
            do mp = 0, l
              do ip = 1, ipmax
                isu = isu+1
                b(ibm) = b(ibm) + sm(ip,mp)*dmatl(isu)
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
    b(1:ibmax) = -b(1:ibmax)/fpisq ! scale
    b(1) = fpisq + b(1) ! add constant sqrt(4*pi) in the l=0,m=0 channel
    
    nonzero = nonzero .or. (abs(b(:)) > 1d-6)
    thetas_s(ir,1:ibmax) = b(1:ibmax)
    
  enddo meshloop

  !now rearrange thetas_s array that it contains only non-zero shapefunctions
  !this is done "in-place"

  lmifun_s = 0 ! lm-index of this shape function

  jbm = 0
  do ibm = 1, ibmax
    if (nonzero(ibm)) then

      jbm = jbm + 1
      lmifun_s(jbm) = lofm(ibm)*lofm(ibm) + lofm(ibm) + mofm(ibm) + 1

      if (jbm < ibm) then
        thetas_s(1:meshn,jbm) = thetas_s(1:meshn,ibm)
      endif ! jbm < ibm

    else
      thetas_s(1:meshn,ibm) = 0.d0
    endif ! nonzero
  enddo ! ibm

  nfun = count(nonzero(1:ibmax)) ! count non-zero shape functions

endsubroutine shapeIntegration

endmodule ShapeIntegration_mod
