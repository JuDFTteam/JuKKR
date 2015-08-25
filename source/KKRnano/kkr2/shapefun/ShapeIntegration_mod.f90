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

subroutine shapeintegration(lmax, nface, meshn, xrn, dlt, thetas_s, lmifun_s, nfun, meshnd, ibmaxd)

  use shape_constants_mod, only: pi, lmaxd1, isumd, icd, iced
  use tetrahedra_common, only: rd, ntt, r0, fa, fb, fd, isignu
  use angles_common, only: alpha, beta, gamma
  use shapeintegrationhelpers_mod, only: pintg, ccoef, d_real

  integer :: lmax
  integer :: nface
  integer :: meshn
  real*8 ::  xrn(meshnd)
  real*8 :: dlt
  real*8 ::  thetas_s(meshnd,ibmaxd)
  integer :: lmifun_s(ibmaxd)
  integer, intent(in) :: meshnd
  integer, intent(in) :: ibmaxd
  integer, intent(out) :: nfun

  real*8 ::  dmatl(isumd)

  real*8 ::     cl(icd)
  real*8 ::     c(iced)

  real*8 ::  r
  real*8 ::  rap
  real*8 ::  rdown

  real*8 ::  arg1
  real*8 ::  arg2

  integer :: ivtot
  integer :: iface
  integer :: ntet
  integer :: itet
  integer :: i
  integer :: m
  integer :: ic
  integer :: ice
  integer :: ib
  integer :: isu
  integer :: l
  integer :: mp
  integer :: k0
  integer :: k
  integer :: is
  integer :: imax
  integer :: mo
  integer :: ip
  integer :: ipmax
  integer :: icount

  integer :: ibm
  integer :: ibmax
  integer :: n

  real*8, allocatable ::  rupsq(:) ! dim nvtotd

  real*8 ::     s(-lmaxd1:lmaxd1,0:lmaxd1)
  real*8 ::    s1(-lmaxd1:lmaxd1,0:lmaxd1)
  real*8 ::    s2(-lmaxd1:lmaxd1,0:lmaxd1)
  real*8 ::    s3(-lmaxd1:lmaxd1,0:lmaxd1)
  real*8 ::    sum(0:lmaxd1,2)
  real*8 ::    fk
  real*8 ::    fl
  real*8 ::    fpisq

  ! local automatic arrays
  integer ::   lofm(ibmaxd)
  integer ::   mofm(ibmaxd)
  real*8 ::    b(ibmaxd)
  integer ::   isw(ibmaxd)
  ! isw index array, value is 1 for lm for which shapefunction is non-zero
  ! otherwise the value is 0

  allocate(rupsq(size(rd))) ! dim: nvtotd

  fpisq=dsqrt(4.d0*pi)

  ibmax=(lmax+1)*(lmax+1)

  thetas_s = 0.0d0

  !.......................................................................
  !     e x p a n s i o n    c o e f f i c i e n t s
  !.......................................................................
  call ccoef(lmax,cl,c)
  ivtot=0
  do iface=1,nface
    ntet=ntt(iface)
    do itet=1,ntet
      ivtot=ivtot+1
      rupsq(ivtot)=sqrt((rd(ivtot)-r0(iface))*(rd(ivtot)+r0(iface)))
    enddo
  enddo
  do ibm=1,ibmax
    isw(ibm)=0
  enddo

  !===================== split ??? =======================================

  !.......................................................................
  !     l o o p    o v e r    r a d i a l    m e s h    p o i n t s
  !.......................................................................
  meshloop: do n=1,meshn
    r=xrn(n)
    do ibm=1,ibmax
      b(ibm)=0.d0
    enddo
    ivtot=0
    !.......................................................................
    !     l o o p    o v e r    p y r a m i d s
    !.......................................................................
py: do iface=1,nface
      ntet=ntt(iface)

      if(r > r0(iface))  goto 31
      ivtot=ivtot+ntet
      do i=0,lmax
        s(0,i) =0.d0
      enddo
      do m=1,lmax
        do i=0,lmax-m
          s(-m,i)=0.d0
          s( m,i)=0.d0
        enddo
      enddo
      goto 13
31    continue
      !if(newsch(iface) == 1) goto 35
      !ivtot=ivtot+ntet
      !goto 32
35    arg1=r0(iface)/r
      rdown=sqrt((r-r0(iface))*(r+r0(iface)))
      do i=0,lmax
        s(0,i) =0.d0
      enddo
      do m=1,lmax
        do i=0,lmax-m
          s(-m,i)=0.d0
          s( m,i)=0.d0
        enddo
      enddo
      !.......................................................................
      !     l o o p     o v e r     t e t r a h e d r a
      !.......................................................................
      do itet=1,ntet
        ivtot=ivtot+1
        if(r <= rd(ivtot))      then
          call pintg(fa(ivtot),fb(ivtot),dlt,s1,lmax,isignu(ivtot), &
          arg1,fd(ivtot),0)
          do i=0,lmax
            s(0,i)=s(0,i)+s1(0,i)
          enddo
          do m=1,lmax
            do i=0,lmax-m
              s(-m,i)=s(-m,i)+s1(-m,i)
              s( m,i)=s( m,i)+s1( m,i)
            enddo
          enddo
        else
          rap =rupsq(ivtot)/rdown
          arg2=rupsq(ivtot)/r0(iface)
          fk=fd(ivtot)-acos(rap)
          fl=fd(ivtot)+acos(rap)

          fk=dmax1(fa(ivtot),fk)
          fl=dmax1(fa(ivtot),fl)
          fk=dmin1(fb(ivtot),fk)
          fl=dmin1(fb(ivtot),fl)
          call pintg(fa(ivtot),fk,dlt,s1,lmax,isignu(ivtot), &
          arg1,fd(ivtot),0)
          call pintg(fk       ,fl,dlt,s2,lmax,isignu(ivtot), &
          arg2,fd(ivtot),1)
          call pintg(fl,fb(ivtot),dlt,s3,lmax,isignu(ivtot), &
          arg1,fd(ivtot),0)
          do i=0,lmax
            s(0,i)=s(0,i)+s1(0,i)+s2(0,i)+s3(0,i)
          enddo
          do m=1,lmax
            do i=0,lmax-m
              s(-m,i)=s(-m,i)+s1(-m,i)+s2(-m,i)+s3(-m,i)
              s( m,i)=s( m,i)+s1( m,i)+s2( m,i)+s3( m,i)
            enddo
          enddo
        endif

      enddo  ! tetraeder loop

32    continue

      !.......................................................................
      !     i n t e g r a l   e x p a n s i o n        b a c k - r o t a t i o
      !.......................................................................

      ib=0
      ic=0
      ice=0
      ! calculate transformation matrices for spherical harmonics
      call d_real(lmax,alpha(iface),beta(iface),gamma(iface),dmatl,isumd,lmaxd1)

      isu=0
      do l=0,lmax
        ib=ib+l+1
        do mp=l,1,-1
          sum(mp,1)=0.d0
          sum(mp,2)=0.d0
          ice=ice+1
          k0=(l+mp+1)/2
          do k=l,k0,-1
            is=2*k-l-mp
            ic=ic+1
            sum(mp,2)=sum(mp,2)+cl(ic)*s(-mp,is)
            sum(mp,1)=sum(mp,1)+cl(ic)*s( mp,is)
          enddo
          sum(mp,2)=sum(mp,2)*c(ice)
          sum(mp,1)=sum(mp,1)*c(ice)
        enddo
        sum(0,1)=0.d0
        ice=ice+1
        k0=(l+1)/2
        do k=l,k0,-1
          is=2*k-l
          ic=ic+1
          sum(0,1)=sum(0,1)+cl(ic)*s(0,is)
        enddo
        sum(0,1)=sum(0,1)*c(ice)
        imax=1
        m=0
  8     continue
        do i=1,imax
          mo=(3-2*i)*m
          ibm=ib+mo
          lofm(ibm)=l
          mofm(ibm)=mo
          ipmax=1
          mp=0
    16    continue
          do ip=1,ipmax
            isu=isu+1
            b(ibm)=b(ibm)+sum(mp,ip)*dmatl(isu)
          enddo
          ipmax=2
          mp=mp+1
          if(mp <= l) goto 16
        enddo
        imax=2
        m=m+1
        if(m <= l) goto 8
        ib=ib+l
      enddo ! loop over l

  13 continue
    enddo py
    !.......................................................................
    !     d e f i n e s   a n d    s a v e s   s h a p e    f u n c t i o n ???
    !.......................................................................
    b(1)=fpisq-b(1)/fpisq
    do ibm=2,ibmax
      b(ibm)=-b(ibm)/fpisq
    enddo
    do ibm=1,ibmax
      !     write(6,*) ibm,b(ibm)
      if(abs(b(ibm)) > 1.d-6) isw(ibm)=1
      !irec=(ibm-1)*meshn+n
      !write(11,rec=irec) b(ibm)
      thetas_s(n, ibm) = b(ibm)
    enddo
  enddo meshloop

  !now rearrange thetas_s array that it contains only non-zero shapefunctions
  !this is done "in-place"

  lmifun_s = 0

  icount = 1
  do ibm=1,ibmax
    if (isw(ibm) == 1) then

      lmifun_s(icount) = lofm(ibm)*lofm(ibm)+lofm(ibm)+mofm(ibm)+1

      if (icount .ne. ibm) then
        do n = 1, meshn
          thetas_s(n, icount) = thetas_s(n, ibm)
        enddo
      endif
      icount = icount + 1

    else

      do n = 1, meshn
        thetas_s(n, ibm) = 0.0d0
      enddo
    endif
  enddo

  ! count non-zero shape functions
  nfun=0
  do ibm=1,ibmax
    if(isw(ibm) == 1)  nfun=nfun+1
  enddo

endsubroutine

endmodule ShapeIntegration_mod
