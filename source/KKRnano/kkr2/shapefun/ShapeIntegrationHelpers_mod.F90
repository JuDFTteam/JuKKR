!>    Auxillary module needed for shape function calculation.

      module ShapeIntegrationHelpers_mod
      implicit none
      private
      public :: pintg, ccoef, d_real

      contains

!-----------------------------------------------------------------------
!>    this routine  accomplishes the fi-integration of real spherical
!>    harmonics by the repeated simpson's method , or analytically accor
!>    ding to the value of itype. the obtained results have to be multi-
!>    plied by the appropriate expansion coefficients.
!-----------------------------------------------------------------------
      subroutine pintg(x1, x2, dlt, s, lmax, isi, arg, fd, itype)
      use shape_constants_mod, only: lmaxd1, ndim
!     include 'inc.geometry':   integer, parameter (lmaxd=25,ndim=1000)

      real*8, intent(in) :: x1, x2, dlt, arg, fd
      integer, intent(in) :: lmax
      integer(kind=1), intent(in) :: isi ! sign
      integer, intent(in) :: itype ! itype in [0, 1]
      real*8, intent(out) :: s(-lmaxd1:lmaxd1,0:lmaxd1) ! uses more memory than needed, about 55% used
      
      integer :: i, m, n, k
      real*8 :: x, theta, w
      real*8 :: xx(ndim), ww(ndim)

!-----------------------------------------------------------------------
      if (lmax > lmaxd1) then
        write(*,fmt="(3x,'from pintg: lmax=',i4,' greater than dimensioned',i4)")lmax,lmaxd1
        stop
      endif
    
      s = 0.d0
      
      if (itype == 0) then
        theta = acos(arg)
        call recur0(lmax,x1,theta,-dble(isi),s)
        call recur0(lmax,x2,theta, dble(isi),s)
      else
        n = (x2 - x1)/dlt + 3
        if (n > ndim) stop 'increase ndim'
        call gauleg(x1, x2, xx, ww, n)
        do k = 1, n
          x = xx(k)
          w = dble(isi)*ww(k)
          theta = atan(arg/cos(x - fd))
          call recur(lmax, x, theta, w, s)
        enddo ! k
      endif ! itype == 0
      
      endsubroutine pintg

!======================================================================

!-----------------------------------------------------------------------
!>    this routine is used to perform the fi-integration of real sphe-
!>    rical harmonics .the theta-integration is performed analytically
!>    using recurrence relations.
!-----------------------------------------------------------------------
      subroutine recur(lmax,x,theta,fac,s)
      use shape_constants_mod, only: lmaxd1
      integer, intent(in) :: lmax
      real*8, intent(in) :: x, theta, fac
      real*8, intent(inout) :: s(-lmaxd1:lmaxd1,0:lmaxd1)

      integer :: m,i
      real*8 :: ol0,ol,el0,el,c1,c2,ss,cc
      real*8 :: c01(lmaxd1),c02(lmaxd1),ssa(lmaxd1+2),cca(lmaxd1+2)
!-----------------------------------------------------------------------
      ss=sin(theta)
      cc=cos(theta)
      do 13 i=1,lmax
      c01(i)=fac*sin(dble(i)*x)
   13 c02(i)=fac*cos(dble(i)*x)
      do 11 i=1,lmax+2
      ssa(i)=ss**i
      cca(i)=cc**i
   11 continue
      ol0=(theta-ss*cc)/2.d0
      el0=0d0
      do 1 m=1,lmax,2
      ol=ol0
      el=el0
      c1=c01(m)
      c2=c02(m)
      do 2 i=0,lmax-m,2
      s(-m,i)=s(-m,i)+c1*(ol-el)
      s( m,i)=s( m,i)+c2*(ol-el)
      el= dble(i+1)*el/dble(i+m+3)
    2 ol=(dble(i+1)*ol+(ssa(m+2))*(cca(i+1)))/dble(i+m+3)
      el0= dble(m+2)*el0/dble(m+3)
      ol0=(dble(m+2)*ol0-(ssa(m+2))*cc)/dble(m+3)
    1 continue
      ol0=ssa(3)/3.d0
      el0=0d0
      do 3 m=1,lmax,2
      ol=ol0
      el=el0
      c1=c01(m)
      c2=c02(m)
      do 4 i=1,lmax-m,2
      s(-m,i)=s(-m,i)+c1*(ol-el)
      s( m,i)=s( m,i)+c2*(ol-el)
      ol=(dble(i+1)*ol+ssa(m+2)*cca(i+1))/dble(i+m+3)
    4 el= dble(i+1)*el/dble(i+m+3)
      ol0=(dble(m+2)*ol0-(ssa(m+2))*cca(2))/dble(m+4)
      el0= dble(m+2)*el0/dble(m+4)
    3 continue
      ol0=-cc
      el0=-1d0
      ol=ol0
      el=el0
      c2=fac
      do 9 i=0,lmax,2
      s(0,i)=s(0,i)+c2*(ol-el)
      ol=(dble(i+1)*ol+ssa(2)*cca(i+1))/dble(i+3)
    9 el= dble(i+1)*el/dble(i+3)
      ol0=(2.d0*ol0-ssa(2)*cc)/3.d0
      el0= 2.d0*el0/3.d0
      do 5 m=2,lmax,2
      ol=ol0
      el=el0
      c1=c01(m)
      c2=c02(m)
      do 6 i=0,lmax-m,2
      s(-m,i)=s(-m,i)+c1*(ol-el)
      s( m,i)=s( m,i)+c2*(ol-el)
      ol=(dble(i+1)*ol+ssa(m+2)*cca(i+1))/dble(i+m+3)
    6 el= dble(i+1)*el/dble(i+m+3)
      ol0=(dble(m+2)*ol0-ssa(m+2)*cc)/dble(m+3)
      el0= dble(m+2)*el0/dble(m+3)
    5 continue
      ol0=-cca(2)/2d0
      el0=-0.5d0
      ol=ol0
      el=el0
      c2=fac
      do 10 i=1,lmax,2
      s(0,i)=s(0,i)+c2*(ol-el)
      ol=(dble(i+1)*ol+ss*ss*cca(i+1))/dble(i+3)
   10 el= dble(i+1)*el/dble(i+3)
      ol0=(2.d0*ol0-ssa(2)*cca(2))/4.d0
      el0= 2.d0*el0/4.d0
      do 7 m=2,lmax,2
      ol=ol0
      el=el0
      c1=c01(m)
      c2=c02(m)
      do 8 i=1,lmax-m,2
      s(-m,i)=s(-m,i)+c1*(ol-el)
      s( m,i)=s( m,i)+c2*(ol-el)
      ol=(dble(i+1)*ol+ssa(m+2)*cca(i+1))/dble(i+m+3)
    8 el= dble(i+1)*el/dble(i+m+3)
      ol0=(dble(m+2)*ol0-ssa(m+2)*cca(2))/dble(m+4)
      el0= dble(m+2)*el0/dble(m+4)
    7 continue
      endsubroutine recur

!======================================================================

!-----------------------------------------------------------------------
!>    this routine is used to perform  the  fi - integration  of real sp
!>    rical harmonics analytically.  the  theta-integration  is   perfor
!>    also analytically using recurrence relations.(theta is fi-independ
!-----------------------------------------------------------------------
      subroutine recur0(lmax,x,theta,fac,s)
      use shape_constants_mod, only: lmaxd1
      integer, intent(in) :: lmax
      real*8, intent(in) :: x, theta, fac
      real*8, intent(inout) :: s(-lmaxd1:lmaxd1,0:lmaxd1)

      integer :: m, i
      real*8 :: ol0,ol,el0,el,c1,c2,ss,cc
!-----------------------------------------------------------------------
      ss=sin(theta)
      cc=cos(theta)
      ol0=(theta-ss*cc)/2d0
      el0=0d0
      do 1 m=1,lmax,2
      ol=ol0
      el=el0
      c1=-fac*cos(dble(m)*x)/dble(m)
      c2= fac*sin(dble(m)*x)/dble(m)
      do 2 i=0,lmax-m,2
      s(-m,i)=s(-m,i)+c1*(ol-el)
      s( m,i)=s( m,i)+c2*(ol-el)
      el= dble(i+1)*el/dble(i+m+3)
    2 ol=(dble(i+1)*ol+(ss**(m+2))*(cc**(i+1)))/dble(i+m+3)
      el0= dble(m+2)*el0/dble(m+3)
      ol0=(dble(m+2)*ol0-(ss**(m+2))*cc)/dble(m+3)
    1 continue
      ol0=ss*ss*ss/3d0
      el0=0d0
      do 3 m=1,lmax,2
      ol=ol0
      el=el0
      c1=-fac*cos(dble(m)*x)/dble(m)
      c2= fac*sin(dble(m)*x)/dble(m)
      do 4 i=1,lmax-m,2
      s(-m,i)=s(-m,i)+c1*(ol-el)
      s( m,i)=s( m,i)+c2*(ol-el)
      ol=(dble(i+1)*ol+(ss**(m+2))*(cc**(i+1)))/dble(i+m+3)
    4 el= dble(i+1)*el/dble(i+m+3)
      ol0=(dble(m+2)*ol0-(ss**(m+2))*cc*cc)/dble(m+4)
      el0= dble(m+2)*el0/dble(m+4)
    3 continue
      ol0=-cc
      el0=-1d0
      ol=ol0
      el=el0
      c2=fac*x
      do 9 i=0,lmax,2
      s(0,i)=s(0,i)+c2*(ol-el)
      ol=(dble(i+1)*ol+ss*ss*(cc**(i+1)))/dble(i+3)
    9 el= dble(i+1)*el/dble(i+3)
      ol0=(2.d0*ol0-ss*ss*cc)/3.d0
      el0= 2.d0*el0/3.d0
      do 5 m=2,lmax,2
      ol=ol0
      el=el0
      c1=-fac*cos(dble(m)*x)/dble(m)
      c2 =fac*sin(dble(m)*x)/dble(m)
      do 6 i=0,lmax-m,2
      s(-m,i)=s(-m,i)+c1*(ol-el)
      s( m,i)=s( m,i)+c2*(ol-el)
      ol=(dble(i+1)*ol+(ss**(m+2))*(cc**(i+1)))/dble(i+m+3)
    6 el= dble(i+1)*el/dble(i+m+3)
      ol0=(dble(m+2)*ol0-(ss**(m+2))*cc)/dble(m+3)
      el0= dble(m+2)*el0/dble(m+3)
    5 continue
      ol0=-cc*cc/2d0
      el0=-0.5d0
      ol=ol0
      el=el0
      c2=fac*x
      do 10 i=1,lmax,2
      s(0,i)=s(0,i)+c2*(ol-el)
      ol=(dble(i+1)*ol+ss*ss*(cc**(i+1)))/dble(i+3)
   10 el= dble(i+1)*el/dble(i+3)
      ol0=(2.d0*ol0-ss*ss*cc*cc)/4.d0
      el0= 2.d0*el0/4.d0
      do 7 m=2,lmax,2
      ol=ol0
      el=el0
      c1=-fac*cos(dble(m)*x)/dble(m)
      c2= fac*sin(dble(m)*x)/dble(m)
      do 8 i=1,lmax-m,2
      s(-m,i)=s(-m,i)+c1*(ol-el)
      s( m,i)=s( m,i)+c2*(ol-el)
      ol=(dble(i+1)*ol+(ss**(m+2))*(cc**(i+1)))/dble(i+m+3)
    8 el= dble(i+1)*el/dble(i+m+3)
      ol0=(dble(m+2)*ol0-(ss**(m+2))*cc*cc)/dble(m+4)
      el0= dble(m+2)*el0/dble(m+4)
    7 continue
      endsubroutine recur0

!======================================================================

!     ----------------------------------------------------------------
!>    ginen the lower and upper limits of integration  x1 and x2, and
!>    given n, this subroutine returns the  arrays x(1:n) and  w(1:n)
!>    of length n, containing the abscissas and weights of the  gauss
!>    legendre n-point quadrature formula (numerical recipes,2nd ed.).
!     ----------------------------------------------------------------
  subroutine gauleg(x1,x2,x,w,n)
    use shape_constants_mod, only: pi
    
    integer, intent(in) :: n
    real*8, intent(in) :: x1, x2
    real*8, intent(out) :: x(1:), w(1:)

!     ----------------------------------------------------------------
    integer :: i, j, m
    real*8 :: p1,p2,p3,pp,xl,xm,z,z1

    m = (n+1)/2
    xm = 0.5d0*(x2 + x1)
    xl = 0.5d0*(x2 - x1)
    do i = 1, m
      z = cos(pi*(i - .25d0)/(n + .5d0))
      z1 = z + 99.
      do while (abs(z - z1) > 3.d-14)
        p1 = 1.d0
        p2 = 0.d0
        do j = 1, n
          p3 = p2
          p2 = p1
          p1 = ((2*j-1)*z*p2 - (j-1)*p3)/dble(j)
        enddo ! n
        pp = n*(z*p1 - p2)/(z*z - 1.d0)
        z1 = z
        z = z1 - p1/pp
      enddo ! while (abs(z - z1) > 3.d-14)
      x(i) = xm - xl*z
      x(n+1-i) = xm + xl*z
      w(i) = 2.d0*xl/((1.d0 - z*z)*pp*pp)
      w(n+1-i) = w(i)
    enddo ! i
  endsubroutine gauleg

!-----------------------------------------------------------------------
!>    this routine calculates the coefficients of a polynomial expansion
!>    of renormalized legendre functions in powers of cosines.
!>    the possibility of overflow (high lmax) is avoided by using facto-
!>    rized forms for the numbers.
!-----------------------------------------------------------------------
      subroutine ccoef(lmax, cl, coe)
!     changed: get constants from module shape_constants_mod instead of inc.geometry
      use shape_constants_mod, only: icd, iced, lmaxd1
      integer, intent(in) :: lmax
      real*8, intent(out) :: cl(icd) ! warning: old m-ordering: ???
      real*8, intent(out) :: coe(iced) ! warning: old m-ordering: ((m, m=l...0), l=0,lmax)
     
      integer, parameter :: ifmx=25, lma2d=lmaxd1/2+1
      integer :: icmax,l,ice,ic,i,m,k,k0,isi,ire,ir,ic1,ic2, la, lb, ieupsq, ieint
      real*8 :: up,down,upsq
      integer :: ie(ifmx,lma2d), ied(ifmx)
      integer, dimension(ifmx) :: l1st, l2st, l1, l2, jm0, iea, ieb, il2p
      integer, parameter :: primes(ifmx) = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97]
!-----------------------------------------------------------------------
      icmax = sum([( (lmax+1-l)*(l/2+1), l=0,lmax)])
      
      if (lmax > lmaxd1 .or. icmax > icd) then
        write(*,fmt="(13x,'from ccoef: inconsistency data-dimension'/14x,'lmax:',2i5/13x,'icmax:',2i5)") lmax,lmaxd1,icmax,icd
        stop 'lmax > lmaxd1 .or. icmax > icd'
      endif ! lmax > lmaxd1 .or. icmax > icd
!     if (test('shape   ')) write(6,fmt="(13x,'there are',i5,'  coefficients'/)") icmax
      ice = 0
      ic = 1
      l = 0
      l1st(:) = 0
      l2st(:) = 0
      
    1 continue
      call factorize(2*l+1, primes, il2p)
      il2p(1) = il2p(1)+1 ! increase the exponent of the lowest prime factor (which is 2)
      
      m = l
      l1(:) = l1st(:)
      l2(:) = l2st(:)
      jm0(:) = 0
      
    2 continue
      ice = ice+1
      isi = 1
!  this is changed: old: isi = 1-2*mod(m,2)
      k0 = (l+m+1)/2
      k = l
      ire = 1
      ic1 = ic
      ie(:,ire) = jm0(:)
      
    3 continue
      if ((k-1) < k0) goto 30
      ire = ire+1
      ic = ic+1
      la = (2*k-l-m)*(2*k-l-m-1)
      lb = 2*(2*k-1)*(l-k+1)
      call factorize(la, primes, iea)
      call factorize(lb, primes, ieb)
      ie(:,ire) = ie(:,ire-1) + iea(:) - ieb(:)
      k = k-1
      goto 3
   30 continue
   
      ic2 = ic
      do i = 1, ifmx
        ied(i) = ie(i,1)
        do ir = 2, ire
          if(ie(i,ir) < ied(i)) ied(i) = ie(i,ir)
        enddo ! ir
        do ir = 1, ire
          ie(i,ir) = ie(i,ir) - ied(i)
        enddo ! ir
      enddo ! i
      ir = 0
      do ic = ic1, ic2
        ir = ir+1
        cl(ic) = 1.d0
        cl(ic) = cl(ic)*product(primes(:)**ie(:,ir))
        cl(ic) = isi*cl(ic)
        isi = -isi
      enddo ! ic
      if (m == 0) il2p(1) = il2p(1)-1 ! undo the increasing done above for m=0
      
      up   = 1.d0
      upsq = 1.d0
      down = 1.d0
      do i = 1, ifmx
        ieupsq = 2*ied(i) + il2p(i) + l1(i)
        ieint = ieupsq/2 - l2(i)
        upsq  = upsq*primes(i)**mod(ieupsq, 2)
        if (ieint > 0) then 
          up = up*primes(i)**ieint
        elseif (ieint < 0) then
          down = down*primes(i)**(-ieint)
        endif ! for the ==0 case, we need to do neither of both operations
      enddo ! i
      coe(ice) = sqrt(upsq)*up/down
      
!     if (test('shape   ')) write(6,fmt="(2x,'l=',i2,' m=',i2,f10.3,' *sqrt(',f16.2,')/',f10.3/2x,'cl  :',6f14.2)") l,m,up,upsq,down,(cl(ic),ic=ic1,ic2)
      if (m == 0) goto 20
      la = l+m
      lb = l-m+1
      call factorize(la, primes, iea)
      call factorize(lb, primes, ieb)
      jm0(:) = jm0(:) + iea(:) - ieb(:)
      l1(:)  = l1(:)  - iea(:) + ieb(:)
      m = m-1
      goto 2
   20 continue
   
!     if (test('shape   ')) write(6,fmt="(80('*'))")
      if (l == lmax) return
      
      la = (2*l+1)*(2*l+2)
      lb = (l+1)*2
      call factorize(la, primes, iea)
      call factorize(lb, primes, ieb)
      l1st(:) = l1st(:) + iea(:)
      l2st(:) = l2st(:) + ieb(:)
      l = l+1
      goto 1
      
      endsubroutine ccoef


!-----------------------------------------------------------------------
!>    this routine reduces a positive integer   input number 'nmbr'
!>    to a product  of  first  numbers 'primes' , at powers  'iexp'.
!-----------------------------------------------------------------------
    subroutine factorize(nmbr, primes, iexp)
      integer, intent(in)  :: nmbr
      integer, intent(in)  :: primes(1:) ! set of lowest prime numbers
      integer, intent(out) :: iexp(1:)
!-----------------------------------------------------------------------

      integer :: i, nmb, ifmx
      
      if (nmbr <= 0) then
        write(*,fmt="(3x,i15,'  non positive number')") nmbr
        stop
      endif ! nmbr <= 0
      
      ifmx = min(size(primes), size(iexp))
       
      iexp(:) = 0
      if (nmbr == 1) return ! all iexp zero
      
      nmb = nmbr ! copy
      do i = 1, ifmx
!       iexp(i) = 0 ! redundant
        do while (mod(nmb, primes(i)) == 0)
          nmb = nmb/primes(i) ! integer division, reduce by the prime factor primes(i)
          iexp(i) = iexp(i)+1 ! and count up the number of exponents
        enddo ! while
        if (nmb == 1) return ! default exit point of this routine
      enddo ! i ! loop over all prime factors
      
      write(*,fmt="(3x,i15,'  cannot be reduced in the basis of first numbers given'/20x,'increase the basis of first numbers')") nmbr
      stop
      
    endsubroutine factorize

!------------------------------------------------------------------
!>    this routine computes transformation matrices associated to
!>    the rotation through the euler angles alpha,beta,gamma  for
!>    real spherical harmonics up to quantum number lmax. the re-
!>    sults are stored in dmatl(isumd)
!------------------------------------------------------------------
      subroutine d_real(lmax,euler, dmatl, isumd, lmaxd1)
      integer, intent(in) :: lmaxd1, isumd, lmax
      real*8, intent(in) :: euler(1:3) ! alpha, beta, gamma
      real*8, intent(out) :: dmatl(isumd)
      
!-----------------------------------------------------------------------
      integer :: l,m,mp,i,imax,ip,ipmax,isu!,isum
      real*8 :: fac,fac1,fac2,d,d1,d2
      real*8 :: dmn(lmaxd1+1,lmaxd1+1), dpl(lmaxd1+1,lmaxd1+1)
      real*8, parameter :: sqr2 = sqrt(2.d0)
      isu = 0 ! outer order is s,p,d,f,... and inner order for m and mp is 0,1,-1,2,-2,...,l,-l
      ! therefore as a function of lmax, isumd can be as low as sum_l=0...lmax (2*l+1)^2
      ! = (lmax*(3+4*(lmax+1)*(lmax+2)))/3+1
      do l = 0, lmax
      
        fac2 = 1.d0
        do m = 0, l
          fac1 = 1.d0
          do mp = 0, m
            fac = 0.5d0*fac1*fac2
#define beta euler(2)
            d1 = drot(l,mp, m,beta)
            d2 = drot(l,mp,-m,beta)
#undef  beta            
            if (mod(m, 2) /= 0) d2 = -d2
            dpl(mp+1,m+1) = (d1 + d2)*fac
            dmn(mp+1,m+1) = (d1 - d2)*fac
            if (mod(m+mp, 2) /= 0) then
              dmn(m+1,mp+1) = -dmn(mp+1,m+1)
              dpl(m+1,mp+1) = -dpl(mp+1,m+1)
            else
              dmn(m+1,mp+1) =  dmn(mp+1,m+1)
              dpl(m+1,mp+1) =  dpl(mp+1,m+1)
            endif
            fac1 = sqr2
          enddo ! mp
          fac2 = sqr2
        enddo ! m
        
        imax = 0
        do m = 0, l
          do i = 0, imax
            ipmax = 0
            do mp = 0, l
              do ip = 0, ipmax
#define alpha euler(1)
#define gamma euler(3)
                if (ip == 1) then
                  if (i == 1) then
                    d = -sin(mp*alpha)*sin(m*gamma)*dpl(mp+1,m+1) + cos(mp*alpha)*cos(m*gamma)*dmn(mp+1,m+1) ! i==1, ip==1
                  else  ! i == 1
                    d =  sin(mp*alpha)*cos(m*gamma)*dpl(mp+1,m+1) + cos(mp*alpha)*sin(m*gamma)*dmn(mp+1,m+1) ! i==0, ip==1
                  endif ! i == 1
                else  ! ip == 1
                  if (i == 1) then
                    d = -cos(mp*alpha)*sin(m*gamma)*dpl(mp+1,m+1) - sin(mp*alpha)*cos(m*gamma)*dmn(mp+1,m+1) ! i==1, ip==0
                  else  ! i == 1
                    d =  cos(mp*alpha)*cos(m*gamma)*dpl(mp+1,m+1) - sin(mp*alpha)*sin(m*gamma)*dmn(mp+1,m+1) ! i==0, ip==0
                  endif ! i == 1
                endif ! ip == 1
#undef  gamma
#undef  alpha
                if (mod(m+mp, 2) /= 0) d = -d
                isu = isu+1
                dmatl(isu) = d
              enddo ! ip
              ipmax = 1
            enddo ! mp
          enddo ! i
          imax = 1
        enddo ! m
      enddo ! l
!     isum = isu ! unused: minimal value of isumd
      endsubroutine d_real

!-----------------------------------------------------------------------
!>    calculation of d coefficient according to rose, elementary theory angular momentum,j.wiley & sons ,1957 , eq. (4.13).
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!>    calculation of d coefficient according to rose, elementary theory
!>    angular momentum,j.wiley & sons ,1957 , eq. (4.13).
!-----------------------------------------------------------------------
      real*8 function drot(l, mp, m, beta)
      
      integer, intent(in) :: l,m,mp
      real*8, intent(in) :: beta

      integer :: i,kmin,kmax,ltrm,n,k
      real*8  :: cosb,sinb,term,ff
      integer :: nf(4)
!-----------------------------------------------------------------------
      drot = 0.d0 ! init result for early return
      
      if (iabs(m) > l .or. iabs(mp) > l) then
        write(*,fmt="('     l=',i5,'    m=',i5,'    mp =',i5)") l,m,mp
        stop 'drot: |m| or |m''| is larger than l' ! error
      endif
      
      nf(1:4) = l + [m, -m, mp, -mp]
      ff = 1.d0
      do n = 1, 4
        if (nf(n) == 0) cycle
        do i = 1, nf(n)
          ff = ff*i
        enddo ! i
      enddo ! n
      ff = sqrt(ff)
      cosb =  cos(beta*0.5d0)
      sinb = -sin(beta*0.5d0)

      if (abs(cosb) < 1d-4) then
        ltrm = l
        term = sinb
      elseif(abs(sinb) < 1d-4) then
        ltrm = 0
        term = cosb
      else         
        kmax = min(l-mp, l+m)
        kmin = max(m-mp, 0)
        term = cosb**(2*l+m-mp-2*kmin) * sinb**(mp-m+2*kmin) * ff
        ltrm = -1 ! -1: do not execute the following if branch
      endif
      
      if (ltrm >= 0) then
        ! either |cos| or |sin| is small
        if (mod(m-mp, 2) /= 0) return ! 0.d0
        kmax = ltrm + (m-mp)/2
        if (kmax < max(m-mp, 0) .or. kmax > min(l-mp, l+m)) return ! 0.d0
!       if (term < 0.d0 .and. mod(mp-m, 2) /= 0) then; term = -ff; else; term = ff; endif
!       ! or term = sign(ff, term) ! uses the sign of term, i.e. sinb or cosb and the absolute of ff
        term = ff ! since the negative case is never reached because of the earlier return statement
        kmin = kmax
      endif ! ltrm >= 0
      
      if (mod(kmin, 2) /= 0) term = -term
   
      nf(1) = l-kmin-mp
      nf(2) = l-kmin+m
      nf(3) =   kmin+mp-m
      nf(4) =   kmin
      ff = 1.d0
      do n = 1, 4
        if (nf(n) == 0) cycle
        do i = 1, nf(n)
          ff = ff*i
        enddo ! i
      enddo ! n
      term = term/ff
      
      drot = term
      if (kmin == kmax) return
      kmin = kmin+1
      cosb = cosb**2
      sinb = sinb**2
      nf(3) = nf(3)+1
      do k = kmin, kmax
        term = -term*nf(1)*nf(2)*sinb/(cosb*k*nf(3))
        drot = drot + term
        nf(1:4) = nf(1:4) + [-1,-1,1,0]
      enddo ! k
      endfunction drot

      endmodule ShapeIntegrationHelpers_mod
