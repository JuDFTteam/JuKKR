module Harmonics_mod
!-------------------------------------------------------------------------------
!> Summary: Spherical harmonic functions and Wigner 3j symbols (aka Gaunt coefficents)
!> Author: Bernhard H Drittler, Marc Weinert, Elias Rabel, Paul F Baumeister
!> Category: KKRnano, special-functions
!> ToDo: Check the simlarity of gaunt and gaunt2 with those in GauntCoefficients_mod
!-------------------------------------------------------------------------------
  implicit none
  private
  
  public :: ymy, gaunt, gaunt2, shapeg, shapeg_count
  
  contains

! ************************************************************************
  subroutine gaunt(lmax, lpot, w, yr, cleb, loflm, icleb, iend, jend, ncleb)
! ************************************************************************
!
!   - fills the array cleb with the gaunt coeffients ,i.e.
!      the integral of y(l1,m1)*y(l2,m2)*y(l3,m3)
!      but only for lm2 <= lm1 and lm3>1
!   - calculate the pointer array jend  to project the indices
!      array cleb with the same lm3,l1,l2 values - because of
!      the special ordering of array cleb only the last index
!      has to be determined .
!     (the parameter n has to be chosen that l1+l2+l3  <  2*n)
!     using gaussian quadrature as given by
!     m. abramowitz and i.a. stegun, handbook of mathematical functions,
!     nbs applied mathematics series 55 (1968), pages 887 and 916
!     m. weinert and e. wimmer
!     northwestern university march 1980
!
!     an index array -icleb- is used to save storage place .
!     fills the array loflm which is used to determine the
!     l-value of a given lm-value .
!     this subroutine has to be called only once !
!
!                               b.drittler   november 1987
!
!     modified gaunt coefficients are als calculated defined by
!     the integral of y(l1,m1)*y(l2,m2)*y(l3,m3)*i**(l2-l1+l3)
!-----------------------------------------------------------------------
!
!---> attention : ncleb is an empirical factor - it has to be optimized
    integer, intent(out) :: iend
    integer, intent(in) :: lmax
    integer, intent(in) :: lpot
    integer, intent(in) :: ncleb
    double precision, intent(out) :: cleb(ncleb,2)
    double precision, intent(in) :: w(1:)
    double precision, intent(in) :: yr(4*lmax,0:4*lmax,0:4*lmax) ! (1:n,0:n,0:n)
    integer, intent(out) :: icleb(ncleb,3)
    integer, intent(out) :: jend((lpot+1)**2,0:lmax,0:lmax) ! (lmpotd,0:lmaxd,0:lmaxd)
    integer, intent(out) :: loflm(1:)

    double complex, parameter :: ci=(0.d0,1.d0)
    double precision clecg,factor,fci,s
    integer :: n,i,j,l,l1,l1p,l2,l2p,l3,lm1,lm2,lm3,lm3p,m,m1,m1a,m1s,m2,m2a,m2s,m3,m3a,m3s

    n = 4*lmax

    i = 0
    do l = 0, 2*lmax
      do m = -l, l
        i = i + 1
        loflm(i) = l
      enddo ! m
    enddo ! l
!
    if (lpot == 0) then
      iend = 1
      icleb(1,1) = (lmax+1)**2
      icleb(1,2) = 1
      icleb(1,3) = 1
      return
    endif
!
!---> set up of the gaunt coefficients with an index field
!
    i = 0
    do l3 = 1, lpot
      do m3 = -l3, l3
        m3s = sign(1, m3)
        m3a = abs(m3)
        do l1 = 0, lmax
          do l2 = 0, l1
            if (mod(l1+l2+l3, 2) /= 1 .and. l1+l2 >= l3 .and. l1+l3 >= l2 .and. l2+l3 >= l1) then
              fci = dble(ci**(l2-l1+l3)) ! real part of i^(l2-l1+l3)
              do m1 = -l1, l1
                lm1 = l1*l1 + l1 + m1 + 1
                m1s = sign(1, m1)
                m1a = abs(m1)
                do m2 = -l2, l2
                  lm2 = l2*l2 + l2 + m2 + 1

!---> store only gaunt coeffients for lm2 <= lm1
                  if (lm2 <= lm1) then
                    m2s = sign(1, m2)
                    if (m1s*m2s*m3s >= 0) then
                      factor = 0.d0
                      m2a = abs(m2)

                      if (m1a+m2a == m3a) factor = factor + 0.125d0*(3*m3s + sign(1, -m3))
                      if (m1a-m2a == m3a) factor = factor + 0.25d0*m1s
                      if (m2a-m1a == m3a) factor = factor + 0.25d0*m2s

                      if (factor /= 0.d0) then
                        if (m1s*m2s /= 1 .or. m2s*m3s /= 1 .or. m1s*m3s /= 1) factor = -factor

                        s = 0.d0
                        do j = 1, n
                          s = s + w(j)*yr(j,l1,m1a)*yr(j,l2,m2a)*yr(j,l3,m3a)
                        enddo ! j
                        
                        clecg = s*factor
                        if (abs(clecg) > 1.d-10) then
                          i = i + 1
                          cleb(i,1) = clecg
                          cleb(i,2) = fci*clecg
                          icleb(i,1) = lm1
                          icleb(i,2) = lm2
                          icleb(i,3) = l3*l3 + l3 + m3 + 1
                        endif ! |c| > 10^-10
                      endif ! factor /= 0
                    endif ! m1s*m2s*m3s >= 0
                  endif ! lm2 <= lm1
                enddo ! m2
              enddo ! m1
            endif ! mod(l1+l2+l3, 2) /= 1 .and. l1+l2 >= l3 .and. l1+l3 >= l2 .and. l2+l3 >= l1
          enddo ! l2
        enddo ! l1
      enddo ! m3
    enddo ! l3
    iend = i
    if (ncleb < iend) then
      write (6,fmt="(13x,'error stop in gaunt : dimension of NCLEB = ',i10,' too small ',/,13x,'change NCLEB to ',i6)") ncleb,iend
      stop '33' ! call rcstop('33      ')
    endif
!
!---> set up of the pointer array jend,use explicitly the ordering of the gaunt coeffients
!
    jend(:,:,:) = 0
!     lmpot = (lpot+1)**2
!     do l1 = 0, lmax
!       do l2 = 0, l1
!         do lm3 = 2, lmpot
!           jend(lm3,l1,l2) = 0
!         enddo ! lm3
!       enddo ! l2
!     enddo ! l1

      j = 1
    l1 = loflm(icleb(j,1))
    l2 = loflm(icleb(j,2))
    lm3 =      icleb(j,3)

    do j = 2, iend
      l1p = loflm(icleb(j,1))
      l2p = loflm(icleb(j,2))
      lm3p =      icleb(j,3)

      if (lm3 /= lm3p .or. l1 /= l1p .or. l2 /= l2p) then
        jend(lm3,l1,l2) = j - 1
        l1 = l1p
        l2 = l2p
        lm3 = lm3p
      endif
    enddo ! j
    jend(lm3,l1,l2) = iend
      
  endsubroutine ! gaunt

      
      
      
!>    sets up values needed for gaunt.

!>       m. weinert  january 1982
!>
!>    changed for calculating with real spherical harmonics
!>                                          b.drittler  july 1987
!>
!>    W(N)        integration weights on 4*LMAX points in the intervall
!>                (-1,0) (from routine GRULE)
!>
!>    YR(N,L,M)   spherical harmonics on 4*LMAX points to angular
!>                momentum indices (l,m) scaled with a factor
!>                of RF=(4*pi)**(1/3)
!>
!>    LMAX        l-cutoff, added after inc.p remove
!***********************************************************************
  subroutine gaunt2(w, yr, lmax)
! ************************************************************************
!     sets up values needed for gaunt.

!        m. weinert  january 1982
!
!     changed for calculating with real spherical harmonics
!                                           b.drittler  july 1987
!
!     W(N)        integration weights on 4*LMAX points in the intervall
!                 (-1,0) (from routine GRULE)
!
!     YR(N,L,M)   spherical harmonics on 4*LMAX points to angular
!                 momentum indices (l,m) scaled with a factor
!                 of RF=(4*pi)**(1/3)
!
!     LMAX        l-cutoff, added after inc.p remove
!
!-----------------------------------------------------------------------
    integer, intent(in) :: lmax
    double precision, intent(out) :: w(1:)
    double precision, intent(out) :: yr(4*lmax,0:4*lmax,0:4*lmax)

    double precision a,cd,cth,fac,fpi,rf,sth
    integer k,l,lomax,m,n
    double precision :: p(0:4*lmax+1,0:4*lmax), x(4*lmax)
    
    n = 4*lmax
    fpi = 16.d0*atan(1.d0)
    rf = fpi**(1.d0/3.d0)
    lomax = n
!
!--->    obtain gauss-legendre points and weights
!
    call grule(2*n, x, w)
!
!--->    generate associated legendre functions for m >= 0
!
    do k = 1, n
      cth = x(k)
      sth = sqrt(1.d0 - cth*cth)
      fac = 1.d0
!
!--->    loop over m values
!
      do m = 0, lomax
        fac = (1-2*m)*fac
        p(m,m) = fac
        p(m+1,m) = (2*m+1)*cth*fac
!
!--->    recurse upward in l
!
        do l = m+2, lomax
          p(l,m) = ((2*l-1)*cth*p(l-1,m) - (l+m-1)*p(l-2,m))/dble(l-m)
        enddo ! l
        fac = fac*sth
      enddo ! m
!
!--->    multiply in the normalization factors
!
      do l = 0, lomax
        a = rf*sqrt((2*l+1)/fpi)
        cd = 1.d0
        yr(k,l,0) = a*p(l,0)
        do m = 1, l
          cd = cd/dble((l+1-m)*(l+m))
          yr(k,l,m) = a*sqrt(2.d0*cd)*p(l,m)
        enddo ! m
      enddo ! l
    enddo ! k
    
  endsubroutine ! gaunt2
  
  
  
!>    determines the (n+1)/2 nonnegative points x(i) and
!>    the corresponding weights w(i) of the n-point
!>    gauss-legendre integration rule, normalized to the
!>    interval [-1,1]. the x(i) appear in descending order.

!>    this routine is from 'methods of numerical integration',
!>    p.j. davis and p. rabinowitz, page 369.

! **********************************************************************
  subroutine grule(n, x, w)
!
!***********************************************************************
!
!     determines the (n+1)/2 nonnegative points x(i) and
!     the corresponding weights w(i) of the n-point
!     gauss-legendre integration rule, normalized to the
!     interval [-1,1]. the x(i) appear in descending order.
!
!     this routine is from 'methods of numerical integration',
!     p.j. davis and p. rabinowitz, page 369.
!
!***********************************************************************
    integer, intent(in) :: n
    double precision, intent(out) :: w(1:), x(1:)

    double precision :: d1,d2pn,d3pn,d4pn,den,dp,dpn,e1,fx,h,p,pi,pk,pkm1,pkp1,t,t1,u,v,x0
    integer :: i,it,k,m

    pi = 4.d0*atan(1.d0)
    m  = (n+1)/2 ! integer divide
    e1 = dble(n*(n+1))
    do i = 1, m
      t = (4*i-1)*pi/dble(4*n+2)
      x0 = (1.d0 - (1.d0 - 1.d0/dble(n))/(8.d0*n*n))*cos(t)
!
!--->    iterate on the value  (m.w. jan. 1982)
!
      do it = 1, 2
        pkm1 = 1.
        pk = x0
        do k = 2, n
          t1 = x0*pk
          pkp1 = t1 - pkm1 - (t1-pkm1)/k + t1
          pkm1 = pk
          pk = pkp1
        enddo ! k
        den = 1. - x0*x0
        d1 = n* (pkm1 - x0*pk)
        dpn = d1/den
        d2pn = (2.d0*x0*dpn - e1*pk)/den
        d3pn = (4.d0*x0*d2pn + (2.d0-e1)*dpn)/den
        d4pn = (6.d0*x0*d3pn + (6.d0-e1)*d2pn)/den
        u = pk/dpn
        v = d2pn/dpn
        h = -u*(1.d0 + 0.5d0*u*(v + u*(v*v - u*d3pn/(3.d0*dpn))))
        p = pk + h*(dpn + 0.5d0*h*(d2pn + h/3.d0*(d3pn + 0.25d0*h*d4pn)))
        dp = dpn + h*(d2pn + 0.5d0*h*(d3pn + h*d4pn/3.d0))
        h = h - p/dp
        x0 = x0 + h
      enddo ! it
      x(i) = x0
      fx = d1 - h*e1*(pk + 0.5d0*h*(dpn + h/3.d0*(d2pn + 0.25d0*h*(d3pn + 0.2d0*h*d4pn))))
      w(i) = 2.d0*(1. - x(i)*x(i))/(fx*fx)
    enddo ! i
    if (2*m > n) x(m) = 0.d0 ! center if n is odd
  endsubroutine ! grule
      
      
      
!-----------------------------------------------------------------------
!>  - prepares shape corrections.
!>    (the parameter n has to be chosen that l1+l2+l3  <=  2*n)
!>    using gaussian quadrature as given by
!>    m. abramowitz and i.a. stegun, handbook of mathematical functions,
!>    nbs applied mathematics series 55 (1968), pages 887 and 916
!-----------------------------------------------------------------------
!***********************************************************************
  subroutine shapeg_impl(lpot, gsh, ilm, imaxsh, w, yr, lmax, ngshd, ncount)
!***********************************************************************
    integer, intent(in) :: ngshd
    integer, intent(in) :: lmax
    integer, intent(in) :: lpot
    double precision, intent(in) :: w(1:)
    double precision, intent(in) :: yr(4*lmax,0:4*lmax,0:4*lmax)
    
    double precision, intent(out) :: gsh(1:)
    integer, intent(out) :: ilm(ngshd,3)
    integer, intent(out) :: imaxsh(0:(lpot+1)**2) ! (0:lmpotd)
    
    integer, intent(out), optional :: ncount ! if present, arrays gsh and ilm are not written
    
    double precision factor,gaunt,s
    integer i,j,l1,l2,l3,lm1,lm2,lm3,m1,m1a,m1s,m2,m2a,m2s,m3,m3a,m3s
    integer n,lassld,lmpotd,lmxspd

    n = 4*lmax
    lassld = n
    lmpotd = (lpot+1)**2
    lmxspd = (2*lpot+1)**2
    
!---> set up of the gaunt coefficients with an index field so c(lm,lm',lm'') is mapped to c(i)
    i = 0
    do l1 = 0, lpot
      do m1 = -l1, l1
        lm1 = l1*l1 + l1 + m1 + 1
        imaxsh(lm1-1) = i
        do l3 = 0, 2*lpot
          do m3 = -l3, l3
            lm3 = l3*l3 + l3 + m3 + 1
            do l2 = 0, lpot
              if (mod(l1+l2+l3, 2) /= 1 .and. l1+l2 >= l3 .and. l1+l3 >= l2 .and. l2+l3 >= l1) then
                do m2 = -l2, l2
                  lm2 = l2*l2 + l2 + m2 + 1
!---> use the m-conditions for the gaunt coefficients not to be 0
                  m1s = sign(1, m1)
                  m2s = sign(1, m2)
                  m3s = sign(1, m3)

                  if (m1s*m2s*m3s >= 0) then
                    m1a = abs(m1)
                    m2a = abs(m2)
                    m3a = abs(m3)

                    factor = 0.d0

                    if (m1a+m2a == m3a) factor = factor + 0.125d0*(3*m3s + sign(1, -m3))
                    if (m1a-m2a == m3a) factor = factor + 0.25d0*m1s
                    if (m2a-m1a == m3a) factor = factor + 0.25d0*m2s

                    if (factor /= 0.d0) then

                      if (m1s*m2s /= 1 .or. m2s*m3s /= 1 .or. m1s*m3s /= 1) factor = -factor
                      s = 0.d0
                      do j = 1, n
                        s = s + w(j)*yr(j,l1,m1a)*yr(j,l2,m2a)*yr(j,l3,m3a)
                      enddo ! j
                      gaunt = s*factor
                      if (abs(gaunt) > 1.d-10) then

                        i = i + 1
                        
                        if (.not. present(ncount)) then
!                           e. r.: add check to avoid write add of bounds
                          if (i > ngshd) then
                            write(*,*) "shapeg: NGSHD too small.", i, ngshd
                            stop
                          endif

                          gsh(i) = gaunt
                          ilm(i,1) = lm1
                          ilm(i,2) = lm2
                          ilm(i,3) = lm3
                        endif ! present ncount
                        
                      endif
                    endif
                  endif
                enddo ! m2
              endif
            enddo ! l2
          enddo ! m3
        enddo ! l3
      enddo ! m1
    enddo ! l1
    imaxsh(lm1) = i
    
    if (present(ncount)) ncount = i

!     WRITE (*,FMT="(' >>> SHAPE : IMAXSH(',I6,'),NGSHD :',2I6)") IMAXSH(LM1),NGSHD
!     IF (IMAXSH(LM1) > NGSHD) CALL RCSTOP('SHAPE   ')
  endsubroutine ! shapeg_impl

!======================================================================


  subroutine shapeg_count(lpot, w, yr, lmax, ncount)
    integer, intent(in) :: lpot
    integer, intent(in) :: lmax
    double precision, intent(in) :: w(1:)
    double precision, intent(in) :: yr(4*lmax,0:4*lmax,0:4*lmax)
    integer, intent(out) :: ncount

    double precision :: gsh(1:1) ! dummy
    integer :: ilm(1:1,3) ! dummy
    integer :: imaxsh(0:(lpot+1)**2) ! (0:lmpotd)
    
    !< calling shapeg with the optional argument ncount present will
    ! not write to the arrays gsh and ilm and simply return the maximum index
    call shapeg_impl(lpot, gsh, ilm, imaxsh, w, yr, lmax, ngshd=1, ncount=ncount)
    
  endsubroutine ! shapeg_count
  
  
  subroutine shapeg(lpot, gsh, ilm, imaxsh, w, yr, lmax, ngshd)
!***********************************************************************
    integer, intent(in) :: ngshd
    integer, intent(in) :: lmax
    integer, intent(in) :: lpot
    double precision, intent(in) :: w(1:)
    double precision, intent(in) :: yr(4*lmax,0:4*lmax,0:4*lmax)
    
    double precision, intent(out) :: gsh(1:)
    integer, intent(out) :: ilm(ngshd,3)
    integer, intent(out) :: imaxsh(0:(lpot+1)**2) ! (0:lmpotd)
    
    call shapeg_impl(lpot, gsh, ilm, imaxsh, w, yr, lmax, ngshd)
    
  endsubroutine ! shapeg
  
  
! **********************************************************************
!>   this subroutine calculates real spherical harmonics with the
!>    normalization : <y|y> =1
!>   returns also r = length of vector v
!
!>    generate the complex spherical harmonics for the vector v
!>    using a stable upward recursion in l.  (see notes
!>    by m. weinert.)
!>                                 m.weinert  1982
!
!>    converted to real spherical harmonics .
!>                                 b.drittler 1987
! **********************************************************************
  subroutine ymy(v1, v2, v3, r, ylm, lmax)
    use Constants_mod, only: pi
    double precision, intent(in) :: v1, v2, v3
    double precision, intent(out) :: r
    integer, intent(in) :: lmax
    double precision, intent(out) :: ylm(1:)
    
    double precision, parameter :: szero=1.d-20
    
    double precision :: a,cd,cph,cth,fac,fpi,rtwo,sgm,sph,sth,t,xy2,xyz2,xy
    double precision :: c(0:lmax), p(0:lmax,0:lmax), s(0:lmax)
    integer :: i, l, m

    fpi = 4.d0*pi
    rtwo = sqrt(2.d0)

!--->    calculate sin and cos of theta and phi
    xy2 = v1*v1 + v2*v2
    xyz2 = xy2 + v3*v3

    if (xyz2 <= 0.d0) stop 'ERROR Ylm: [x,y,z]=0 !!!' ! call rcstop('ylm=0   ')
    r = sqrt(xyz2)
 
    if (xy2 > szero*xyz2) then
      xy = sqrt(xy2)
      ! xyz = sqrt(xyz2) ! done for r
      cth = v3/r
      sth = xy/r
      cph = v1/xy
      sph = v2/xy
    else
      sth = 0.d0
      cth = 1.d0
      if (v3 < 0) cth = -1.d0
      cph = 1.d0
      sph = 0.d0
    endif

!--->    generate associated legendre functions for m >= 0
!        loop over m values
    fac = 1.d0
    do m = 0, lmax-1
      fac = -fac*(2*m-1)
      p(m,m) = fac
      p(m+1,m) = (2*m+1)*cth*fac

!--->    recurse upward in l
      do l = m+2, lmax
        p(l,m) = ((2*l-1)*cth*p(l-1,m) - (l+m-1)*p(l-2,m))/dble(l-m)
      enddo ! l
      fac = fac*sth
    enddo ! m
    p(lmax,lmax) = -fac*(2*lmax-1)

!--->    determine sin and cos of phi ! todo: can this be done much easier with complex multiplication?
    s(0:1) = [0.d0, sph]
    c(0:1) = [1.d0, cph]
    do m = 2, lmax
      s(m) = 2.d0*cph*s(m-1) - s(m-2)
      c(m) = 2.d0*cph*c(m-1) - c(m-2)
    enddo ! m

!--->    multiply in the normalization factors
    i = 0
    do l = 0, lmax
      i = i + l + 1 ! forward i to the center (m==0) of this l-section
      a = sqrt(dble(2*l+1)/fpi)
      cd = 1.d0
      ylm(i) = a*p(l,0) ! m == 0
      sgm = -rtwo
      do m = 1, l
        cd = cd/dble((l+1-m)*(l+m))
        t = a*sqrt(cd)
        ylm(i+m) = sgm*t*p(l,m)*c(m)
        ylm(i-m) = sgm*t*p(l,m)*s(m)
        sgm = -sgm
      enddo ! m
      i = i + l ! forward i to the end (m==l) of this l-section
    enddo ! l
      
  endsubroutine ! ymy

endmodule ! Harmonics_mod
