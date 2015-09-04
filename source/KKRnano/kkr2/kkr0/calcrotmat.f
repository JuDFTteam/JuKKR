      subroutine calcrotmat(nk, irel, alfdeg, betdeg, gamdeg,
     &                   rot, fact, nkmmax)
!   ********************************************************************
!   *                                                                  *
!   *   SETS UP THE ROTATION-MATRICES FOR THE EULER ANGLES             *
!   *           (ALFDEG, BETDEG, GAMDEG)                             *
!   *                                                                  *
!   *   SEE:     E.M. ROSE  ELEMENTARY THEORY OF ANGULAR MOMENTUM      *
!   *            EQS. (4.8), (4.12) AND (4.13)                         *
!   *                                                                  *
!   *   for IREL=0,1   NK == NL           non-relativistic (l,m_l)     *
!   *       IREL=3     NK == odd          relativistic (kappa,mue)     *
!   *                                                                  *
!   *   12/11/96  HE  deal with beta = 0                               *
!   ********************************************************************

      implicit none
!
      integer, intent(in) :: nk, irel, nkmmax
      real*8, intent(in) :: alfdeg, betdeg, gamdeg
      
      complex*16 ci, c0
      parameter (ci = (0.d0,1.d0), c0 = (0.d0,0.d0))
      real*8 pi
      parameter (pi = 3.141592653589793238462643d0)   
!
      real*8     num, msb05, msb05sq, msb05pw, j,m1,m2, rfac,x
      real*8     fact(0:100)

      integer    s, slow, shigh, off
      complex*16 emim2a, emim1g, rot(nkmmax,nkmmax) 
      integer :: i2, i1, k, l, nmue, im2, im1
      real*8 :: cb05, cb05sq, sm, cb05pw, dom
!                       
! inline function    factorial for real argument
      rfac(x) = fact(nint(x))
!
      if(irel == 2) call errortrap('calcrotmat',12,1) 
      if(irel == 3 .and. mod(nk, 2) == 0) 
     &            call errortrap('calcrotmat',13,1) 

      do 20 i2=1,nkmmax 
      do 20 i1=1,nkmmax 
20    rot(i1,i2) = c0
!
       cb05   =   dcos(betdeg*0.5d0*pi/180.d0)
       cb05sq =   cb05 *  cb05
      msb05   = - dsin(betdeg*0.5d0*pi/180.d0)
      msb05sq =  msb05 * msb05
!     
      off = 0
      do 100 k=1,nk 
      if(irel < 2) then
         l = k - 1
         j = l
      else
         l = k/2
         if(l*2 == k) then
            j = l - 0.5d0
         else 
            j = l + 0.5d0
         end if         
      end if

      nmue = nint(2*j + 1)
!
         do 90 im2 = 1, nmue
         m2 = - j + (im2-1.d0)
         emim2a = cdexp(-ci*m2*alfdeg*pi/180.d0)
!
            do 80 im1 = 1, nmue
            m1 = - j + (im1-1.d0)
            emim1g = cdexp(-ci*m1*gamdeg*pi/180.d0)
!
            if(dabs(betdeg) < 1d-8) then
               if(im1 == im2) then 
                  sm = 1.d0
               else
                  sm = 0.d0
               end if
            else
               slow   = max(         0, nint(m1-m2))
               shigh  = min(nint(j-m2), nint(j+m1))
                cb05pw =  cb05**nint(2*j+m1-m2-2*slow    +2)
               msb05pw = msb05**nint(   m2-m1+2*slow    -2)
               dom = (-1.d0)**(slow-1) *
     &          dsqrt(rfac(j+m1)*rfac(j-m1)*rfac(j+m2)*rfac(j-m2))
               sm = 0.d0
!
               do s=slow,shigh
                  dom = -dom
                  num =    fact(s) * rfac(j-m2-s)
     &                             * rfac(j+m1-s) * rfac(m2-m1+s)
                   cb05pw =  cb05pw /  cb05sq
                  msb05pw = msb05pw * msb05sq
                  sm = sm + (dom/num) * cb05pw * msb05pw
               end do
            end if
!
80          rot(off+im2,off+im1) = emim1g * sm * emim2a
!
90       continue

      off = off + nmue
100   continue
!
      return
      end
