       SUBROUTINE FPLANEG(lamda,g,pref,lmax,ga,vol)
c **************************************************************************
c This sub calculates the derivatives of the inverce  
c space contribution to the ewald sum
c 
c  
c     l     (                                                      )
c    d   pi*( exp(gz)*erfc(g/lam/lam + 2*z)*lam/2                  )  /
c           (             +   exp(-gz)*erfc(g/lam/lam - 2*z)*lam/2 ) / (g Vol)
c    ---    (                                                      )
c      l
c    dz
c
c   And the limit z -> 0 is taken (lam is the lamda parameter 
c
c *********************************************************************
      implicit none
      integer lmax
      double precision alpha,g(0:4),pref(0:lmax)
      integer l
      double precision derfc,lamda,er,ex,pi,ga,vol
      double precision sqpi
c
      do l=0,4
        g(l) = 0.d0
      end do
      do l=0,lmax
         pref(l) = 0.d0
      end do
      pi = 4.d0*atan(1.d0)
      sqpi=sqrt(pi)
      alpha = ga/2.d0/lamda
      er = derfc(alpha)
      ex = exp(-alpha*alpha)
c   
      if (abs(ga).gt.1.d-6) then
         g(0) = 2.d0*pi/vol*er/ga
      else
         g(0) = 0.d0
      end if
c
      pref(2) = sqrt(5.d0/pi)/2.d0/2.d0
      g(2) = pref(2)/vol*( 2.d0*pi*ga*er - ex*4.d0*sqpi*lamda)
c
      pref(4) = 3.d0*sqrt(9.d0/pi)/16.d0/9.d0
      g(4) = pref(4)/vol*(2.d0*pi*ga*ga*ga*er - ex*4.d0*sqpi* 
     &                                (lamda*ga*ga - 2.d0*lamda**3))
c
c      pref(6) = sqrt(13.d0/pi)/32.d0/13.d0
c      g(6) = pref(6)/vol*(2.d0*pi*ga**5*er - ex*4.d0*sqpi*
c     &                               (lamda*ga**4
c     &                               -2.d0*ga*ga*lamda*lamda
c     &                               +12.d0*lamda**5))                          
      end 
