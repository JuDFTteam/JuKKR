      SUBROUTINE derivat(x,y,n,yp1,ypn)
      INTEGER n
      parameter (NMAX=1000)
      double precision yp1,ypn,x(NMAX),y(NMAX)
      yp1=(y(2)-y(1))/(x(2)-x(1))
      ypn=(y(n)-y(n-1))/(x(n)-x(n-1))
C      print*,'yp1=',yp1,' ypn=',ypn
      return
      end

C     Subroutines from Numerical Recipes chapter 3.3
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      parameter (NMAX=1000)
      double precision yp1,ypn,x(NMAX),y(NMAX),y2(NMAX)
      integer i,k
      double precision p,qn,sig,un,u(NMAX)
C     
C     calculates the second derivative y2
C
C      print*,'the spline subroutine starts'
      if (yp1.gt.0.99d30) then
          y2(1)=0.0d0
          u(1)=0.0d0
      else
          y2(1)=-0.5d0
          u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))/yp1)
      endif
C      print*,'FIRST boundary'
      do i=2,n-1
          sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
          p=sig*y2(i-1)+2.d0
          y2(i)=(sig-1.d0)/p
          u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     &         /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo
      if (ypn.gt.0.99d30) then
          qn=0.0d0
          un=0.0d0
      else
          qn=0.5d0
          un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
C      print*,'SECOND boundary'
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
C     backsubstitution of the tridiagonal algorithm
      do k=n-1,1,-1
          y2(k)=y2(k)*y2(k+1)+u(k)
      enddo
      return
      end
C
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
C     supply the old xa(1:n),ya(1:n) and second derivative y2a(1:n)
C     calculates y's for new x-mesh
      INTEGER n
      parameter (NMAX=1000)
      double precision x,y,xa(NMAX),ya(NMAX),y2a(NMAX)
      integer k,khi,klo
      double precision a,b,h
C
      klo=1
      khi=n
 1    if (khi-klo.gt.1) then
          k=(khi+klo)/2
          if (xa(k).gt.x) then
              khi=k
          else
              klo=k
          endif
          goto 1
      endif
C      print*,'n=',n
C      print*,'k=',k
C      print*,'xa(',k,')=',xa(k),'x=',x
C      print*,'xa(',khi,')=',xa(khi),'xa(',klo,')=',xa(klo)
      h=xa(khi)-xa(klo)
      if (h.eq.0.d0) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     &     ((a**3.d0-a)*y2a(klo)+(b**3.d0-b)*y2a(khi))*(h**2.d0)/6.d0
      return
      end
      
