      SUBROUTINE CROSSPOIPLANE(x1,y1,z1,x2,y2,z2,a3,b3,c3,d3,xcut,
     &     ycut,zcut,a)
      implicit none
c#@# KKRtags: VORONOI geometry
c       
c      R = Ro + a V      (eq of line )
c      -   -      -
c 
c      C.R = D           (eq of plane)
c      - -
c
c
c  Cross point if we find  a
c               
c       a = (D - C.Ro ) / C.V 
c                - -      - -
c
c   Ro = (x1,y1,z1), V = (x2-x1,y2-y1,z2-z1), C = (A3,B3,C3), D = D3 
c
      real*8 x1,y1,z1,x2,y2,z2,a3,b3,c3,d3,xcut,ycut,zcut,a 
c     
      a = ( d3 - (a3*x1+b3*y1+c3*z1) ) / 
     &                      ( a3*(x2-x1)+b3*(y2-y1)+c3*(z2-z1) )
      xcut = x1 + a*(x2-x1)
      ycut = y1 + a*(y2-y1)
      zcut = z1 + a*(z2-z1)  
      return 
      end 
