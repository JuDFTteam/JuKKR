c***********************************************************************
      SUBROUTINE NORMALPLANE(X1,Y1,Z1,X2,Y2,Z2,TAU,A,B,C,D)
c Given two points in space, r1=(X1,Y1,Z1) and r2=(X2,Y2,Z2), this
c subroutine returns the coefficients defining a plane through the 
c equation A*x+B*y+C*z=D, which is normal to the vector r2-r1 and passes
c through the point (1.-TAU)*r1 + TAU*r2 (TAU thus being a parameter
c defining how close the plane is to each of the two points).
      implicit none
c#@# KKRtags: VORONOI geometry deprecated
c Input:
      REAL*8           X1,Y1,Z1,X2,Y2,Z2,TAU
c Output:
      REAL*8           A,B,C,D
c Inside:
      REAL*8           ONEMTAU
c The plane is defined as 
c (A,B,C)*(X-X1,Y-Y1,Z-Z1)=const=
c                         =(distance from r1 to (1.-TAU)*r1 + TAU*r2)**2
c so A,B,C are the coords. of a vector connecting the point r1 to
c the point (1.-TAU)*r1 + TAU*r2.
      ONEMTAU = 1.D0 - TAU

      A = ONEMTAU * X1 + TAU * X2
      B = ONEMTAU * Y1 + TAU * Y2
      C = ONEMTAU * Z1 + TAU * Z2
      D = A*(A+X1) + B*(B+Y1) + C*(C+Z1)

      RETURN
      END
