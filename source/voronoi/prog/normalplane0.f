c***********************************************************************
      SUBROUTINE NORMALPLANE0(
     >     X1,Y1,Z1,TAU,
     <     A,B,C,D)
c Given a point in space, r1=(X1,Y1,Z1), this
c subroutine returns the coefficients defining a plane through the 
c equation A*x+B*y+C*z=D, which is normal to the vector r1 and passes
c through the point TAU*r1 (TAU thus being a parameter
c defining how close the plane is to the point).
      implicit none
c#@# KKRtags: VORONOI geometry
c Input:
      REAL*8           X1,Y1,Z1,TAU
c Output:
      REAL*8           A,B,C,D
c Inside:
c The plane is defined by 
c (A,B,C)*(X,Y,Z) = D = (tau * r1)**2
c so A,B,C are the coords. of the vector tau * r1.
c If tau=0 (plane passes through the origin), then D=0.

      IF (TAU.NE.0.D0) THEN
         A = TAU * X1
         B = TAU * Y1
         C = TAU * Z1
         D = A*A + B*B + C*C
      ELSE
         A = X1
         B = Y1
         C = Z1
         D = 0.D0
      ENDIF

      RETURN
      END
