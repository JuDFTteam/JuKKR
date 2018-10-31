c***********************************************************************
      REAL*8           FUNCTION DISTPLANE(A,B,C,D)
c Returns the distance of a plane A*x+B*y+C*z=D to the origin.
c#@# KKRtags: VORONOI geometry
      implicit none
      REAL*8           A,B,C,D
      REAL*8           ABCSQ

      ABCSQ = A*A + B*B + C*C

      IF (ABCSQ.LT.1.D-100) STOP 'DISTPLANE'

      DISTPLANE = DABS(D)/DSQRT(ABCSQ)  

      RETURN
      END
