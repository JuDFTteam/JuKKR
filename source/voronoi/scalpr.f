C ************************************************************************
      Subroutine SCALPR(x,y,z)
c#@# KKRtags: VORONOI
C ************************************************************************
C     SCALSP COMPUTES THE scalar PRODUCT OF X AND Y RETURNING
C     IT INTO Z.
C ------------------------------------------------------------------------
      REAL*8           X(*), Y(*), Z
      z= x(1)*y(1)+x(2)*y(2)+x(3)*y(3)
      return
      END
