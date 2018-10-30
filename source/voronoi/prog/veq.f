C ************************************************************************
      SUBROUTINE VEQ(a,b)
c#@# KKRtags: VORONOI
c#@# KKRmerge: remove this, only used in rrgen2000.f
C ************************************************************************
      REAL*8          a(*),b(*)
      integer i
      do 1 i=1,3
        b(i)=a(i)
 1    continue
      return
      END
