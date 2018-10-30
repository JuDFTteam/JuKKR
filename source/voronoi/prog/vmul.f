C ************************************************************************
      SUBROUTINE VMUL(a,b,c)
c#@# KKRtags: VORONOI
c#@# KKRmerge: remove this, only used in rrgen2000.f
C ************************************************************************
      REAL*8          a(*),b,c(*)
      integer i
      do 1 i=1,3
        c(i)=b*a(i)
 1    continue
      return
      END
