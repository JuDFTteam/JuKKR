c ************************************************************************
      INTEGER FUNCTION LENGTH(S,MAX)
c#@# KKRtags: VORONOI
c#@# KKRmerge: can be replaced by Fortran intrinsic trim_len
C ************************************************************************
      CHARACTER*1 S(*)
      INTEGER MAX

      INTEGER I
c ------------------------------------------------------------------------
      I = MAX

      DO 10 WHILE (S(I).EQ.' ')
        I = I - 1
 10   END DO

      LENGTH = I

      RETURN
      END
