      SUBROUTINE MINMAX(N,A,AMIN,AMAX)
c Given an array A(N), this subroutine returns the minimum and maximum
c value of A.
      implicit none
c#@# KKRtags: VORONOI deprecated
c#@# KKRmerge: can be replaced by Fortran intrinsics
c Input:
      INTEGER N
      double precision A(*)
c Output:
      double precision AMIN,AMAX
c Inside:
      INTEGER I

      AMIN = A(1)
      AMAX = A(1)
      DO I = 2,N
         IF (A(I).LT.AMIN) THEN
            AMIN = A(I)
         ELSE
            IF (A(I).GT.AMAX) AMAX = A(I)
         ENDIF
      ENDDO

      RETURN
      END
