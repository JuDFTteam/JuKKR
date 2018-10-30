c***********************************************************************
      LOGICAL FUNCTION HALFSPACE(A,B,C,D,X,Y,Z,TOLHS)
c#@# KKRtags: VORONOI geometry
c Given a plane A*x+B*y+C*z=D, and a point (X,Y,Z) in space, this 
c func takes the value TRUE if (X,Y,Z) lies in the half-space 
c defined by the plane and the origin (0,0,0) (including the plane 
c itself). Else, the value FALSE is returned.
c
c The criterion used is that the inner product of the vector (X,Y,Z) 
c with the vector d connecting the origin to the plane vertically be 
c less than or equal to d**2:  (d_x,d_y,d_z)*(X,Y,Z) =< d**2.
c
c Input:
      REAL*8           A,B,C,D,X,Y,Z,TEST,D1,D2,TOLHS

      IF (DABS(A)+DABS(B)+DABS(C).LT.1.D-80) 
     &                                STOP 'HALFSPACE: A,B,C too small.'

      HALFSPACE = .FALSE.

!      IF (D*(A*X+B*Y+C*Z).LE.D*D) HALFSPACE = .TRUE.  
      IF (D*(A*X+B*Y+C*Z)-D*D.LT.TOLHS) HALFSPACE = .TRUE.  
c (re-checked 31May2008 FM)

      RETURN
      END
