c***********************************************************************
      SUBROUTINE VERTEX3D(
     >                    NPLANE,A3,B3,C3,D3,NVERTMAX,TOLVDIST,TOLHS,
     <                    NFACE,NVERT,XVERT,YVERT,ZVERT)
c Given a set of planes, defined by A3*x+B3*y+C3*z=D3 and defining
c a convex part of space (the minimal one containing the origin, 
c usually a WS-polyhedron), this subroutine returns the vertices
c of this polyhedron in cartesian coordinates. For the planes that
c are not faces of the polyhedron, a value NVERT(IPLANE)=0 is returned.
c The total no. of faces found is returned as NFACE.
c
c Uses logical function HALFSPACE
      implicit none
c Input:
      INTEGER NPLANE                ! Number of planes.
      INTEGER NVERTMAX              ! Max. number of vertices per plane.
      REAL*8           A3(*),B3(*),C3(*),D3(*)  ! Coefs. defining the planes, 
c                                     ! dimensioned >= NPLANE.
      REAL*8 TOLVDIST               ! Min. distance between vertices
      REAL*8 TOLHS                  ! Tolerance for halfspace routine
c Output:
      INTEGER NVERT(*)  ! Number of vertices found for each face
      INTEGER NFACE     ! Number of faces found (with nvert>0).
      REAL*8  XVERT(NVERTMAX,*),YVERT(NVERTMAX,*),ZVERT(NVERTMAX,*)
c                            ! Cartesian coords. of vertices for each plane
c                            ! (2nd index is for planes).
c Inside:
      INTEGER IPLANE1,IPLANE2,IPLANE3,IPLANE,KPLANE ! Plane indices
      INTEGER IVERT                    ! Vertex index
      REAL*8           XCUT,YCUT,ZCUT ! Cut point of three planes.
      REAL*8           DET,DETX,DETY,DETZ ! Determinants of 3x3 system for 
c                               ! XCUT,YCUT,ZCUT.
      REAL*8           DISTANCE   ! A distance criterium of two points in space
c The following are for sorting the vertices of each face:
      REAL*8           V1(3),V2(3),V3(3)   ! Auxiliary vectors...
      REAL*8           SINFIV1V2,COSFIV1V2 ! ...and their inner and outer products
      REAL*8           FI(NVERTMAX)        ! ...and also their relative angles.
      REAL*8           UV(3),VL,SN         ! Unit vector, length, sign of sin(fi)

      LOGICAL HALFSPACE    ! Function used, see function itself.
      LOGICAL LACCEPT      ! Determining whether a cut point is inside
c                          !                            the polyhedron.
c---------------------------------------------------------------
c Check & initialize
      IF (NPLANE.LT.4) WRITE(*,*) 'VERT3D: Error:NPLANE was only',NPLANE
      DO IPLANE = 1,NPLANE
         NVERT(IPLANE) = 0
      ENDDO
c===============================================================
c Start loop over all planes that can be cut:
      DO 120 IPLANE1 = 1,NPLANE
c Start loop over all other planes:
      DO 110 IPLANE2 = 1,NPLANE
      IF (IPLANE2.EQ.IPLANE1) GOTO 110 
c Start loop over all other-other (!) planes. Do from IPLANE2+1 to 
c NPLANE so that no pair is considered twice.
      DO 100 IPLANE3 = IPLANE2+1,NPLANE        ! nikos  IPLANE2+1,NPLANE
      IF (IPLANE3.EQ.IPLANE1) GOTO 100
c     IF (IPLANE3.EQ.IPLANE2) GOTO 100 ! added by nikos
c Solve the 3x3 system to find the cut point.
      DET= A3(IPLANE1)*(B3(IPLANE2)*C3(IPLANE3)-B3(IPLANE3)*C3(IPLANE2))
     &   + A3(IPLANE2)*(B3(IPLANE3)*C3(IPLANE1)-B3(IPLANE1)*C3(IPLANE3))
     &   + A3(IPLANE3)*(B3(IPLANE1)*C3(IPLANE2)-B3(IPLANE2)*C3(IPLANE1))

c---------------------------------------------------------------
      IF (DABS(DET).GT.1.D-12) THEN ! there is a cut point 

      DETX=D3(IPLANE1)*(B3(IPLANE2)*C3(IPLANE3)-B3(IPLANE3)*C3(IPLANE2))
     &   + D3(IPLANE2)*(B3(IPLANE3)*C3(IPLANE1)-B3(IPLANE1)*C3(IPLANE3))
     &   + D3(IPLANE3)*(B3(IPLANE1)*C3(IPLANE2)-B3(IPLANE2)*C3(IPLANE1))

      DETY=A3(IPLANE1)*(D3(IPLANE2)*C3(IPLANE3)-D3(IPLANE3)*C3(IPLANE2))
     &   + A3(IPLANE2)*(D3(IPLANE3)*C3(IPLANE1)-D3(IPLANE1)*C3(IPLANE3))
     &   + A3(IPLANE3)*(D3(IPLANE1)*C3(IPLANE2)-D3(IPLANE2)*C3(IPLANE1))

      DETZ=A3(IPLANE1)*(B3(IPLANE2)*D3(IPLANE3)-B3(IPLANE3)*D3(IPLANE2))
     &   + A3(IPLANE2)*(B3(IPLANE3)*D3(IPLANE1)-B3(IPLANE1)*D3(IPLANE3))
     &   + A3(IPLANE3)*(B3(IPLANE1)*D3(IPLANE2)-B3(IPLANE2)*D3(IPLANE1))

      XCUT = DETX/DET
      YCUT = DETY/DET
      ZCUT = DETZ/DET
c     write(6,333) IPLANE1,IPLANE2,IPLANE3,XCUT,YCUT,ZCUT
c333  format('Cutting point of planes ',3I5,':',3D15.7)
c-----------------------------------
c Accept this cut point as a vertex, if it belongs to the polyhedron. So,
c make a loop over all other (than IPLANE1,2,3) planes:
      LACCEPT = .TRUE.
      DO 50 KPLANE = 1,NPLANE
         IF (KPLANE.EQ.IPLANE1.OR.KPLANE.EQ.IPLANE2
     &                        .OR.KPLANE.EQ.IPLANE3) GOTO 50
         LACCEPT = LACCEPT.AND.HALFSPACE(A3(KPLANE),B3(KPLANE),
     &                       C3(KPLANE),D3(KPLANE),XCUT,YCUT,ZCUT,TOLHS)
 50   CONTINUE
c-----------------------------------
      IF (LACCEPT) THEN
c If the cut point found belongs to the cell, we accept it unless it has
c occured before for this face (IPLANE1). Such a situation is possible
c when 4 or more planes pass through the same point (e.g. for the vertices
c of the fcc WS-cell). So...
         DO IVERT = 1,NVERT(IPLANE1)
            DISTANCE = (XVERT(IVERT,IPLANE1)-XCUT)**2
     &               + (YVERT(IVERT,IPLANE1)-YCUT)**2
     &               + (ZVERT(IVERT,IPLANE1)-ZCUT)**2
            DISTANCE = DSQRT(DISTANCE)
            IF (DISTANCE.LT.TOLVDIST ) THEN 
               LACCEPT = .FALSE. ! vertex is too close to a previous one.
               EXIT              ! Jump loop, no need to continue.
            ENDIF
         ENDDO
      ENDIF
c Now we're ready to add the point to the vertex list.
      IF (LACCEPT) THEN
         NVERT(IPLANE1) = NVERT(IPLANE1) + 1
         XVERT(NVERT(IPLANE1),IPLANE1) = XCUT
         YVERT(NVERT(IPLANE1),IPLANE1) = YCUT
         ZVERT(NVERT(IPLANE1),IPLANE1) = ZCUT
      ENDIF

      ENDIF ! (DABS(DET).GT.1.D-12)
c---------------------------------------------------------------

 100  CONTINUE       ! IPLANE3
 110  CONTINUE       ! IPLANE2
c     write(6,*) 
c    & 'Number of vertices for plane ',iplane1,'  :',nvert(iplane1)
 120  CONTINUE       ! IPLANE1

c===============================================================
c Each plane should finally have either at least 3 vertices, if it is a
c face of the polyhedron, or none at all. Check this:
      DO IPLANE = 1,NPLANE
      IF (NVERT(IPLANE).EQ.1.OR.NVERT(IPLANE).EQ.2) THEN
      WRITE(*,*) 'VERTEX3D: Error:There is a problem with the vertices.'
      WRITE(*,*) 'For plane',IPLANE,
     &  ' ,only ',NVERT(IPLANE),' vertices were found.'
      ENDIF
      ENDDO

c===============================================================
c For each face of the polyhedron, sort the vertices in a consecutive order
c as vertices of the polygon. The order is not necessarily mathematically
c positive.
      NFACE = 0
      DO IPLANE = 1,NPLANE
      IF (NVERT(IPLANE).GE.3) THEN
         NFACE = NFACE + 1      ! Count the faces
         FI(1) = -4.D0          ! Just a number smaller than -pi.
c Unit vector in the direction of first vertex:
         VL = DSQRT( XVERT(1,IPLANE)**2 + 
     &               YVERT(1,IPLANE)**2 + ZVERT(1,IPLANE)**2 )
         UV(1) = XVERT(1,IPLANE) / VL
         UV(2) = YVERT(1,IPLANE) / VL
         UV(3) = ZVERT(1,IPLANE) / VL

c Define the vector connecting the first vertex to the (now-) second:
         V1(1) = XVERT(2,IPLANE) - XVERT(1,IPLANE)
         V1(2) = YVERT(2,IPLANE) - YVERT(1,IPLANE)
         V1(3) = ZVERT(2,IPLANE) - ZVERT(1,IPLANE)

         DO IVERT = 2,NVERT(IPLANE)
c Define the vector connecting the first vertex to the current one:
            V2(1) = XVERT(IVERT,IPLANE) - XVERT(1,IPLANE)
            V2(2) = YVERT(IVERT,IPLANE) - YVERT(1,IPLANE)
            V2(3) = ZVERT(IVERT,IPLANE) - ZVERT(1,IPLANE)
c Find the angle fi between v1 and v2 
c ( always, -pi < fi < pi from the definition of DATAN2 )
            COSFIV1V2 = V1(1)*V2(1) + V1(2)*V2(2) + V1(3)*V2(3)
            CALL CROSPR(V1,V2,V3) ! Cross product = |v1|*|v2|*sinfi
            SINFIV1V2 = DSQRT(V3(1)*V3(1) + V3(2)*V3(2) + V3(3)*V3(3))
c Sign of sinfi is defined with respect to unit vector uv (see above)
            SN = UV(1)*V3(1) + UV(2)*V3(2) + UV(3)*V3(3)
            IF (SN.LT.0) SINFIV1V2 = -SINFIV1V2
            
            IF (SINFIV1V2.EQ.0.D0.AND.COSFIV1V2.EQ.0.D0) THEN
c Point falls exactly on 1st vertex...
               FI(IVERT) = -4.D0
               WRITE(*,*) 
     &         'VERTEX3D: Error: Found two identical vertex points'
c ...while it shouldn't ! (this was checked earlier)
            ELSE
               FI(IVERT) = DATAN2(SINFIV1V2,COSFIV1V2)
            ENDIF

         ENDDO                  ! IVERT = 3,NVERT(IPLANE)

c Store with respect to the angle found:
         CALL SORTVERTICES(NVERT(IPLANE),FI,XVERT(1,IPLANE),
     &                   YVERT(1,IPLANE),ZVERT(1,IPLANE))



      ENDIF                     ! (NVERT(IPLANE).GE.3)
      ENDDO                     ! IPLANE = 1,NPLANE

      RETURN
      END












