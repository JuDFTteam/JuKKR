      SUBROUTINE POLYHEDRON08(
     >                        NPLANE,NVERTMAX,NFACED,TOLVDIST,TOLAREA,
     X                        TOLHS,A3,B3,C3,D3,
     <                        NFACE,NVERT,XVERT,YVERT,ZVERT)
c Given a set of planes, defined by A3*x+B3*y+C3*z=D3 and defining
c a convex part of space (the minimal one containing the origin, 
c usually a WS-polyhedron), this subroutine returns the actual faces of
c the polyhedron, discarding the planes that do not contain faces. Also,
c the coordinates of the verticess of the faces XVERT,YVERT,ZVERT and their 
c number NVERT per face are returned. The coefficients of the actual
c faces are returned in the same arrays A3,B3,C3, and D3.
c
c Uses subroutine VERTEX3D. 
      implicit none
c Input:
      INTEGER NPLANE                ! Initial number of planes.
      INTEGER NVERTMAX              ! Max. number of vertices per plane.
      INTEGER NFACED                ! Max. number of faces.
      REAL*8  TOLVDIST              ! Max. tolerance for distance of two vertices
      REAL*8  TOLAREA               ! Max. tolerance for area of polygon face
      REAL*8  TOLHS                 ! Tolerance for halfspace routine
c Input and Output
      REAL*8            A3(*),B3(*),C3(*),D3(*)  ! Coefs. defining the planes, 
c                                     ! dimensioned >= NPLANE.
c Output:
      INTEGER NVERT(*)  ! Number of vertices found for each face
      INTEGER NFACE     ! Number of faces found (with nvert>0).
      INTEGER NVERTTOT  ! Total number of vertices
      REAL*8            XVERT(NVERTMAX,NFACED),YVERT(NVERTMAX,NFACED),
     &                  ZVERT(NVERTMAX,NFACED)
c                            ! Cartesian coords. of vertices for each plane
c                            ! (2nd index is for planes).
c Inside:
      INTEGER IPLANEHIGH,INDEXHIGH,INDEXLOW,IVERT,id
      LOGICAL TEST ! test options
c---------------------------------------------------------------
c Find all faces and vertices of the polyhedron. On output, the vertices of 
c each face are sorted (clockwise or anticlockwise).
      IF (TEST('verb1   ')) WRITE(*,*) 'Entering VERTEX3D'
      CALL VERTEX3D(
     >               NPLANE,A3,B3,C3,D3,NVERTMAX,TOLVDIST,TOLHS,
     <               NFACE,NVERT,XVERT,YVERT,ZVERT)
      IF (TEST('verb1   ')) 
     & WRITE(*,*) 'VERTEX3D found',NFACE,' faces with >3 vertices.'
c---------------------------------------------------------------
c Analyze the faces and vertices of the polyhedron.
c Use criteria for rejecting faces that are too small 
c or vertices that are too close to each other.
c On output, number of faces and vertices may be reduced 
c after some rejections have taken place.
      IF (TEST('verb1   ')) WRITE(*,*) 'Entering ANALYZEVERT3D'
      CALL ANALYZEVERT3D(
     >                    NVERTMAX,NFACED,TOLVDIST,TOLAREA,NPLANE,
     X                    NFACE,NVERT,XVERT,YVERT,ZVERT,
     X                    A3,B3,C3,D3)
      IF (TEST('verb1   '))
     &     WRITE(*,*) 'ANALYZEVERT3D accepted',NFACE,' faces.'

c---------------------------------------------------------------
c Pack the planes that contain faces at the beginning of the arrays
c A3, B3, C3, D3, and do the same for NVERT,XVERT,YVERT,ZVERT. The
c order is changed.

      RETURN


C **********************                                               FROM HERE ON NOT USED
C **********************                                               FROM HERE ON NOT USED
C **********************                                               FROM HERE ON NOT USED
C **********************                                               FROM HERE ON NOT USED
C **********************                                               FROM HERE ON NOT USED
C **********************                                               FROM HERE ON NOT USED

c You have to fill up the arrays up to NFACE, so...
         indexlow = 1
      DO while (INDEXLOW.le.nface)
          
         IF (NVERT(INDEXLOW).EQ.0) THEN
c     promote all planes by one
            do id = indexlow+1,nplane
               A3(Id-1) = A3(id)
               B3(Id-1) = B3(Id)
               C3(Id-1) = C3(Id)
               D3(Id-1) = D3(Id)
               NVERT(Id-1) = NVERT(Id) 
               DO IVERT = 1,NVERT(Id-1) 
                  XVERT(IVERT,Id-1) = XVERT(IVERT,Id)
                  YVERT(IVERT,Id-1) = YVERT(IVERT,Id)
                  ZVERT(IVERT,Id-1) = ZVERT(IVERT,Id)
               ENDDO
            end do
         else
            indexlow = indexlow + 1
         end if
      end do

      END







