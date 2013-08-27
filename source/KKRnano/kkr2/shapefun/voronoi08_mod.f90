MODULE VORONOI08_MOD

CONTAINS

SUBROUTINE POLYHEDRON08( &
                        NPLANE,NVERTMAX,NFACED,TOLVDIST,TOLAREA, &
                        A3,B3,C3,D3, &
                        NFACE,NVERT,XVERT,YVERT,ZVERT,output)
! Given a set of planes, defined by A3*x+B3*y+C3*z=D3 and defining
! a convex part of space (the minimal one containing the origin, 
! usually a WS-polyhedron), this subroutine returns the actual faces of
! the polyhedron, discarding the planes that do not contain faces. Also,
! the coordinates of the verticess of the faces XVERT,YVERT,ZVERT and their 
! number NVERT per face are returned. The coefficients of the actual
! faces are returned in the same arrays A3,B3,C3, and D3.
!
! Uses subroutine VERTEX3D. 
implicit none

logical :: output

! Input:
INTEGER NPLANE                ! Initial number of planes.
INTEGER NVERTMAX              ! Max. number of vertices per plane.
INTEGER NFACED                ! Max. number of faces.
REAL*8  TOLVDIST              ! Max. tolerance for distance of two vertices
REAL*8  TOLAREA               ! Max. tolerance for area of polygon face
! Input and Output
REAL*8            A3(*),B3(*),C3(*),D3(*)  ! Coefs. defining the planes,
!                                     ! dimensioned >= NPLANE.
! Output:
INTEGER NVERT(*)  ! Number of vertices found for each face
INTEGER NFACE     ! Number of faces found (with nvert>0).
INTEGER NVERTTOT  ! Total number of vertices
REAL*8            XVERT(NVERTMAX,NFACED),YVERT(NVERTMAX,NFACED), &
                  ZVERT(NVERTMAX,NFACED)
!                            ! Cartesian coords. of vertices for each plane
!                            ! (2nd index is for planes).
! Inside:
INTEGER IPLANEHIGH,INDEXHIGH,INDEXLOW,IVERT,id

!---------------------------------------------------------------
! Find all faces and vertices of the polyhedron. On output, the vertices of 
! each face are sorted (clockwise or anticlockwise).
IF (output) then
  WRITE(*,*) 'Entering VERTEX3D'
END IF

CALL VERTEX3D( &
               NPLANE,A3,B3,C3,D3,NVERTMAX, &
               NFACE,NVERT,XVERT,YVERT,ZVERT,output)
IF (output) WRITE(*,*) 'VERTEX3D found',NFACE, &
            ' faces with >3 vertices.'
!---------------------------------------------------------------
! Analyze the faces and vertices of the polyhedron.
! Use criteria for rejecting faces that are too small 
! or vertices that are too close to each other.
! On output, number of faces and vertices may be reduced 
! after some rejections have taken place.
IF (output) WRITE(*,*) 'Entering ANALYZEVERT3D'
CALL ANALYZEVERT3D( &
                    NVERTMAX,NFACED,TOLVDIST,TOLAREA,NPLANE, &
                    NFACE,NVERT,XVERT,YVERT,ZVERT, &
                    A3,B3,C3,D3, output)
IF (output) WRITE(*,*) 'ANALYZEVERT3D accepted',NFACE,' faces.'

!---------------------------------------------------------------
! Pack the planes that contain faces at the beginning of the arrays
! A3, B3, C3, D3, and do the same for NVERT,XVERT,YVERT,ZVERT. The
! order is changed.

RETURN


! **********************                                               FROM HERE ON NOT USED
! **********************                                               FROM HERE ON NOT USED
! **********************                                               FROM HERE ON NOT USED
! **********************                                               FROM HERE ON NOT USED
! **********************                                               FROM HERE ON NOT USED
! **********************                                               FROM HERE ON NOT USED

! You have to fill up the arrays up to NFACE, so...
   indexlow = 1
DO while (INDEXLOW.le.nface)
    
   IF (NVERT(INDEXLOW).EQ.0) THEN
!     promote all planes by one
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

END SUBROUTINE







!***********************************************************************
SUBROUTINE VERTEX3D( &
                    NPLANE,A3,B3,C3,D3,NVERTMAX, &
                    NFACE,NVERT,XVERT,YVERT,ZVERT, output)
! Given a set of planes, defined by A3*x+B3*y+C3*z=D3 and defining
! a convex part of space (the minimal one containing the origin, 
! usually a WS-polyhedron), this subroutine returns the vertices
! of this polyhedron in cartesian coordinates. For the planes that
! are not faces of the polyhedron, a value NVERT(IPLANE)=0 is returned.
! The total no. of faces found is returned as NFACE.
!
! Uses logical function HALFSPACE
implicit none

logical :: output

! Input:
INTEGER NPLANE                ! Number of planes.
INTEGER NVERTMAX              ! Max. number of vertices per plane.
REAL*8           A3(*),B3(*),C3(*),D3(*)  ! Coefs. defining the planes, 
!                                     ! dimensioned >= NPLANE.
! Output:
INTEGER NVERT(*)  ! Number of vertices found for each face
INTEGER NFACE     ! Number of faces found (with nvert>0).
REAL*8  XVERT(NVERTMAX,*),YVERT(NVERTMAX,*),ZVERT(NVERTMAX,*)
!                            ! Cartesian coords. of vertices for each plane
!                            ! (2nd index is for planes).
! Inside:
INTEGER IPLANE1,IPLANE2,IPLANE3,IPLANE,KPLANE ! Plane indices
INTEGER IVERT                    ! Vertex index
REAL*8           XCUT,YCUT,ZCUT ! Cut point of three planes.
REAL*8           DET,DETX,DETY,DETZ ! Determinants of 3x3 system for 
!                               ! XCUT,YCUT,ZCUT.
REAL*8           DISTANCE   ! A distance criterium of two points in space
! The following are for sorting the vertices of each face:
REAL*8           V1(3),V2(3),V3(3)   ! Auxiliary vectors...
REAL*8           SINFIV1V2,COSFIV1V2 ! ...and their inner and outer products
REAL*8           FI(NVERTMAX)        ! ...and also their relative angles.

!LOGICAL HALFSPACE    ! Function used, see function itself.
LOGICAL LACCEPT      ! Determining whether a cut point is inside
!                          !                            the polyhedron.
!---------------------------------------------------------------
! Check & initialize

IF (NPLANE.LT.4 .and. output) &
WRITE(*,*) 'VERT3D: NPLANE was only',NPLANE

DO IPLANE = 1,NPLANE
   NVERT(IPLANE) = 0
ENDDO
!===============================================================
! Start loop over all planes that can be cut:
DO 120 IPLANE1 = 1,NPLANE
! Start loop over all other planes:
DO 110 IPLANE2 = 1,NPLANE
IF (IPLANE2.EQ.IPLANE1) GOTO 110 
! Start loop over all other-other (!) planes. Do from IPLANE2+1 to 
! NPLANE so that no pair is considered twice.
DO 100 IPLANE3 = IPLANE2+1,NPLANE        ! nikos  IPLANE2+1,NPLANE
IF (IPLANE3.EQ.IPLANE1) GOTO 100
!     IF (IPLANE3.EQ.IPLANE2) GOTO 100 ! added by nikos
! Solve the 3x3 system to find the cut point.
DET= A3(IPLANE1)*(B3(IPLANE2)*C3(IPLANE3)-B3(IPLANE3)*C3(IPLANE2)) &
   + A3(IPLANE2)*(B3(IPLANE3)*C3(IPLANE1)-B3(IPLANE1)*C3(IPLANE3)) &
   + A3(IPLANE3)*(B3(IPLANE1)*C3(IPLANE2)-B3(IPLANE2)*C3(IPLANE1))

!---------------------------------------------------------------
IF (DABS(DET).GT.1.D-12) THEN ! there is a cut point 

DETX=D3(IPLANE1)*(B3(IPLANE2)*C3(IPLANE3)-B3(IPLANE3)*C3(IPLANE2)) &
   + D3(IPLANE2)*(B3(IPLANE3)*C3(IPLANE1)-B3(IPLANE1)*C3(IPLANE3)) &
   + D3(IPLANE3)*(B3(IPLANE1)*C3(IPLANE2)-B3(IPLANE2)*C3(IPLANE1))

DETY=A3(IPLANE1)*(D3(IPLANE2)*C3(IPLANE3)-D3(IPLANE3)*C3(IPLANE2)) &
   + A3(IPLANE2)*(D3(IPLANE3)*C3(IPLANE1)-D3(IPLANE1)*C3(IPLANE3)) &
   + A3(IPLANE3)*(D3(IPLANE1)*C3(IPLANE2)-D3(IPLANE2)*C3(IPLANE1))

DETZ=A3(IPLANE1)*(B3(IPLANE2)*D3(IPLANE3)-B3(IPLANE3)*D3(IPLANE2)) &
   + A3(IPLANE2)*(B3(IPLANE3)*D3(IPLANE1)-B3(IPLANE1)*D3(IPLANE3)) &
   + A3(IPLANE3)*(B3(IPLANE1)*D3(IPLANE2)-B3(IPLANE2)*D3(IPLANE1))

XCUT = DETX/DET
YCUT = DETY/DET
ZCUT = DETZ/DET
!     write(6,333) IPLANE1,IPLANE2,IPLANE3,XCUT,YCUT,ZCUT
!333  format('Cutting point of planes ',3I5,':',3D15.7)
!-----------------------------------
! Accept this cut point as a vertex, if it belongs to the polyhedron. So,
! make a loop over all other (than IPLANE1,2,3) planes:
LACCEPT = .TRUE.
DO 50 KPLANE = 1,NPLANE
   IF (KPLANE.EQ.IPLANE1.OR.KPLANE.EQ.IPLANE2 &
                        .OR.KPLANE.EQ.IPLANE3) GOTO 50
   LACCEPT = LACCEPT.AND.HALFSPACE(A3(KPLANE),B3(KPLANE), &
                         C3(KPLANE),D3(KPLANE),XCUT,YCUT,ZCUT)
50 CONTINUE
!-----------------------------------
IF (LACCEPT) THEN
! If the cut point found belongs to the cell, we accept it unless it has
! occured before for this face (IPLANE1). Such a situation is possible
! when 4 or more planes pass through the same point (e.g. for the vertices
! of the fcc WS-cell). So...
   DO IVERT = 1,NVERT(IPLANE1)
      DISTANCE = DABS(XVERT(IVERT,IPLANE1)-XCUT) &
               + DABS(YVERT(IVERT,IPLANE1)-YCUT) &
               + DABS(ZVERT(IVERT,IPLANE1)-ZCUT)
!           IF (DISTANCE.LT.1.D-40 ) write(6,*) 'vertices plane reject'
      IF (DISTANCE.LT.1.D-10 ) LACCEPT = .FALSE.     ! accuracy problem
   ENDDO
ENDIF
! Now we're ready to add the point to the vertex list.
IF (LACCEPT) THEN
   NVERT(IPLANE1) = NVERT(IPLANE1) + 1
   XVERT(NVERT(IPLANE1),IPLANE1) = XCUT
   YVERT(NVERT(IPLANE1),IPLANE1) = YCUT
   ZVERT(NVERT(IPLANE1),IPLANE1) = ZCUT
ENDIF

ENDIF ! (DABS(DET).GT.1.D-12)
!---------------------------------------------------------------

100 CONTINUE       ! IPLANE3
110 CONTINUE       ! IPLANE2
!     write(6,*) 
!    & 'Number of vertices for plane ',iplane1,'  :',nvert(iplane1)
120 CONTINUE       ! IPLANE1

!===============================================================
! Each plane should finally have either at least 3 vertices, if it is a
! face of the polyhedron, or none at all. Check this:
DO IPLANE = 1,NPLANE
IF (NVERT(IPLANE).EQ.1.OR.NVERT(IPLANE).EQ.2) THEN
IF (output) THEN
WRITE(*,*) 'VERTEX3D: There might be a problem with the vertices.'
WRITE(*,*) 'For plane',IPLANE, &
  ' ,only ',NVERT(IPLANE),' vertices were found.'
ENDIF
ENDIF
ENDDO

!===============================================================
! For each face of the polyhedron, sort the vertices in a consecutive order
! as vertices of the polygon. The order is not necessarily mathematically
! positive.
NFACE = 0
DO IPLANE = 1,NPLANE
IF (NVERT(IPLANE).GE.3) THEN
   NFACE = NFACE + 1      ! Count the faces
   FI(1) = -4.D0          ! Just a number smaller than -pi.

! Define the vector connecting the first vertex to the (now-) second:
   V1(1) = XVERT(2,IPLANE) - XVERT(1,IPLANE)
   V1(2) = YVERT(2,IPLANE) - YVERT(1,IPLANE)
   V1(3) = ZVERT(2,IPLANE) - ZVERT(1,IPLANE)

   DO IVERT = 2,NVERT(IPLANE)
! Define the vector connecting the first vertex to the current one:
      V2(1) = XVERT(IVERT,IPLANE) - XVERT(1,IPLANE)
      V2(2) = YVERT(IVERT,IPLANE) - YVERT(1,IPLANE)
      V2(3) = ZVERT(IVERT,IPLANE) - ZVERT(1,IPLANE)
! Find the angle fi between v1 and v2 
! ( always, -pi < fi < pi from the definition of DATAN2 )
      COSFIV1V2 = V1(1)*V2(1) + V1(2)*V2(2) + V1(3)*V2(3)
      CALL CROSPR(V1,V2,V3) ! Cross product = |v1|*|v2|*sinfi
      SINFIV1V2 = DSQRT(V3(1)*V3(1) + V3(2)*V3(2) + V3(3)*V3(3))

      IF (SINFIV1V2.EQ.0.D0.AND.COSFIV1V2.EQ.0.D0) THEN
! Point falls exactly on 1st vertex...
         FI(IVERT) = -4.D0
         WRITE(*,*) 'VERTEX3D: Found two identical vertex points' 
! ...while it shouldn't ! (this was checked earlier)
      ELSE
         FI(IVERT) = DATAN2(SINFIV1V2,COSFIV1V2)
      ENDIF

   ENDDO                  ! IVERT = 3,NVERT(IPLANE)

! Store with respect to the angle found:
   CALL SORTVERTICES(NVERT(IPLANE),FI,XVERT(1,IPLANE), &
                   YVERT(1,IPLANE),ZVERT(1,IPLANE))

ENDIF                     ! (NVERT(IPLANE).GE.3)
ENDDO                     ! IPLANE = 1,NPLANE

RETURN
END SUBROUTINE












SUBROUTINE ANALYZEVERT3D( &
                          NVERTMAX,NFACED,TOLVDIST,TOLAREA,NPLANE, &
                          NFACE,NVERT,XVERT,YVERT,ZVERT, &
                          A3,B3,C3,D3,output)
! Analyze the faces and vertices of the polyhedron.
! Use criteria for rejecting faces that are too small 
! or vertices that are too close to each other.
! On output, number of faces and vertices may be reduced 
! after some rejections have taken place.
implicit none

logical :: output

! Input:
INTEGER NVERTMAX,NFACED
INTEGER NPLANE                 ! Number of planes
! Input and output:
INTEGER NFACE                  ! Number of faces (usually much less than nplane)
INTEGER NVERT(NFACED)          ! Number of vertices for each face
REAL*8  XVERT(NVERTMAX,NFACED),YVERT(NVERTMAX,NFACED) &
       ,ZVERT(NVERTMAX,NFACED)
REAL*8 TOLVDIST,TOLAREA  ! Max. tolerance for distance of two vertices and area of face.
REAL*8 A3(*),B3(*),C3(*),D3(*)  ! Coefs. defining the planes, to be reordered at end
! Inside
REAL*8 VDIST             ! Distance between consecutive vertices

LOGICAL LACCEPTVERT(NVERTMAX),LACCTOT,LFOUNDNEXT,LTHISISTHELAST
LOGICAL LACCEPTFACE(NFACED)
INTEGER NEWINDEXFACE(NFACED)

INTEGER IFACE,IVERT,IVERT2,INEXT,IPLANE
INTEGER IFACENEWCOUNT,IVERTNEWCOUNT
REAL*8 DX,DY,DZ,X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3,TRIANGLEAREA
REAL*8 FACEAREA(NFACED)


REAL*8 X4,Y4,Z4,DET,DETSUM

! First analyze vertices.
! Reject doubles (also vertices which fall almost on the previous vertex).
DO 100 IPLANE = 1,NPLANE
   IF (NVERT(IPLANE).EQ.0) GOTO 100  ! no need checking this one, proceed to next!

   LACCTOT = .TRUE.
! First vertices 1st with 2nd, 2nd with 3rd, etc...
   LACCEPTVERT(1) = .TRUE.                 ! First vertex is always accepted.
   IVERT = 1
   LTHISISTHELAST = .FALSE.   ! This flag will go up at the last acceptable vertex.

! Double loop: First loop is over all vertices; 
! but if during the loop vertices are found that have to be rejected, they are jumped over.
   DO WHILE (IVERT.LT.NVERT(IPLANE).AND..NOT.LTHISISTHELAST)
      LFOUNDNEXT = .FALSE.    ! This flag will become true when the next acceptable vertex is found.
      IVERT2 = IVERT + 1
! Second loop is over subsequent vertices (i.e., vertices ivert2 > ivert).
! Look for the first acceptable-as-next vertex, but do not go beyond last vertex.
! Stop loop as soon as the acceptable next vertex is found.
      DO WHILE (.NOT.LFOUNDNEXT.AND.IVERT2.LE.NVERT(IPLANE))
         DX = XVERT(IVERT2,IPLANE) -  XVERT(IVERT,IPLANE)
         DY = YVERT(IVERT2,IPLANE) -  YVERT(IVERT,IPLANE)
         DZ = ZVERT(IVERT2,IPLANE) -  ZVERT(IVERT,IPLANE)
         VDIST = DSQRT(DX*DX + DY*DY + DZ*DZ)

         IF (VDIST.GE.TOLVDIST) THEN
            LACCEPTVERT(IVERT2) = .TRUE.   ! Vertex is to be accepted
            INEXT = IVERT2                 ! Set this as the next vertex
            LFOUNDNEXT = .TRUE.            ! and we have a winner, exit loop.
         ELSE
            LACCEPTVERT(IVERT2) = .FALSE.  ! Remember that vertex is to be rejected later
            LACCTOT = .FALSE.              ! Remember that at least one vertex has to be rejected.
            IVERT2 = IVERT2 + 1            ! Now compare to the next vertex
        ENDIF  
      ENDDO

      IF (.NOT.LFOUNDNEXT) LTHISISTHELAST = .TRUE. ! If there is no next acceptable vertex,
                                                        ! then this was the last one. Jump out.
      IVERT = INEXT
      
   ENDDO

! ...and now 1st with last to close the cycle:
   IVERT = 1
   IVERT2 = NVERT(IPLANE)
   DX = XVERT(IVERT2,IPLANE) -  XVERT(IVERT,IPLANE)
   DY = YVERT(IVERT2,IPLANE) -  YVERT(IVERT,IPLANE)
   DZ = ZVERT(IVERT2,IPLANE) -  ZVERT(IVERT,IPLANE)
   VDIST = DSQRT(DX*DX + DY*DY + DZ*DZ)
   IF (VDIST.GE.TOLVDIST) THEN
      LACCEPTVERT(IVERT2) = .TRUE.        ! Vertex is to be accepted
   ELSE
      LACCEPTVERT(IVERT2) = .FALSE.       ! Remember that vertex is to be rejected later
      LACCTOT = .FALSE.                   ! Remember that at least one vertex has to be rejected.
   ENDIF


! Reject vertices which were found inappropriate and re-index vertices in each plane:
   IF (.NOT.LACCTOT) THEN
      IVERTNEWCOUNT = 0
      DO IVERT = 1,NVERT(IPLANE)
         IF (LACCEPTVERT(IVERT)) THEN
            IVERTNEWCOUNT = IVERTNEWCOUNT + 1    ! One more vertex to accept 
            IF (IVERTNEWCOUNT.NE.IVERT) THEN     ! Otherwise the correct value is already at the correct place
               XVERT(IVERTNEWCOUNT,IPLANE) = XVERT(IVERT,IPLANE) ! Re-index vertex
               YVERT(IVERTNEWCOUNT,IPLANE) = YVERT(IVERT,IPLANE)
               ZVERT(IVERTNEWCOUNT,IPLANE) = ZVERT(IVERT,IPLANE)
            ENDIF
         ENDIF
      ENDDO
      NVERT(IPLANE) = IVERTNEWCOUNT
   ENDIF

100 ENDDO



! Now analyze faces, reject faces with less than three vertices and faces of very small area.
DO 200 IPLANE = 1,NPLANE
   IF (NVERT(IPLANE).GE.3) THEN  ! calculate area
      X1 = XVERT(1,IPLANE)
      Y1 = YVERT(1,IPLANE)
      Z1 = ZVERT(1,IPLANE)
      FACEAREA(IPLANE) = 0.d0
      DO IVERT = 2,NVERT(IPLANE)-1
         X2 = XVERT(IVERT,IPLANE)
         Y2 = YVERT(IVERT,IPLANE)
         Z2 = ZVERT(IVERT,IPLANE)
         X3 = XVERT(IVERT+1,IPLANE)
         Y3 = YVERT(IVERT+1,IPLANE)
         Z3 = ZVERT(IVERT+1,IPLANE)
         TRIANGLEAREA = 0.5d0 * DABS( &
              (X2-X1)*(X3-X1)+(Y2-Y1)*(Y3-Y1)+(Z2-Z1)*(Z3-Z1) )
         FACEAREA(IPLANE) = FACEAREA(IPLANE)+ TRIANGLEAREA
      ENDDO

      IF (output) WRITE(*,8000) IPLANE,FACEAREA(IPLANE)

      IF (FACEAREA(IPLANE).GE.TOLAREA) THEN
         LACCEPTFACE(IPLANE) = .TRUE.
      ELSE
         LACCEPTFACE(IPLANE) = .FALSE.  ! Reject facees with small area
         IF (output) WRITE(*,8010) TOLAREA 
      ENDIF
      
   ELSE
      LACCEPTFACE(IPLANE) = .FALSE.  ! Reject planes with less than 3 vertices
      IF (output) WRITE(*,8020) IPLANE,NVERT(IPLANE)
   ENDIF

200 ENDDO


! Re-order the faces so that the accepted ones are in the first NFACE positions (NFACE is recalculated);
! The rest of the array entries are not taken care of, and can contain garbage.
IFACENEWCOUNT = 0
DO IPLANE = 1,NPLANE
   IF (LACCEPTFACE(IPLANE)) THEN
      IFACENEWCOUNT = IFACENEWCOUNT + 1  ! One more face to accept
      IF (IFACENEWCOUNT.NE.IPLANE) THEN   ! Otherwise the correct value is already at the correct place
         NVERT(IFACENEWCOUNT) = NVERT(IPLANE) ! Re-index face vertex number
         DO IVERT = 1,NVERT(IPLANE)
            XVERT(IVERT,IFACENEWCOUNT) = XVERT(IVERT,IPLANE) ! Re-index face vertices
            YVERT(IVERT,IFACENEWCOUNT) = YVERT(IVERT,IPLANE) ! Re-index face vertices
            ZVERT(IVERT,IFACENEWCOUNT) = ZVERT(IVERT,IPLANE) ! Re-index face vertices
         ENDDO
         A3(IFACENEWCOUNT) = A3(IPLANE) ! Re-index face equation parameters
         B3(IFACENEWCOUNT) = B3(IPLANE) ! Re-index face equation parameters
         C3(IFACENEWCOUNT) = C3(IPLANE) ! Re-index face equation parameters
         D3(IFACENEWCOUNT) = D3(IPLANE) ! Re-index face equation parameters
      ENDIF
   ENDIF
ENDDO
NFACE = IFACENEWCOUNT

! Check for every face that all veritces lie on the same plane
! by checking linear dependence
DO IFACE = 1,NFACE
   X2 = XVERT(2,IFACE) - XVERT(1,IFACE)
   Y2 = XVERT(2,IFACE) - YVERT(1,IFACE)
   Z2 = XVERT(2,IFACE) - ZVERT(1,IFACE)
   X3 = XVERT(3,IFACE) - XVERT(1,IFACE)
   Y3 = XVERT(3,IFACE) - YVERT(1,IFACE)
   Z3 = XVERT(3,IFACE) - ZVERT(1,IFACE)
   DETSUM = 0.d0
   DO IVERT = 4,NVERT(IFACE)
      X4 = XVERT(IVERT,IFACE) - XVERT(1,IFACE)
      Y4 = XVERT(IVERT,IFACE) - YVERT(1,IFACE)
      Z4 = XVERT(IVERT,IFACE) - ZVERT(1,IFACE)
      DET = X2*(Y3*Z4-Y4*Z3)+Y2*(Z3*X4-Z4*X3)+Z2*(X3*Y4-X4*Y3)
      DETSUM = DETSUM + DABS(DET)
      IF (DABS(DET).GT.1.D-16 .AND. output) THEN
        WRITE(*,9000) IFACE,IVERT,DET
      END IF
   ENDDO
   IF (output) WRITE(*,9010) IFACE,DETSUM
ENDDO


8000 FORMAT('ANALYZEVERT3D: Face',I5,' has area',E12.4)
8010 FORMAT('Face will be rejected ; Max. area tolerance=',E12.4)
8020 FORMAT('Plane',I5,' has only',I3,' vertices and is rejected')
9000 FORMAT('Error from ANALYZEVERT3D: Vertices not on single plane.', &
        ' IFACE=',I5,' IVERT=',I5,' DETERMINANT=',E12.4)
9010 FORMAT('ANALYZEVERT3D: Checking that vertices lie on plane.', &
     ' IFACE=',I5,' ; Determinants sum to be zero=',E12.4)

END SUBROUTINE

!***********************************************************************
LOGICAL FUNCTION HALFSPACE(A,B,C,D,X,Y,Z)
! Given a plane A*x+B*y+C*z=D, and a point (X,Y,Z) in space, this 
! function takes the value TRUE if (X,Y,Z) lies in the half-space 
! defined by the plane and the origin (0,0,0) (including the plane 
! itself). Else, the value FALSE is returned.
!
! The criterion used is that the inner product of the vector (X,Y,Z) 
! with the vector d connecting the origin to the plane vertically be 
! less than or equal to d**2:  (d_x,d_y,d_z)*(X,Y,Z) =< d**2.
!
! Input:
REAL*8           A,B,C,D,X,Y,Z,TEST,D1,D2

IF (DABS(A)+DABS(B)+DABS(C).LT.1.D-80) &
                                STOP 'HALFSPACE: A,B,C too small.'

HALFSPACE = .FALSE.

IF (D*(A*X+B*Y+C*Z).LE.D*D) HALFSPACE = .TRUE.
! (re-checked 31May2008 FM)

RETURN
END FUNCTION

!***********************************************************************
SUBROUTINE NORMALPLANE(X1,Y1,Z1,X2,Y2,Z2,TAU,A,B,C,D)
! Given two points in space, r1=(X1,Y1,Z1) and r2=(X2,Y2,Z2), this
! subroutine returns the coefficients defining a plane through the 
! equation A*x+B*y+C*z=D, which is normal to the vector r2-r1 and passes
! through the point (1.-TAU)*r1 + TAU*r2 (TAU thus being a parameter
! defining how close the plane is to each of the two points).
implicit none
! Input:
REAL*8           X1,Y1,Z1,X2,Y2,Z2,TAU
! Output:
REAL*8           A,B,C,D
! Inside:
REAL*8           ONEMTAU
! The plane is defined as 
! (A,B,C)*(X-X1,Y-Y1,Z-Z1)=const=
!                         =(distance from r1 to (1.-TAU)*r1 + TAU*r2)**2
! so A,B,C are the coords. of a vector connecting the point r1 to
! the point (1.-TAU)*r1 + TAU*r2.
ONEMTAU = 1.D0 - TAU

A = ONEMTAU * X1 + TAU * X2
B = ONEMTAU * Y1 + TAU * Y2
C = ONEMTAU * Z1 + TAU * Z2
D = A*(A+X1) + B*(B+Y1) + C*(C+Z1)

RETURN
END SUBROUTINE
!***********************************************************************
SUBROUTINE NORMALPLANE0( &
     X1,Y1,Z1,TAU, &
     A,B,C,D)
! Given a point in space, r1=(X1,Y1,Z1), this
! subroutine returns the coefficients defining a plane through the 
! equation A*x+B*y+C*z=D, which is normal to the vector r1 and passes
! through the point TAU*r1 (TAU thus being a parameter
! defining how close the plane is to the point).
implicit none
! Input:
REAL*8           X1,Y1,Z1,TAU
! Output:
REAL*8           A,B,C,D
! Inside:
! The plane is defined by 
! (A,B,C)*(X,Y,Z) = D = (tau * r1)**2
! so A,B,C are the coords. of the vector tau * r1.
! If tau=0 (plane passes through the origin), then D=0.

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
END SUBROUTINE

!***********************************************************************
SUBROUTINE SORTVERTICES(N,SINFI,X,Y,Z)
! Sorts the array SINFI(N) in ascending order using straight insertion. 
! The arrays Z(N), Y(N), and Z(N) follow.
! On output, arrays SINFI, X, Y, and Z return sorted.
implicit none
INTEGER N,I,J
REAL*8           SINFI(*),X(*),Y(*),Z(*),TMPS,TMPX,TMPY,TMPZ

DO J = 2,N
   TMPS = SINFI(J)
   TMPX = X(J)
   TMPY = Y(J)
   TMPZ = Z(J)
   DO I = J-1,1,-1
      IF (SINFI(I).LE.TMPS) GOTO 10
      SINFI(I+1) = SINFI(I)
      X(I+1) = X(I)
      Y(I+1) = Y(I)
      Z(I+1) = Z(I)
   ENDDO
   I = 0
10    SINFI(I+1) = TMPS
   X(I+1) = TMPX
   Y(I+1) = TMPY
   Z(I+1) = TMPZ
ENDDO

RETURN
END SUBROUTINE

! ************************************************************************
SUBROUTINE CROSPR(X,Y,Z)
! ************************************************************************
!     CROSP COMPUTES THE CROSS PRODUCT OF X AND Y RETURNING
!     IT INTO Z.
! ------------------------------------------------------------------------
REAL*8           X(*), Y(*), Z(*)
Z(1)=X(2)*Y(3)-X(3)*Y(2)
Z(2)=X(3)*Y(1)-X(1)*Y(3)
Z(3)=X(1)*Y(2)-X(2)*Y(1)
RETURN
END SUBROUTINE

!***********************************************************************
REAL*8           FUNCTION DISTPLANE(A,B,C,D)
! Returns the distance of a plane A*x+B*y+C*z=D to the origin.
implicit none
REAL*8           A,B,C,D
REAL*8           ABCSQ

ABCSQ = A*A + B*B + C*C

IF (ABCSQ.LT.1.D-100) STOP 'DISTPLANE'

DISTPLANE = DABS(D)/DSQRT(ABCSQ)  

RETURN
END FUNCTION

! ************************************************************************
SUBROUTINE DSORT (W,IND,MAX,POS)
implicit none
! ************************************************************************
!     p.zahn, april 96
!     W   is the original array returned unchanged
!     IND is an array that holds the new positions
!     max number of ellements to be sorted
!     pos the position where the first element is found
! ------------------------------------------------------------------------
INTEGER MAX,POS
REAL*8            W(*)
REAL*8            BOUND, DIFF
INTEGER IND(*)

INTEGER I,II,J,JJ,K
DATA BOUND /1.0D-12/
! ------------------------------------------------------------------------
DO 10 I = 1,MAX
  IND(I) = I
10 END DO

J = MAX
J = 1
DO 60 WHILE (J.LT.MAX/3)
  J = 3*J+1
60 END DO

DO 20 WHILE (J.GT.1)
  J = J/3
  JJ = 1
  DO 30 WHILE (JJ.EQ.1)
    JJ = 0
    DO 40 K=1,MAX-J
      DIFF = ABS( W(IND(K)) - W(IND(K+J)) )
      IF ( W(IND(K)) .GT. W(IND(K+J)) .AND. &
           DIFF.GT.BOUND ) THEN
        II       = IND(K)
        IND(K)   = IND(K+J)
        IND(K+J) = II
        JJ = 1
      END IF
40     END DO                    ! K=1,MAX-J
30   END DO                      ! WHILE (JJ.EQ.1)
20 END DO

DO 50 I=1,MAX
  IF (IND(I) .EQ. 1) POS=I
50 END DO

RETURN
END SUBROUTINE

END MODULE
