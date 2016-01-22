      SUBROUTINE ANALYZEVERT3D(
     >                          NVERTMAX,NFACED,TOLVDIST,TOLAREA,NPLANE,
     X                          NFACE,NVERT,XVERT,YVERT,ZVERT,
     X                          A3,B3,C3,D3)
c Analyze the faces and vertices of the polyhedron.
c Use criteria for rejecting faces that are too small 
c or vertices that are too close to each other.
c On output, number of faces and vertices may be reduced 
c after some rejections have taken place.
      implicit none
c Input:
      INTEGER NVERTMAX,NFACED
      INTEGER NPLANE                 ! Number of planes
c Input and output:
      INTEGER NFACE                  ! Number of faces (usually much less than nplane)
      INTEGER NVERT(NFACED)          ! Number of vertices for each face
      REAL*8  XVERT(NVERTMAX,NFACED),YVERT(NVERTMAX,NFACED)
     &       ,ZVERT(NVERTMAX,NFACED)
      REAL*8 TOLVDIST,TOLAREA  ! Max. tolerance for distance of two vertices and area of face.
      REAL*8 A3(*),B3(*),C3(*),D3(*)  ! Coefs. defining the planes, to be reordered at end
c Inside
      REAL*8 VDIST             ! Distance between consecutive vertices

      LOGICAL LACCEPTVERT(NVERTMAX),LACCTOT,LFOUNDNEXT,LTHISISTHELAST
      LOGICAL LACCEPTFACE(NFACED)
      INTEGER NEWINDEXFACE(NFACED)

      INTEGER IFACE,IVERT,IVERT2,INEXT,IPLANE
      INTEGER IFACENEWCOUNT,IVERTNEWCOUNT
      REAL*8 DX,DY,DZ,X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3,TRIANGLEAREA
      REAL*8 FACEAREA(NFACED)
      LOGICAL TEST ! test options


      REAL*8 X4,Y4,Z4,DET,DETSUM


c First analyze vertices.
c Reject doubles (also vertices which fall almost on the previous vertex).
      DO 100 IPLANE = 1,NPLANE
         IF (NVERT(IPLANE).EQ.0) GOTO 100  ! no need checking this one, proceed to next!

         LACCTOT = .TRUE.
c First vertices 1st with 2nd, 2nd with 3rd, etc...
         LACCEPTVERT(1) = .TRUE.                 ! First vertex is always accepted.
         IVERT = 1
         LTHISISTHELAST = .FALSE.   ! This flag will go up at the last acceptable vertex.

c Double loop: First loop is over all vertices; 
c but if during the loop vertices are found that have to be rejected, they are jumped over.
         DO WHILE (IVERT.LT.NVERT(IPLANE).AND..NOT.LTHISISTHELAST)
            LFOUNDNEXT = .FALSE.    ! This flag will become true when the next acceptable vertex is found.
            IVERT2 = IVERT + 1
c Second loop is over subsequent vertices (i.e., vertices ivert2 > ivert).
c Look for the first acceptable-as-next vertex, but do not go beyond last vertex.
c Stop loop as soon as the acceptable next vertex is found.
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

c ...and now 1st with last to close the cycle:
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


c Reject vertices which were found inappropriate and re-index vertices in each plane:
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

 100  ENDDO



c Now analyze faces, reject faces with less than three vertices and faces of very small area.
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
               TRIANGLEAREA = 0.5d0 * DABS(
     &              (X2-X1)*(X3-X1)+(Y2-Y1)*(Y3-Y1)+(Z2-Z1)*(Z3-Z1) )
               FACEAREA(IPLANE) = FACEAREA(IPLANE)+ TRIANGLEAREA
            ENDDO

            IF (TEST('verb0   ')) WRITE(*,8000) IPLANE,FACEAREA(IPLANE)

            IF (FACEAREA(IPLANE).GE.TOLAREA) THEN
               LACCEPTFACE(IPLANE) = .TRUE.
            ELSE
               LACCEPTFACE(IPLANE) = .FALSE.  ! Reject faces with small area
               WRITE(*,8010) TOLAREA 
            ENDIF
            
         ELSE
            LACCEPTFACE(IPLANE) = .FALSE.  ! Reject planes with less than 3 vertices
            IF (TEST('verb0   ')) WRITE(*,8020) IPLANE,NVERT(IPLANE)
         ENDIF

 200  ENDDO


c Re-order the faces so that the accepted ones are in the first NFACE positions (NFACE is recalculated);
c The rest of the array entries are not taken care of, and can contain garbage.
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
 
c Check for every face that all veritces lie on the same plane
c by checking linear dependence
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
            IF (DABS(DET).GT.1.D-16) WRITE(*,9000) IFACE,IVERT,DET
         ENDDO
         IF (TEST('verb1   ')) WRITE(*,9010) IFACE,DETSUM
      ENDDO
      

 8000 FORMAT('ANALYZEVERT3D: Face',I5,' has area',E12.4)
 8010 FORMAT('Face will be rejected ; Max. area tolerance=',E12.4)
 8020 FORMAT('Plane',I5,' has only',I3,' vertices and is rejected')
 9000 FORMAT('Error from ANALYZEVERT3D: Vertices not on single plane.',
     &        ' IFACE=',I5,' IVERT=',I5,' DETERMINANT=',E12.4,
     &        ' should be zero.')
 9010 FORMAT('ANALYZEVERT3D: Checking that vertices lie on plane.',
     &     ' IFACE=',I5,' ; Determinants sum to be zero=',E12.4)

      END

