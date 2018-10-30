c***********************************************************************
      SUBROUTINE VORONOI12(
     >  NVEC,RVEC,NVERTMAX,NFACED,WEIGHT0,WEIGHT,TOLVDIST,TOLAREA,TOLHS,
     <  RMT,ROUT,VOLUME,NFACE,A3,B3,C3,D3,NVERT,XVERT,YVERT,ZVERT)
c Given a cluster of atomic positions at RVEC(3,NVEC), this subroutine
c returns information about the Voronoi cell around the origin. It is
c supposed, of course, that the origin corresponds to an atomic position
c which is not included in the RVEC(3,N). The information returned is:
c
c RMT: Muffin-tin radius (radius of inscribed sphere centered at the 
c      origin).
c
c ROUT: Radius of circumscribed sphere centered at the origin.
c
c VOLUME: Volume of the Voronoi cell. 
c
c NFACE: Number of faces of the cell.
c
c A3(NFACE),B3(NFACE),C3(NFACE),D3(NFACE): Coefficients defining the
c faces via A3*x+B3*y+C3*z=D3. The arrays are filled in the first
c NFACE positions with the information, the rest can be garbage.
c
c NVERT(NFACE): Number of vertices corresponting to each face.
c
c XVERT(NVERTMAX,NFACE),YVERT(NVERTMAX,NFACE),ZVERT(NVERTMAX,NFACE): 
c Coordinates of theese vertices.
c
c The Voronoi construction performed here allows for different than
c 50%/50% bisections. For this, each atomic site, positioned say at
c vector r(i), is assigned a weight w(i). Then a point r in space
c belongs to the cell i if, for all j, 
c |r-r(i)|**2 - w(i) < |r-r(j)|**2 - w(j).    (1)
c The points for which the inequality becomes an equality is the 
c bisector, and it can be easily shown that it is a plane perpendicular
c to the vector r(j)-r(i).
c One can parametrize the segment connecting r(i) and r(j) using a
c parameter t, 0 < t < 1. For t=0 we are at r(i), for t=1 at r(j), and
c in-between we have a linear dependence of the distance from r(i) with
c t. I.e., the position vector is r = r(i)*(1-t) + r(j)*t.
c The point where the bisector cuts this segment will then be at
c t=(1/2)*( (w(i)-w(j))/dist**2 + 1 )         (2)
c where dist is the distance of the two points. As a special case we see
c that, for w(i)=w(j), the bisector will be in the middle, else the
c atomic site with the bigger weight gains more space.
c
c The above procedure is guaranteed to divide space into tesselating 
c (:=space-filling) convex polyhedra (one of their sides could be at 
c infinity in special cases). The convexity of the polyhedra is 
c guaranteed by the constuction, i.e. as mutual cut of half-spaces, and
c the property of tesselation is guaranteed by the fact that every point
c in space is assigned to some atom.
c
c However, this procedure might place the polyhedron in such a way that
c its own atomic site is not contained in it. This can happen if there
c is a big enough difference in the weights, whence, from eq.(2) above,
c t can become <0 or >1 (then the bisector has passed one of the two
c points). This cannot be allowed in a KKR calculation, and therefore
c each time it is checked that this is not the case. If such a case
c occurs, a different set of weights must be chosen.
c
c
c Uses subroutines NORMALPLANE, POLYHEDRON, and function DISTPLANE.
      implicit none
c#@# KKRtags: VORONOI geometry
c Input:
      INTEGER  NVEC,NVERTMAX,NFACED
      REAL*8   RVEC(3,NFACED)
      REAL*8   WEIGHT0,WEIGHT(NFACED)  ! Weight of central atom, 
c                          !   and of all others (dimensioned as RVEC). 
      REAL*8  TOLVDIST              ! Max. tolerance for distance of two vertices
      REAL*8  TOLAREA               ! Max. tolerance for area of polygon face
      REAL*8  TOLHS                 ! Tolerance for halfspace routine
c Output:
      INTEGER NFACE
      INTEGER NVERT(*)
      REAL*8              VOLUME
      REAL*8              A3(NFACED),B3(NFACED),C3(NFACED),D3(NFACED)
      REAL*8              XVERT(NVERTMAX,NFACED),YVERT(NVERTMAX,NFACED),
     &                    ZVERT(NVERTMAX,NFACED)
      REAL*8              RMT,ROUT
c Inside:
      INTEGER IVEC,IFACE,IVERT,I,NVERTTOT
      INTEGER NPLANE
      REAL*8           X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,TETRVOL,RSQ,TAU
      REAL*8           FACEAREA(NFACED),TRIANGLEAREA   
      REAL*8           TEMP
      REAL*8           DISTPLANE          ! Function used.
      INTEGER ISORT(NFACED),POS,INEW    ! Index for sorting
      REAL*8 RSORT(4,NFACED)              ! Aux. function for sorting
      REAL*8 ATMP(NFACED),BTMP(NFACED),CTMP(NFACED),DTMP(NFACED) ! For sorting
      REAL*8 XVTMP(NVERTMAX,NFACED),YVTMP(NVERTMAX,NFACED)       ! For sorting
      REAL*8 ZVTMP(NVERTMAX,NFACED)                              ! For sorting
      INTEGER NVTMP(NFACED)                                      ! For sorting
      INTEGER NPANEL
      LOGICAL TEST  ! test options

c---------------------------------------------------------------
c Check that the origin is not included in RVEC.
      DO IVEC = 1,NVEC
         IF (DABS(RVEC(1,IVEC))+DABS(RVEC(2,IVEC))+DABS(RVEC(3,IVEC))
     &                                               .LT.1.D-12) THEN
            WRITE(*,*) 'VORONOI: Vector',IVEC,'is zero.'
            STOP 'VORONOI'
         ENDIF
      ENDDO
c---------------------------------------------------------------
c Define the planes as normal to the vectors RVEC, passing from t as
c in eq. (2) above:
      DO IVEC = 1,NVEC
         RSQ = RVEC(1,IVEC) * RVEC(1,IVEC)
     &       + RVEC(2,IVEC) * RVEC(2,IVEC)
     &       + RVEC(3,IVEC) * RVEC(3,IVEC)
         TAU = 0.5D0*((WEIGHT0-WEIGHT(IVEC))/RSQ + 1.D0)
         CALL NORMALPLANE0(
     >                    RVEC(1,IVEC),RVEC(2,IVEC),RVEC(3,IVEC),TAU,
     <                    A3(IVEC),B3(IVEC),C3(IVEC),D3(IVEC))
      ENDDO
c---------------------------------------------------------------
c Find the Voronoi polyhedron.
      NPLANE = NVEC
      IF (TEST('verb1   ')) WRITE(*,*) 'Entering POLYHEDRON08'
      CALL POLYHEDRON08(
     >                  NPLANE,NVERTMAX,NFACED,TOLVDIST,TOLAREA,
     X                  TOLHS,A3,B3,C3,D3,
     <                  NFACE,NVERT,XVERT,YVERT,ZVERT)
      IF (TEST('verb1   ')) WRITE(*,*) 'Exited POLYHEDRON08'

c---------------------------------------------------------------
c Calculate the volume as sum of the volumes of all tetrahedra 
c connecting the origin to the faces. Use for each tetrahedron
c volume = det((r0-r1),(r0-r2),(r0-r3))/6, where r0 is here the
c origin and r1,r2,r3 the vectors of the 3 other vertices.
c Algorithm requires that the face is a convex polygon with the 
c vertices ordered (doesn't matter if they are clock- or anticlockwise).
      VOLUME = 0.d0
      DO IFACE = 1,NFACE
         X1 = XVERT(1,IFACE)
         Y1 = YVERT(1,IFACE)
         Z1 = ZVERT(1,IFACE)
         FACEAREA(IFACE) = 0.d0
         DO IVERT = 2,NVERT(IFACE)-1
            X2 = XVERT(IVERT,IFACE)
            Y2 = YVERT(IVERT,IFACE)
            Z2 = ZVERT(IVERT,IFACE)
            X3 = XVERT(IVERT+1,IFACE)
            Y3 = YVERT(IVERT+1,IFACE)
            Z3 = ZVERT(IVERT+1,IFACE)
            TETRVOL = X1*(Y2*Z3-Y3*Z2)+X2*(Y3*Z1-Y1*Z3)+X3*(Y1*Z2-Y2*Z1)
            VOLUME = VOLUME + DABS(TETRVOL)
            TRIANGLEAREA = 0.5d0 * DSQRT(
     &           ( X1*Y2 + X2*Y3 + X3*Y1 - Y2*X3 - Y3*X1 - Y1*X2)**2
     &         + ( Y1*Z2 + Y2*Z3 + Y3*Z1 - Z2*Y3 - Z3*Y1 - Z1*Y2)**2
     &         + ( Z1*X2 + Z2*X3 + Z3*X1 - X2*Z3 - X3*Z1 - X1*Z2)**2 )
! wrong           TRIANGLEAREA = 0.5d0 * DABS(
!     &           (X2-X1)*(X3-X1)+(Y2-Y1)*(Y3-Y1)+(Z2-Z1)*(Z3-Z1) )
            FACEAREA(IFACE) = FACEAREA(IFACE)+ TRIANGLEAREA
         ENDDO
      ENDDO
      VOLUME = VOLUME/6.D0

      WRITE(6,*) ' Polyhedron properties '
      WRITE(6,*) ' Number of faces : ',nface

      IF (TEST('verb0   ')) THEN
      DO IFACE=1,NFACE
         WRITE(6,201) IFACE,NVERT(IFACE),FACEAREA(IFACE)
 201     FORMAT(' Face ',I4,'   has ',I4,' vertices ','; Area= ',E12.4)  
         DO IVERT=1,NVERT(IFACE)
            WRITE(*,9000) IVERT,XVERT(IVERT,IFACE),YVERT(IVERT,IFACE)
     &                ,ZVERT(IVERT,IFACE)
         ENDDO
         WRITE(*,9010) A3(IFACE),B3(IFACE),C3(IFACE),D3(IFACE)
      END DO 
      ENDIF

      WRITE(6,*) 'The Volume is : ',VOLUME

 9000 FORMAT(I5,4E16.8)
 9010 FORMAT(' Face coefficients:',4E16.8)
c---------------------------------------------------------------
c Find RMT:
      RMT = DISTPLANE(A3(1),B3(1),C3(1),D3(1))
      DO IFACE = 2,NFACE
         TEMP = DISTPLANE(A3(IFACE),B3(IFACE),C3(IFACE),D3(IFACE))
         IF (TEMP.LT.RMT) RMT = TEMP
      ENDDO
c Fint ROUT:
      ROUT = 0.D0
      DO IFACE = 1,NFACE
         DO IVERT = 1,NVERT(IFACE)
            TEMP = XVERT(IVERT,IFACE)*XVERT(IVERT,IFACE) 
     &           + YVERT(IVERT,IFACE)*YVERT(IVERT,IFACE)
     &           + ZVERT(IVERT,IFACE)*ZVERT(IVERT,IFACE)
            IF (TEMP.GT.ROUT) ROUT = TEMP
         ENDDO
      ENDDO
      ROUT = DSQRT(ROUT)
      WRITE(*,9020) RMT,ROUT,RMT*100/ROUT
 9020 FORMAT('Voronoi subroutine: RMT=',E16.8,'; ROUT=',E16.8,
     &       '; RATIO=',F12.2,' %')
c---------------------------------------------------------------
c Sort faces:
      DO IFACE = 1,NFACE
         RSORT(1,IFACE) = 
     &        DISTPLANE(A3(IFACE),B3(IFACE),C3(IFACE),D3(IFACE))
         RSORT(2,IFACE) = C3(IFACE) 
         RSORT(3,IFACE) = B3(IFACE) 
         RSORT(4,IFACE) = A3(IFACE) 
      ENDDO
!      CALL DSORT(RSORT,ISORT,NFACE,POS)

      CALL DSORT_NCOMP(RSORT,4,.FALSE.,ISORT,NFACE,POS)
c     Rearrange using a temporary array
      ATMP(:) = A3(:)
      BTMP(:) = B3(:)
      CTMP(:) = C3(:)
      DTMP(:) = D3(:)
      XVTMP(:,:) = XVERT(:,:)
      YVTMP(:,:) = YVERT(:,:)
      ZVTMP(:,:) = ZVERT(:,:)
      NVTMP(1:NFACE) = NVERT(1:NFACE)
      DO IFACE = 1,NFACE
         INEW = ISORT(IFACE)
         A3(INEW) = ATMP(IFACE)
         B3(INEW) = BTMP(IFACE)
         C3(INEW) = CTMP(IFACE)
         D3(INEW) = DTMP(IFACE)
         XVERT(:,INEW) = XVTMP(:,IFACE)
         YVERT(:,INEW) = YVTMP(:,IFACE)
         ZVERT(:,INEW) = ZVTMP(:,IFACE)
         NVERT(INEW) = NVTMP(IFACE)
      ENDDO
c---------------------------------------------------------------


      RETURN
      END











