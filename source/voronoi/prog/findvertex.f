      SUBROUTINE findvertex(NPLANE,A3,B3,C3,D3,NEDGEMAX,
     &                              NFACE,NEDGE,XEDGE,YEDGE,ZEDGE,
     &                              NVERTEX)
c Given a set of planes, defined by A3*x+B3*y+C3*z=D3 and defining
c a convex part of space (the minimal one containing the origin, 
c usually a WS-polyhedron), this subroutine returns the edges
c of this polyhedron in cartesian coordinates. For the planes that
c are not faces of the polyhedron, a value NEDGE(IPLANE)=0 is returned.
c The total no. of faces found is returned as NFACE.
c
c Updated on 3.9.2001 Algorithm is changed
c Present Algorithm:
c   Define a cube around the  point.
c   Define: POI2PLANE  (Point belongs to planes... )
c           POI2POI    (Point is connected with other points...)
c           
c   Now suppose a new canditate plane, first check if all points are
c   on the correct side of space (same side as origin)
c   find how many are in wrong side, and look at all connections
c   between points in wrong side to points in correct side (use POI2POI)
c
c   Define new points and also arrays POI2POI,PO2PLANE, for these points
c   Bookkeeping...
c   If you found 3 new points then accept the plane and the new points  
c   a few points may already exist, this is taken care later!
c   goto next plane..
c
c Uses logical function HALFSPACE
      implicit none
      integer npoimax,nneimax,nplanemax
      parameter (npoimax=300,nneimax=100,nplanemax=100)
c -------------------------------------------------------
c Input:
      INTEGER NPLANE                ! Number of planes. (changed on out)
      INTEGER NEDGEMAX              ! Max. number of edges per plane.
      REAL*8           A3(*),B3(*),C3(*),D3(*)  ! Coefs. defining the planes, 
c                                               ! dimensioned >= NPLANE.
c                                               ! changed on output  
c Output:
      INTEGER NEDGE(*)  ! Number of edges found for each face
      INTEGER NFACE     ! Number of faces found (with nedge>0).
      integer NVERTEX   ! number of vertices 
      REAL*8           XEDGE(NEDGEMAX,*),YEDGE(NEDGEMAX,*),
     &                 ZEDGE(NEDGEMAX,*)
      

c Local 
      real*8 alfa
      integer npolypoi,npolyplan
      integer poi2poi(0:nneimax,npoimax),poi2plane(0:nplanemax,npoimax)
      integer poiofplane(nplanemax,npoimax)
      real*8 polyplaneA3(nplanemax),polyplaneB3(nplanemax),
     &     polyplaneC3(nplanemax),polyplaneD3(nplanemax) 
      real*8 polypoints(3,npoimax)
      integer iv,i1,i,ii,ntrueedge,itest,k,npoi
      real*8 r,dot,r1,r2,cosfi
c -------------------------------------------------------

      integer index(nedgemax,nplanemax),order(nedgemax)
c                            ! Cartesian coords. of edges for each plane
c                            ! (2nd index is for planes).
c Inside:
      INTEGER IPLANE1,IPLANE2,IPLANE3,IPLANE,KPLANE ! Plane indices
      INTEGER IEDGE           ! Edge index
      REAL*8           XCUT,YCUT,ZCUT     ! Cut point of three planes.
      REAL*8           DET,DETX,DETY,DETZ ! Determinants of 3x3 system for 
c                               ! XCUT,YCUT,ZCUT.
      REAL*8           DISTANCE   ! A distance criterium of two points in space
c The following are for sorting the edges of each face:
      REAL*8           V1(3),V2(3),V3(3)   ! Auxiliary vectors...
      REAL*8           SINFIV1V2,COSFIV1V2 ! ...and their inner and outer products
      REAL*8           FI(NEDGEMAX)        ! ...and also their relative angles.
      real*8 small
      LOGICAL HALFSPACE    ! Function used, see function itself.
      LOGICAL LACCEPT      ! Determining whether a cut point is inside
c                          !                            the polyhedron.
      logical lkeep(nedgemax),ltaken(npoimax)
      data small/1.d-5/
c---------------------------------------------------------------9
c Check & initialize
      
      IF (NPLANE.LT.4) WRITE(*,*) 'EDGE3D: NPLANE was only',NPLANE
      DO IPLANE = 1,NPLANE
         NEDGE(IPLANE) = 0
      ENDDO 
      do iplane=1,nplanemax
         do i1=1,npoimax          
            poiofplane(iplane,i1) = 0
         end do
      end do
c
c Define bounding cube this is the initial polyhedron
c     
      alfa = 30.d0
      CALL DEFCUBE(alfa,npolypoi,npolyplan,polypoints,poi2poi,poi2plane,
     &     polyplaneA3,polyplaneB3,polyplaneC3,polyplaneD3)
      
c
c now cut original cube with all planes availiable
c       
      do iplane = 1, nplane
         
         CALL CUTPLANE(npolypoi,npolyplan,polypoints,poi2poi,poi2plane,
     &        polyplaneA3,polyplaneB3,polyplaneC3,polyplaneD3,
     &        A3(iplane),B3(iplane),C3(iplane),D3(iplane),NEDGE,
     &        POIOFPLANE)
      end do
c
c canditate plane a3,b3,c3,d3 is checked if it belongs to polyhedron :
c  1. Is copied in polyplaneA3, polyplaneB3 ...
c     this includes also the original cube planes this is checked
c     later. 
c  2. The arrays poi2poi, poi2plane are updated this is used to find
c     all points in a face.
c
c Now all info is availiable, just prepare arrays neaded.

      do i=1,6
        if (nedge(i).gt.0) write(6,*) 
     &        ' C A U T I O N : Polyhedron maybe open!'
      end do
c the first 6 planes are the original cube...

      do iplane =1,nplane
         do iv = 1,nedge(iplane)
            
            index(iv,iplane) = poiofplane(iplane,iv)
            XEDGE(iv,iplane) = polypoints(1,poiofplane(iplane,iv))  
            YEDGE(iv,iplane) = polypoints(2,poiofplane(iplane,iv)) 
            ZEDGE(iv,iplane) = polypoints(3,poiofplane(iplane,iv))
            
c     write(6,fmt='(2I5,3F10.5)') iplane,iv,XEDGE(iv,iplane),
c     &           YEDGE(iv,iplane),
c     &           ZEDGE(iv,iplane) 
         end do
      end do    
c
c Now go on to find vertices of each face
c

      nface = 0
      do iplane=1,nplane 
         IF (nedge(iplane).ge.3) THEN
            nface = nface + 1
c     
c     now order using previous info about neighbours
c     
            do i=1,npoimax             ! flag all atoms exept atoms in plane
               ltaken(i) = .true.
            end do 
            do i=1,nedge(iplane)
               ltaken(index(i,iplane)) = .false.
               order(i) = 0 
            end do
            ! write(6,*) 'next plane',nedge(iplane)
            k = 1 
            order(k) = index(1,iplane)
            ltaken(order(k)) = .true. 
            npoi = order(k)
            k = k + 1

            
            do while (k.le.nedge(iplane))
               
               if (k.ge.4) ltaken(order(1)) = .false. ! allow to return back
               i=0 
               do while (i.lt.poi2poi(0,npoi))
                  i = i + 1
                  itest = poi2poi(i,npoi)
                  if (.not.ltaken(itest)) then
                     ltaken(itest) = .true.
                     order(k) = itest
                     i = poi2poi(0,npoi) ! go out of loop
                  end if
               end do
               if (order(k).eq.0) then
                  write(6,*) ' Could not find neighbours to connect '
                  write(6,*) ' Edge not closing  STOPING ...        '
                  write(6,*) ' atom :', itest,k,iplane
                  stop
               end if
               !  write(6,*) 'ORDERING ',k,order(k)
               npoi = order(k)
               k = k + 1  
            end do              ! do while (k.le.nedge(iplane))

            do iedge = 1,nedge(iplane)
               iv = order(iedge)
               XEDGE(iedge,iplane) = polypoints(1,iv)  
               YEDGE(iedge,iplane) = polypoints(2,iv) 
               ZEDGE(iedge,iplane) = polypoints(3,iv)
            end do
         ENDIF                  ! (NEDGE(IPLANE).GE.3)
c
c Now ordering is done 
c


c         write(76,*),iplane,a3(iplane),b3(iplane),c3(iplane),d3(iplane)
c          do iedge =1,nedge(iplane)
c          write(76,*) xedge(iedge,iplane),yedge(iedge,iplane),
c     &            zedge(iedge,iplane)
c          end do
      ENDDO                     ! IPLANE = 1,NPLANE
c
cc Now reject vertices that occur more than once
c
      do 20 iplane =1,nplane
         do iedge=1,nedge(iplane)
            lkeep(iedge) = .true.
         end do
         
         do iedge=1,nedge(iplane)
            do itest=iedge+1,nedge(iplane)
               V1(1) = XEDGE(itest,IPLANE) - XEDGE(iedge,IPLANE)
               V1(2) = YEDGE(itest,IPLANE) - YEDGE(iedge,IPLANE)
               V1(3) = ZEDGE(itest,IPLANE) - ZEDGE(iedge,IPLANE)
               r = sqrt(v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3))
               if (r.lt.small) then
                  lkeep(itest) = .false.
                  write(6,*) 'Vertex ',itest,' in plane ',iplane,
     &                 ' found again'
               end if
            end do
         end do
c    now get rid of double points
         ntrueedge= 0
         do iedge=1,nedge(iplane)
            if (lkeep(iedge))  ntrueedge = ntrueedge + 1  
         end do
         if (ntrueedge.eq.nedge(iplane)) goto 20 ! next plane
 
         iedge = 1
         do while (iedge.le.ntrueedge) 
         ! do iedge = 1,ntrueedge
         !   write(6,*) iedge,lkeep(iedge)
            if (.not.lkeep(iedge)) then ! promote by one all verices
               do itest=iedge+1,nedge(iplane)
                  XEDGE(itest-1,IPLANE) = XEDGE(itest,IPLANE)
                  YEDGE(itest-1,IPLANE) = YEDGE(itest,IPLANE)
                  ZEDGE(itest-1,IPLANE) = ZEDGE(itest,IPLANE)
                  index(itest-1,iplane) = index(itest,iplane)
                  lkeep(itest-1) = lkeep(itest)
               end do
            else
               iedge = iedge + 1
            end if 
         end do   
 
          nedge(iplane) = ntrueedge
 
 20   CONTINUE

c
c Now update the plane coordinates
c
      nplane = npolyplan
      do iplane =1,nplane
         A3(iplane) = polyplaneA3(iplane)
         B3(iplane) = polyplaneB3(iplane)
         C3(iplane) = polyplaneC3(iplane)
         D3(iplane) = polyplaneD3(iplane)
      end do
      NVERTEX = npolypoi

      END 


















