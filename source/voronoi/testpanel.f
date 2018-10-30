      subroutine testpanel(
     >     nface,a3,b3,c3,d3,nvert,xvert,yvert,zvert,
     <     npanel)
      implicit none
c#@# KKRtags: VORONOI radial-grid
      include 'inc.geometry'
c Find radial panels of polyhedron
c
c Not perfected yet
      
c Input:
      integer nface,nvert(nfaced)
      real*8 a3(nfaced),b3(nfaced),c3(nfaced),d3(nfaced)
      real*8 xvert(nvertd,nfaced),yvert(nvertd,nfaced),
     &       zvert(nvertd,nfaced)
c Output:
      integer npanel
c Local:
      integer maxcritd
      parameter (maxcritd=1000)
      real*8 rcrit(maxcritd)
      real*8 rsort(maxcritd) ! sorting
      integer isort(maxcritd),ipos ! sorting
      real*8 vec(3,nvertd),vert(3,nvertd),distsq,r1,r2,dist
      integer iface,ivert,iedge,ncrit,nvert0,icrit,ncritfinal

      real*8 tol
      real*8 distplane ! function used

      tol = 1.d-6

      ncrit = 0

! Faces
      do iface = 1,nface
         ncrit = ncrit + 1
         rsort(ncrit) = distplane(
     &                     a3(iface),b3(iface),c3(iface),d3(iface) )
      enddo

! Edges
      do iface = 1,nface
         nvert0 = nvert(iface)
         vec(1,1) = xvert(1,iface) - xvert(nvert0,iface)
         vec(2,1) = yvert(1,iface) - yvert(nvert0,iface)
         vec(3,1) = zvert(1,iface) - zvert(nvert0,iface)

         do ivert = 2,nvert0
            vec(1,ivert) = xvert(ivert,iface) - xvert(ivert-1,iface)
            vec(2,ivert) = yvert(ivert,iface) - yvert(ivert-1,iface)
            vec(3,ivert) = zvert(ivert,iface) - zvert(ivert-1,iface)
         enddo

         do ivert = 1,nvert0
            distsq = xvert(ivert,iface)**2 + 
     &               yvert(ivert,iface)**2 +
     &               zvert(ivert,iface)**2 -
     &             ( xvert(ivert,iface) * vec(1,ivert) +
     &               yvert(ivert,iface) * vec(2,ivert) +
     &               zvert(ivert,iface) * vec(3,ivert)  )**2
            ncrit = ncrit + 1
            rsort(ncrit) = dsqrt(distsq)
         enddo
      enddo

! Vertices
      do iface = 1,nface
         do ivert = 1,nvert(iface)
            distsq = xvert(ivert,iface)**2 + 
     &               yvert(ivert,iface)**2 + 
     &               zvert(ivert,iface)**2 
            ncrit = ncrit + 1
            rsort(ncrit) = dsqrt(distsq)
         enddo
      enddo

      call dsort(rsort,isort,ncrit,ipos)

      ncritfinal = 1
      rcrit(1) = rsort(isort(1))
      do icrit = 2,ncrit
         r1 = rsort(isort(icrit))
         r2 = rcrit(ncritfinal)
         dist = r1 - r2
         if (dist.ge.tol) then
            ncritfinal = ncritfinal + 1
            rcrit(ncritfinal) = r1
         endif
      enddo

      do icrit = 1,ncritfinal
         write(*,*) rcrit(icrit)
      enddo
      npanel = ncritfinal - 1

      end
