      SUBROUTINE FINDPANELS(
     >     NFACE,A3,B3,C3,D3,NVERT,XVERT,YVERT,ZVERT,TOLEULER,TOLVDIST,
     <     NPAN,CRT)
      implicit none
c This subroutine finds the critical points and panels for a
c suitable radial mesh for a Voronoi polyhedron. It calls the 
c same subroutine CRIT that is called by subr. shape.
      INCLUDE 'inc.geometry'

c Input:
      INTEGER NFACE,NVERT(NFACED) ! Number of faces and vertices per face
      REAL*8 A3(NFACED),B3(NFACED),C3(NFACED),D3(NFACED) 
      ! Parameters for eq. defining face: A3*x+B3*y+C3*z=D3
      REAL*8 XVERT(NVERTD,NFACED),YVERT(NVERTD,NFACED)
      REAL*8 ZVERT(NVERTD,NFACED)  ! xyz coordinates of vertices per face
      REAL*8 TOLEULER,TOLVDIST
c Output:
      INTEGER NPAN ! Number of panels for polyhedron
      REAL*8 CRT(NPAND) ! Critical points (panel end-points)
c Local:
      REAL*8 FACE(3) ! Face eq: FACE(1)*x + FACE(2)*y + FACE(3)*z = 1
      REAL*8 VRT(3,NVERTD)  ! Vertices
      REAL*8 C1 ! For sorting
      INTEGER NVRT ! Number of vertices
      INTEGER IVTOT ! Total number of vertices
      INTEGER IORD ! For sorting
      INTEGER IFACE,IP,IVERT
      

      NPAN = 0
      IVTOT = 0
      DO IFACE = 1,NFACE
         FACE(1) = A3(IFACE)/D3(IFACE)
         FACE(2) = B3(IFACE)/D3(IFACE)
         FACE(3) = C3(IFACE)/D3(IFACE)
         NVRT = NVERT(IFACE)
         DO IVERT = 1,NVRT
            VRT(1,IVERT) = XVERT(IVERT,IFACE)
            VRT(2,IVERT) = YVERT(IVERT,IFACE)
            VRT(3,IVERT) = ZVERT(IVERT,IFACE)
         ENDDO
         ! NPAN and CRT are updated at every call of CRIT (new points added)
         CALL CRIT(IFACE,NVRT,VRT,FACE,NPAN,IVTOT,TOLEULER,TOLVDIST,CRT)
      ENDDO

! Sort critical points
      DO 10 IORD = 1,NPAN
         C1 = CRT(IORD)
         DO 20 IP = NPAN,IORD,-1
            IF(CRT(IP).GT.C1)   GO TO 20
            C1 = CRT(IP)
            CRT(IP) = CRT(IORD)
            CRT(IORD) = C1
 20      CONTINUE
 10   CONTINUE


      END
