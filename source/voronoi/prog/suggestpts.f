      SUBROUTINE SUGGESTPTS(
     >     NPAN,CRIT,DENPT,ALAT,NMIN,
     <     NMESH)
      implicit none
c#@# KKRtags: VORONOI radial-grid geometry
      INCLUDE 'inc.geometry'
c Suggests number of radial points between MT radius and outer
c radius based on the number of panels
c
c Input:
      INTEGER NPAN ! Number of panels
      REAL*8 CRIT(NPAND) ! Critical points (end-radius of each panel)
      REAL*8 DENPT  ! "Optimal" density of radial points, e.g. 100/a_Bohr
      REAL*8 ALAT   ! Lattice constant
      INTEGER NMIN  ! Minimum number of points per panel
c Output:
      INTEGER NMESH ! Suggested number of points
c Local:
      INTEGER IP,N1
      REAL*8 WIDTH ! Panel width

      NMESH = 0

c First panel starts at 0 and ends at MT.
      DO IP = 2,NPAN
         WIDTH = (CRIT(IP) - CRIT(IP-1)) * ALAT
         N1 = NINT(DENPT*WIDTH)
         IF (N1.LT.NMIN) N1 = NMIN
         NMESH = NMESH + N1
      ENDDO

      ! The following avoids problems in routine mesh0
      ! in case that all panels have NMIN points.
      ! (Use NPAN-1 instead of NPAN because panel 1 is always 
      !  the inner part from the nucleus to the muffin tin)
      IF (NMESH.LE.(NPAN-1)*NMIN) NMESH = (NPAN-1)*NMIN + 1

      END
