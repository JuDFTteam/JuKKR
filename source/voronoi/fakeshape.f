      SUBROUTINE FAKESHAPE(
     >     VOLUME,RMT,ALAT,NPOI,
     <     NPAN,MESHN,NM,NFUN,XRN,DRN,LMIFUN,THETAS)
      implicit none
c#@# KKRtags: VORONOI atomic-sphere-approx shape-functions  
! Produces fake shape functions for usage in ASA-mode

! Input:
      INTEGER NPOI      ! Number of points in the shape function
      REAL*8 ALAT       ! Latt. constant
      REAL*8 RMT,VOLUME ! Muffin-tin radius in a.u.,volume in units of alat**3
! Output:
      INTEGER NPAN              ! Number of panels (here: one panel only)
      INTEGER MESHN             ! Number of mesh points (here: equal to npoi)
      INTEGER NM(*)             ! Number of mesh points per panel (here: equal to npoi)
      INTEGER NFUN              ! Number of different lm components for lm > 1 (here:zero)
      INTEGER LMIFUN(*)         ! lm-component of each function (here: lmifun(1)=1)
      REAL*8 XRN(*),DRN(*)      ! Mesh points and weights (uniform mesh)
      REAL*8 THETAS(*)          ! Shape functions (only lm=1 component, others are zero)
! Local:
      INTEGER IR
      REAL*8  FPI,RFPI,RWS,WT,RMTALAT

      FPI = 16.D0*DATAN(1.D0)
      RFPI = DSQRT(FPI) 

      RMTALAT = RMT / ALAT                    ! Muffin-tin radius in units of alat
      RWS = (3.D0*VOLUME/FPI)**(1.D0/3.D0)    ! Wigner-Seitz radius in units of alat

      MESHN = INT( 100.D0 * (RWS*ALAT - RMT) )  ! Rule of thumb: 100 points / a.u.
      MESHN = MIN(MESHN,NPOI)
      IF (MOD(MESHN,2).EQ.0) MESHN = MESHN - 1


      DO IR = 1,MESHN
         XRN(IR) = (RMTALAT * DFLOAT(MESHN-IR) + RWS * DFLOAT(IR-1)) / 
     &               DFLOAT(MESHN - 1)
      ENDDO

      WT = (RWS - RMTALAT) / DFLOAT(MESHN - 1)         ! Integration weight

      DO IR = 1,MESHN
         DRN(IR) = WT
      ENDDO

      NPAN = 1
      NM(1) = MESHN
      NFUN = 1
      LMIFUN(1) = 1
      THETAS(1:MESHN) = RFPI


      END
