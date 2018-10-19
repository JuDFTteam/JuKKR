  !-------------------------------------------------------------------------------
  !> Summary: Generate an angular mesh and spherical harmonics at those mesh points. For an angular integration the weights are generated .
  !> Author: R. Zeller
  !> Date: February 1996
  !> Category: special-functions, radial-mesh, KKRimp
  !> Deprecated: False 
  !> Generate an angular mesh and spherical harmonics at those
  !> mesh points. For an angular integration the weights are generated  
  !-------------------------------------------------------------------------------
      SUBROUTINE SPHERE_NOGGA(LMAX,YR,WTYR,RIJ,IJD)
C     ..
C     .. Scalar Arguments ..
      USE MOD_YMY
      IMPLICIT NONE
      INTEGER IJD,LMAX
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION PI,R,R1,R2,R3
      INTEGER IJ,LM1
C     ..
C     .. External Subroutines ..
!       EXTERNAL YMY
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION RIJ(IJD,3),WTYR(IJD,(LMAX+1)**2), 
     +                 YR(IJD,(LMAX+1)**2)
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION WGHT,Y(1000)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN
C     ..
      PI = 4.D0*ATAN(1.D0)
      WRITE (1337,*) ' SPHERE : read LEBEDEV mesh'
      IF (IJD.GT.1000) STOP ' SPHERE '
c
c
      DO 30 IJ = 1,IJD
        CALL LEBEDEV(IJ,R1,R2,R3,WGHT)
        RIJ(IJ,1) = R1
        RIJ(IJ,2) = R2
        RIJ(IJ,3) = R3
        CALL YMY(R1,R2,R3,R,Y,LMAX)
        DO 10 LM1 = 1, (LMAX+1)**2
          YR(IJ,LM1) = Y(LM1)
   10   CONTINUE
c
c--->   multiply the spherical harmonics with the weights
c
        DO 20 LM1 = 1, (LMAX+1)**2
          WTYR(IJ,LM1) = YR(IJ,LM1)*WGHT*PI*4.D0
   20   CONTINUE
   30 CONTINUE

      END SUBROUTINE
