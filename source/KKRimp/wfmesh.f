!------------------------------------------------------------------------------------
!> Summary: Creation of the radial mesh for the old-solver 
!> Author: 
!> Creation of the radial mesh for the old-solver 
!------------------------------------------------------------------------------------
      MODULE mod_WFMESH
      CONTAINS
  !-------------------------------------------------------------------------------
  !> Summary: Creation of the radial mesh for the old-solver 
  !> Author: 
  !> Category: radial-grid, KKRimp
  !> Deprecated: False 
  !> Creation of the radial mesh for the old-solver 
  !-------------------------------------------------------------------------------      
      SUBROUTINE WFMESH(E,EK,CVLIGHT,NSRA,Z,R,S,RS,NRMAX,LMAXD)
        USE NRTYPE
        IMPLICIT NONE
!     Interface arguments
      COMPLEX(kind=DPC),intent(in)       ::  E
      COMPLEX(kind=DPC),intent(out)      ::  EK
      INTEGER,intent(in)                 ::  NSRA
      REAL(kind=DP),intent(in)           ::  CVLIGHT
      REAL(kind=DP),intent(in)           ::  Z
      REAL(kind=DP),intent(in)           ::  R(NRMAX)
      REAL(kind=DP),intent(out)          ::  RS(NRMAX,0:LMAXD)
      REAL(kind=DP),intent(out)          ::  S(0:LMAXD)
      INTEGER,intent(in)                 ::  NRMAX
      INTEGER,intent(in)                 ::  LMAXD
!     Local arguments
      DOUBLE PRECISION S1
      INTEGER IR,L

!     .. Intrinsic Functions ..
      INTRINSIC DBLE,SQRT

!     ..
      IF (NSRA.EQ.1) EK = SQRT(E)
      IF (NSRA.EQ.2) EK = SQRT(E+E*E/ (CVLIGHT*CVLIGHT))
      DO L = 0,LMAXD

        IF (NSRA.EQ.2) THEN
          S1 = SQRT(DBLE(L*L+L+1)-4.0D0*Z*Z/ (CVLIGHT*CVLIGHT))
          IF (Z.EQ.0.0D0) S1 = DBLE(L)
        ELSE
          S1 = DBLE(L)
        END IF
        S(L) = S1
        RS(1,L) = 0.0D0
        DO IR = 2,NRMAX
          RS(IR,L) = R(IR)**S1
        END DO
!       DO IR = NRMAX+1, NRMAXD
!          RS(IR,L) = 0.0D0
!       END DO

      END DO
      RETURN
      END SUBROUTINE
      END MODULE mod_WFMESH
