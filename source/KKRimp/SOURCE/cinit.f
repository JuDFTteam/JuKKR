!-------------------------------------------------------------------------------
!> Summary: Initialize array elements to complex zero
!> Author: 
!> Deprecated: False 
!-------------------------------------------------------------------------------
!> Setting the first N values of a double complex array A to zero
!>
!> @note
!> This routine could be easily replaced by: a(:n) = czero
!> @endnote
!-------------------------------------------------------------------------------
      MODULE mod_cinit
      CONTAINS
!-------------------------------------------------------------------------------
!> Summary: Initialize array elements to complex zero
!> Author: 
!> Category: KKRimp, numerical-tools 
!> Deprecated: False 
!-------------------------------------------------------------------------------
      SUBROUTINE CINIT(N,A)
        USE NRTYPE
        IMPLICIT NONE
C **********************************************************************
C * Setting the first N values of a double complex array A to zero     *
C **********************************************************************
C     ..
C     .. Arguments ..
      INTEGER,intent(in)              ::  N
      complex(kind=DPC),intent(inout) ::  A(:)
C     .. Locals
      INTEGER                         ::  I,M,MP1
      complex(kind=DPC),parameter     ::  CZERO=(0.0D0,0.0D0)
C     ..
      M = MOD(N,5)
      IF ( M.NE.0 ) THEN
         DO I = 1,M
            A(I) = CZERO
         END DO
         IF ( N.LT.5 ) RETURN
      END IF
      MP1 = M + 1
      DO I = MP1,N,5
        A(I  ) = CZERO
        A(I+1) = CZERO
        A(I+2) = CZERO
        A(I+3) = CZERO
        A(I+4) = CZERO
      END DO
C
      END SUBROUTINE CINIT
      END MODULE mod_cinit
