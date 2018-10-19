  !-------------------------------------------------------------------------------
  !> Summary: From complex to real (differenciated spherical harmonics)
  !> Author: 
  !> Category: special-functions, numerical-tools, KKRimp
  !> Deprecated: True 
  !> From complex to real (differenciated spherical harmonics)
  !-------------------------------------------------------------------------------
      SUBROUTINE TRAREA(A,B,LMAX)
c
C     .. Parameters ..
      DOUBLE PRECISION RTWO
      DOUBLE COMPLEX CI
      PARAMETER (RTWO=1.414213562373d0,CI= (0.d0,1.d0))
C     ..
C     .. Scalar Arguments ..
      INTEGER LMAX
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX A(*)
      DOUBLE PRECISION B(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION SGM
      INTEGER I,L,M
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC CONJG,DBLE
C     ..
C
C    calculate real the spherical harmonics derivetived
C
      I = 0
      DO 20 L = 0,LMAX
        I = I + L + 1
        B(I) = DBLE(A(I))
        SGM = -1.D0
        DO 10 M = 1,L
          B(I-M) = DBLE(CI* (A(I-M)-CONJG(A(I-M))))/RTWO
          B(I+M) = SGM*DBLE((A(I+M)+CONJG(A(I+M))))/RTWO
          SGM = -SGM
   10   CONTINUE
        I = I + L
   20 CONTINUE
      RETURN
      END SUBROUTINE
