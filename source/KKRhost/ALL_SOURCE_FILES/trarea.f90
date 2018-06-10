SUBROUTINE trarea(a,b,lmax)
! from complex to real  (differenciated spherical harmonics)

!.. Parameters ..
DOUBLE PRECISION RTWO
DOUBLE COMPLEX CI
PARAMETER (RTWO=1.414213562373d0,CI= (0.d0,1.d0))
!..
!.. Scalar Arguments ..
INTEGER LMAX
!..
!.. Array Arguments ..
DOUBLE COMPLEX A(*)
DOUBLE PRECISION B(*)
!..
!.. Local Scalars ..
DOUBLE PRECISION SGM
INTEGER I,L,M
!..
!.. Intrinsic Functions ..
INTRINSIC CONJG,DBLE

!    calculate real the spherical harmonics derivetived
i = 0
DO  l = 0,lmax
  i = i + l + 1
  b(i) = DBLE(a(i))
  sgm = -1.d0
  DO  m = 1,l
    b(i-m) = DBLE(ci* (a(i-m)-CONJG(a(i-m))))/rtwo
    b(i+m) = sgm*DBLE((a(i+m)+CONJG(a(i+m))))/rtwo
    sgm = -sgm
  END DO
  i = i + l
END DO
RETURN
END SUBROUTINE trarea
