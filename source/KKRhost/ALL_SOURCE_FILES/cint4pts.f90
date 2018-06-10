SUBROUTINE cint4pts(y,jtop,z)
!   ********************************************************************
!   *                                                                  *
!   *      perform the integral  Z(i)   =  INT   Y(i') di'             *
!   *                                    R=0..R(i)                     *
!   *                                                                  *
!   *      via a 4-point integration formula                           *
!   *                                                                  *
!   *      JTOP:     Y is tabulated form 1 .. JTOP                     *
!   *      Y(i):     function to be integrated                         *
!   *                                                                  *
!   *                       COMPLEX - VERSION                          *
!   *                                                                  *
!   ********************************************************************
IMPLICIT NONE

! Dummy arguments
INTEGER JTOP
COMPLEX*16 Y(JTOP),Z(JTOP)

! Local variables
INTEGER I,IG,J,K,M,N1,N2
REAL*8 Q(5,5),Q5(5,5)
COMPLEX*16 S,SVN

DATA q5/0.d0,251.d0,232.d0,243.d0,224.d0,0.d0,646.d0,992.d0,  &
    918.d0,1024.d0,0.d0, - 264.d0,192.d0,648.d0,384.d0,0.d0,  &
    106.d0,32.d0,378.d0,1024.d0,0.d0, - 19.d0, - 8.d0, - 27.d0, 224.d0/

DO i = 1,5
  DO j = 1,5
    q(i,j) = q5(i,j)/720.0D0
  END DO
END DO

z(1) = DCMPLX(0.d0,0.d0)
svn = z(1)

DO ig = 1,jtop - 4,4
  n1 = ig
  n2 = ig + 4
  DO m = n1 + 1,n2
    i = m - n1 + 1
    s = svn
    DO k = n1,n2
      j = k - n1 + 1
      s = s + q(i,j)*y(k)
    END DO
    z(m) = s
  END DO
  svn = z(n2)
END DO

IF ( n2 /= jtop ) THEN
  n1 = jtop - 4
  n2 = jtop
  svn = z(n1)
  DO m = n1 + 1,n2
    i = m - n1 + 1
    s = svn
    DO k = n1,n2
      j = k - n1 + 1
      s = s + q(i,j)*y(k)
    END DO
    z(m) = s
  END DO
END IF

END SUBROUTINE cint4pts
